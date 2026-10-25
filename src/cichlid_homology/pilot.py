"""The pilot: rooting, relationships, dating and the comparison of TreeFam and Compara for one
focal species. Genes are sampled once and the same genes are called in both methods. Reads only;
writes to a scratch folder."""
from __future__ import annotations

import collections
import csv
import json
import random
import time
import traceback
from dataclasses import asdict
from pathlib import Path

from . import __version__
from .calls import REFS, call_gene
from .consensus import label_agreement
from .dating import Params
from .master import load_species_tree
from .names import attach_names, identity_agreement, load_names
from .orthology import RTree, relationship
from .read import load_family, read_tsv
from .rooting import root_tree
from .species import SHADOW_RULE, drop_shadowed, load_leaf_species
from .symbols import attach_symbol_check, matched, own_symbols
from .tree import leaves, parse_newick, prune

ROOT_RULE = ("reconcile with the master species tree (fewest duplications, then fewest losses); "
             "never midpoint")
PARTS = ("TreeFam", "Compara")

BASE_COLS = ["gene_id", "method", "family", "family_mzeb_genes", "tree_mode", "sampled", "placement_status", "lwr",
             "assign_flags", "fusion_partner", "rooting", "root_cost", "root_losses", "n_tied_roots",
             "n_tied_dup_only", "roots_checked", "leaves_after_prune", "shadowed_dropped",
             "parent_children", "species_in_tree"]
REF_COLS = ["relationship", "n_ref", "n_focal", "closest", "closest_is_ortholog",
            "n_ref_in_family", "closest_tie", "stable_across_roots", "paralog_level", "paralog_level_overlap",
            "dup_level", "dup_level_overlap", "dup_ref_levels_sides", "dup_class", "dup_resolved", "dup_levels_sides", "dup_outgroup_single",
            "dup_outgroup_found"]
TEL_COLS = ["teleost_tel_species", "teleost_tel_one", "teleost_tel_two_plus", "teleost_tel_zero",
            "teleost_tel_frac_one_of_present", "teleost_tel_frac_one_of_all", "teleost_tel_coverage"]
NAME_COLS = ["closest_symbol", "closest_id", "closest_name_status", "closest_name_reliable",
             "name_confidence", "orthologs_ids", "orthologs_symbols", "orthologs_unreliable"]
B_COLS = (["b_anchor", "b_label", "b_label_raw", "expanded_family", "b_pattern", "b_other_state",
           "b_cichlid_state", "b_clade_leaves", "b_teleost_species_in_tree",
           "b_other_n", "b_other_with2p", "b_other_with3p", "b_other_share2p", "b_other_share3p",
           "b_cichlid_with2p", "b_cichlid_share2p", "b_cichlid_share3p", "b_cichlid_extra_species"]
          + [f"b_{lv}_{k}" for lv in ("teleost", "percomorph", "ovalentaria", "cichlid")
             for k in ("n", "1", "2", "3p")])
SYMBOL_COLS = ["own_symbol", "human_symbol_match", "zebrafish_symbol_match"]
COLUMNS = BASE_COLS + B_COLS + SYMBOL_COLS + [f"{k}_{c}" for k in REFS for c in NAME_COLS] + [f"{k}_{c}" for k in REFS for c in REF_COLS] + TEL_COLS
SIZE_BINS = {"1": (1, 1), "2-3": (2, 3), "4-9": (4, 9), "10+": (10, 10 ** 9)}
EXPANDED = 10     # n_ref or n_focal at least this large = an expanded family, reported separately


def prepare(fam, ls, shadowed, st, prm):
    """Parse, drop shadowed/unknown leaves, prune, root. Returns (rooting, rtrees, dropped, unknown)
    or None when nothing is left."""
    root = parse_newick(fam.tree)
    sp_of = {x: ls.species(x) for x in fam.leaves}
    _, dropped = drop_shadowed(sp_of, shadowed, ls.is_project_gene)
    unknown = {x for x, sp in sp_of.items() if sp is None}
    root = prune(root, set(dropped) | unknown)
    if root is None:
        return None
    for x in leaves(root):
        x.sp = sp_of[x.name]
    rooting = root_tree(root, st, fam.rooted, k=prm.k_roots)
    return rooting, [RTree(r) for r in rooting.roots], dropped, unknown


def process_family(fam, ls, shadowed, st, clades, prm, own_prefix, prefixes, alt_t, meta):
    t0 = time.time()
    prep = prepare(fam, ls, shadowed, st, prm)
    if prep is None:
        return [], {"error": "empty after pruning"}
    rooting, rts, dropped, unknown = prep
    by_name = [{x.name: x for x in rt.leaves} for rt in rts]
    rows, sens = [], collections.Counter()
    own_names = sorted(n for n in by_name[0] if n.startswith(own_prefix + "_"))
    for name in own_names:
        g = by_name[0][name]
        own = g.sp
        row = call_gene(rts[0], g, own, prm, clades, fam.tree_mode, prefixes)
        for key in REFS:
            rels = {row[f"{key}_relationship"]}
            for rt, bn in zip(rts[1:], by_name[1:]):
                rels.add(relationship(rt, bn[name], REFS[key], own, prm.overlap_t)["rel"])
            row[f"{key}_stable_across_roots"] = len(rels) == 1
        for t in alt_t:
            for key, ref in REFS.items():
                if relationship(rts[0], g, ref, own, t)["rel"] != row[f"{key}_relationship"]:
                    sens[(t, key)] += 1
        row["expanded_family"] = len(own_names) >= prm.expanded_min
        row["b_label_raw"] = row.get("b_label", "")
        if row["expanded_family"] and row.get("b_label") in ("partial_retention", "multi_copy"):
            row["b_label"] = "expanded_family"
        row.update({"gene_id": name, "family": fam.family_id, "tree_mode": fam.tree_mode,
                    "family_mzeb_genes": len(own_names),
                    "rooting": rooting.how, "root_cost": rooting.cost, "root_losses": rooting.losses,
                    "n_tied_roots": rooting.n_tied, "n_tied_dup_only": rooting.n_tied_dup,
                    "roots_checked": len(rts), "leaves_after_prune": len(rts[0].leaves),
                    "shadowed_dropped": len(dropped)})
        rows.append(row)
    info = {"leaves": len(rts[0].leaves), "n_tied": rooting.n_tied, "n_tied_dup": rooting.n_tied_dup,
            "cost": rooting.cost, "how": rooting.how, "dropped": len(dropped),
            "unknown": len(unknown), "genes": len(rows),
            "seconds": round(time.time() - t0, 3), "sens": sens}
    return rows, info


def choose_genes(candidates, n, seed=1):
    return sorted(random.Random(seed).sample(sorted(candidates), min(n, len(candidates))))


def run_pilot(layout, cat, clades, species, n, prm: Params, alt_t, out_dir: Path,
              n_denovo: int = 100) -> dict:
    prefixes = {s.prefix for s in cat.species(release="v1")}
    st = load_species_tree()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = {"cichlid_homology": __version__, "species": species,
               "parameters": {**asdict(prm), "alt_overlap_t": list(alt_t), "n_genes": n,
                              "n_denovo_genes_per_method": n_denovo, "shadow_rule": SHADOW_RULE,
                              "rooting_rule": ROOT_RULE, "mouse": "out of scope for release 1"}}
    res = {p: layout.results("v1", p) for p in PARTS}
    assign = {p: {r["gene_id"]: r for r in read_tsv(res[p].assign(species).local)} for p in PARTS}
    froot = {p: res[p].family("X").local.parent for p in PARTS}
    gp = {p: {r["gene_id"]: r for r in read_tsv(Path(res[p].native().local) / "summary" /
                                                 "gene_placement.tsv")
              if r["gene_id"].startswith(species + "_")} for p in PARTS}

    def folder_ok(p, fam):
        return fam not in ("", "NoHit") and (froot[p] / fam / "done.json").exists()

    # the same genes in both methods: genes that have a finished family in each
    both = [g for g in assign["TreeFam"] if g in assign["Compara"]
            and folder_ok("TreeFam", assign["TreeFam"][g]["best_family"])
            and folder_ok("Compara", assign["Compara"][g]["best_family"])]
    sampled = set(choose_genes(both, n))
    summary["genes_with_a_family_in_both_methods"] = len(both)
    summary["genes_sampled"] = len(sampled)
    all_rows = {}
    try:
        own = own_symbols(layout, species)
    except Exception:
        own = {}
        summary["own_symbols_error"] = traceback.format_exc()[-800:]
    for part in PARTS:
        t0 = time.time()
        gset = set(sampled)
        dn_all = [g for g, r in gp[part].items() if r["mode"] == "denovo"
                  and folder_ok(part, r["family"])]
        gset |= set(choose_genes(dn_all, n_denovo, seed=2))
        fam_ids = sorted({assign[part][g]["best_family"] for g in gset if g in assign[part]
                          and folder_ok(part, assign[part][g]["best_family"])})
        fams = {}
        for f in fam_ids:
            try:
                fams[f] = load_family(froot[part] / f, prefixes, f)
            except Exception:
                fams[f] = None
        ls = load_leaf_species(layout, part, only={x for fm in fams.values() if fm for x in fm.leaves})
        shadowed = set(ls.shadowed())
        rows_part, errors, per_family = [], [], []
        sens = collections.Counter()
        for f, fam in fams.items():
            if fam is None:
                errors.append({"family": f, "error": "unreadable family folder"})
                continue
            try:
                rows, info = process_family(fam, ls, shadowed, st, clades, prm, species,
                                            prefixes, alt_t, None)
            except Exception:
                errors.append({"family": f, "error": traceback.format_exc()[-1500:]})
                continue
            for r in rows:
                a, p = assign[part].get(r["gene_id"], {}), gp[part].get(r["gene_id"], {})
                r.update(method=part, sampled=r["gene_id"] in gset,
                         placement_status=p.get("status", ""), lwr=p.get("lwr", ""),
                         assign_flags=a.get("flags", ""), fusion_partner=a.get("fusion_partner", ""))
            rows_part.extend(rows)
            sens.update(info.pop("sens", {}))
            per_family.append({"family": f, **info})
        names_error = None
        try:
            attach_names(rows_part, load_names(layout, part))
            attach_symbol_check(rows_part, own)
        except Exception:
            names_error = traceback.format_exc()[-1200:]
        all_rows[part] = rows_part
        with open(out_dir / f"pilot_genes.{part}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore", restval="")
            w.writeheader()
            w.writerows(rows_part)
        summary[part] = _summarize(rows_part, per_family, errors, sens, shadowed, time.time() - t0)
        summary[part]["reference_names"] = names_error or _names_summary(rows_part)
        summary[part]["symbol_check"] = names_error or _symbol_summary(rows_part)
    summary["agreement"] = _agreement(all_rows, sampled, out_dir)
    (out_dir / "pilot_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    return summary


def _summarize(rows, fams, errors, sens, shadowed, seconds):
    C = collections.Counter
    out = {"families_ok": len(fams), "families_errors": len(errors), "first_errors": errors[:3],
           "genes_called": len(rows), "genes_sampled_flag": sum(1 for r in rows if r["sampled"]),
           "tree_mode_genes": dict(C(r["tree_mode"] for r in rows)),
           "shadowed": sorted(shadowed), "seconds": round(seconds, 1),
           "rooting_how": dict(C(f["how"] for f in fams)),
           "root_dups_min_med_max": _mmm([f["cost"] for f in fams if f["cost"] is not None]),
           "families_by_n_tied_roots_(dup,loss)": dict(C(min(f["n_tied"], 5) for f in fams)),
           "families_by_n_tied_roots_(dup_only)": dict(C(min(f["n_tied_dup"], 5) for f in fams)),
           "family_seconds_max": max((f["seconds"] for f in fams), default=None),
           "unknown_species_leaves": sum(f["unknown"] for f in fams)}
    for key in REFS:
        exp = [r for r in rows if _int(r[f"{key}_n_ref"]) >= EXPANDED or _int(r[f"{key}_n_focal"]) >= EXPANDED]
        norm = [r for r in rows if r not in exp]
        out[key] = {
            "relationship_all": dict(C(r[f"{key}_relationship"] for r in rows)),
            "relationship_excluding_expanded": dict(C(r[f"{key}_relationship"] for r in norm)),
            "genes_in_expanded_families": len(exp),
            "stable_across_roots": _frac(sum(1 for r in rows if r[f"{key}_stable_across_roots"]), len(rows)),
            "closest_tie": sum(1 for r in rows if r[f"{key}_closest_tie"]),
            "paralog_level": dict(C(r.get(f"{key}_paralog_level") for r in rows if r.get(f"{key}_paralog_level"))),
            "dup_class": dict(C(r.get(f"{key}_dup_class") for r in rows if r.get(f"{key}_dup_class"))),
            "dup_unresolved": sum(1 for r in rows if r.get(f"{key}_dup_resolved") is False),
            "threshold_changes": {str(t): v for (t, k), v in sorted(sens.items()) if k == key},
            "relationship_by_family_size": {
                lab: dict(C(r[f"{key}_relationship"] for r in rows if lo <= _int(r["family_mzeb_genes"]) <= hi))
                for lab, (lo, hi) in SIZE_BINS.items()},
            "paralog_level_overlap": dict(C(r.get(f"{key}_paralog_level_overlap") or "none_shared"
                                            for r in rows if r.get(f"{key}_paralog_level"))),
            "paralog_level_vs_overlap_level": [[list(k), v] for k, v in C(
                (r[f"{key}_paralog_level"], r.get(f"{key}_paralog_level_overlap") or "none_shared")
                for r in rows if r.get(f"{key}_paralog_level")).most_common(12)],
            "ohnolog_ref_levels_weaker_side": dict(C(
                min(map(int, r[f"{key}_dup_ref_levels_sides"].split("/")))
                for r in rows if r.get(f"{key}_dup_class", "").startswith("teleost_ohnolog")))}
    one = [r for r in rows if r["human_relationship"] == "1:1" and r.get("teleost_tel_frac_one_of_all") is not None]
    allf = [r["teleost_tel_frac_one_of_all"] for r in one]
    cov = [r for r in one if (r["teleost_tel_coverage"] or 0) >= 0.5]
    out["teleost_pattern_human_1to1"] = {
        "n": len(one), "n_with_coverage_ge_0.5": len(cov),
        "frac_one_of_all_deciles": _deciles(allf),
        **{f"share_ge_{t}_of_all": _frac(sum(1 for x in allf if x >= t), len(allf)) for t in (0.7, 0.8, 0.9)},
        **{f"share_ge_{t}_of_all_given_coverage": _frac(
            sum(1 for r in cov if r["teleost_tel_frac_one_of_all"] >= t), len(cov)) for t in (0.7, 0.8, 0.9)},
        "share_any_species_with_2plus": _frac(sum(1 for r in one if int(r["teleost_tel_two_plus"]) > 0), len(one))}
    out["goal_B"] = _goal_b(rows)
    out["low_support_(lwr<0.5)_human_relationships"] = dict(C(
        r["human_relationship"] for r in rows if _float(r["lwr"]) is not None and _float(r["lwr"]) < 0.5))
    return out


NO_B = (None, "", "no_orthologous_anchor", "insufficient_species", "expanded_family")


def _goal_b(rows):
    C = collections.Counter
    scored = [r for r in rows if r.get("b_label_raw") not in NO_B and r.get("b_other_share2p") not in (None, "")]
    small = [r for r in scored if _int(r["family_mzeb_genes"]) < 10]
    carriers = C()
    for r in rows:
        if r.get("b_label") == "cichlid_specific_duplication":
            carriers[f"{_int(r['b_cichlid_with2p'])} of {_int(r['b_cichlid_n'])} cichlid species"] += 1
    return {"label_counts": dict(C(r.get("b_label") for r in rows)),
            "anchor_counts": dict(C(r.get("b_anchor") for r in rows)),
            "label_by_family_size": {lab: dict(C(r.get("b_label") for r in rows
                                                 if lo <= _int(r["family_mzeb_genes"]) <= hi))
                                     for lab, (lo, hi) in SIZE_BINS.items()},
            "label_by_human_relationship": {rel: dict(C(r.get("b_label") for r in rows
                                                        if r["human_relationship"] == rel))
                                            for rel in ("1:1", "1:many", "many:1", "many:many", "paralog", "none")},
            "histogram_share_with_2plus_copies_non_cichlid": _hist([_float(r["b_other_share2p"]) for r in scored]),
            "histogram_share_with_2plus_copies_cichlid": _hist([_float(r["b_cichlid_share2p"]) for r in scored]),
            "patterns_small_families_(teleost,percomorph,ovalentaria,cichlid)": dict(
                C(r["b_pattern"] for r in small).most_common(15)),
            "b2_cichlid_carriers_for_cichlid_specific_duplications": dict(carriers.most_common(15)),
            "expanded_family_genes": sum(1 for r in rows if r.get("expanded_family")),
            "thresholds_note": "kept/lost/many thresholds are in parameters; read them off the histograms"}


def _names_summary(rows):
    C = collections.Counter
    out = {}
    for key in REFS:
        have = [r for r in rows if r.get(f"{key}_closest")]
        out[key] = {
            "closest_name_status": dict(C(r.get(f"{key}_closest_name_status", "") for r in have)),
            "share_low_name_confidence": _frac(sum(1 for r in have if r.get(f"{key}_name_confidence") == "low"),
                                               len(have)),
            "genes_with_unreliable_orthologs": sum(1 for r in rows if _int(r.get(f"{key}_orthologs_unreliable")) > 0)}
    return out


def _symbol_summary(rows):
    """How often the closest human / zebrafish gene has the gene's own NCBI symbol."""
    C = collections.Counter
    out = {"own_symbol_named": sum(1 for r in rows if r.get("own_symbol") and not r["own_symbol"].upper().startswith("LOC")),
           "genes": len(rows)}
    for key in REFS:
        g = f"{key}_symbol_match"
        out[key] = {"grades": dict(C(r.get(g) for r in rows))}
        for rel in ("1:1", "1:many", "paralog"):
            sub = [r for r in rows if r[f"{key}_relationship"] == rel and r.get(g) in ("exact", "stem", "none")]
            out[key][f"checkable_{rel}"] = len(sub)
            out[key][f"exact_rate_{rel}"] = _frac(sum(1 for r in sub if r[g] == "exact"), len(sub))
            out[key][f"exact_or_stem_rate_{rel}"] = _frac(sum(1 for r in sub if matched(r[g])), len(sub))
    return out


def _identity(by, genes, out_dir):
    """Do TreeFam and Compara name the same reference gene (official ID, reliable names only)?"""
    ortho = ("1:1", "1:many", "many:1", "many:many")
    res, diffs = {}, []
    adj = {key: collections.Counter() for key in REFS}
    for key in REFS:
        c = collections.Counter()
        for g in genes:
            a, b = by["TreeFam"][g], by["Compara"][g]
            if a[f"{key}_relationship"] not in ortho or b[f"{key}_relationship"] not in ortho:
                c["not_applicable_no_ortholog_in_one_method"] += 1
                continue
            v = identity_agreement(a.get(f"_{key}_ids", set()), b.get(f"_{key}_ids", set()))
            c[v] += 1
            if v == "different_gene":
                ga, gb = a.get(f"{key}_symbol_match", ""), b.get(f"{key}_symbol_match", "")
                if ga in ("unnamed", "no_reference_symbol", "") and gb in ("unnamed", "no_reference_symbol", ""):
                    adj[key]["cannot_adjudicate_gene_unnamed"] += 1
                elif matched(ga) and matched(gb):
                    adj[key]["both_methods_match_own_symbol"] += 1
                elif matched(ga):
                    adj[key]["TreeFam_matches_own_symbol"] += 1
                elif matched(gb):
                    adj[key]["Compara_matches_own_symbol"] += 1
                else:
                    adj[key]["neither_matches_own_symbol"] += 1
            if v == "different_gene" and len(diffs) < 2000:
                diffs.append({"gene_id": g, "question": key,
                              "own_symbol": a.get("own_symbol", ""),
                              "TreeFam_symbol_match": a.get(f"{key}_symbol_match", ""),
                              "Compara_symbol_match": b.get(f"{key}_symbol_match", ""),
                              "TreeFam_symbols": a.get(f"{key}_orthologs_symbols", ""),
                              "Compara_symbols": b.get(f"{key}_orthologs_symbols", ""),
                              "TreeFam_ids": a.get(f"{key}_orthologs_ids", ""),
                              "Compara_ids": b.get(f"{key}_orthologs_ids", ""),
                              "TreeFam_relationship": a[f"{key}_relationship"],
                              "Compara_relationship": b[f"{key}_relationship"],
                              "TreeFam_family": a["family"], "Compara_family": b["family"]})
        ans = c["same_gene"] + c["different_gene"]
        res[key] = {**dict(c), "same_gene_rate_when_comparable": _frac(c["same_gene"], ans),
                    "adjudication_by_own_ncbi_symbol": dict(adj[key])}
    if diffs:
        with open(out_dir / "pilot_identity_differences.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(diffs[0]), delimiter="\t")
            w.writeheader()
            w.writerows(diffs)
    return res


def _hist(v, bins=10):
    h = [0] * bins
    for x in v:
        h[min(int(x * bins), bins - 1)] += 1
    return h


def _agreement(all_rows, sampled, out_dir):
    by = {m: {r["gene_id"]: r for r in rows} for m, rows in all_rows.items()}
    genes = sorted(set(by["TreeFam"]) & set(by["Compara"]))
    res = {"genes_called_in_both_methods": len(genes), "sampled_genes_in_both": len(sampled & set(genes))}
    review = []
    small = [g for g in genes if _int(by["TreeFam"][g]["family_mzeb_genes"]) < 10
             and _int(by["Compara"][g]["family_mzeb_genes"]) < 10]
    small_set = set(small)
    res["genes_in_small_families_both"] = len(small)
    for key in REFS:
        labels, pairs = collections.Counter(), collections.Counter()
        small_labels = collections.Counter()
        for g in genes:
            vals = {}
            for m in by:
                v = by[m][g][f"{key}_relationship"]
                vals[m] = None if v == "none" else v      # none = no such gene in this family
            lab, _, _ = label_agreement(vals)
            labels[lab] += 1
            if g in small_set:
                small_labels[lab] += 1
            pairs[(vals["TreeFam"], vals["Compara"])] += 1
            if lab == "split":
                review.append({"gene_id": g, "question": f"{key}_relationship",
                               "TreeFam": vals["TreeFam"], "Compara": vals["Compara"],
                               "TreeFam_family": by["TreeFam"][g]["family"],
                               "Compara_family": by["Compara"][g]["family"],
                               "TreeFam_closest": by["TreeFam"][g][f"{key}_closest"],
                               "Compara_closest": by["Compara"][g][f"{key}_closest"],
                               "TreeFam_tied_roots": by["TreeFam"][g]["n_tied_roots"],
                               "Compara_tied_roots": by["Compara"][g]["n_tied_roots"]})
        ans = labels["agree"] + labels["split"]
        sans = small_labels["agree"] + small_labels["split"]
        res[key] = {"labels": dict(labels), "labels_small_families": dict(small_labels),
                    "agree_rate_when_both_answer": _frac(labels["agree"], ans),
                    "agree_rate_when_both_answer_small_families": _frac(small_labels["agree"], sans),
                    "pairs_TreeFam_Compara": [[list(map(str, k)), v] for k, v in pairs.most_common(20)]}
    labs = collections.Counter()
    pairs = collections.Counter()
    for g in genes:
        vals = {}
        for m in by:
            v = by[m][g].get("b_label")
            vals[m] = None if v in NO_B else v
        lab, _, _ = label_agreement(vals)
        labs[lab] += 1
        pairs[(vals["TreeFam"], vals["Compara"])] += 1
    res["goal_B_label"] = {"labels": dict(labs),
                           "pairs_TreeFam_Compara": [[list(map(str, k)), v] for k, v in pairs.most_common(15)]}
    if review:
        with open(out_dir / "pilot_review.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(review[0]), delimiter="\t")
            w.writeheader()
            w.writerows(review)
    res["identity_of_the_reference_gene"] = _identity(by, genes, out_dir)
    res["review_rows"] = len(review)
    res["note"] = ("relationships are compared by category; the reference gene is compared by official ID "
                   "(HGNC / ZFIN) where both methods have reliable names, see identity_of_the_reference_gene")
    return res


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _frac(k, n):
    return round(k / n, 4) if n else None


def _mmm(v):
    v = sorted(v)
    return [v[0], v[len(v) // 2], v[-1]] if v else None


def _deciles(v):
    v = sorted(v)
    return [round(v[int(q / 10 * (len(v) - 1))], 3) for q in range(11)] if v else None
