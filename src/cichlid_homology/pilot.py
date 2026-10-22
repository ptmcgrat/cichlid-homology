"""The pilot: rooting, relationships, dating and the comparison of TreeFam and Compara for one
focal species, on a sample of families. Reads only; writes to a scratch folder."""
from __future__ import annotations

import collections
import csv
import json
import time
import traceback
from dataclasses import asdict
from pathlib import Path

from . import __version__
from .calls import REFS, call_gene
from .consensus import label_agreement
from .dating import Params
from .master import load_species_tree
from .orthology import RTree
from .read import load_family, read_tsv, sample_families
from .rooting import root_tree
from .species import SHADOW_RULE, drop_shadowed, load_leaf_species
from .tree import leaves, parse_newick, prune

ROOT_RULE = "reconcile with the master species tree (fewest duplications); never midpoint"

BASE_COLS = ["gene_id", "method", "family", "tree_mode", "placement_status", "lwr", "assign_flags",
             "fusion_partner", "rooting", "root_cost", "n_tied_roots", "roots_checked",
             "leaves_after_prune", "shadowed_dropped", "parent_children", "species_in_tree"]
REF_COLS = ["relationship", "n_ref", "n_focal", "closest", "closest_is_ortholog",
            "n_ref_in_family", "closest_tie", "stable_across_roots", "paralog_level",
            "dup_level", "dup_class", "dup_resolved", "dup_levels_sides", "dup_outgroup_single"]
TEL_COLS = ["teleost_tel_species", "teleost_tel_one", "teleost_tel_two_plus", "teleost_tel_zero",
            "teleost_tel_frac_one_of_present"]
COLUMNS = BASE_COLS + [f"{k}_{c}" for k in REFS for c in REF_COLS] + TEL_COLS


def process_family(fam, ls, shadowed, st, clades, prm, own_prefix, prefixes, alt_t, meta):
    root = parse_newick(fam.tree)
    sp_of = {x: ls.species(x) for x in fam.leaves}
    kept, dropped = drop_shadowed(sp_of, shadowed, ls.is_project_gene)
    unknown = {x for x, sp in sp_of.items() if sp is None}
    root = prune(root, set(dropped) | unknown)
    if root is None:
        return [], {"error": "empty after pruning"}
    for x in leaves(root):
        x.sp = sp_of[x.name]
    t0 = time.time()
    rooting = root_tree(root, st, fam.rooted, k=prm.k_roots)
    rts = [RTree(r) for r in rooting.roots]
    by_name = [{x.name: x for x in rt.leaves} for rt in rts]
    rows, sens = [], collections.Counter()
    for name in sorted(n for n in by_name[0] if n.startswith(own_prefix + "_")):
        g = by_name[0][name]
        own = g.sp
        row = call_gene(rts[0], g, own, prm, clades, fam.tree_mode, prefixes)
        # stability across tied roots
        for key in REFS:
            rels = {row[f"{key}_relationship"]}
            for rt, bn in zip(rts[1:], by_name[1:]):
                from .orthology import relationship
                rels.add(relationship(rt, bn[name], REFS[key], own, prm.overlap_t)["rel"])
            row[f"{key}_stable_across_roots"] = len(rels) == 1
        # sensitivity to the overlap threshold
        from .orthology import relationship
        for t in alt_t:
            for key, ref in REFS.items():
                if relationship(rts[0], g, ref, own, t)["rel"] != row[f"{key}_relationship"]:
                    sens[(t, key)] += 1
        row.update({"gene_id": name, "family": fam.family_id, "tree_mode": fam.tree_mode,
                    "rooting": rooting.how, "root_cost": rooting.cost,
                    "n_tied_roots": rooting.n_tied, "roots_checked": len(rts),
                    "leaves_after_prune": len(rts[0].leaves), "shadowed_dropped": len(dropped)})
        rows.append(row)
    info = {"leaves": len(rts[0].leaves), "n_tied": rooting.n_tied, "cost": rooting.cost,
            "how": rooting.how, "dropped": len(dropped), "unknown": len(unknown),
            "genes": len(rows), "seconds": round(time.time() - t0, 3), "sens": sens}
    return rows, info


def run_pilot(layout, cat, clades, species, n, prm: Params, alt_t, out_dir: Path) -> dict:
    prefixes = {s.prefix for s in cat.species(release="v1")}
    st = load_species_tree()
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = {"cichlid_homology": __version__, "species": species,
               "parameters": {**asdict(prm), "alt_overlap_t": list(alt_t), "n_families": n,
                              "shadow_rule": SHADOW_RULE, "rooting_rule": ROOT_RULE,
                              "mouse": "out of scope for release 1"}}
    all_rows = {}
    for part in ("TreeFam", "Compara"):
        t0 = time.time()
        res = layout.results("v1", part)
        gp = {r["gene_id"]: r for r in read_tsv(
            Path(res.native().local) / "summary" / "gene_placement.tsv") if r["gene_id"].startswith(species + "_")}
        root, chosen, assign_rows = sample_families(layout, part, species, n)
        assign = {r["gene_id"]: r for r in assign_rows}
        fams = {}
        for f in chosen:
            try:
                fams[f] = load_family(root / f, prefixes, f)
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
                a, p = assign.get(r["gene_id"], {}), gp.get(r["gene_id"], {})
                r.update(method=part, placement_status=p.get("status", ""), lwr=p.get("lwr", ""),
                         assign_flags=a.get("flags", ""), fusion_partner=a.get("fusion_partner", ""))
            rows_part.extend(rows)
            sens.update(info.pop("sens", {}))
            per_family.append({"family": f, **info})
        all_rows[part] = rows_part
        with open(out_dir / f"pilot_genes.{part}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore", restval="")
            w.writeheader()
            w.writerows(rows_part)
        summary[part] = _summarize(rows_part, per_family, errors, sens, prm, shadowed, time.time() - t0)
    summary["agreement"] = _agreement(all_rows, out_dir)
    (out_dir / "pilot_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    return summary


def _summarize(rows, fams, errors, sens, prm, shadowed, seconds):
    C = collections.Counter
    out = {"families_ok": len(fams), "families_errors": len(errors),
           "first_errors": errors[:3], "genes_called": len(rows), "shadowed": sorted(shadowed),
           "seconds": round(seconds, 1),
           "rooting_how": dict(C(f["how"] for f in fams)),
           "root_cost_min_med_max": _mmm([f["cost"] for f in fams if f["cost"] is not None]),
           "n_tied_roots": dict(C(min(f["n_tied"], 5) for f in fams)),
           "family_seconds_max": max((f["seconds"] for f in fams), default=None),
           "unknown_species_leaves": sum(f["unknown"] for f in fams)}
    for key in REFS:
        out[key] = {
            "relationship": dict(C(r[f"{key}_relationship"] for r in rows)),
            "closest_is_ortholog_false_but_closest_present": sum(
                1 for r in rows if r[f"{key}_closest"] and r[f"{key}_closest_is_ortholog"] is False),
            "stable_across_roots": _frac(sum(1 for r in rows if r[f"{key}_stable_across_roots"]), len(rows)),
            "closest_tie": sum(1 for r in rows if r[f"{key}_closest_tie"]),
            "paralog_level": dict(C(r.get(f"{key}_paralog_level") for r in rows if r.get(f"{key}_paralog_level"))),
            "dup_class": dict(C(r.get(f"{key}_dup_class") for r in rows if r.get(f"{key}_dup_class"))),
            "dup_unresolved": sum(1 for r in rows if r.get(f"{key}_dup_resolved") is False),
            "threshold_changes": {str(t): v for (t, k), v in sens.items() if k == key},
        }
    # goal 5: teleost pattern of one-to-one human orthologs
    one = [r["teleost_tel_frac_one_of_present"] for r in rows
           if r["human_relationship"] == "1:1" and r.get("teleost_tel_frac_one_of_present") is not None]
    out["teleost_pattern_human_1to1"] = {
        "n": len(one), "frac_one_deciles": _deciles(one),
        "share_ge_0.7": _frac(sum(1 for x in one if x >= 0.7), len(one)),
        "share_ge_0.8": _frac(sum(1 for x in one if x >= 0.8), len(one)),
        "share_ge_0.9": _frac(sum(1 for x in one if x >= 0.9), len(one))}
    lw = []
    for r in rows:
        try:
            lw.append((float(r["lwr"]) < 0.5, r["human_relationship"]))
        except (TypeError, ValueError):
            pass
    out["low_support_(lwr<0.5)_relationships"] = dict(C(rel for low, rel in lw if low))
    return out


def _agreement(all_rows, out_dir):
    by = {m: {r["gene_id"]: r for r in rows} for m, rows in all_rows.items()}
    genes = sorted(set().union(*[set(d) for d in by.values()]))
    res, review = {}, []
    for key in REFS:
        labels, pairs = collections.Counter(), collections.Counter()
        for g in genes:
            vals = {}
            for m in by:
                r = by[m].get(g)
                v = r[f"{key}_relationship"] if r else None
                vals[m] = None if v == "none" else v          # none = no such gene in this family
            lab, val, _ = label_agreement(vals)
            labels[lab] += 1
            pairs[(vals.get("TreeFam"), vals.get("Compara"))] += 1
            if lab == "split":
                review.append({"gene_id": g, "question": f"{key}_relationship",
                               **{f"{m}": vals[m] for m in vals},
                               **{f"{m}_family": (by[m].get(g) or {}).get("family", "") for m in by}})
        res[key] = {"labels": dict(labels),
                    "top_pairs_TreeFam_Compara": [[list(map(str, k)), v] for k, v in pairs.most_common(15)]}
    if review:
        with open(out_dir / "pilot_review.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(review[0]), delimiter="\t")
            w.writeheader()
            w.writerows(review)
    res["note"] = ("relationship category only; identity of the closest human/zebrafish gene is not "
                   "compared until the identifier mappings exist (symbols only)")
    res["review_rows"] = len(review)
    return res


def _frac(k, n):
    return round(k / n, 4) if n else None


def _mmm(v):
    v = sorted(v)
    return [v[0], v[len(v) // 2], v[-1]] if v else None


def _deciles(v):
    v = sorted(v)
    return [round(v[int(q / 10 * (len(v) - 1))], 3) for q in range(11)] if v else None
