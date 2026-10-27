"""Summaries of the 0.7.0 counts: overall numbers per method, agreement between TreeFam and Compara,
the best reference gene, the NCBI-symbol check, and the teleost distribution (goal B) by category."""
from __future__ import annotations

import collections
import csv
from pathlib import Path

from .clade_counts import CATS
from .consensus import label_agreement
from .symbols import matched

REFS = ("human", "zebrafish")
A, B = "TreeFam", "Compara"
SIZE_BINS = {"1": (1, 1), "2-3": (2, 3), "4-9": (4, 9), "10+": (10, 10 ** 9)}
NO_B = (None, "", "no_anchor", "no_orthologous_anchor", "insufficient_species", "expanded_family")


def _i(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _pct(k, n):
    return round(100 * k / n, 1) if n else None


def _median(v):
    v = sorted(v)
    return v[len(v) // 2] if v else None


def _hist(v, bins=10):
    h = [0] * bins
    for x in v:
        if x is not None:
            h[min(int(x * bins), bins - 1)] += 1
    return h


def overall(rows, key):
    n = len(rows)
    c = collections.Counter(r[f"a_{key}_cat"] for r in rows)
    have = [r for r in rows if r[f"a_{key}_cat"] != "None"]
    return {
        "genes": n,
        "categories": {cat: {"n": c.get(cat, 0), "pct": _pct(c.get(cat, 0), n)} for cat in CATS},
        "by_family_size": {lab: dict(collections.Counter(
            r[f"a_{key}_cat"] for r in rows if lo <= _i(r["family_mzeb_genes"]) <= hi))
            for lab, (lo, hi) in SIZE_BINS.items()},
        "genes_with_a_reference_gene": len(have),
        "best_gene_tied_on_path_length": sum(1 for r in have if r[f"a_{key}_tie"]),
        "best_gene_name_confidence_low": sum(1 for r in have if r.get(f"a_{key}_name_confidence") == "low"),
        "median_clade_leaves": {cat: _median([_i(r[f"a_{key}_clade_leaves"]) for r in have
                                              if r[f"a_{key}_cat"] == cat]) for cat in CATS[:-1]},
    }


def _side(r, key):
    n_ref, n_foc = _i(r[f"a_{key}_ref_n"]), _i(r[f"a_{key}_focal_n"])
    return ("1" if n_ref == 1 else "2+"), ("1" if n_foc == 1 else "2" if n_foc == 2 else "3+")


def agreement(by, genes, key):
    labels, pairs = collections.Counter(), collections.Counter()
    both = same = ref_same = foc_same = one_none = both_none = 0
    for g in genes:
        ra, rb = by[A][g], by[B][g]
        a, b = ra[f"a_{key}_cat"], rb[f"a_{key}_cat"]
        lab, _, _ = label_agreement({A: None if a == "None" else a, B: None if b == "None" else b})
        labels[lab] += 1
        pairs[(a, b)] += 1
        if a == "None" and b == "None":
            both_none += 1
        elif a == "None" or b == "None":
            one_none += 1
        else:
            both += 1
            same += a == b
            sa, sb = _side(ra, key), _side(rb, key)
            ref_same += sa[0] == sb[0]
            foc_same += sa[1] == sb[1]
    return {"genes_in_both_methods": len(genes), "labels": dict(labels),
            "both_gave_a_category": both, "same_category": same, "same_category_pct": _pct(same, both),
            "reference_side_same_(1_vs_2plus)": _pct(ref_same, both),
            "focal_side_same_(1_vs_2_vs_3plus)": _pct(foc_same, both),
            "only_one_method_has_a_reference_gene": one_none, "neither_has_one": both_none,
            "category_pairs_TreeFam_Compara": [[list(k), v] for k, v in pairs.most_common(25)]}


def _adjudicate(ra, rb, key):
    ga, gb = ra.get(f"a_{key}_symbol_match", ""), rb.get(f"a_{key}_symbol_match", "")
    skip = ("unnamed", "no_reference_symbol", "")
    if ga in skip and gb in skip:
        return "cannot_adjudicate_gene_unnamed"
    if matched(ga) and matched(gb):
        return "both_match_own_symbol"
    if matched(ga):
        return "TreeFam_matches_own_symbol"
    if matched(gb):
        return "Compara_matches_own_symbol"
    return "neither_matches_own_symbol"


def best_gene_agreement(by, genes, key, out_dir: Path | None = None):
    states, adj, diffs = collections.Counter(), collections.Counter(), []
    for g in genes:
        ra, rb = by[A][g], by[B][g]
        if ra[f"a_{key}_cat"] == "None" or rb[f"a_{key}_cat"] == "None":
            continue
        ia = ra.get(f"a_{key}_best_id", "") if ra.get(f"a_{key}_best_name_reliable") else ""
        ib = rb.get(f"a_{key}_best_id", "") if rb.get(f"a_{key}_best_name_reliable") else ""
        if not ia or not ib:
            state = "not_comparable_no_reliable_id"
        elif ia == ib:
            state = "same_best_gene"
        elif ia in rb.get(f"_a_{key}_ids", set()) or ib in ra.get(f"_a_{key}_ids", set()):
            state = "different_best_gene_but_inside_the_other_clade"
        else:
            state = "different_gene"
        states[state] += 1
        if state in ("different_best_gene_but_inside_the_other_clade", "different_gene"):
            v = _adjudicate(ra, rb, key)
            adj[(state, v)] += 1
            if len(diffs) < 3000:
                diffs.append({"gene_id": g, "question": key, "state": state, "own_symbol": ra.get("own_symbol", ""),
                              "TreeFam_best": ra.get(f"a_{key}_best_symbol", ""), "Compara_best": rb.get(f"a_{key}_best_symbol", ""),
                              "TreeFam_cat": ra[f"a_{key}_cat"], "Compara_cat": rb[f"a_{key}_cat"],
                              "TreeFam_symbol_match": ra.get(f"a_{key}_symbol_match", ""),
                              "Compara_symbol_match": rb.get(f"a_{key}_symbol_match", ""),
                              "TreeFam_family": ra["family"], "Compara_family": rb["family"]})
    comparable = sum(v for k, v in states.items() if not k.startswith("not_comparable"))
    out = {**dict(states), "same_best_gene_pct_of_comparable": _pct(states["same_best_gene"], comparable),
           "adjudication_by_own_ncbi_symbol": {f"{s} | {v}": n for (s, v), n in sorted(adj.items())}}
    if diffs and out_dir is not None:
        with open(out_dir / f"pilot_best_gene_differences.{key}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(diffs[0]), delimiter="\t")
            w.writeheader()
            w.writerows(diffs)
    return out


def symbol_summary(rows):
    out = {"genes": len(rows),
           "genes_with_a_real_own_symbol": sum(1 for r in rows if r.get("own_symbol")
                                                and not r["own_symbol"].upper().startswith("LOC"))}
    for key in REFS:
        g = f"a_{key}_symbol_match"
        per = {}
        for cat in CATS[:-1]:
            sub = [r for r in rows if r[f"a_{key}_cat"] == cat and r.get(g) in ("exact", "stem", "none")]
            per[cat] = {"checkable": len(sub), "exact_pct": _pct(sum(1 for r in sub if r[g] == "exact"), len(sub)),
                        "exact_or_stem_pct": _pct(sum(1 for r in sub if matched(r[g])), len(sub))}
        out[key] = per
    return out


def goal_b(rows):
    C = collections.Counter
    out = {}
    for cat in CATS:
        sub = [r for r in rows if r["a_human_cat"] == cat]
        out[f"human_{cat}"] = {"genes": len(sub), "labels": dict(C(r.get("b_label") for r in sub)),
                               "anchor": dict(C(r.get("b_anchor") for r in sub))}
    for cat in ("OneToOne", "OneToTwo"):
        sc = [r for r in rows if r["a_human_cat"] == cat
              and r.get("b_label_raw", r.get("b_label")) not in NO_B and _f(r.get("b_other_share2p")) is not None]
        carriers = C(f"{_i(r['b_cichlid_with2p'])} of {_i(r['b_cichlid_n'])} cichlid species" for r in sc
                     if r.get("b_label") == "cichlid_specific_duplication")
        out[f"distribution_human_{cat}"] = {
            "genes_scored": len(sc),
            "histogram_share_of_non_cichlid_teleost_species_with_2plus_copies": _hist([_f(r["b_other_share2p"]) for r in sc]),
            "histogram_share_of_cichlid_species_with_2plus_copies": _hist([_f(r["b_cichlid_share2p"]) for r in sc]),
            "patterns_(teleost,percomorph,ovalentaria,cichlid)": dict(C(r["b_pattern"] for r in sc).most_common(12)),
            "cichlid_carriers_when_cichlid_specific": dict(carriers.most_common(12))}
    return out


def goal_b_agreement(by, genes):
    labs, pairs = collections.Counter(), collections.Counter()
    focus = ("OneToOne", "OneToTwo")
    n = 0
    for g in genes:
        ra, rb = by[A][g], by[B][g]
        if ra["a_human_cat"] not in focus or rb["a_human_cat"] not in focus:
            continue
        n += 1
        va = None if ra.get("b_label") in NO_B else ra["b_label"]
        vb = None if rb.get("b_label") in NO_B else rb["b_label"]
        lab, _, _ = label_agreement({A: va, B: vb})
        labs[lab] += 1
        pairs[(va, vb)] += 1
    both = labs["agree"] + labs["split"]
    return {"genes_OneToOne_or_OneToTwo_in_both_methods": n, "labels": dict(labs),
            "agree_pct_when_both_labelled": _pct(labs["agree"], both),
            "pairs_TreeFam_Compara": [[list(map(str, k)), v] for k, v in pairs.most_common(15)]}


def write_splits(by, genes, out_dir: Path):
    rows = []
    for key in REFS:
        for g in genes:
            ra, rb = by[A][g], by[B][g]
            a, b = ra[f"a_{key}_cat"], rb[f"a_{key}_cat"]
            if a == "None" or b == "None" or a == b:
                continue
            rows.append({"gene_id": g, "question": key, "own_symbol": ra.get("own_symbol", ""),
                         "TreeFam_cat": a, "Compara_cat": b,
                         "TreeFam_ref_n": ra[f"a_{key}_ref_n"], "TreeFam_focal_n": ra[f"a_{key}_focal_n"],
                         "Compara_ref_n": rb[f"a_{key}_ref_n"], "Compara_focal_n": rb[f"a_{key}_focal_n"],
                         "TreeFam_clade_leaves": ra[f"a_{key}_clade_leaves"], "Compara_clade_leaves": rb[f"a_{key}_clade_leaves"],
                         "TreeFam_best": ra.get(f"a_{key}_best_symbol", ""), "Compara_best": rb.get(f"a_{key}_best_symbol", ""),
                         "TreeFam_symbol_match": ra.get(f"a_{key}_symbol_match", ""),
                         "Compara_symbol_match": rb.get(f"a_{key}_symbol_match", ""),
                         "TreeFam_family": ra["family"], "Compara_family": rb["family"]})
    if rows:
        with open(out_dir / "pilot_category_splits.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    return len(rows)


def build(all_rows, out_dir: Path) -> dict:
    by = {m: {r["gene_id"]: r for r in rows} for m, rows in all_rows.items()}
    genes = sorted(set(by[A]) & set(by[B]))
    out = {"genes_called": {m: len(rows) for m, rows in all_rows.items()}, "genes_in_both_methods": len(genes),
           "overall": {m: {k: overall(rows, k) for k in REFS} for m, rows in all_rows.items()},
           "agreement": {k: agreement(by, genes, k) for k in REFS},
           "best_gene": {k: best_gene_agreement(by, genes, k, out_dir) for k in REFS},
           "symbol_check": {m: symbol_summary(rows) for m, rows in all_rows.items()},
           "goal_B": {m: goal_b(rows) for m, rows in all_rows.items()},
           "goal_B_agreement": goal_b_agreement(by, genes)}
    out["category_splits_written"] = write_splits(by, genes, out_dir)
    return out
