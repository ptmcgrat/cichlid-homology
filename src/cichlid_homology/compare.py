"""Compare two gene tables of the same species (0.20.1): how many calls change after an input changed.

Made for the change of the master species tree (coordinator v0.14.2 resolved its fish polytomies), which moves the parsimony and the gene-tree dating of
the pilot and the cuts of the clade split, so every call built on them can move. Run the pilot, the clade split and the gene table again, then compare the
new gene table with the old one: which genes changed which call, the transitions between categories, and which ohnolog pairs were kept, lost, gained or
re-paired. Genes are matched by `gene_id`; a gene in only one table is counted apart. No call is recomputed here."""
from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

csv.field_size_limit(sys.maxsize)

GROUPS = {
    "closest genes": ["human_gene", "zebrafish_gene"],
    "B1 call": ["b1_category", "b1_confidence", "vote_thresholds", "vote_parsimony", "vote_tree", "tree_call", "later_dup_level"],
    "clade and trend": ["clade_id", "clade_type", "teleost_trend", "rel_human", "rel_zebrafish", "copy_vs_trend", "out_paralog", "human_clean_tier"],
    "ohnolog": ["ohnolog_status", "ohnolog_gene", "ohnolog_confidence", "ev_lineages", "ev_synteny"],
}
HEADLINE = ["b1_category", "tree_call", "ohnolog_status", "ohnolog_gene", "teleost_trend", "rel_human", "rel_zebrafish", "clade_type", "out_paralog", "later_dup_level"]
PAIRISH = ("paired", "paired (resolved by synteny)", "pair retained in few teleosts (synteny-supported)", "possible pair (the other teleosts lost a copy)")


def read(path) -> list:
    with open(Path(path).expanduser(), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _pairs(rows) -> dict:
    out = {}
    for r in rows:
        if r.get("ohnolog_status") in PAIRISH and r.get("ohnolog_gene"):
            out[frozenset((r["gene_id"], r["ohnolog_gene"]))] = r["ohnolog_status"]
    return out


def compare_gene_tables(old, new) -> dict:
    o, n = {r["gene_id"]: r for r in old}, {r["gene_id"]: r for r in new}
    both = sorted(set(o) & set(n))
    res = {"genes_old": len(o), "genes_new": len(n), "genes_compared": len(both), "only_in_old": len(set(o) - set(n)), "only_in_new": len(set(n) - set(o))}
    cols = {}
    changed = collections.defaultdict(dict)
    for g, cl in GROUPS.items():
        for c in cl:
            if c not in (old[0] if old else {}) or c not in (new[0] if new else {}):
                continue
            tr = collections.Counter()
            for gid in both:
                a, b = o[gid].get(c, ""), n[gid].get(c, "")
                if a != b:
                    tr[(a, b)] += 1
                    changed[gid][c] = (a, b)
            cols[c] = {"group": g, "changed": sum(tr.values()), "transitions": [[a, b, k] for (a, b), k in tr.most_common(12)]}
    res["columns"] = cols
    heads = [c for c in HEADLINE if c in cols]
    any_call = {gid for gid, ch in changed.items() if any(c in ch for c in heads)}
    res["headline_columns"] = heads
    res["genes_with_a_changed_call"] = len(any_call)
    res["genes_with_a_changed_call_share"] = round(len(any_call) / len(both), 4) if both else 0
    by_group = {g: len({gid for gid, ch in changed.items() if any(c in ch for c in cl)}) for g, cl in GROUPS.items()}
    res["genes_changed_by_group"] = by_group
    po, pn = _pairs(old), _pairs(new)
    kept, lost, gained = set(po) & set(pn), set(po) - set(pn), set(pn) - set(po)
    res["pairs"] = {"old": len(po), "new": len(pn), "kept": len(kept), "lost": len(lost), "gained": len(gained),
                    "kept_with_a_different_status": sum(1 for k in kept if po[k] != pn[k]),
                    "genes_with_a_different_partner": sum(1 for gid in both if o[gid].get("ohnolog_gene") and n[gid].get("ohnolog_gene") and o[gid]["ohnolog_gene"] != n[gid]["ohnolog_gene"])}
    res["changed_genes"] = {gid: {c: list(v) for c, v in ch.items()} for gid, ch in changed.items()}
    return res


def compare_clades(old, new) -> dict:
    key = lambda r: f"{r['family']}:{r['group']}"
    o, n = {key(r): r for r in old}, {key(r): r for r in new}
    both = set(o) & set(n)
    fo = collections.Counter(r["family"] for r in old)
    fn = collections.Counter(r["family"] for r in new)
    res = {"clades_old": len(o), "clades_new": len(n), "clades_in_both": len(both), "families_with_a_different_number_of_clades": sum(1 for f in set(fo) | set(fn) if fo[f] != fn[f])}
    for c in ("teleost_trend", "n_human", "tetrapod_single_copy_share", "zebrafish_genes", "teleost_median_genes"):
        if old and new and c in old[0] and c in new[0]:
            res[f"changed_{c}"] = sum(1 for k in both if o[k][c] != n[k][c])
    return res


def report_text(res, clades=None) -> str:
    L = ["# How many calls change", "",
         f"{res['genes_compared']} genes compared ({res['only_in_old']} only in the old table, {res['only_in_new']} only in the new).", "",
         f"**{res['genes_with_a_changed_call']} genes ({100 * res['genes_with_a_changed_call_share']:.1f}%) change at least one call** among: {', '.join(res['headline_columns'])}.", "",
         "## Genes changed, by group of columns", "", "| group | genes with a change |", "|---|---|"]
    L += [f"| {g} | {v} |" for g, v in res["genes_changed_by_group"].items()]
    p = res["pairs"]
    L += ["", "## Ohnolog pairs", "",
          f"Old {p['old']}, new {p['new']}: {p['kept']} kept ({p['kept_with_a_different_status']} with a different status), {p['lost']} lost, {p['gained']} gained; "
          f"{p['genes_with_a_different_partner']} genes have a different partner.", "", "## Column by column", ""]
    for c, d in sorted(res["columns"].items(), key=lambda kv: -kv[1]["changed"]):
        L += [f"### {c} ({d['group']}): {d['changed']} genes change", ""]
        L += [f"- {a or '(empty)'} to {b or '(empty)'}: {k}" for a, b, k in d["transitions"][:8]] + [""]
    if clades:
        L += ["## Clade tables", "", json.dumps(clades, indent=1), ""]
    return "\n".join(L)


def write_outputs(out, res, clades=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    ch = res.pop("changed_genes")
    (out / "compare.md").write_text(report_text(res, clades))
    (out / "compare.json").write_text(json.dumps({**res, "clades": clades}, indent=1))
    cols = sorted({c for v in ch.values() for c in v})
    with open(out / "changed_genes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id"] + [f"{c}_{s}" for c in cols for s in ("old", "new")])
        for gid in sorted(ch):
            row = [gid]
            for c in cols:
                a, b = ch[gid].get(c, ("", ""))
                row += [a, b]
            w.writerow(row)
    res["changed_genes"] = ch
