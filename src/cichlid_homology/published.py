"""Our pairs against the published teleost-duplication pairs (OHNOLOGS v2), zebrafish (0.23.0). Read-only; symbols are matched, as `ohnologs.py` does.

Recall: of the published pairs whose two zebrafish genes each have a gene of ours as closest match (the represented ones), how many are our pair. The rest are counted by where our two genes are: in one clade but not
partners, in one family with different clades, in different families, or labelled older (gar or reedfish carry both copies). A pair of ours whose two genes are both closest to one of the two zebrafish genes
is recovered too, and is counted apart (`recovered, both genes closest to one zebrafish gene`). Precision: of our pairs whose two genes match two different reliable zebrafish names, how many are published,
and of the rest how many have the a and b names of one stem (plausible: the published lists are incomplete) and how many have unrelated names."""
from __future__ import annotations

import collections
import csv
import re
import sys
from pathlib import Path

from .ohnologs import find_files

csv.field_size_limit(sys.maxsize)
PAIRS = ("paired", "paired (resolved by synteny)", "pair retained in few teleosts (synteny-supported)", "possible pair (the other teleosts lost a copy)")
OLDER = "older duplicate pair (gar or reedfish carry both copies)"
LEVELS = ("Strict", "Intermediate", "Relaxed")


def load_zebrafish_pairs(folder) -> dict:
    files = find_files(folder)
    out = {}
    for lvl in LEVELS:
        f = files.get(("drerio", lvl))
        if not f:
            continue
        s = set()
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("Gene type") != "protein_coding" or r.get("Duplication time") != "FishWGD":
                    continue
                a, b = (r.get("Symbol1") or "").strip().lower(), (r.get("Symbol2") or "").strip().lower()
                if a and b and a != b and "(" not in a + b and not a.startswith("ens") and not b.startswith("ens"):
                    s.add(frozenset((a, b)))
        out[lvl] = s
    return out


def reliable(s) -> str:
    s = (s or "").strip().lower()
    return s if s and not s.startswith("loc") and ":" not in s and "(" not in s and not s.startswith("ens") else ""


def recall(rows: list, pairs: set) -> dict:
    by = collections.defaultdict(list)
    t = {r["gene_id"]: r for r in rows}
    for r in rows:
        z = reliable(r.get("zebrafish_symbol"))
        if z:
            by[z].append(r["gene_id"])
    c = collections.Counter()
    for pr in pairs:
        a, b = sorted(pr)
        ga, gb = by.get(a, []), by.get(b, [])
        if not ga or not gb:
            side = ga or gb
            if len(side) >= 2 and any(t[x]["ohnolog_gene"] == y and t[x]["ohnolog_status"] in PAIRS + (OLDER,) for x in side for y in side if x != y):
                c["recovered, both genes closest to one zebrafish gene"] += 1
            else:
                c["not testable (a zebrafish gene has no gene of ours)"] += 1
            continue
        best = None
        rank = {"recovered": 0, "recovered, labelled older duplicate": 0, "in one clade, not partners": 1, "in one family, different clades": 2, "in different families": 3}
        for x in ga:
            for y in gb:
                rx, ry = t[x], t[y]
                if rx["ohnolog_gene"] == y and rx["ohnolog_status"] in PAIRS:
                    k = "recovered"
                elif rx["ohnolog_gene"] == y and rx["ohnolog_status"] == OLDER:
                    k = "recovered, labelled older duplicate"
                elif rx["clade_id"] and rx["clade_id"] == ry["clade_id"]:
                    k = "in one clade, not partners"
                elif rx["family"] and rx["family"] == ry["family"] and rx["family"] != "NoHit":
                    k = "in one family, different clades"
                else:
                    k = "in different families"
                if best is None or rank[k] < rank[best]:
                    best = k
        c[best] += 1
    n = len(pairs)
    hidden = c["recovered, both genes closest to one zebrafish gene"]
    represented = n - c["not testable (a zebrafish gene has no gene of ours)"] - hidden
    return {"published": n, "represented": represented, "recovered": c["recovered"], "recovered_share_of_represented": round(c["recovered"] / represented, 4) if represented else None,
            "hidden_recoveries": hidden, "recovered_incl_hidden_share": round((c["recovered"] + hidden) / (represented + hidden), 4) if represented + hidden else None, "kinds": dict(c)}


def precision(rows: list, pubs: dict) -> dict:
    t = {r["gene_id"]: r for r in rows}
    seen, c = set(), collections.Counter()
    for r in rows:
        if r["ohnolog_status"] not in PAIRS or r["ohnolog_gene"] not in t:
            continue
        key = frozenset((r["gene_id"], r["ohnolog_gene"]))
        if key in seen:
            continue
        seen.add(key)
        z1, z2 = reliable(r["zebrafish_symbol"]), reliable(t[r["ohnolog_gene"]]["zebrafish_symbol"])
        if not z1 or not z2:
            c["no zebrafish name for one gene"] += 1
        elif z1 == z2:
            c["both closest to the same zebrafish gene"] += 1
        else:
            pk = frozenset((z1, z2))
            lvl = next((l for l in LEVELS if pk in pubs.get(l, ())), None)
            if lvl:
                c[f"published, {lvl.lower()}"] += 1
            elif z1[-1] in "ab" and z2[-1] in "ab" and z1[:-1] == z2[:-1]:
                c["not published; a/b names of one stem"] += 1
            else:
                c["not published; unrelated names"] += 1
    ev = sum(v for k, v in c.items() if k.startswith(("published", "not published")))
    pub = sum(v for k, v in c.items() if k.startswith("published"))
    plausible = pub + c["not published; a/b names of one stem"]
    return {"pairs": len(seen), "evaluable": ev, "published": pub, "published_or_a_b_stem": plausible, "share_published_or_a_b_stem": round(plausible / ev, 4) if ev else None, "kinds": dict(c)}


def read(path) -> list:
    with open(Path(path).expanduser(), newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def report_text(res: dict) -> str:
    L = ["# Our pairs against the published zebrafish pairs (OHNOLOGS v2)", ""]
    for name, d in res.items():
        L += [f"## {name}", "", "| level | published | represented | recovered | share of represented | hidden recoveries | share incl. hidden | where the others are |", "|---|---|---|---|---|---|---|---|"]
        for lvl, r in d["recall"].items():
            others = {k: v for k, v in r["kinds"].items() if k not in ("recovered",) and not k.startswith("not testable")}
            L.append(f"| {lvl} | {r['published']} | {r['represented']} | {r['recovered']} | {100 * (r['recovered_share_of_represented'] or 0):.1f}% | {r['hidden_recoveries']} | {100 * (r['recovered_incl_hidden_share'] or 0):.1f}% | "
                     f"{ {k: v for k, v in others.items() if not k.startswith('recovered, both')} } |")
        p = d["precision"]
        L += ["", f"Precision: {p['pairs']} pairs of ours; {p['evaluable']} have two different reliable zebrafish names; {p['published']} are published; "
              f"{p['published_or_a_b_stem']} ({100 * (p['share_published_or_a_b_stem'] or 0):.0f}%) are published or have the a and b names of one stem.", "", f"{p['kinds']}", ""]
    return "\n".join(L)
