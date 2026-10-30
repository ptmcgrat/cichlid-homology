"""Where the non-teleost genes fall relative to the two teleost lineages of our ohnolog pairs (0.22.1).

For every clade with a teleost split (`lineages.py`) that holds one of our pairs, each outgroup species (gar, reedfish, shark, coelacanth, and the tetrapods together) is put in one of five
cases: absent from the clade; outside both lineages (the expected place: the gene is the ancestral copy, sister to the split); inside one lineage; in both lineages (a gene in each: the species
carries both copies, which is the evidence for a duplication older than that species' split); or mixed (some of its genes outside, some inside). A gene inside one lineage is also placed as the
slower lineage (the shorter stem) or the faster one, because the tree pulls an outgroup toward the slower copy. The unit is the clade, not the gene: the genes of a clade are not independent, and
a clade that holds several pairs is counted once (the number of pairs is given beside it)."""
from __future__ import annotations

import collections
import csv
import json
from pathlib import Path

OUTGROUPS = (("gar", "lepisosteus_oculatus"), ("reedfish", "erpetoichthys_calabaricus"), ("shark", "callorhinchus_milii"), ("coelacanth", "latimeria_chalumnae"))
PAIR_STATUSES = ("paired", "paired (resolved by synteny)", "pair retained in few teleosts (synteny-supported)", "possible pair (the other teleosts lost a copy)")
CASES = ("outside both lineages", "inside one lineage", "in both lineages", "mixed (some outside, some inside)")


def _case(inside: list, outside: list):
    if not inside and not outside:
        return "absent", None
    if not inside:
        return "outside both lineages", None
    lins = sorted({o["lineage"] for o in inside})
    if outside:
        return "mixed (some outside, some inside)", lins
    if len(lins) > 1:
        return "in both lineages", lins
    return "inside one lineage", lins


def species_case(rec: dict, species: str):
    return _case([o for o in rec["outgroup"]["inside"] if o["species"] == species], [o for o in rec["outgroup"]["outside"] if o["species"] == species])


def tetrapod_case(rec: dict):
    tet = {s for s, g in zip(rec["strip"]["species"], rec["strip"]["groups"]) if g == "tetrapod"}
    return _case([o for o in rec["outgroup"]["inside"] if o["species"] in tet], [o for o in rec["outgroup"]["outside"] if o["species"] in tet])


def eligible(rec: dict, pairs_only: bool = True) -> bool:
    if rec.get("mode") != "split" or (rec.get("split") or {}).get("kind") != "teleost":
        return False
    return (not pairs_only) or any((g.get("calls") or {}).get("ohnolog_status") in PAIR_STATUSES for g in rec["genes"])


def _pair_genes(rec: dict) -> list:
    return [g for g in rec["genes"] if (g.get("calls") or {}).get("ohnolog_status") in PAIR_STATUSES]


def _slower(rec: dict, lineage: str):
    stems = {l["id"]: l["stem"] for l in rec["lineages"] if l.get("stem")}
    if lineage not in stems or len(stems) < 2 or len(set(stems.values())) == 1:
        return None
    return stems[lineage] == min(stems.values())


def clade_row(rec: dict) -> dict:
    pg = _pair_genes(rec)
    calls = [g.get("calls") or {} for g in pg]
    row = {"clade": rec["clade"], "family": rec["family"], "pair_genes": len(pg), "level": (rec.get("split") or {}).get("level", ""), "overlap": (rec.get("split") or {}).get("overlap", ""),
           "stem_ratio": rec["triage"].get("stem_ratio"), "deep_order_score": rec["triage"].get("deep_order_score"),
           "published": any((c.get("published_pair") or "").startswith("the two") for c in calls),
           "confidence": ("high" if any(c.get("ohnolog_confidence") == "high" for c in calls) else ("medium" if any(c.get("ohnolog_confidence") == "medium" for c in calls) else "low"))}
    inside_n = 0
    for name, sp in OUTGROUPS:
        case, lins = species_case(rec, sp)
        row[name] = case + (f" ({','.join(lins)})" if lins else "")
        inside_n += case != "absent" and case != "outside both lineages"
    case, lins = tetrapod_case(rec)
    row["tetrapods"] = case + (f" ({','.join(lins)})" if lins else "")
    row["outgroup_species_inside"] = inside_n
    return row


def summarize(records: list, pairs_only: bool = True) -> dict:
    recs = [r for r in records if eligible(r, pairs_only)]
    out = {"clades": len(recs), "pair_genes": sum(len(_pair_genes(r)) for r in recs), "pairs_only": pairs_only, "species": {}}
    for name, sp in OUTGROUPS + (("tetrapods", None),):
        c = collections.Counter()
        slower = collections.Counter()
        byletter = collections.Counter()
        for r in recs:
            case, lins = tetrapod_case(r) if sp is None else species_case(r, sp)
            c[case] += 1
            if case == "inside one lineage":
                s = _slower(r, lins[0])
                slower["slower lineage" if s else ("faster lineage" if s is False else "stems equal or missing")] += 1
                byletter[lins[0]] += 1
        present = len(recs) - c["absent"]
        out["species"][name] = {"clades_with_the_species": present, "absent": c["absent"], **{k: c[k] for k in CASES}, "inside_any_lineage": c["inside one lineage"] + c["in both lineages"] + c["mixed (some outside, some inside)"],
                               "inside_one_lineage_which": dict(slower), "inside_one_lineage_by_letter": dict(byletter)}
    n_inside = collections.Counter(clade_row(r)["outgroup_species_inside"] for r in recs)
    out["clades_by_number_of_outgroup_species_inside"] = dict(sorted(n_inside.items()))
    same = diff = 0
    for r in recs:
        lins = []
        for _n, sp in OUTGROUPS:
            case, l = species_case(r, sp)
            if case == "inside one lineage":
                lins.append(l[0])
        if len(lins) >= 2:
            same += len(set(lins)) == 1
            diff += len(set(lins)) > 1
    out["clades_with_two_or_more_species_inside_one_lineage"] = {"all_in_the_same_lineage": same, "in_different_lineages": diff}
    # does it go with the quality of the pair?
    q = collections.defaultdict(collections.Counter)
    for r in recs:
        row = clade_row(r)
        for key in (("published" if row["published"] else "not published"), f"confidence {row['confidence']}"):
            q[key]["clades"] += 1
            q[key]["any outgroup species inside a lineage"] += row["outgroup_species_inside"] > 0                      # gar, reedfish, shark or coelacanth
            q[key]["tetrapods inside a lineage"] += row["tetrapods"].startswith(("inside", "in both", "mixed"))
    out["by_pair_quality"] = {k: dict(v) for k, v in q.items()}
    stems = collections.defaultdict(list)
    for r in recs:
        sr = r["triage"].get("stem_ratio")
        if sr:
            stems["none inside" if clade_row(r)["outgroup_species_inside"] == 0 else "some inside"].append(sr)
    out["median_stem_ratio"] = {k: sorted(v)[len(v) // 2] for k, v in stems.items()}
    return out


def report_text(s: dict) -> str:
    L = ["# Where the outgroups fall relative to the two teleost lineages", "",
         f"{s['clades']} clades with a teleost split{' that hold one of our pairs' if s['pairs_only'] else ''} ({s['pair_genes']} genes in pairs). The unit is the clade.", "",
         "| outgroup | clades that have it | outside both lineages | inside one lineage | in both lineages | mixed | inside any lineage |", "|---|---|---|---|---|---|---|"]
    for name, d in s["species"].items():
        n = d["clades_with_the_species"] or 1
        pc = lambda k: f"{d[k]} ({100 * d[k] / n:.0f}%)"
        L.append(f"| {name} | {d['clades_with_the_species']} | {pc('outside both lineages')} | {pc('inside one lineage')} | {pc('in both lineages')} | {pc('mixed (some outside, some inside)')} | {pc('inside_any_lineage')} |")
    L += ["", "Percentages are of the clades that have the species. \"In both lineages\" is a gene in each: the species carries both copies.", "", "## Inside one lineage: the slower or the faster one", ""]
    for name, d in s["species"].items():
        if d["inside one lineage"]:
            L.append(f"- {name}: {d['inside_one_lineage_which']}; by letter {d['inside_one_lineage_by_letter']}")
    L += ["", "## How many of the four outgroup species are inside a lineage in one clade", "", json.dumps(s["clades_by_number_of_outgroup_species_inside"]),
          "", f"With two or more inside one lineage each: {s['clades_with_two_or_more_species_inside_one_lineage']}", "", "## Does it go with the quality of the pair?", ""]
    for k, v in s["by_pair_quality"].items():
        f, t = v.get("any outgroup species inside a lineage", 0), v.get("tetrapods inside a lineage", 0)
        L.append(f"- {k}: {v['clades']} clades; gar, reedfish, shark or coelacanth inside a lineage in {f} ({100 * f / v['clades']:.0f}%); tetrapods inside a lineage in {t} ({100 * t / v['clades']:.0f}%)")
    L += ["", f"Median stem ratio of the two lineages: {s['median_stem_ratio']}", ""]
    return "\n".join(L)


def write_outputs(out, records: list, pairs_only: bool = True) -> dict:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    s = summarize(records, pairs_only)
    (out / "outgroup_summary.md").write_text(report_text(s))
    (out / "outgroup_summary.json").write_text(json.dumps(s, indent=1))
    rows = [clade_row(r) for r in records if eligible(r, pairs_only)]
    if rows:
        with open(out / "outgroup_by_clade.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    return s
