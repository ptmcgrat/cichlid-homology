"""Goals B1 and B2: the copy-number profile of a family across teleost species.

For a cichlid gene g, take the smallest clade that contains g and an orthologous non-teleost anchor
(human if the tree has it, otherwise any other non-teleost jawed vertebrate). That clade is rooted
outside the teleost genome duplication, so it holds both duplicates if both were kept. Count the genes
of each teleost species in it. Species are grouped into the non-cichlid teleosts (levels teleost,
percomorph, ovalentaria, merged because TreeFam has only two to four species per level) and the
cichlids. For each group the share of species with two or more copies decides its state:

    kept     share >= kept_frac   (the pilot's histograms have a valley near 0.5)
    lost     share <  lost_frac
    partial  in between

B1 label (non-cichlid state, cichlid state):
    lost, lost                 single                       a duplicate was lost early
    kept, kept or partial      retained_pair                both copies kept across teleosts
    lost, kept or partial      cichlid_specific_duplication a later duplicate (B2 counts the carriers)
    kept, lost                 cichlid_lost_one_copy
    partial, any               partial_retention            the middle ground (read `b_pattern`)
    3+ copies widely           multi_copy                   an expansion, not a pair
B2 (counts only; OrthoFinder is the authority for where carriers meet):
    b_cichlid_n, b_cichlid_with2p, b_cichlid_extra_species.
`b_pattern` has one character per level (teleost, percomorph, ovalentaria, cichlid):
    + extra copy kept   . absent   ~ in between   - no species.
"""
from __future__ import annotations

from .clade_counts import smallest_clade
from .dating import Params, level_of_species
from .orthology import TELEOST_LEVELS, RTree
from .tree import Node

LEVELS = ("teleost", "percomorph", "ovalentaria", "cichlid")
NON_CICHLID = ("teleost", "percomorph", "ovalentaria")


def anchor_clade(rt: RTree, g: Node, anchor_sets, t: float):
    """First (nearest) ancestor of g whose other children hold an *orthologous* anchor leaf.
    anchor_sets: ordered [(name, set of species)]. Returns (name, ancestor) or (None, None)."""
    walk = rt.walk(g, t)
    for name, species in anchor_sets:
        if not species:
            continue
        for (a, _, c, spec) in walk:
            if spec and any(x.sp in species for x in rt.leaves_in(c)):
                return name, a
    return None, None


def copy_counts(rt: RTree, clade: Node, level_of) -> dict:
    """{level: {species: copies}} for the teleost species in the clade."""
    out = {lv: {} for lv in LEVELS}
    for x in rt.leaves_in(clade):
        lv = level_of(x.sp)
        if lv in out:
            out[lv][x.sp] = out[lv].get(x.sp, 0) + 1
    return out


def state(share, prm: Params) -> str:
    if share == "" or share is None:
        return "none"
    if share >= prm.kept_frac:
        return "kept"
    if share < prm.lost_frac:
        return "lost"
    return "partial"


def _group(prefix: str, d: dict) -> dict:
    n = len(d)
    two = sum(1 for v in d.values() if v >= 2)
    three = sum(1 for v in d.values() if v >= 3)
    return {f"b_{prefix}_n": n, f"b_{prefix}_with2p": two, f"b_{prefix}_with3p": three,
            f"b_{prefix}_share2p": round(two / n, 4) if n else "",
            f"b_{prefix}_share3p": round(three / n, 4) if n else ""}


def summarize(counts: dict, prm: Params) -> dict:
    row = {}
    for lv in LEVELS:
        c = counts[lv]
        row[f"b_{lv}_n"] = len(c)
        row[f"b_{lv}_1"] = sum(1 for v in c.values() if v == 1)
        row[f"b_{lv}_2"] = sum(1 for v in c.values() if v == 2)
        row[f"b_{lv}_3p"] = sum(1 for v in c.values() if v >= 3)
    other = {}
    for lv in NON_CICHLID:
        other.update(counts[lv])
    row.update(_group("other", other))
    row.update(_group("cichlid", counts["cichlid"]))
    code = {"kept": "+", "lost": ".", "partial": "~", "none": "-"}
    row["b_pattern"] = "".join(
        code[state(_share(counts[lv]), prm)] for lv in LEVELS)
    row["b_cichlid_extra_species"] = ";".join(sorted(s for s, v in counts["cichlid"].items() if v >= 2))
    row["b_other_state"] = state(row["b_other_share2p"], prm)
    row["b_cichlid_state"] = state(row["b_cichlid_share2p"], prm)
    return row


def _share(d: dict):
    return (sum(1 for v in d.values() if v >= 2) / len(d)) if d else ""


def label(row: dict, prm: Params) -> str:
    if row["b_other_n"] < prm.min_other or row["b_cichlid_n"] < 1:
        return "insufficient_species"
    if row["b_other_share3p"] != "" and row["b_other_share3p"] >= prm.many_frac:
        return "multi_copy"
    so, sc = row["b_other_state"], row["b_cichlid_state"]
    if so == "partial":
        return "partial_retention"
    if so == "lost":
        return "single" if sc == "lost" else "cichlid_specific_duplication"
    return "cichlid_lost_one_copy" if sc == "lost" else "retained_pair"


def clade_profile(rt: RTree, g: Node, clades, t: float, prm: Params | None = None) -> dict:
    prm = prm or Params()
    lv = lambda s: level_of_species(clades, s)
    human = {"homo_sapiens"} & set(rt.species)
    other = {s for s in rt.species if s != "homo_sapiens" and lv(s) == "gnathostome"}
    name, a = anchor_clade(rt, g, [("human", human), ("other_gnathostome", other)], t)
    row = {"b_anchor": name or "none"}
    if a is None:
        row["b_label"] = "no_orthologous_anchor"
        return row
    counts = copy_counts(rt, a, lv)
    row.update(summarize(counts, prm))
    row["b_clade_leaves"] = a.hi - a.lo
    row["b_teleost_species_in_tree"] = sum(1 for s in rt.species if lv(s) in TELEOST_LEVELS)
    row["b_label"] = label(row, prm)
    return row


def clade_profile_v2(rt: RTree, g: Node, clades, prm: Params | None, human_clade: Node | None) -> dict:
    """Goal B on the clade of the best human gene (the A1 clade), with no orthology test. If the tree has
    no human gene, the smallest clade with any other non-teleost jawed vertebrate is used instead."""
    prm = prm or Params()
    lv = lambda s: level_of_species(clades, s)
    if human_clade is not None:
        name, a = "human", human_clade
    else:
        other = {s for s in rt.species if lv(s) == "gnathostome"}
        a = smallest_clade(rt, g, other)
        name = "other_gnathostome" if a is not None else None
    row = {"b_anchor": name or "none"}
    if a is None:
        row["b_label"] = "no_anchor"
        return row
    counts = copy_counts(rt, a, lv)
    row.update(summarize(counts, prm))
    row["b_clade_leaves"] = a.hi - a.lo
    row["b_teleost_species_in_tree"] = sum(1 for s in rt.species if lv(s) in TELEOST_LEVELS)
    row["b_label"] = label(row, prm)
    return row
