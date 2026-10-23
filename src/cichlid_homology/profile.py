"""Goal B: the copy-number profile of a family across teleost species.

For a cichlid gene g, take the smallest clade that contains g and an orthologous non-teleost anchor
(human if the tree has it, otherwise any other non-teleost jawed vertebrate). That clade is rooted
outside the teleost genome duplication, so it holds both duplicates if both were kept. Then count the
genes of each teleost species in it. A species with one copy: the duplicate was lost (or never
existed); two copies: both were kept; the species are grouped by ladder level so a duplication that
only the cichlids carry shows up as cichlids with 2 and the other teleosts with 1.

No node dating and no weight on any one outgroup branch. The labels use provisional thresholds
(`frac`), and the pilot prints the distributions so they can be set from the data.
"""
from __future__ import annotations

from .dating import level_of_species
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


def summarize(counts: dict) -> dict:
    row = {}
    for lv in LEVELS:
        c = counts[lv]
        row[f"b_{lv}_n"] = len(c)
        row[f"b_{lv}_1"] = sum(1 for v in c.values() if v == 1)
        row[f"b_{lv}_2"] = sum(1 for v in c.values() if v == 2)
        row[f"b_{lv}_3p"] = sum(1 for v in c.values() if v >= 3)
    return row


def _group(row, levels, k):
    return sum(row[f"b_{lv}_{k}"] for lv in levels)


def label(row: dict, frac: float = 0.7, min_other: int = 3, min_cichlid: int = 1) -> str:
    """Provisional label from the per-level counts."""
    n_o, n_c = _group(row, NON_CICHLID, "n"), row["b_cichlid_n"]
    if n_o < min_other or n_c < min_cichlid:
        return "insufficient_species"
    o1, o2 = _group(row, NON_CICHLID, "1") / n_o, _group(row, NON_CICHLID, "2") / n_o
    c1, c2 = row["b_cichlid_1"] / n_c, row["b_cichlid_2"] / n_c
    if o1 >= frac and c1 >= frac:
        return "single"
    if o2 >= frac and c2 >= frac:
        return "retained_pair"
    if o1 >= frac and c2 >= frac:
        return "cichlid_specific_duplication"
    if o2 >= frac and c1 >= frac:
        return "cichlid_lost_one_copy"
    return "mixed"


def clade_profile(rt: RTree, g: Node, clades, t: float, frac: float = 0.7) -> dict:
    lv = lambda s: level_of_species(clades, s)
    human = {"homo_sapiens"} & set(rt.species)
    other = {s for s in rt.species if s != "homo_sapiens" and lv(s) == "gnathostome"}
    name, a = anchor_clade(rt, g, [("human", human), ("other_gnathostome", other)], t)
    row = {"b_anchor": name or "none"}
    if a is None:
        row["b_label"] = "no_orthologous_anchor"
        return row
    counts = copy_counts(rt, a, lv)
    row.update(summarize(counts))
    row["b_clade_leaves"] = a.hi - a.lo
    row["b_teleost_species_in_tree"] = sum(1 for s in rt.species if lv(s) in TELEOST_LEVELS)
    row["b_label"] = label(row, frac)
    return row
