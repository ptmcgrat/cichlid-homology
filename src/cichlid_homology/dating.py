"""Dating a duplication on the clade ladder, and telling retained teleost ohnologs from recent
duplications. The ladder level of a duplication node is `common_level` of the species below it."""
from __future__ import annotations

from dataclasses import dataclass

from .orthology import NON_TELEOST_FISH, TELEOST_LEVELS, RTree
from .tree import Node


@dataclass
class Params:
    overlap_t: float = 0.0          # species-overlap threshold (settle in the pilot)
    min_levels_each_side: int = 2   # "several teleost levels" on each side of a retained ohnolog pair
    mostly_one_to_one: float = 0.8  # fraction of teleost species with exactly one copy (pilot)
    k_roots: int = 3                # tied rootings checked for stability


def level_of_species(clades, sp):
    try:
        return clades.level(sp)
    except Exception:
        return None


def level_of_set(clades, species: list[str]):
    try:
        return clades.common_level(species)
    except Exception:
        return None


def date_node(rt: RTree, node: Node, clades) -> str | None:
    """Ladder level of the last common ancestor of the species below a node.
    The disputed level is reported as 'cyclostome_disputed' and never as a date."""
    lv = level_of_set(clades, rt.names_of_mask(node.mask))
    return "cyclostome_disputed" if lv == "cyclostome" else lv


def _is_focal(leaf: Node, prefixes) -> bool:
    head, sep, _ = (leaf.name or "").partition("_")
    return bool(sep) and head in prefixes


def birth_of_copy(rt: RTree, g: Node, co_orthologs: set, clades, prm: Params,
                  tree_mode: str, prefixes) -> dict | None:
    """Nearest duplication separating g from another co-ortholog (a leaf whose id is in
    `co_orthologs`, g excluded). Returns None if g has no co-ortholog."""
    for (a, c0, c, _) in rt.walk(g, prm.overlap_t):
        hit = [x for x in rt.leaves_in(c) if id(x) in co_orthologs and x is not g]
        if not hit:
            continue
        level = date_node(rt, a, clades)
        side_a = rt.names_of_mask(c0.mask)
        side_b = rt.names_of_mask(c.mask)
        lv_a = {level_of_species(clades, s) for s in side_a} & TELEOST_LEVELS
        lv_b = {level_of_species(clades, s) for s in side_b} & TELEOST_LEVELS
        # place mode cannot order focal genes among themselves: if both sides hold only our own
        # genes, the node is a multifurcation artefact, not a resolved duplication
        only_focal = all(_is_focal(x, prefixes) for x in rt.leaves_in(c0)) and \
            all(_is_focal(x, prefixes) for x in rt.leaves_in(c))
        resolved = not (tree_mode == "place" and only_focal)
        out_sp = {}
        for s in rt.species:
            try:
                if clades.lineage(s) == NON_TELEOST_FISH:
                    out_sp[s] = len(rt.by_species[s])
            except Exception:
                pass
        outgroup_single = bool(out_sp) and max(out_sp.values()) <= 1
        if not resolved:
            cls = "unresolved_place_mode"
        elif level == "teleost" and len(lv_a) >= prm.min_levels_each_side \
                and len(lv_b) >= prm.min_levels_each_side:
            cls = "teleost_ohnolog" if outgroup_single else (
                "teleost_ohnolog_no_outgroup" if not out_sp else "teleost_level_outgroup_multi")
        elif level == "teleost":
            cls = "teleost_level_not_spread"
        elif level in ("percomorph", "ovalentaria", "cichlid"):
            cls = f"lineage_specific:{level}"
        elif level == "cyclostome_disputed":
            cls = "cyclostome_disputed"
        elif level in ("gnathostome", "outside_vertebrates"):
            cls = f"older_than_teleost:{level}"
        else:
            cls = "unknown_level"
        return {"level": level, "class": cls, "resolved": resolved,
                "levels_side_a": len(lv_a), "levels_side_b": len(lv_b),
                "outgroup_species": len(out_sp), "outgroup_single": outgroup_single,
                "n_species_below": bin(a.mask).count("1")}
    return None
