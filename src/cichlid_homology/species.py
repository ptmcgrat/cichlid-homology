"""Species of the leaves of a family, with the decided rule for reference leaves that are the same
genome as one of our placed species ("shadowed" species).

Rule (coordinator v0.13.0, decided; the owner may overrule): ignore. Before species overlap is
computed, drop the reference leaves of shadowed species and let our placed genes stand for the
species. The Ensembl copy and our RefSeq gene are two annotations of one genome. Record the rule
in run.json (SHADOW_RULE).
"""
from __future__ import annotations

SHADOW_RULE = "drop reference leaves of shadowed species; placed genes stand for the species"


def drop_shadowed(species_of: dict, shadowed: set, is_project_gene) -> tuple[dict, list]:
    """species_of: {leaf: canonical species or None}. Returns (kept, dropped).
    A leaf is dropped when it is NOT one of our genes and its species is shadowed.
    Leaves with unknown species (None) are kept; the caller decides what to do with them."""
    kept, dropped = {}, []
    for leaf, sp in species_of.items():
        if sp is not None and sp in shadowed and not is_project_gene(leaf):
            dropped.append(leaf)
        else:
            kept[leaf] = sp
    return kept, dropped


def leaf_species_map(ls, leaves) -> dict:
    """{leaf: species or None} using LeafSpecies.species(leaf) (None for an unknown leaf)."""
    return {leaf: ls.species(leaf) for leaf in leaves}


def load_leaf_species(layout, part: str, only=None):
    """LeafSpecies for a part. `only` (a set of leaves) limits the large mapping file; if the
    installed coordinator rejects the argument, load without it."""
    from cichlid_data import LeafSpecies

    if only is not None:
        try:
            return LeafSpecies.load(layout, part, only=set(only))
        except TypeError:
            pass
    return LeafSpecies.load(layout, part)
