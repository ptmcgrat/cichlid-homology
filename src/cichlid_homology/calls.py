"""All per-gene calls on one rooted tree (pure; no coordinator needed)."""
from __future__ import annotations

from .dating import Params, birth_of_copy, date_mask, date_node, level_of_species
from .orthology import RTree, relationship, teleost_pattern
from .clade_counts import clade_counts, group
from .profile import clade_profile_v2
from .tree import Node

REFS = {"human": "homo_sapiens", "zebrafish": "danio_rerio"}


def call_gene(rt: RTree, g: Node, own: str, prm: Params, clades, tree_mode: str, prefixes) -> dict:
    row = {"parent_children": len(g.parent.children) if g.parent is not None else 0,
           "species_in_tree": len(rt.species)}
    human_node = None
    for key, ref in REFS.items():
        c = clade_counts(rt, g, ref, own)
        row[f"a_{key}_cat"] = c["cat"]
        row[f"a_{key}_group"] = group(key, c["cat"])
        row[f"_a_{key}_focal_leaves"] = c["focal_leaves"]
        row[f"a_{key}_ref_n"] = c["ref_n"]
        row[f"a_{key}_focal_n"] = c["focal_n"]
        row[f"a_{key}_clade_leaves"] = c["clade_leaves"]
        row[f"a_{key}_best"] = c["best"]
        row[f"a_{key}_best_dist"] = c["best_dist"]
        row[f"a_{key}_tie"] = c["tie"]
        row[f"_a_{key}_ref_leaves"] = c["ref_leaves"]
        if key == "human":
            human_node = c["node"]
        r = relationship(rt, g, ref, own, prm.overlap_t)
        row[f"{key}_relationship"] = r["rel"]
        row[f"{key}_n_ref"] = r["n_ref"]
        row[f"{key}_n_focal"] = r["n_focal"]
        row[f"{key}_closest"] = r["closest"]
        row[f"{key}_closest_is_ortholog"] = r["closest_is_ortholog"]
        row[f"{key}_n_ref_in_family"] = r["n_ref_in_family"]
        row[f"_{key}_orth_leaves"] = r["ref_orthologs"]
        row[f"{key}_closest_tie"] = r["tie"]
        if r["rel"] == "paralog" and r["meet"] is not None:
            row[f"{key}_paralog_level"] = date_node(rt, r["meet"], clades)
            c0, c = r["meet_sides"]
            row[f"{key}_paralog_level_overlap"] = date_mask(rt, c0.mask & c.mask, clades) or ""
        if r["n_focal"] > 1:
            b = birth_of_copy(rt, g, r["focal_ids"], clades, prm, tree_mode, prefixes)
            if b:
                row[f"{key}_dup_level"] = b["level"]
                row[f"{key}_dup_class"] = b["class"]
                row[f"{key}_dup_resolved"] = b["resolved"]
                row[f"{key}_dup_levels_sides"] = f"{b['levels_side_a']}/{b['levels_side_b']}"
                row[f"{key}_dup_outgroup_single"] = b["outgroup_single"]
                row[f"{key}_dup_outgroup_found"] = b["outgroup_found"]
                row[f"{key}_dup_level_overlap"] = b["level_overlap"] or ""
                row[f"{key}_dup_ref_levels_sides"] = f"{b['ref_levels_a']}/{b['ref_levels_b']}"
    row.update(clade_profile_v2(rt, g, clades, prm, human_node))
    tp = teleost_pattern(rt, g, own, prm.overlap_t, lambda s: level_of_species(clades, s))
    row.update({f"teleost_{k}": v for k, v in tp.items()})
    return row
