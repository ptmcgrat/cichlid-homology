"""A1 and A2 as counts (version 0.7.0 pilot).

For a focal-species gene g and a reference species R (human, then zebrafish):
  1. the clade is the smallest clade of the rooted tree that contains g and at least one R gene;
  2. count the R genes in it, and the genes of g's own species in it (other cichlids are not counted);
  3. the best R gene is the one with the shortest path from g (ties are flagged and broken by name).
No duplication/speciation call is made, so one stray species cannot demote an ortholog.

Category (reference genes : focal genes):
    OneToOne 1:1   OneToTwo 1:2   OneToMany 1:3 or more
    ManyToOne 2+:1   ManyToMany 2+:2+      None = no gene of R in this method's tree.
"""
from __future__ import annotations

import bisect

from .orthology import RTree, _path_len
from .tree import Node

CATS = ("OneToOne", "OneToTwo", "OneToMany", "ManyToOne", "ManyToMany", "None")
TOL = 1e-9


def category(n_ref: int, n_focal: int) -> str:
    if n_ref == 0:
        return "None"
    if n_ref == 1:
        return "OneToOne" if n_focal == 1 else "OneToTwo" if n_focal == 2 else "OneToMany"
    return "ManyToOne" if n_focal == 1 else "ManyToMany"


def smallest_clade(rt: RTree, g: Node, species: set) -> Node | None:
    """Nearest ancestor of g whose subtree holds a leaf of any of `species`."""
    pos = sorted(x.lo for sp in species for x in rt.by_species.get(sp, []))
    if not pos:
        return None
    a = g.parent
    while a is not None:
        i = bisect.bisect_left(pos, a.lo)
        if i < len(pos) and pos[i] < a.hi:
            return a
        a = a.parent
    return None


def clade_counts(rt: RTree, g: Node, ref: str, own: str) -> dict:
    out = {"cat": "None", "ref_n": 0, "focal_n": 0, "clade_leaves": 0, "best": "", "best_dist": "",
           "tie": False, "ref_leaves": [], "node": None}
    a = smallest_clade(rt, g, {ref})
    if a is None:
        return out
    ref_in = [x for x in rt.by_species.get(ref, []) if a.lo <= x.lo < a.hi]
    focal_in = [x for x in rt.by_species.get(own, []) if a.lo <= x.lo < a.hi]
    dg = _path_len(g, a)
    scored = sorted((dg + _path_len(x, a), x.name) for x in ref_in)
    best_d, best = scored[0]
    out.update(cat=category(len(ref_in), len(focal_in)), ref_n=len(ref_in), focal_n=len(focal_in),
               clade_leaves=a.hi - a.lo, best=best, best_dist=round(best_d, 6),
               tie=sum(1 for d, _ in scored if abs(d - best_d) <= TOL) > 1,
               ref_leaves=sorted(x.name for x in ref_in), node=a)
    return out
