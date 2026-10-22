"""Orthology calls on one rooted gene tree.

Species overlap rule: at a node, two child subtrees give a *duplication* when their species sets
overlap by more than `t` (overlap = |A and B| / min(|A|, |B|); t = 0 means any shared species),
otherwise a *speciation*. The rule is pairwise between the two children that hold the two genes, so a
multifurcation is handled child against child. Shadowed reference leaves must be dropped before the
tree is built (see species.py).
"""
from __future__ import annotations

from .tree import Node, leaves, postorder

TELEOST_LEVELS = {"teleost", "percomorph", "ovalentaria", "cichlid"}
NON_TELEOST_FISH = "non-teleost actinopterygian"


class RTree:
    """A rooted tree annotated for orthology queries. Every leaf must carry `.sp`."""

    def __init__(self, root: Node):
        self.root = root
        self.leaves = leaves(root)
        sp = sorted({x.sp for x in self.leaves})
        self.species = sp
        self.bit = {s: i for i, s in enumerate(sp)}
        self.by_species: dict[str, list[Node]] = {}
        for i, x in enumerate(self.leaves):
            x.lo, x.hi = i, i + 1
            x.mask = 1 << self.bit[x.sp]
            self.by_species.setdefault(x.sp, []).append(x)
        for n in postorder(root):
            if n.children:
                n.mask = 0
                n.lo = min(c.lo for c in n.children)
                n.hi = max(c.hi for c in n.children)
                for c in n.children:
                    n.mask |= c.mask
        # leaves() is a left-to-right DFS, so every node covers a contiguous interval

    def names_of_mask(self, mask: int) -> list[str]:
        return [s for s, b in self.bit.items() if mask >> b & 1]

    def inside(self, leaf: Node, node: Node) -> bool:
        return node.lo <= leaf.lo < node.hi

    def walk(self, g: Node, t: float):
        """Yield (ancestor, own_child, sibling_child, is_speciation), nearest ancestor first."""
        out = []
        c0 = g
        a = g.parent
        while a is not None:
            for c in a.children:
                if c is c0:
                    continue
                out.append((a, c0, c, is_speciation(c0.mask, c.mask, t)))
            c0 = a
            a = a.parent
        return out

    def leaves_in(self, node: Node):
        return self.leaves[node.lo:node.hi]


def is_speciation(m1: int, m2: int, t: float) -> bool:
    inter = m1 & m2
    if not inter:
        return True
    if t <= 0:
        return False
    return inter.bit_count() / min(m1.bit_count(), m2.bit_count()) <= t


def _path_len(node: Node, ancestor: Node) -> float:
    s = 0.0
    while node is not None and node is not ancestor:
        s += node.length or 0.0
        node = node.parent
    return s


def orthologs(rt: RTree, g: Node, t: float):
    """(orthologous leaves, paralogous leaves) of g in the tree."""
    orth, para = [], []
    for _, _, c, spec in rt.walk(g, t):
        (orth if spec else para).extend(rt.leaves_in(c))
    return orth, para


def relationship(rt: RTree, g: Node, ref: str, own: str, t: float) -> dict:
    """Closest homolog of g in reference species `ref`, and the relationship.
    own = the species of g (to count focal co-orthologs)."""
    r_all = rt.by_species.get(ref, [])
    res = {"rel": "none", "n_ref": 0, "n_focal": 0, "closest": "", "closest_is_ortholog": "",
           "n_ref_in_family": len(r_all), "ref_orthologs": [], "meet": None, "tie": False,
           "focal_ids": set()}
    if not r_all:
        return res
    w = rt.walk(g, t)
    r_orth = [x for (_, _, c, spec) in w if spec for x in rt.leaves_in(c) if x.sp == ref]
    # closest by position: R leaves under the nearest ancestor that has any
    first_a, cand = None, []
    for (a, _, c, spec) in w:
        if first_a is not None and a is not first_a:
            break
        hit = [x for x in rt.leaves_in(c) if x.sp == ref]
        if hit:
            first_a = a
            cand.extend((x, spec) for x in hit)
    orth_cand = [x for (x, spec) in cand if spec]
    pool = orth_cand or [x for (x, _) in cand]
    if len(pool) > 1:
        res["tie"] = True
    pool.sort(key=lambda x: (_path_len(x, first_a), x.name))
    best = pool[0]
    res.update(closest=best.name, closest_is_ortholog=bool(orth_cand), meet=first_a)
    if not r_orth:
        res.update(rel="paralog", n_ref=0, n_focal=0)
        return res
    focal = set()
    for r in r_orth:
        for _, _, c, spec in rt.walk(r, t):
            if spec:
                focal.update(id(x) for x in rt.leaves_in(c) if x.sp == own)
    res["focal_ids"] = focal
    res["ref_orthologs"] = sorted(x.name for x in r_orth)
    n_r, n_f = len(r_orth), len(focal)
    res["n_ref"], res["n_focal"] = n_r, n_f
    res["rel"] = f"{'1' if n_r == 1 else 'many'}:{'1' if n_f == 1 else 'many'}"
    return res


def teleost_pattern(rt: RTree, g: Node, own: str, t: float, level_of) -> dict:
    """Copies of the gene's orthologs per teleost species (species in teleost-clade levels)."""
    orth, _ = orthologs(rt, g, t)
    count = {}
    for x in orth:
        count[x.sp] = count.get(x.sp, 0) + 1
    count[own] = count.get(own, 0) + 1          # g itself
    tel = [s for s in rt.species if level_of(s) in TELEOST_LEVELS]
    n1 = sum(1 for s in tel if count.get(s, 0) == 1)
    n2 = sum(1 for s in tel if count.get(s, 0) >= 2)
    n0 = sum(1 for s in tel if count.get(s, 0) == 0)
    return {"tel_species": len(tel), "tel_one": n1, "tel_two_plus": n2, "tel_zero": n0,
            "tel_frac_one_of_present": round(n1 / (n1 + n2), 4) if (n1 + n2) else None}
