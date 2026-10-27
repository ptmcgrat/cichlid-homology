"""Tree-based dating of the duplications of a family (0.11.0): the third line of evidence for B1 and B2.

Counts per species cannot tell "a pair at the teleost root, with early losses" from "a new duplicate in the percomorphs"
(rows 3 and 4 of docs/NOTES.md), nor reciprocal losses from a lost duplicate. The shape of the tree can. Here a rooted
gene tree is reconciled with the master species tree: a node is a *duplication* when the species-tree node that is the
lowest common ancestor of its species (M) equals that of one of its children, otherwise a speciation. No share
thresholds are used. Only jawed-vertebrate species count (the same scope as the parsimony). In a placement tree (TreeFam,
Compara) the placed project genes are not used to find duplications, because genes on one branch cannot be ordered:
the structure comes from the reference leaves, and each project gene is assigned to the lineage it was placed in.
In an OrthoFinder tree (denovo) every leaf counts.

For a gene g in its anchor clade (the same clade as B1):
  tr    the lowest node above every teleost leaf of the clade (the teleost root of this family)
  d_root  `pair` if tr is a duplication (the teleost genes form two or more lineages), `single` if a speciation
  d_root_level  the ladder level of the species at tr: `teleost` when early teleosts are there (a duplication
        at that level is a teleost-genome-duplication pair), `percomorph` or lower when they are not (then the
        tree cannot tell a teleost pair lost in the early teleosts from a later duplicate)
  d_support  for a pair, `both_early` (two lineages hold early teleosts), `one_early` or `none_early`
  d_overlap  for a pair, the species the two largest lineages have in common as a share of all the species in either
        (shared / union): near 1 for a fully retained pair, 0 when no species has both (reciprocal losses), small when
        one lineage is a few genes (copies lost in most species, or one stray gene in the wrong clade that makes a
        duplication that is not there): read a pair with a small overlap with care
  d_gene_lineage  A, B, ... the lineage of tr that holds g (the letters follow the order in the rooted tree and are arbitrary); d_focal_lineages the number of lineages with g's species
  d_later_levels  the duplications below tr on the path to g, nearest first, as `level:species` (species = number
        of species below it); d_first_later_level the deepest one that involves two or more species
  d_carriers  (denovo only) the cichlid species that hold both copies below the deepest later duplication
  d_call  pair | lost_duplicate | later_duplicate_in_cichlids | later_duplicate_before_cichlids | pair_unresolved_depth
        (in a placement tree duplications among project genes are invisible, so `lost_duplicate` means none among
        the reference species below the teleost root)"""
from __future__ import annotations

from collections import Counter

from .parsimony import JAWED
from .tree import postorder

RANK = {lv: i for i, lv in enumerate(JAWED)}            # gnathostome 0, teleost 1, percomorph 2, ovalentaria 3, cichlid 4
D_COLS = ["d_status", "d_call", "d_root", "d_root_level", "d_support", "d_lineages", "d_overlap", "d_gene_lineage",
          "d_focal_lineages", "d_later_levels", "d_first_later_level", "d_species_dups", "d_carriers", "d_carrier_n"]
LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


def _bits(m: int) -> int:
    return bin(m).count("1")


class _Tree:
    """Per-tree facts, computed once: the species-tree node M of every node, its informative species and levels."""

    def __init__(self, rt, ctx, mode, prefixes):
        self.rt, self.ctx = rt, ctx
        st = ctx.st
        rank = {sp: RANK[ctx.level(sp)] for sp in rt.species if ctx.in_scope(sp)}
        self.rank = rank
        self.early = 0
        for sp, r in rank.items():
            if r == 1:
                self.early |= 1 << rt.bit[sp]
        self.M, self.smask, self.minrank = {}, {}, {}
        for n in postorder(rt.root):
            k = id(n)
            if not n.children:
                project = prefixes and n.name.split("_")[0] in prefixes
                ok = n.sp in rank and (mode == "denovo" or not project)
                self.M[k] = st.id_of(n.sp) if ok else None
                self.smask[k] = (1 << rt.bit[n.sp]) if ok else 0
                self.minrank[k] = rank[n.sp] if ok else None
                continue
            m, mask, mr = None, 0, None
            for c in n.children:
                mc = self.M[id(c)]
                if mc is None:
                    continue
                m = mc if m is None else st.lca(m, mc)
                mask |= self.smask[id(c)]
                r = self.minrank[id(c)]
                mr = r if mr is None else min(mr, r)
            self.M[k], self.smask[k], self.minrank[k] = m, mask, mr

    def dup(self, n) -> bool:
        if not n.children or self.M[id(n)] is None:
            return False
        inf = [c for c in n.children if self.M[id(c)] is not None]
        return len(inf) >= 2 and any(self.M[id(c)] == self.M[id(n)] for c in inf)

    def level(self, n):
        r = self.minrank[id(n)]
        return None if r is None else JAWED[r]

    def n_species(self, n) -> int:
        return _bits(self.smask[id(n)])


def _tree_for(rt, ctx, mode, prefixes) -> _Tree:
    key = (mode, tuple(sorted(prefixes or ())))
    cache = getattr(rt, "_treedate", None)
    if cache is None or cache[0] != key:
        cache = (key, _Tree(rt, ctx, mode, prefixes))
        rt._treedate = cache
    return cache[1]


def treedate_profile(rt, g, anchor, ctx, mode: str, prefixes=None) -> dict:
    name, a = anchor
    out = {c: "" for c in D_COLS}
    if a is None:
        out["d_status"] = "no_anchor"
        return out
    T = _tree_for(rt, ctx, mode, prefixes)
    tl = [x for x in rt.leaves_in(a) if T.M[id(x)] is not None and T.minrank[id(x)] >= 1]
    if not tl:
        out["d_status"] = "no_teleost"
        return out
    lo, hi = min(x.lo for x in tl), max(x.hi for x in tl)
    tr = min(tl, key=lambda x: x.lo)
    while tr.parent is not None and not (tr.lo <= lo and tr.hi >= hi):
        tr = tr.parent
    if not tr.children:
        out["d_status"] = "one_teleost_gene"
        return out
    if T.minrank[id(tr)] == 0:
        out["d_status"] = "teleosts_not_monophyletic"
        return out
    if T.n_species(tr) < 2:
        out["d_status"] = "one_teleost_species"
        return out
    out["d_status"] = "ok"
    is_dup = T.dup(tr)
    lineages = [c for c in tr.children if T.M[id(c)] is not None and T.minrank[id(c)] >= 1]
    out["d_root"] = "pair" if is_dup else "single"
    out["d_root_level"] = T.level(tr)
    if is_dup:
        n_early = sum(1 for c in lineages if T.smask[id(c)] & T.early)
        out["d_lineages"] = len(lineages)
        out["d_support"] = "both_early" if n_early >= 2 else ("one_early" if n_early == 1 else "none_early")
        big = sorted(lineages, key=lambda c: -T.n_species(c))[:2]
        union = _bits(T.smask[id(big[0])] | T.smask[id(big[1])])
        out["d_overlap"] = round(_bits(T.smask[id(big[0])] & T.smask[id(big[1])]) / union, 2) if union else ""
    c0 = g
    while c0.parent is not None and c0.parent is not tr:
        c0 = c0.parent
    inside = c0.parent is tr
    if is_dup and inside and c0 in lineages:
        out["d_gene_lineage"] = LETTERS[lineages.index(c0)]
    elif is_dup:
        out["d_gene_lineage"] = "unassigned"
    out["d_focal_lineages"] = sum(1 for c in tr.children if any(x.sp == g.sp for x in rt.leaves_in(c)))
    later, v, deepest = [], g.parent, None
    while v is not None and v is not tr and inside:
        if T.dup(v):
            later.append((T.level(v), T.n_species(v)))
            if T.n_species(v) >= 2:
                deepest = v
        v = v.parent
    out["d_later_levels"] = ";".join(f"{lv}:{n}" for lv, n in later)
    out["d_species_dups"] = sum(1 for _, n in later if n == 1)
    clade_level = [lv for lv, n in later if n >= 2]
    first = T.level(deepest) if deepest is not None else ""
    out["d_first_later_level"] = first
    if mode == "denovo" and deepest is not None:
        cnt = Counter()
        for c in deepest.children:
            for sp in {x.sp for x in rt.leaves_in(c)}:
                if sp in T.rank and T.rank[sp] == 4:
                    cnt[sp] += 1
        car = sorted(sp for sp, k in cnt.items() if k >= 2)
        out["d_carriers"], out["d_carrier_n"] = ";".join(car), len(car)
    if is_dup:
        out["d_call"] = "pair" if out["d_root_level"] == "teleost" else "pair_unresolved_depth"
    elif clade_level:
        out["d_call"] = "later_duplicate_in_cichlids" if first == "cichlid" else "later_duplicate_before_cichlids"
    else:
        out["d_call"] = "lost_duplicate"
    return out


def treedate_summary(rows) -> dict:
    C = Counter
    ok = [r for r in rows if r.get("d_status") == "ok"]
    out = {"genes": len(rows), "status": dict(C(r.get("d_status", "") for r in rows)), "scored": len(ok)}
    if not ok:
        return out
    out["calls"] = dict(C(r["d_call"] for r in ok).most_common())
    out["pair_support"] = dict(C(r["d_support"] for r in ok if r["d_root"] == "pair"))
    out["pair_root_level"] = dict(C(r["d_root_level"] for r in ok if r["d_root"] == "pair"))
    pairs = [r for r in ok if r["d_call"] == "pair"]
    out["pairs_with_focal_species_in_two_or_more_lineages"] = sum(1 for r in pairs if int(r["d_focal_lineages"] or 0) >= 2)
    out["by_category"] = {}
    for r in ok:
        k = r.get("t_category")
        if k:
            out["by_category"].setdefault(k, C())[r["d_call"]] += 1
    out["by_category"] = {k: dict(v.most_common()) for k, v in out["by_category"].items()}
    return out
