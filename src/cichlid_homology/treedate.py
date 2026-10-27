"""Tree-based dating of the duplications of a family (0.12): the third line of evidence for B1 and B2.

Counts per species cannot tell "a pair at the teleost root, with early losses" from "a new duplicate in the percomorphs"
(rows 3 and 4 of docs/NOTES.md), nor reciprocal losses from a lost duplicate. The shape of the tree can, in part. A rooted
gene tree is reconciled with the master species tree: a node is a duplication when the species-tree node that is the
lowest common ancestor of its species (M) equals that of one of its children, otherwise a speciation. Only jawed-vertebrate
species count (the same scope as the parsimony). In a placement tree (TreeFam, Compara) the placed project genes are not
used to find duplications, because genes on one branch cannot be ordered: the structure comes from the reference leaves,
and each project gene is assigned to the lineage it was placed in. In an OrthoFinder tree (denovo) every leaf counts.

What the runs on Compara taught (the numbers are on the 2,637 published-pair positives and 8,192 negatives):
  - A duplication by reconciliation alone is not evidence: the gene tree often disagrees with the master species tree (which
    has a polytomy among the early teleosts), and every node of that kind is called a duplication. Two lineages must also
    share species: `d_overlap` (shared / all, of the two largest children) is bimodal, with a valley between 0.3 and 0.5, and
    0.5 (`pair_overlap`) sits in it. With the overlap rule at the root only, 56.5% of positives and 0.3% of negatives were pairs.
  - The top of the gene's teleost clade is often where the tree disagrees with the species tree, and the real duplication is a
    node below it. So the pair node is the DEEPEST duplication with enough overlap on the path from the gene up through its
    teleost clade, not the top of the clade: 82.4% of positives and 0.5% of negatives (thresholds on counts: 84.0% and 0.8%;
    thresholds or tree: 91.0% and 0.9%).
  - A stray early-teleost gene inside one percomorph lineage must not turn a percomorph duplicate into a teleost pair, so the
    level of a duplication is `teleost` only if early teleosts are in both lineages, or an early teleost has genes on both
    sides, or every early-teleost gene of the clade is inside the node (a pair whose second copy was lost in the early teleosts:
    a later duplicate would leave the early genes outside); otherwise it is the level of the species the two sides share.

For a gene g (the same anchor clade as B1):
  tr    the highest ancestor of g, inside the anchor clade, that holds no non-teleost jawed vertebrate: g's own teleost clade
  d_root  tr itself: `pair` (a supported duplication), `thin` (the lineages share some species, less than half), `disjoint`
        (reconciliation says duplication, no species in common: the gene tree disagrees with the species tree, or reciprocal
        losses; counts cannot tell these apart), `single` (a speciation)
  d_root_level  the ladder level of the species at tr
  d_dup_level  the level of the deepest supported duplication on the path from g to tr (clade-level: two or more species):
        teleost | percomorph | ovalentaria | cichlid; d_dup_at_root yes if that node is tr itself
  d_support, d_lineages, d_overlap  at that node (or at tr when there is no supported duplication): `both_early` (early teleosts
        in two lineages), `shared_early` (an early teleost has genes on both sides), `early_within` (all the early-teleost
        genes of the clade are inside the node: a pair whose second copy was lost in the early teleosts), `one_early`, `none_early`; the number of
        lineages; shared species / all species of the two largest
  d_gene_lineage  A, B, ... the lineage of that node that holds g (letters follow the order in the rooted tree and are
        arbitrary); d_focal_lineages the number of lineages with g's species
  d_later_levels  every supported duplication on the path, nearest first, as `level:species`; d_species_dups the number of
        supported duplications inside one species; d_first_later_level the level of the deepest one when it is a later duplicate
  d_carriers  (denovo only) the cichlid species that hold both copies below the deepest later duplication
  d_call  pair (the deepest supported duplication is at the teleost level) | pair_unresolved_depth (the gene's whole teleost
        clade is a duplication with no early teleosts: a pair kept only in the percomorphs or a later duplicate) |
        later_duplicate_before_cichlids | later_duplicate_in_cichlids | thin_pair | lost_duplicate (in a placement tree
        duplications among project genes are invisible, so `lost_duplicate` means none among the reference species)"""
from __future__ import annotations

import bisect
from collections import Counter

from .parsimony import JAWED
from .tree import postorder

RANK = {lv: i for i, lv in enumerate(JAWED)}            # gnathostome 0, teleost 1, percomorph 2, ovalentaria 3, cichlid 4
D_COLS = ["d_status", "d_call", "d_root", "d_root_level", "d_dup_level", "d_dup_at_root", "d_support", "d_lineages",
          "d_overlap", "d_gene_lineage", "d_focal_lineages", "d_later_levels", "d_first_later_level", "d_species_dups",
          "d_carriers", "d_carrier_n"]
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
        self.M, self.smask, self.minrank, self._ov = {}, {}, {}, {}
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

    def overlap(self, n) -> float:
        """Shared species / all species of the two largest informative children of n (0 when it has fewer than two)."""
        k = id(n)
        if k not in self._ov:
            inf = [c for c in n.children if self.M[id(c)] is not None]
            v = 0.0
            if len(inf) >= 2:
                big = sorted(inf, key=lambda c: -self.n_species(c))[:2]
                union = _bits(self.smask[id(big[0])] | self.smask[id(big[1])])
                v = _bits(self.smask[id(big[0])] & self.smask[id(big[1])]) / union if union else 0.0
            self._ov[k] = v
        return self._ov[k]

    def _early_positions(self, tr):
        k = ("early", id(tr))
        if k not in self._ov:
            self._ov[k] = sorted(x.lo for x in self.rt.leaves_in(tr)
                                 if self.M[id(x)] is not None and self.minrank[id(x)] == 1)
        return self._ov[k]

    def evaluate(self, n, tr):
        """Overlap, level and support of the duplication at n (the two largest informative children), within teleost clade tr.

        The level is `teleost` when early teleosts are in two lineages, or an early teleost has genes on both sides, or
        every early-teleost gene of the clade is inside the node (a pair whose second copy was lost in all the early
        teleosts: a later duplicate would leave the early genes outside); a stray early gene inside one lineage while others
        sit outside does not make it one. Otherwise it is the level of the species the two sides share."""
        k = ("ev", id(n), id(tr))
        if k not in self._ov:
            inf = [c for c in n.children if self.M[id(c)] is not None and self.minrank[id(c)] >= 1]
            big = sorted(inf, key=lambda c: -self.n_species(c))[:2]
            ov, shared_rank, n_early = 0.0, None, sum(1 for c in inf if self.smask[id(c)] & self.early)
            if len(big) == 2:
                inter = self.smask[id(big[0])] & self.smask[id(big[1])]
                union = _bits(self.smask[id(big[0])] | self.smask[id(big[1])])
                ov = _bits(inter) / union if union else 0.0
                ranks = [r for sp, r in self.rank.items() if inter >> self.rt.bit[sp] & 1]
                shared_rank = min(ranks) if ranks else None
            pos = self._early_positions(tr)
            inside = bisect.bisect_left(pos, n.hi) - bisect.bisect_left(pos, n.lo)
            outside = len(pos) - inside
            if n_early >= 2:
                level, support = "teleost", "both_early"
            elif shared_rank == 1:
                level, support = "teleost", "shared_early"
            elif inside >= 1 and outside == 0:
                level, support = "teleost", "early_within"
            elif shared_rank is not None:
                level, support = JAWED[shared_rank], ("one_early" if inside else "none_early")
            else:
                level, support = None, ("one_early" if inside else "none_early")
            self._ov[k] = {"ov": ov, "level": level, "support": support, "lineages": inf, "nsp": self.n_species(n)}
        return self._ov[k]

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


def treedate_profile(rt, g, anchor, ctx, mode: str, prefixes=None, pair_overlap: float = 0.5) -> dict:
    name, a = anchor
    out = {c: "" for c in D_COLS}
    if a is None:
        out["d_status"] = "no_anchor"
        return out
    T = _tree_for(rt, ctx, mode, prefixes)
    tr = g
    while tr.parent is not None and tr is not a and T.minrank[id(tr.parent)] != 0:
        tr = tr.parent
    if not tr.children:
        out["d_status"] = "one_teleost_gene"
        return out
    if T.minrank[id(tr)] is None:
        out["d_status"] = "no_teleost"
        return out
    if T.n_species(tr) < 2:
        out["d_status"] = "one_teleost_species"
        return out
    out["d_status"] = "ok"
    out["d_root_level"] = T.level(tr)
    root = "single"
    if T.dup(tr):
        ev = T.evaluate(tr, tr)
        root = "pair" if ev["ov"] >= pair_overlap else ("thin" if ev["ov"] > 0 else "disjoint")
    out["d_root"] = root
    supported, v = [], g.parent
    while v is not None:
        if T.dup(v):
            ev = T.evaluate(v, tr)
            if ev["ov"] >= pair_overlap and ev["level"] is not None:
                supported.append((v, ev))
        if v is tr:
            break
        v = v.parent
    out["d_later_levels"] = ";".join(f"{e['level']}:{e['nsp']}" for _, e in supported)
    out["d_species_dups"] = sum(1 for _, e in supported if e["nsp"] == 1)
    clade = [(n, e) for n, e in supported if e["nsp"] >= 2]
    top, tev = clade[-1] if clade else (None, None)
    node = top if top is not None else (tr if root != "single" else None)
    if node is not None:
        e = T.evaluate(node, tr)
        out["d_support"], out["d_lineages"], out["d_overlap"] = e["support"], len(e["lineages"]), round(e["ov"], 2)
        c0 = g
        while c0.parent is not node:
            c0 = c0.parent
        out["d_gene_lineage"] = LETTERS[e["lineages"].index(c0)] if c0 in e["lineages"] else "unassigned"
        out["d_focal_lineages"] = sum(1 for c in node.children if any(x.sp == g.sp for x in rt.leaves_in(c)))
    else:
        out["d_focal_lineages"] = sum(1 for c in tr.children if any(x.sp == g.sp for x in rt.leaves_in(c)))
    if top is None:
        out["d_call"] = "thin_pair" if root == "thin" else "lost_duplicate"
        return out
    out["d_dup_level"], out["d_dup_at_root"] = tev["level"], "yes" if top is tr else "no"
    if tev["level"] == "teleost":
        out["d_call"] = "pair"
    elif top is tr:
        out["d_call"] = "pair_unresolved_depth"
    else:
        out["d_call"] = "later_duplicate_in_cichlids" if tev["level"] == "cichlid" else "later_duplicate_before_cichlids"
        out["d_first_later_level"] = tev["level"]
        if mode == "denovo":
            cnt = Counter()
            for c in top.children:
                for sp in {x.sp for x in rt.leaves_in(c)}:
                    if sp in T.rank and T.rank[sp] == 4:
                        cnt[sp] += 1
            car = sorted(sp for sp, k in cnt.items() if k >= 2)
            out["d_carriers"], out["d_carrier_n"] = ";".join(car), len(car)
    return out


def treedate_summary(rows) -> dict:
    C = Counter
    ok = [r for r in rows if r.get("d_status") == "ok"]
    out = {"genes": len(rows), "status": dict(C(r.get("d_status", "") for r in rows)), "scored": len(ok)}
    if not ok:
        return out
    out["calls"] = dict(C(r["d_call"] for r in ok).most_common())
    out["pair_support"] = dict(C(r["d_support"] for r in ok if r["d_call"] == "pair"))
    out["pair_root_level"] = dict(C(r["d_dup_at_root"] for r in ok if r["d_call"] == "pair"))
    pairs = [r for r in ok if r["d_call"] == "pair"]
    ovs = [float(r["d_overlap"]) for r in ok if r["d_overlap"] not in ("", None)]
    out["root_classes"] = dict(C(r["d_root"] for r in ok).most_common())
    out["overlap_histogram_tenths"] = [sum(1 for x in ovs if a <= x < (b if b < 1 else 1.01)) for a, b in
                                       zip([0, 0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], [0.001, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])]
    out["pairs_with_focal_species_in_two_or_more_lineages"] = sum(1 for r in pairs if int(r["d_focal_lineages"] or 0) >= 2)
    out["by_category"] = {}
    for r in ok:
        k = r.get("t_category")
        if k:
            out["by_category"].setdefault(k, C())[r["d_call"]] += 1
    out["by_category"] = {k: dict(v.most_common()) for k, v in out["by_category"].items()}
    return out
