"""Lost duplicate or retained ohnolog by parsimony on the species tree (0.9.0 pilot; no share thresholds).

For a clade, count the genes of each jawed-vertebrate species. Explain those counts with copy-number
changes on the branches of the species tree: a gain of one copy costs `dup`, a loss of one copy costs
`loss`; a lineage with no copy cannot regain one (no gene is born from nothing). The deepest ancestor is
assumed to have one gene (a conserved single-copy clade). Then ask which number of copies the teleost
ancestor had, 1, 2 or 3, and take the history that needs the fewest events:

    1 copy   lost_duplicate     the teleost duplicate was lost early (later gains are separate events)
    2 copies retained_ohnolog   both teleost copies were kept
    3 copies multi_copy

The margin between the best and the second-best history is the confidence: a tie is `ambiguous`, and a call
that changes when the relative cost of a loss and a gain changes (1:1, 1:2, 2:1) is `ambiguous` as well.
Species with no gene in the clade are counted as losses (`zero_mode` loss) or ignored (`missing`); both
are scored. Scope: species that are jawed vertebrates in the clade table (the lamprey and the invertebrates
are left out). The events of the best history are reported: where gains and losses fall on the tree."""
from __future__ import annotations

from collections import Counter, defaultdict

from .dating import level_of_set, level_of_species
from .orthology import TELEOST_LEVELS

INF = 10 ** 6
JAWED = ("gnathostome", "teleost", "percomorph", "ovalentaria", "cichlid")
LABEL = {1: "lost_duplicate", 2: "retained_ohnolog", 3: "multi_copy"}
COST_SETS = ((1, 1), (1, 2), (2, 1))          # (cost of a gain, cost of a loss)
# calibration grid (optional): (cost of a gain, cost of a loss, teleost duplication free)
GRID = tuple((g, 1, free) for free in (False, True) for g in (1, 2, 3, 4, 6, 8))


def grid_name(g: int, l: int, free: bool) -> str:
    return f"g{g}l{l}" + ("f" if free else "")


GRID_COLS = [f"pg_{grid_name(*x)}" for x in GRID] + [f"pgm_{grid_name(*x)}" for x in GRID]


def trans_matrix(K: int, dup: int, loss: int):
    t = [[INF] * (K + 1) for _ in range(K + 1)]
    for p in range(K + 1):
        for c in range(K + 1):
            if p == 0:
                t[p][c] = 0 if c == 0 else INF          # nothing is born from nothing
            elif c >= p:
                t[p][c] = dup * (c - p)
            else:
                t[p][c] = loss * (p - c)
    return t


def trans_matrix_wgd(K: int, dup: int, loss: int):
    """The branch above the teleost ancestor: one extra copy is free (the genome duplication doubles every
    gene); further gains and all losses cost as usual."""
    t = trans_matrix(K, dup, loss)
    for p in range(1, K + 1):
        for c in range(p + 1, K + 1):
            t[p][c] = dup * (c - p - 1)
    return t


class PTree:
    """The species tree pruned to a set of species (unary nodes removed), ready for the DP."""
    __slots__ = ("children", "parent", "species", "order", "root", "T", "under", "node_of_species")


def _build_ptree(st, st_children, species, teleost):
    id2sp = {st.id_of(s): s for s in species}
    if len(id2sp) < 2:
        return None
    cnt = defaultdict(int)
    for i in id2sp:
        v = i
        while v != -1:
            cnt[v] += 1
            v = st.parent[v]
    inc = lambda v: [c for c in st_children[v] if cnt.get(c, 0) > 0]

    def collapse(v):
        while v not in id2sp and len(inc(v)) == 1:
            v = inc(v)[0]
        return v

    nodes, children, parent = [], [], []
    stack = [(collapse(0), -1)]
    while stack:
        v, p = stack.pop()
        i = len(nodes)
        nodes.append(v)
        parent.append(p)
        children.append([])
        if p >= 0:
            children[p].append(i)
        if v not in id2sp:
            for c in inc(v):
                stack.append((collapse(c), i))
    n = len(nodes)
    pt = PTree()
    pt.children, pt.parent, pt.root = children, parent, 0
    pt.order = list(range(n - 1, -1, -1))                 # every child has a larger index than its parent
    pt.species = [id2sp.get(v) for v in nodes]
    pt.node_of_species = {sp: i for i, sp in enumerate(pt.species) if sp}
    under = [None] * n
    for i in pt.order:
        under[i] = [pt.species[i]] if pt.species[i] else [s for c in children[i] for s in under[c]]
    pt.under = under
    tel = [pt.node_of_species[s] for s in species if s in teleost]
    pt.T = None
    if len(tel) >= 2:
        t = tel[0]
        while not all(_is_ancestor(t, x, parent) for x in tel):
            t = parent[t]
        pt.T = t                                  # the lowest node above every teleost species
    return pt


def _is_ancestor(a, v, parent):
    while v != -1:
        if v == a:
            return True
        v = parent[v]
    return False


def _tables(pt, obs, K, trans, forced=None, want_arg=False, tT=None):
    R = range(K + 1)
    n = len(pt.children)
    f, g = [None] * n, [None] * n
    arg = [None] * n if want_arg else None
    for v in pt.order:
        ch = pt.children[v]
        if not ch:
            o = obs[v]
            f[v] = [0] * (K + 1) if o is None else [0 if s == o else INF for s in R]
        else:
            f[v] = [min(INF, sum(g[c][p] for c in ch)) for p in R]
        if forced is not None and v == pt.T:
            f[v] = [f[v][x] if x == forced else INF for x in R]
        if v != pt.root:
            fv = f[v]
            tv = tT if (tT is not None and v == pt.T) else trans
            best = [min(range(K + 1), key=lambda c, p=p: (tv[p][c] + fv[c], c)) for p in R]
            g[v] = [min(INF, tv[p][best[p]] + fv[best[p]]) for p in R]
            if want_arg:
                arg[v] = best
    return f, g, arg


def _cost_forced(pt, f, g, s, K, trans, root_state=1, tT=None):
    R = range(K + 1)
    T = pt.T
    cur_f = [f[T][x] if x == s else INF for x in R]
    if T == pt.root:
        return cur_f[root_state]
    t0 = tT if tT is not None else trans
    cur_g = [min(INF, min(t0[p][c] + cur_f[c] for c in R)) for p in R]
    node = T
    while True:
        parent = pt.parent[node]
        fp = [min(INF, sum((cur_g[p] if c == node else g[c][p]) for c in pt.children[parent])) for p in R]
        if parent == pt.root:
            return fp[root_state]
        cur_g = [min(INF, min(trans[p][c] + fp[c] for c in R)) for p in R]
        node = parent


def classify(costs: dict):
    srt = sorted(costs.values())
    if srt[0] >= INF:
        return "no_solution", 0
    winners = [s for s, c in costs.items() if c == srt[0]]
    return (LABEL[winners[0]] if len(winners) == 1 else "ambiguous"), srt[1] - srt[0]


class ParsimonyContext:
    def __init__(self, st, clades, kcap: int = 6, grid=()):
        self.st, self.clades, self.K = st, clades, kcap
        self.children = [[] for _ in st.parent]
        for i, p in enumerate(st.parent):
            if p >= 0:
                self.children[p].append(i)
        self._level, self._pt, self._cache = {}, {}, {}
        self.grid = tuple(grid)
        sets = set(COST_SETS) | {(g, l) for g, l, _ in self.grid}
        self.trans = {cs: trans_matrix(kcap, *cs) for cs in sets}
        self.trans_free = {cs: trans_matrix_wgd(kcap, *cs) for cs in sets}

    def level(self, sp):
        if sp not in self._level:
            self._level[sp] = level_of_species(self.clades, sp)
        return self._level[sp]

    def in_scope(self, sp) -> bool:
        return sp in self.st.leaf_id and self.level(sp) in JAWED

    def ptree(self, scope: frozenset):
        if scope not in self._pt:
            teleost = {s for s in scope if self.level(s) in TELEOST_LEVELS}
            self._pt[scope] = _build_ptree(self.st, self.children, sorted(scope), teleost)
        return self._pt[scope]

    def score(self, counts: dict, scope: frozenset) -> dict:
        pt = self.ptree(scope)
        if pt is None or pt.T is None:
            return {"p_status": "insufficient_teleosts"}
        if not any(self.level(s) == "gnathostome" for s in scope) or pt.T == pt.root:
            return {"p_status": "no_outgroup"}
        K = self.K
        clip = lambda c: min(c, K)
        n = len(pt.children)
        obs_loss = [None if pt.species[i] is None else clip(counts.get(pt.species[i], 0)) for i in range(n)]
        obs_miss = [None if (o is None or o == 0) else o for o in obs_loss]
        key = (id(pt), tuple(obs_loss))
        if key not in self._cache:
            self._cache[key] = self._score(pt, obs_loss, obs_miss)
        res = dict(self._cache[key])
        out_sp = [s for s in scope if self.level(s) == "gnathostome"]
        res.update(p_status="ok", p_scope_species=len(scope), p_outgroup_species=len(out_sp),
                   p_teleost_species=len(scope) - len(out_sp),
                   p_outgroup_with_gene=sum(1 for s in out_sp if counts.get(s, 0) > 0),
                   p_outgroup_single=sum(1 for s in out_sp if counts.get(s, 0) == 1))
        return res

    def _costs(self, pt, obs, cs, free=False):
        tT = self.trans_free[cs] if free else None
        f, g, _ = _tables(pt, obs, self.K, self.trans[cs], tT=tT)
        return {s: _cost_forced(pt, f, g, s, self.K, self.trans[cs], tT=tT) for s in (1, 2, 3)}

    def _score(self, pt, obs_loss, obs_miss):
        per_cs = {cs: self._costs(pt, obs_loss, cs) for cs in COST_SETS}
        labels = {cs: classify(c)[0] for cs, c in per_cs.items()}
        unit = per_cs[(1, 1)]
        label_unit, margin = classify(unit)
        robust = len(set(labels.values())) == 1 and label_unit != "ambiguous"
        miss_costs = self._costs(pt, obs_miss, (1, 1))
        label_miss, margin_miss = classify(miss_costs)
        res = {"p_label": label_unit if robust else "ambiguous", "p_label_unit": label_unit, "p_robust": robust,
               "p_margin": margin, "p_cost1": unit[1], "p_cost2": unit[2], "p_cost3": unit[3],
               "p_label_missing": label_miss, "p_margin_missing": margin_miss,
               "p_zero_modes_agree": (label_miss == label_unit) if label_unit != "ambiguous" and label_miss != "ambiguous" else ""}
        for g, l, free in self.grid:
            lab, mar = classify(self._costs(pt, obs_loss, (g, l), free))
            res[f"pg_{grid_name(g, l, free)}"] = lab
            res[f"pgm_{grid_name(g, l, free)}"] = mar
        res.update(self._events(pt, obs_loss, label_unit, unit))
        return res

    def _events(self, pt, obs, label, unit):
        s = min(unit, key=lambda x: (unit[x], x))              # one optimal teleost-ancestor state
        K, trans = self.K, self.trans[(1, 1)]
        f, g, arg = _tables(pt, obs, K, trans, forced=s, want_arg=True)
        state = [None] * len(pt.children)
        state[pt.root] = 1
        for v in range(len(pt.children)):
            for c in pt.children[v]:
                state[c] = arg[c][state[v]]
        in_T = [False] * len(pt.children)
        for v in range(len(pt.children)):
            in_T[v] = v != pt.T and _is_ancestor(pt.T, v, pt.parent)
        gains, losses, gl, lsp, lcl = 0, 0, Counter(), [], Counter()
        g_below = l_below = 0
        for v in range(1, len(pt.children)):
            d = state[v] - state[pt.parent[v]]
            if d == 0:
                continue
            leaf = pt.species[v]
            if d > 0:
                gains += d
                g_below += d if in_T[v] else 0
                lvl = "teleost_ancestor" if v == pt.T else (self.level(leaf) if leaf else level_of_set(self.clades, pt.under[v]))
                gl[f"{lvl}{'(species)' if leaf else ''}"] += d
            else:
                losses += -d
                l_below += -d if in_T[v] else 0
                if leaf:
                    lsp.append(leaf)
                else:
                    lcl[level_of_set(self.clades, pt.under[v]) or "?"] += -d
        return {"p_gains_n": gains, "p_losses_n": losses, "p_gains_below_T": g_below, "p_losses_below_T": l_below,
                "p_gain_levels": ";".join(f"{k}:{v}" for k, v in sorted(gl.items())),
                "p_loss_species": ";".join(sorted(lsp)),
                "p_loss_clades": ";".join(f"{k}:{v}" for k, v in sorted(lcl.items()))}


def parsimony_profile(rt, g, anchor, ctx: "ParsimonyContext") -> dict:
    """Score the copy numbers in the gene's anchor clade (the same clade as the B1 profile)."""
    name, a = anchor
    if a is None:
        return {"p_status": "no_anchor"}
    scope = frozenset(s for s in rt.species if ctx.in_scope(s))
    counts: dict = {}
    for x in rt.leaves_in(a):
        if x.sp in scope:
            counts[x.sp] = counts.get(x.sp, 0) + 1
    res = ctx.score(counts, scope)
    res["p_anchor"] = name
    return res
