"""Rooting an unrooted gene tree by reconciliation with a species tree (never midpoint).

Cost of a rooting = number of duplication nodes under LCA reconciliation: a node is a duplication
when its species-tree image equals the image of one of its children. Computed for every edge in one
pass over directed edges, so all rootings are scored without rebuilding the tree.
"""
from __future__ import annotations

from dataclasses import dataclass

from .speciestree import SpeciesTree
from .tree import Node, leaves, postorder


@dataclass
class Rooting:
    roots: list          # rooted trees (Node), best first, at most k
    cost: int | None     # minimum duplication count
    n_tied: int          # how many edges reach the minimum
    how: str             # "recorded" or "reconciliation"


def _unrooted(root: Node):
    nodes, stack = [], [root]
    while stack:
        n = stack.pop()
        n.idx = len(nodes)
        nodes.append(n)
        stack.extend(reversed(n.children))
    adj = [[] for _ in nodes]
    elen = {}
    for n in nodes:
        for c in n.children:
            adj[n.idx].append(c.idx)
            adj[c.idx].append(n.idx)
            elen[(min(n.idx, c.idx), max(n.idx, c.idx))] = c.length
    if len(root.children) == 2:     # suppress the degree-2 root
        a, b = root.children
        adj[a.idx].remove(root.idx)
        adj[b.idx].remove(root.idx)
        adj[a.idx].append(b.idx)
        adj[b.idx].append(a.idx)
        adj[root.idx] = []
        la = elen.pop((min(root.idx, a.idx), max(root.idx, a.idx)), None)
        lb = elen.pop((min(root.idx, b.idx), max(root.idx, b.idx)), None)
        elen[(min(a.idx, b.idx), max(a.idx, b.idx))] = (
            None if la is None and lb is None else (la or 0.0) + (lb or 0.0))
    return nodes, adj, elen


def _edge_costs(nodes, adj, st: SpeciesTree, sp_id: dict):
    """Return {(u, v): total duplications if the root is placed on edge (u, v)}."""
    M, D = {}, {}

    def solve(p, v):
        stack = [(p, v)]
        while stack:
            key = stack[-1]
            if key in M:
                stack.pop()
                continue
            p_, v_ = key
            kids = [c for c in adj[v_] if c != p_]
            if not kids:
                M[key] = sp_id[v_]
                D[key] = 0
                stack.pop()
                continue
            missing = [(v_, c) for c in kids if (v_, c) not in M]
            if missing:
                stack.extend(missing)
                continue
            ms = [M[(v_, c)] for c in kids]
            m = ms[0]
            for x in ms[1:]:
                m = st.lca(m, x)
            d = sum(D[(v_, c)] for c in kids)
            if len(kids) >= 2 and m in ms:
                d += 1
            M[key], D[key] = m, d
            stack.pop()

    cost = {}
    for u in range(len(nodes)):
        for v in adj[u]:
            if u < v:
                solve(v, u)     # subtree at u, away from v
                solve(u, v)
                mu, mv = M[(v, u)], M[(u, v)]
                m = st.lca(mu, mv)
                cost[(u, v)] = D[(v, u)] + D[(u, v)] + (1 if m in (mu, mv) else 0)
    return cost


def _reroot(nodes, adj, elen, u, v) -> Node:
    root = Node()
    stack = [(root, v, u), (root, u, v)]
    while stack:
        newp, p, x = stack.pop()
        o = nodes[x]
        nn = Node(o.name, elen.get((min(p, x), max(p, x))))
        nn.label, nn.sp = o.label, o.sp
        newp.add(nn)
        for y in adj[x]:
            if y != p:
                stack.append((nn, x, y))
    for c in root.children:
        if c.length is not None:
            c.length = c.length / 2.0
    return root


def root_tree(tree: Node, st: SpeciesTree, already_rooted: bool, k: int = 1) -> Rooting:
    """`tree` leaves must carry `.sp` (a species-tree leaf name). Returns up to k rooted trees,
    best first; ties are ordered by edge index so the result is deterministic."""
    lv = leaves(tree)
    if already_rooted and len(tree.children) == 2:
        return Rooting([tree], None, 1, "recorded")
    if len(lv) < 2:
        return Rooting([tree], 0, 1, "reconciliation")
    nodes, adj, elen = _unrooted(tree)
    sp_id = {n.idx: st.id_of(n.sp) for n in nodes if n.is_leaf}
    cost = _edge_costs(nodes, adj, st, sp_id)
    best = min(cost.values())
    tied = sorted(e for e, c in cost.items() if c == best)
    roots = [_reroot(nodes, adj, elen, u, v) for (u, v) in tied[:max(1, k)]]
    return Rooting(roots, best, len(tied), "reconciliation")
