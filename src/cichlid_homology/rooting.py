"""Rooting an unrooted gene tree by reconciliation with a species tree (never midpoint).

Cost of a rooting = (duplications, losses) under LCA reconciliation, compared in that order: the
fewest duplications decide, and ties are broken by the fewest losses (the 0.2.x pilot found that most
genes had several rootings with equal duplication counts). A node is a duplication when its
species-tree image equals the image of one of its children. All rootings are scored in one pass over
directed edges, so none needs the tree rebuilt.
"""
from __future__ import annotations

from dataclasses import dataclass

from .speciestree import SpeciesTree
from .tree import Node, leaves

W = 10_000_000   # a duplication outweighs any number of losses (losses stay far below W)


@dataclass
class Rooting:
    roots: list          # rooted trees (Node), best first, at most k
    cost: int | None     # duplications at the best root
    losses: int | None   # losses at the best root
    n_tied: int          # rootings tied on (duplications, losses)
    n_tied_dup: int      # rootings tied on duplications alone
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


def _node_cost(st: SpeciesTree, m: int, ms: list):
    """(is_duplication, losses) of a node with image m whose children have images ms."""
    dup = len(ms) >= 2 and m in ms
    dep = st.depth
    loss = sum(dep[x] - dep[m] for x in ms) if dup else sum(dep[x] - dep[m] - 1 for x in ms)
    return dup, max(loss, 0)


def _edge_costs(nodes, adj, st: SpeciesTree, sp_id: dict):
    """{(u, v): W * duplications + losses if the root is placed on edge (u, v)}."""
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
            dup, loss = _node_cost(st, m, ms)
            M[key] = m
            D[key] = sum(D[(v_, c)] for c in kids) + (W if dup else 0) + loss
            stack.pop()

    cost = {}
    for u in range(len(nodes)):
        for v in adj[u]:
            if u < v:
                solve(v, u)     # subtree at u, away from v
                solve(u, v)
                mu, mv = M[(v, u)], M[(u, v)]
                m = st.lca(mu, mv)
                dup, loss = _node_cost(st, m, [mu, mv])
                cost[(u, v)] = D[(v, u)] + D[(u, v)] + (W if dup else 0) + loss + st.depth[m]
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
        return Rooting([tree], None, None, 1, 1, "recorded")
    if len(lv) < 2:
        return Rooting([tree], 0, 0, 1, 1, "reconciliation")
    nodes, adj, elen = _unrooted(tree)
    sp_id = {n.idx: st.id_of(n.sp) for n in nodes if n.is_leaf}
    cost = _edge_costs(nodes, adj, st, sp_id)
    best = min(cost.values())
    tied = sorted(e for e, c in cost.items() if c == best)
    best_dup = min(c // W for c in cost.values())
    n_tied_dup = sum(1 for c in cost.values() if c // W == best_dup)
    roots = [_reroot(nodes, adj, elen, u, v) for (u, v) in tied[:max(1, k)]]
    return Rooting(roots, best // W, best % W, len(tied), n_tied_dup, "reconciliation")
