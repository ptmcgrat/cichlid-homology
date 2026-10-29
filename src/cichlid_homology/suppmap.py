"""Putting Ensembl's support values on our trees, and rooting a placed tree where the reference tree is rooted (0.21.0).

`attach_support` maps the dict of `nhx.block_supports` onto a tree of `Node`s: a node gets the support of the raw clade that holds exactly the same reference proteins. Our own genes (the
ones that are not in the reference) are ignored, so a node that holds placed genes still finds the support of its reference clade.

`reroot_like` roots an unrooted placed tree on the edge that splits the reference leaves the way the reference tree's own root does. The clade cuts and the supports were made on that rooting,
so a figure drawn on another rooting would not line up with them. It returns None when no edge fits (then the caller roots by reconciliation, as the pilot does)."""
from __future__ import annotations

from .rooting import _reroot, _unrooted
from .tree import Node, leaves, postorder


def attach_support(root: Node, support_map: dict, universe: set) -> dict:
    """{id(node): (bootstrap, duplication flag, dubious flag)} for every internal node whose reference leaves have a raw support."""
    sets, out = {}, {}
    for n in postorder(root):
        if not n.children:
            s = frozenset([n.name]) if n.name in universe else frozenset()
        else:
            s = frozenset().union(*(sets[id(c)] for c in n.children))
            if len(s) >= 2 and s in support_map:
                out[id(n)] = support_map[s]
        sets[id(n)] = s
    return out


def _child_sets(ref: Node, universe: set):
    a = frozenset(x.name for x in leaves(ref.children[0]) if x.name in universe)
    return a, frozenset(universe) - a


def reroot_like(placed: Node, ref: Node):
    """Root `placed` (any rooting; its own root is ignored) on the edge that gives the reference leaves the same two groups as the root of `ref`."""
    if not ref.children or len(ref.children) < 2:
        return None
    pnames = {x.name for x in leaves(placed)}
    universe = {x.name for x in leaves(ref)} & pnames
    a, b = _child_sets(ref, universe)
    if not a or not b:
        return None
    nodes, adj, elen = _unrooted(placed)
    start = next((n.idx for n in nodes if not n.children and adj[n.idx]), None)
    if start is None:
        return None
    # depth-first from a leaf: the reference leaves on the far side of every directed edge
    below, order, stack, seen = {}, [], [(start, -1)], {start}
    while stack:
        x, p = stack.pop()
        order.append((x, p))
        for y in adj[x]:
            if y != p and y not in seen:
                seen.add(y)
                stack.append((y, x))
    for x, p in reversed(order):
        n = nodes[x]
        s = {n.name} & universe if not n.children else set()
        for y in adj[x]:
            if y != p and y in below:
                s |= below[y]
        below[x] = frozenset(s)
    for x, p in order:
        if p >= 0 and below[x] in (a, b):
            return _reroot(nodes, adj, elen, p, x)
    return None
