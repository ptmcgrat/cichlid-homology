"""A species tree with LCA queries (the master tree, or a test tree)."""
from __future__ import annotations

from .tree import Node, parse_newick


class SpeciesTree:
    def __init__(self, newick: str):
        root = parse_newick(newick)
        self.parent: list[int] = []
        self.depth: list[int] = []
        self.leaf_id: dict[str, int] = {}
        stack = [(root, -1, 0)]
        while stack:
            n, p, d = stack.pop()
            i = len(self.parent)
            self.parent.append(p)
            self.depth.append(d)
            if n.children:
                for c in reversed(n.children):
                    stack.append((c, i, d + 1))
            else:
                if n.name in self.leaf_id:
                    raise ValueError(f"species {n.name} twice in the species tree")
                self.leaf_id[n.name] = i

    def id_of(self, species: str) -> int:
        try:
            return self.leaf_id[species]
        except KeyError:
            raise KeyError(f"species {species!r} is not a leaf of the species tree") from None

    def lca(self, a: int, b: int) -> int:
        da, db, par = self.depth, self.depth, self.parent
        while da[a] > db[b]:
            a = par[a]
        while db[b] > da[a]:
            b = par[b]
        while a != b:
            a, b = par[a], par[b]
        return a
