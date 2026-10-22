"""Adapter for the coordinator's master species tree.

The sheet names `cichlid_data.master_tree.load_master_tree()` but not the type it returns, so this
accepts a Newick string, a path, or an object that can give Newick, and fails with a clear message
(run `cichlid-homology probe`) otherwise.
"""
from __future__ import annotations

from pathlib import Path

from .speciestree import SpeciesTree


def node_newick(root) -> str:
    """Newick (names only) of any object with `children` and `name`, such as the coordinator's
    cichlid_data.newick.Node, which is what load_master_tree() returns. Iterative."""
    out, stack = [], [(root, 0)]
    while stack:
        n, state = stack.pop()
        if n is None:
            out.append(",")
            continue
        kids = list(n.children)
        if not kids:
            out.append(str(n.name))
        elif state == 0:
            out.append("(")
            stack.append((n, 1))
            for i, c in enumerate(reversed(kids)):
                stack.append((c, 0))
                if i < len(kids) - 1:
                    stack.append((None, 2))
        else:
            out.append(")")
    return "".join(out) + ";"


def newick_of(obj) -> str:
    if isinstance(obj, str):
        s = obj.strip()
        if s.startswith("("):
            return s
        if Path(s).exists():
            return Path(s).read_text()
    if isinstance(obj, Path) and obj.exists():
        return obj.read_text()
    if hasattr(obj, "children") and hasattr(obj, "name"):
        return node_newick(obj)
    for name in ("to_newick", "newick", "write_newick", "as_newick"):
        f = getattr(obj, name, None)
        if f is not None:
            r = f() if callable(f) else f
            if isinstance(r, str) and r.strip().startswith("("):
                return r.strip()
    s = str(obj).strip()
    if s.startswith("("):
        return s
    raise TypeError(f"cannot get Newick from {type(obj).__name__}; run `cichlid-homology probe` "
                    "and tell the coordinator what it prints")


def load_species_tree() -> SpeciesTree:
    from cichlid_data.master_tree import load_master_tree

    return SpeciesTree(newick_of(load_master_tree()))
