"""Small Newick tree toolkit (iterative, so a 1,500-leaf caterpillar cannot hit the recursion limit).

Our own parser keeps this repository independent of the coordinator's internal modules.
"""
from __future__ import annotations

import re


class Node:
    __slots__ = ("name", "length", "label", "children", "parent", "sp", "mask", "lo", "hi", "idx")

    def __init__(self, name=None, length=None):
        self.name = name
        self.length = length
        self.label = None
        self.children: list[Node] = []
        self.parent: Node | None = None
        self.sp = None      # species (leaves)
        self.mask = 0       # bit mask of the species below
        self.lo = self.hi = 0   # leaf-order interval below this node
        self.idx = -1

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def add(self, c: "Node") -> None:
        c.parent = self
        self.children.append(c)


_TOK = re.compile(r"\(|\)|,|;|:[^,();\s]*|[^,();:\s]+")


def parse_newick(text: str) -> Node:
    text = re.sub(r"\[[^\]]*\]", "", text).strip()
    stack: list[Node] = []
    last: Node | None = None
    for tok in _TOK.findall(text):
        if tok == "(":
            n = Node()
            if stack:
                stack[-1].add(n)
            stack.append(n)
            last = None
        elif tok == ",":
            last = None
        elif tok == ")":
            if not stack:
                raise ValueError("unbalanced parentheses")
            last = stack.pop()
        elif tok == ";":
            break
        elif tok[0] == ":":
            if last is None:
                raise ValueError("branch length without a node")
            last.length = float(tok[1:]) if len(tok) > 1 else None
        else:
            if last is None:
                leaf = Node(tok)
                if stack:
                    stack[-1].add(leaf)
                last = leaf
            else:
                last.label = tok
    if stack:
        raise ValueError("unbalanced parentheses")
    if last is None:
        raise ValueError("empty tree")
    return last


def postorder(root: Node) -> list[Node]:
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    out.reverse()
    return out


def leaves(root: Node) -> list[Node]:
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        if n.children:
            stack.extend(reversed(n.children))
        else:
            out.append(n)
    return out


def squeeze(c: Node) -> Node:
    while len(c.children) == 1:
        k = c.children[0]
        if c.length is not None or k.length is not None:
            k.length = (c.length or 0.0) + (k.length or 0.0)
        c = k
    return c


def prune(root: Node, drop: set) -> Node | None:
    """Remove the leaves whose names are in `drop`, and the unary nodes that results. In place.
    Returns the new root, or None if nothing is left."""
    dead = set()
    for n in postorder(root):
        if n.children:
            alive = [c for c in n.children if id(c) not in dead]
            if not alive:
                dead.add(id(n))
                n.children = []
            else:
                n.children = [squeeze(c) for c in alive]
                for c in n.children:
                    c.parent = n
        elif n.name in drop:
            dead.add(id(n))
    if id(root) in dead:
        return None
    root = squeeze(root)
    root.parent = None
    return root


def to_newick(root: Node) -> str:
    # iterative serialisation
    out, stack = [], [(root, 0)]
    while stack:
        n, state = stack.pop()
        if n.children:
            if state == 0:
                out.append("(")
                stack.append((n, 1))
                for i, c in enumerate(reversed(n.children)):
                    stack.append((c, 0))
                    if i < len(n.children) - 1:
                        stack.append((Node(name="\0,"), 0))
            else:
                out.append(")" + (n.label or "") + (f":{n.length}" if n.length is not None else ""))
        else:
            if n.name == "\0,":
                out.append(",")
            else:
                out.append((n.name or "") + (f":{n.length}" if n.length is not None else ""))
    return "".join(out) + ";"
