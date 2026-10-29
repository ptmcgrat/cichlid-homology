"""Lineages of a clade from the teleost genes alone (0.21.0).

Why: in the reference tree of IRX1 the duplication that gave irx1a and irx1b is clear (35 and 34 teleost species, 33 of 37 in both), but the gar and reedfish genes hang inside one
of the two lineages, with a bootstrap of 0 and 25. A rule that wants both sides to be teleost-only threw the node away. Here the lineages are found on the teleost genes only: the
highest node whose two biggest children share at least half of their teleost species (and both hold an early teleost) is the split. Every gene that is not a teleost gene is the
ancestral copy by definition, wherever the tree attaches it. Where it attaches, and how well supported that is, is reported and never used to label anything.

A leaf takes part when `participates(leaf)` is true (for a placed tree: the reference leaves, not our own genes; for a de novo tree: every leaf) and its species is a teleost.
`support` (see `suppmap.attach_support`) is a dict of id(node) to (bootstrap, duplication flag, dubious flag), or None when there are no support values."""
from __future__ import annotations

import collections
from dataclasses import dataclass, field

from .tree import Node, leaves, postorder

TELEOST = ("teleost", "percomorph", "ovalentaria", "cichlid")
LETTERS = "ABCDEFGH"
MIN_SHARED = 2          # two genes of one species are a species-level duplication, not a lineage split: the lineages must share at least two species


@dataclass
class Split:
    node: Node
    kind: str                 # "teleost" (early teleosts on both sides) or "later" (a duplication among the teleosts that has no early teleost on both sides)
    level: str                # the ladder level of the species below the node
    overlap: float            # shared / all species of the two biggest lineages
    lineages: list            # children of the node that are lineages, biggest first
    stray: list               # the other children


@dataclass
class Analysis:
    split: Split | None
    S: dict                                       # id(node) -> frozenset of participating teleost species below it
    member: dict                                  # id(leaf) -> lineage index, "stray" or "out"
    letters: list                                 # one letter per lineage
    counts: dict                                  # species (our own genes: "focal:<species>") -> [outside, lineage 0, lineage 1, ..., stray]
    lineages: list = field(default_factory=list)  # one dict per lineage: id, species, genes, stem, depth, crown_support, anchors
    outgroup_inside: list = field(default_factory=list)
    outgroup_outside: list = field(default_factory=list)
    mode: str = "none"                            # "split", "single" (teleost genes, no split) or "none" (no teleost gene)


def _jacc(a, b) -> float:
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def _sets(g: Node, participates, is_tel, early):
    S, E = {}, {}
    for n in postorder(g):
        if not n.children:
            ok = n.sp is not None and participates(n) and is_tel(n.sp)
            S[id(n)] = frozenset([n.sp]) if ok else frozenset()
            E[id(n)] = frozenset([n.sp]) if ok and early(n.sp) else frozenset()
        else:
            S[id(n)] = frozenset().union(*(S[id(c)] for c in n.children))
            E[id(n)] = frozenset().union(*(E[id(c)] for c in n.children))
    return S, E


def find_split(g: Node, S: dict, E: dict, clades, overlap_t: float = 0.5, order_key=None):
    """The highest node that splits the teleost genes into lineages. A teleost-level split (early teleosts in both) wins over a later one, wherever they are.
    The lineages are returned in the order of `order_key(child)` (default: the bigger first, then by the smallest gene name), so that their letters do not depend on how the tree was drawn."""
    default_key = lambda c: (-len(S[id(c)]), min(x.name for x in leaves(c)))
    queue, later = collections.deque([g]), None
    while queue:
        v = queue.popleft()
        kids = [c for c in v.children if S[id(c)]]
        if len(kids) >= 2:
            big = max(kids, key=lambda c: (len(S[id(c)]), [-ord(ch) for ch in min(x.name for x in leaves(c))]))
            lin = [big] + [c for c in kids if c is not big and _jacc(S[id(big)], S[id(c)]) >= overlap_t and len(S[id(big)] & S[id(c)]) >= MIN_SHARED]
            if len(lin) >= 2:
                by_size = sorted(lin, key=default_key)
                stray = [c for c in v.children if all(c is not x for x in lin)]
                ov = _jacc(S[id(by_size[0])], S[id(by_size[1])])
                sp = frozenset().union(*(S[id(c)] for c in lin))
                ordered = sorted(lin, key=order_key) if order_key else by_size
                if E[id(by_size[0])] and E[id(by_size[1])]:
                    return Split(v, "teleost", clades.common_level(sorted(sp)), ov, ordered, stray)
                if later is None and min(len(S[id(by_size[0])]), len(S[id(by_size[1])])) >= 2:
                    later = Split(v, "later", clades.common_level(sorted(sp)), ov, ordered, stray)
        queue.extend(v.children)
    return later


def crown_of(c: Node, S: dict):
    """Descend from a lineage's root while only one child holds teleost genes: the crown is where they first split. Returns (crown, stem length)."""
    n, stem = c, c.length or 0.0
    while n.children:
        kids = [k for k in n.children if S[id(k)]]
        if len(kids) != 1:
            break
        n = kids[0]
        stem += n.length or 0.0
    return n, stem


def _depth(crown: Node, part_ids: set) -> float | None:
    out, stack = [], [(crown, 0.0)]
    while stack:
        n, d = stack.pop()
        if not n.children:
            if id(n) in part_ids:
                out.append(d)
        else:
            stack.extend((c, d + (c.length or 0.0)) for c in n.children)
    return sum(out) / len(out) if out else None


def _lowest_with_teleosts(leaf: Node, S: dict):
    n = leaf.parent
    while n is not None and not S[id(n)]:
        n = n.parent
    return n


def analyse(g: Node, is_focal, clades, support=None, overlap_t: float = 0.5, participates=None, order_key=None) -> Analysis:
    """Lineages, membership of every leaf, copies per species, and where the non-teleost genes attach. `is_focal(name)` tells our own genes."""
    part = participates or (lambda x: not is_focal(x.name))
    is_tel = lambda s: clades.level(s) in TELEOST
    early = lambda s: clades.level(s) == "teleost"
    S, E = _sets(g, part, is_tel, early)
    split = find_split(g, S, E, clades, overlap_t, order_key)
    lv = leaves(g)
    member, letters = {}, []
    if split:
        letters = list(LETTERS[:len(split.lineages)])
        for i, c in enumerate(split.lineages):
            for x in leaves(c):
                member[id(x)] = i
        for c in split.stray:
            for x in leaves(c):
                member[id(x)] = "stray"
        for x in lv:
            member.setdefault(id(x), "out")
        mode = "split"
    else:
        tel = [x for x in lv if x.sp is not None and is_tel(x.sp)]
        letters = ["A"] if tel else []
        for x in lv:
            member[id(x)] = 0 if (x.sp is not None and is_tel(x.sp)) else "out"
        mode = "single" if tel else "none"
    k = len(letters)
    counts = collections.defaultdict(lambda: [0] * (k + 2))
    for x in lv:
        if x.sp is None:
            continue
        m = member[id(x)]
        outgroup = not is_tel(x.sp) and not is_focal(x.name)       # gar, shark, tetrapods...: the ancestral copy, wherever the tree attaches them
        key = ("focal:" + x.sp) if is_focal(x.name) else x.sp           # a reference species can also be one of ours (stickleback, betta): keep them apart
        counts[key][0 if (m == "out" or outgroup) else (k + 1 if m == "stray" else 1 + m)] += 1
    part_ids = {id(x) for x in lv if S.get(id(x))}
    lins = []
    if split:
        for i, c in enumerate(split.lineages):
            crown, stem = crown_of(c, S)
            lins.append({"id": letters[i], "species": len(S[id(c)]), "genes": sum(1 for x in leaves(c) if id(x) in part_ids), "stem": round(stem, 5),
                         "depth": (round(_depth(crown, part_ids), 5) if _depth(crown, part_ids) is not None else None),
                         "crown_support": (support.get(id(crown), (None,))[0] if support else None), "crown_node": crown})
    elif mode == "single":
        lins.append({"id": "A", "species": len(S[id(g)]), "genes": len(part_ids & {id(x) for x in lv}), "stem": None, "depth": None, "crown_support": None, "crown_node": g})
    inside, outside = [], []
    for x in lv:
        if x.sp is None or (x.sp is not None and (is_tel(x.sp) or is_focal(x.name))):
            continue                                    # the outgroup genes are the ones that are neither teleost nor ours
        m = member[id(x)]
        if split and isinstance(m, int):
            u = _lowest_with_teleosts(x, S)
            inside.append({"gene": x.name, "species": x.sp, "lineage": letters[m],
                           "attach_support": (support.get(id(u), (None,))[0] if support and u is not None else None)})
        else:
            outside.append({"gene": x.name, "species": x.sp})
    return Analysis(split, S, member, letters, dict(counts), lins, inside, outside, mode)
