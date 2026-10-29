"""Reading the raw Compara dump (0.21.0): NHX trees with their support tags, and the EMF file that holds them.

The prepared trees of the Compara part were pruned and lost the tags Ensembl writes on every node: `D=Y` or `D=N` (a duplication), `DD` (a dubious duplication)
and `B=` (a bootstrap value, 0 to 100). The dump `Compara.<release>.protein_default.nhx.emf.gz` still has them. A block of the EMF file is a run of `SEQ` lines (species,
protein, chromosome, start, end, strand, gene, symbol), a line `DATA`, the tree on one line, and a line `//`.

`block_supports` keys the support by the set of reference proteins below a node, so that it can be found again on a pruned or rooted copy of the same tree. When several raw nodes
hold the same reference proteins (the others carry only species we pruned away), the lowest one is used: the smallest clade that holds them."""
from __future__ import annotations

import gzip
import re
from pathlib import Path

_TOK = re.compile(r"\(|\)|,|;|\[[^\]]*\]|:[0-9eE+.\-]+|[^,();:\[\]\s]+")


class NHXNode:
    __slots__ = ("name", "length", "tags", "children", "parent")

    def __init__(self, name=None):
        self.name, self.length, self.tags, self.children, self.parent = name, None, {}, [], None


def _tags(text: str) -> dict:
    out = {}
    for kv in text.strip("[]").replace("&&NHX:", "").replace("&&NHX", "").split(":"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def parse_nhx(text: str) -> NHXNode:
    """A Newick tree with `[&&NHX:k=v:...]` tags. Iterative, so a long caterpillar cannot hit the recursion limit."""
    stack, last, root = [], None, None
    for tok in _TOK.findall(text.strip()):
        if tok == "(":
            n = NHXNode()
            if stack:
                n.parent = stack[-1]
                stack[-1].children.append(n)
            stack.append(n)
            last = None
        elif tok == ",":
            last = None
        elif tok == ")":
            last = stack.pop()
            root = last
        elif tok == ";":
            break
        elif tok[0] == ":":
            if last is not None:
                last.length = float(tok[1:])
        elif tok[0] == "[":
            if last is not None:
                last.tags.update(_tags(tok))
        elif last is not None and last.children and last.name is None:
            last.name = tok                                  # the label of an internal node
        else:
            n = NHXNode(tok)
            if stack:
                n.parent = stack[-1]
                stack[-1].children.append(n)
            else:
                root = n
            last = n
    if root is None:
        raise ValueError("no tree in the text")
    return root


def postorder(root: NHXNode) -> list:
    out, stack = [], [(root, False)]
    while stack:
        n, done = stack.pop()
        if done:
            out.append(n)
        else:
            stack.append((n, True))
            stack.extend((c, False) for c in reversed(n.children))
    return out


def iter_emf(path):
    """Yield (seqs, tree_text) for every block of an EMF file (plain or .gz). `seqs` maps protein to {species, chrom, start, end, strand, gene, symbol}."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, tree, want = {}, None, False
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("SEQ "):
                f = line.split()
                if len(f) >= 8:
                    seqs[f[2]] = {"species": f[1], "chrom": f[3], "start": int(f[4]) if f[4].isdigit() else 0, "end": int(f[5]) if f[5].isdigit() else 0,
                                  "strand": f[6], "gene": f[7], "symbol": f[8] if len(f) > 8 else ""}
            elif line == "DATA":
                want = True
            elif line.startswith("//"):
                if seqs or tree:
                    yield seqs, tree
                seqs, tree, want = {}, None, False
            elif want and line:
                tree, want = line, False
        if seqs or tree:
            yield seqs, tree


def block_supports(tree_text: str, universe: set) -> dict:
    """{frozenset of reference proteins below a node: (bootstrap or None, duplication flag 'Y'/'N' or None, dubious-duplication flag or None)}.
    Only nodes with at least two reference proteins are kept."""
    root = parse_nhx(tree_text)
    sets, out = {}, {}
    for n in postorder(root):
        if not n.children:
            s = frozenset([n.name]) if n.name in universe else frozenset()
        else:
            s = frozenset().union(*(sets[id(c)] for c in n.children))
            if len(s) >= 2 and s not in out:                 # postorder: the lowest node with this set is met first
                b = n.tags.get("B")
                out[s] = (int(b) if b not in (None, "") else None, n.tags.get("D"), n.tags.get("DD"))
        sets[id(n)] = s
    return out


def collect_support(path, wanted: dict, progress=None, stop_when_done=True) -> dict:
    """Stream the EMF file once. `wanted` maps a family id to the set of its reference proteins. Returns {family: {"support": block_supports(...), "symbols": {protein: symbol}}}
    for the families found; a family is matched by any one of its proteins."""
    prot = {p: f for f, ps in wanted.items() for p in ps}
    found, n = {}, 0
    for seqs, tree in iter_emf(path):
        n += 1
        if progress and n % 2000 == 0:
            progress(n, len(found))
        if tree is None:
            continue
        fid = next((prot[p] for p in seqs if p in prot), None)
        if fid is None or fid in found:
            continue
        uni = set(wanted[fid])
        found[fid] = {"support": block_supports(tree, uni), "symbols": {p: s["symbol"] for p, s in seqs.items() if p in uni and s["symbol"]}}
        if stop_when_done and len(found) == len(wanted):
            break
    return found
