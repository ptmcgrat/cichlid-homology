"""Read the families of any method into one common in-memory form.

Pure functions work on a family folder you give them (testable offline). Only `sample_families`
touches the coordinator's layout object, which the caller passes in.
"""
from __future__ import annotations

import csv
import json
import random
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

csv.field_size_limit(sys.maxsize)

_LEAF = re.compile(r"[(,]\s*([^(),:;\s]+)(?=\s*[:),])")


def leaf_names(nwk: str) -> list[str]:
    """Leaf names of a Newick string, in order. Internal-node labels (after ')') are ignored.
    Assumes no quoted names containing commas or parentheses."""
    return _LEAF.findall(nwk)


def root_children(nwk: str) -> int | None:
    """Number of children at the root, or None if the string does not start with '('."""
    s = nwk.strip().rstrip(";")
    if not s.startswith("("):
        return None
    depth, n = 0, 1
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                break
        elif ch == "," and depth == 1:
            n += 1
    return n


def fasta_names(path: Path) -> list[str]:
    with open(path) as fh:
        return [ln[1:].split()[0] for ln in fh if ln.startswith(">")]


def is_focal(leaf: str, prefixes: set[str]) -> bool:
    """A leaf is one of our genes if the text before its first '_' is a catalog prefix
    (the family contract's rule)."""
    head, sep, _ = leaf.partition("_")
    return bool(sep) and head in prefixes


@dataclass
class Family:
    family_id: str
    folder: Path
    tree: str
    leaves: list[str]
    focal: list[str]
    reference: list[str]
    aln_rows: list[str]
    done: dict
    tree_mode: str | None
    rooted: bool
    rooting: str | None
    root_children: int | None
    problems: list[str] = field(default_factory=list)


def load_family(folder: Path, prefixes: set[str], family_id: str | None = None) -> Family:
    """Read one family folder (tree.nwk, combined.aln, done.json). Never writes.
    Rooted only if done.json says so with a rooting text and the root of tree.nwk has two
    children (the contract); anything else is unrooted."""
    folder = Path(folder)
    fid = family_id or folder.name
    done = json.loads((folder / "done.json").read_text())
    tree = (folder / "tree.nwk").read_text()
    lv = leaf_names(tree)
    focal = [x for x in lv if is_focal(x, prefixes)]
    focalset = set(focal)
    ref = [x for x in lv if x not in focalset]
    aln = fasta_names(folder / "combined.aln")
    rc = root_children(tree)

    problems: list[str] = []
    if len(set(lv)) != len(lv):
        problems.append("duplicate_leaf")
    if set(aln) - set(lv):
        problems.append("aln_row_without_leaf")
    if set(lv) - set(aln):
        problems.append("leaf_without_aln_row")

    rooted, rooting = False, None
    claimed = done.get("rooted")
    if claimed is True and done.get("rooting") and rc == 2:
        rooted, rooting = True, str(done["rooting"])
    elif claimed is True:
        problems.append("rooted_claim_not_usable")
    elif claimed not in (False, None):
        problems.append("rooted_malformed")

    return Family(fid, folder, tree, lv, focal, ref, aln, done, done.get("tree_mode"),
                  rooted, rooting, rc, problems)


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sample_families(layout, part: str, species: str, n: int, seed: int = 1, release: str = "v1"):
    """Choose up to n families that hold a gene of `species` and have a finished folder
    (done.json present). Returns (families root, sorted ids chosen, assignment rows)."""
    res = layout.results(release, part)
    rows = read_tsv(res.assign(species).local)
    fams = sorted({r["best_family"] for r in rows if r["best_family"] not in ("", "NoHit")})
    root = res.family("X").local.parent
    done = [f for f in fams if (root / f / "done.json").exists()]
    rng = random.Random(seed)
    chosen = sorted(rng.sample(done, min(n, len(done))))
    return root, chosen, rows
