"""Species profiles of the clades that hold a teleost-genome-duplication node (0.20.0): which species carry which copy.

For each clade of `clade-split` that has a supported duplication of teleost level (early teleosts in both lineages, or an early teleost with genes on
both sides, or all the early-teleost genes inside the node: the same rule as the dating of the gene trees), the highest such node is the pair node. Its two
largest lineages are A and B (a third or later lineage is counted with B). For every species the profile counts the genes outside the node (the ancestral
copy: tetrapods, gar, reedfish, shark and the teleosts outside), in lineage A and in lineage B. A species with no gene in A or in B lost that copy, or the
annotation lacks it; a count of two or more in one lineage is a later duplicate of that ohnolog. Later duplications that cover two or more species are listed
with their lineage, level and species, so a bracket can be drawn over a clade (all salmonids, say). The zebrafish genes of each lineage are listed so that
a focal-species gene, whose closest zebrafish gene sits in one of them, can be assigned to a lineage.

The files written: a species table in species-tree order (with the later whole-genome duplications of the salmonids and carps marked), the species tree pruned
to those species as nested JSON, and one row per clade with the counts as a string `outside.A.B` per species, in the order of the species table."""
from __future__ import annotations

import collections
import json
from pathlib import Path

from .tree import parse_newick

WGD_MARK = {"salmo_salar": "salmonid-specific duplication", "salmo_trutta": "salmonid-specific duplication", "oncorhynchus_mykiss": "salmonid-specific duplication",
            "oncorhynchus_kisutch": "salmonid-specific duplication", "oncorhynchus_tshawytscha": "salmonid-specific duplication",
            "cyprinus_carpio": "carp-specific duplication", "carassius_auratus": "carp-specific duplication", "sinocyclocheilus_grahami": "cyprinid-specific duplication"}
MAX_LATER = 12


def group_label(sp, clades, coelacanth="latimeria_chalumnae") -> str:
    lv = clades.level(sp)
    try:
        lin = clades.lineage(sp)
    except Exception:
        lin = None
    if lv in ("percomorph", "ovalentaria", "cichlid"):
        return {"percomorph": "percomorph", "ovalentaria": "ovalentarian", "cichlid": "cichlid"}[lv]
    if lv == "teleost":
        return "early teleost"
    if lv == "gnathostome":
        if lin == "sarcopterygian":
            return "coelacanth" if sp == coelacanth else "tetrapod"
        if lin == "chondrichthyan":
            return "cartilaginous fish"
        return "gar or reedfish"
    return "other"


def profile_of(fid, gname, g, rt, T, supported, names):
    """The profile dict of clade `g` (a node of rooted tree `rt`) or None when it has no teleost-level duplication node."""
    queue, node = collections.deque([g]), None
    while queue:                                                       # breadth first: the highest such node
        v = queue.popleft()
        if v.children and supported(v) and (T.evaluate(v, g) or {}).get("level") == "teleost":
            node = v
            break
        queue.extend(v.children)
    if node is None:
        return None
    lin = sorted(T.evaluate(node, g)["lineages"], key=lambda c: -T.n_species(c))
    if len(lin) < 2:
        return None
    groups = {"A": [lin[0]], "B": lin[1:]}
    inside = {}
    counts = collections.defaultdict(lambda: [0, 0, 0])
    zf = {"A": [], "B": []}
    for k, (letter, nodes) in enumerate(groups.items(), start=1):
        for c in nodes:
            for x in rt.leaves_in(c):
                if T.M[id(x)] is None:
                    continue
                inside[id(x)] = letter
                counts[x.sp][k] += 1
                if x.sp == names["zebrafish"]:
                    zf[letter].append(x.name)
    for x in rt.leaves_in(g):
        if T.M[id(x)] is not None and id(x) not in inside:
            counts[x.sp][0] += 1
    later, seen = [], set()
    for letter, nodes in groups.items():
        for c in nodes:
            stack = [c]
            while stack:
                u = stack.pop()
                if u.children:
                    stack.extend(u.children)
                    if T.n_species(u) >= 2 and supported(u):
                        sp = tuple(sorted({x.sp for x in rt.leaves_in(u) if T.M[id(x)] is not None}))
                        if (letter, sp) not in seen:
                            seen.add((letter, sp))
                            later.append({"lineage": letter, "level": T.level(u), "species": list(sp)})
    later.sort(key=lambda d: -len(d["species"]))
    return {"clade": f"{fid}:{gname}", "family": fid, "group": gname, "n_lineages": len(lin), "counts": {s: v for s, v in counts.items()},
            "later": later[:MAX_LATER], "later_omitted": max(0, len(later) - MAX_LATER), "zebrafish_a": zf["A"], "zebrafish_b": zf["B"]}


def _pruned(node, keep):
    if not node.children:
        return {"name": node.name} if node.name in keep else None
    kids = [k for k in (_pruned(c, keep) for c in node.children) if k]
    if not kids:
        return None
    if len(kids) == 1:
        return kids[0]
    label = node.name if node.name and not node.name.replace(".", "").replace("-", "").isdigit() else ""
    return {"name": label, "children": kids}


def species_tree_json(newick, keep) -> dict:
    """The species tree pruned to the species in `keep` as nested {name, children}; unary nodes are collapsed and numeric labels dropped."""
    return _pruned(parse_newick(newick), set(keep))


def _leaf_order(t, out):
    if "children" in t:
        for c in t["children"]:
            _leaf_order(c, out)
    else:
        out.append(t["name"])
    return out


def build_files(profiles, clades, newick, coelacanth="latimeria_chalumnae"):
    """(species rows, tree, clade rows) from the profile dicts, in species-tree order."""
    used = sorted({s for p in profiles for s in p["counts"]})
    tree = species_tree_json(newick, used)
    order = _leaf_order(tree, []) if tree else used
    order = order + [s for s in used if s not in order]
    species = [{"order": i + 1, "name": s, "label": s.replace("_", " ").capitalize(), "group": group_label(s, clades, coelacanth), "level": clades.level(s),
                "wgd": WGD_MARK.get(s, "")} for i, s in enumerate(order)]
    rows = []
    for p in profiles:
        z = ",".join(".".join(map(str, p["counts"].get(s, [0, 0, 0]))) for s in order)
        rows.append({"clade": p["clade"], "family": p["family"], "group": p["group"], "counts": z, "n_lineages": p["n_lineages"],
                     "later": json.dumps(p["later"], separators=(",", ":")), "later_omitted": p["later_omitted"],
                     "zebrafish_a": ";".join(p["zebrafish_a"]), "zebrafish_b": ";".join(p["zebrafish_b"])})
    return species, tree, rows


def write_files(out, part, profiles, clades, newick) -> dict:
    species, tree, rows = build_files(profiles, clades, newick)
    out = Path(out)
    for name, obj in (("species", species), ("clades", rows), ("tree", tree)):
        (out / f"clade_species.{part}.{name}.json").write_text(json.dumps(obj, separators=(",", ":"), ensure_ascii=False))
    return {"species": len(species), "clades": len(rows)}
