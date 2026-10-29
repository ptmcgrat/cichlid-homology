"""The summary record of one clade (0.21.0): everything a figure needs, as plain data.

One record per method and clade. The record is the product; the figure is drawn from it (`figrender`), and the site can draw the same record. It holds: the lineages found on the teleost
genes alone (`lineages.analyse`), the collapsed tree with Ensembl's support values where they exist, the copies per species outside and in each lineage, where the non-teleost genes attach
and how well supported that is, our own genes of the focal species with the lineage of their position and the weight of their placement, the calls from the gene table, and a few numbers
to triage by (overlap of the lineages, support of the split, weakest crown, stem ratio). Nothing in it depends on the method: a de novo tree from OrthoFinder gives the same record
(`participates` then takes every leaf; `support` is None until OrthoFinder trees carry support)."""
from __future__ import annotations

import sys

from .cladeprofile import group_label
from .lineages import TELEOST, analyse
from .tree import Node, leaves, postorder

SCHEMA = "cichlid-homology/clade-record/1"


GAR, REEDFISH = "lepisosteus_oculatus", "erpetoichthys_calabaricus"
DEEP_ORDER = ("shark_outside", "coelacanth_with_tetrapods", "reedfish_outside_gar_teleosts", "gar_outside_teleosts")
TOLERANCE = 0.9


def deep_order(g: Node, group_of, is_focal, gar: str = GAR, reedfish: str = REEDFISH, tolerance: float = TOLERANCE) -> dict:
    """Does the tree recover the four relationships we know? Shark outside the bony fishes; coelacanth with the tetrapods; reedfish outside gar and the teleosts; gar outside the teleosts.
    A relationship is tested only when its species are in the clade. Each is tested twice: strictly (the smallest clade that holds all the genes of the group must not hold the outgroup, so one
    stray gene breaks it) and tolerantly (some branch that does not hold the outgroup holds at least `tolerance`, 90%, of the genes of the group). The tolerant test never passes fewer clades than the
    strict one. Gar and reedfish are single genes that branch just before the teleost duplication and are the ones that wander; the coelacanth is the control. The record stores both; the
    score and the flag use the tolerant one, because a single misplaced teleost gene breaks the strict test in half of the clades."""
    post = postorder(g)
    below = {}
    for n in post:
        below[id(n)] = frozenset([id(n)]) if not n.children else frozenset().union(*(below[id(c)] for c in n.children))
    lv = [x for x in leaves(g) if x.sp and not is_focal(x.name)]
    grp = {id(x): group_of(x.sp) for x in lv}
    shark = [x for x in lv if grp[id(x)] == "cartilaginous fish"]
    coel = [x for x in lv if grp[id(x)] == "coelacanth"]
    tet = [x for x in lv if grp[id(x)] == "tetrapod"]
    gars = [x for x in lv if x.sp == gar]
    reeds = [x for x in lv if x.sp == reedfish]
    tel = [x for x in lv if grp[id(x)] in ("early teleost", "percomorph", "ovalentarian", "cichlid")]
    ids = lambda xs: frozenset(id(x) for x in xs)

    def lca(s):
        for n in post:
            if s <= below[id(n)]:
                return n
        return g

    def tolerant(group, forbid, must=None):
        """Some branch (or leaf) holds none of `forbid`, holds `must` if given, and holds at least `tolerance` of the genes of `group`."""
        need = tolerance * len(group)
        for n in post:
            b = below[id(n)]
            if (forbid & b) or (must is not None and not (must & b)):
                continue
            if len(group & b) >= need:
                return True
        return False
    out = {k: None for k in DEEP_ORDER}
    tol = {k: None for k in DEEP_ORDER}
    bony = coel + tet + gars + reeds + tel
    if shark and len(bony) >= 2:
        out["shark_outside"] = not (ids(shark) & below[id(lca(ids(bony)))])
        tol["shark_outside"] = tolerant(ids(bony), ids(shark))
    if coel and tet:
        best = min((lca(frozenset([id(c), id(t)])) for c in coel for t in tet), key=lambda n: len(below[id(n)]))
        out["coelacanth_with_tetrapods"] = tol["coelacanth_with_tetrapods"] = not (ids(gars + reeds + tel) & below[id(best)])
    if reeds and gars and tel:
        out["reedfish_outside_gar_teleosts"] = not (ids(reeds) & below[id(lca(ids(gars + tel)))])
        tol["reedfish_outside_gar_teleosts"] = tolerant(ids(tel), ids(reeds), ids(gars))
    if gars and len({x.sp for x in tel}) >= 2:
        out["gar_outside_teleosts"] = not (ids(gars) & below[id(lca(ids(tel)))])
        tol["gar_outside_teleosts"] = tolerant(ids(tel), ids(gars))

    def score(d):
        tested = [k for k in DEEP_ORDER if d[k] is not None]
        passed = [k for k in tested if d[k]]
        return len(tested), len(passed), (round(len(passed) / len(tested), 2) if tested else None)
    t, p_, sc = score(out)
    tt, tp, tsc = score(tol)
    res = dict(out)
    res.update(tested=t, passed=p_, score=sc)
    res["tolerant"] = dict(tol, tested=tt, passed=tp, score=tsc)
    return res


def _own(x: Node, is_focal, focal_species) -> bool:
    return is_focal(x.name) and x.sp == focal_species


def _labels(g, an, is_focal, focal_species, group_of):
    """The block each leaf belongs to, or None for a leaf that is drawn on its own (the non-teleost fishes). Our own genes sit inside the block of their position."""
    lab = {}
    for x in leaves(g):
        if x.sp is None:
            lab[id(x)] = None
            continue
        grp = group_of(x.sp)
        if grp == "tetrapod":
            lab[id(x)] = "tetrapod"
        elif grp in ("early teleost", "percomorph", "ovalentarian", "cichlid") or is_focal(x.name):
            m = an.member[id(x)]
            lab[id(x)] = an.letters[m] if isinstance(m, int) else ("stray" if m == "stray" else "other teleost")
        else:
            lab[id(x)] = None                       # shark, coelacanth, gar, reedfish: one by one
    return lab


def _collapse(g, an, lab, support, is_focal, focal_species, group_of, symbols, labels):
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    labs = {}
    for n in postorder(g):
        if not n.children:
            l = lab[id(n)]
            labs[id(n)] = frozenset([("leaf", id(n))]) if l is None else frozenset([l])
        else:
            labs[id(n)] = frozenset().union(*(labs[id(c)] for c in n.children))

    def sup(n):
        return support.get(id(n)) if support else None

    def emit(n):
        L = labs[id(n)]
        if len(L) == 1 and not any(isinstance(x, tuple) for x in L):
            ls = leaves(n)
            ref = [x for x in ls if not is_focal(x.name)]               # the counts are of reference genes; our own genes are listed apart
            sp = {x.sp for x in ref if x.sp}
            s = sup(n)
            return {"k": "block", "label": next(iter(L)), "n_genes": len(ref), "n_species": len(sp), "support": s[0] if s else None, "dup": s[1] if s else None,
                    "own": sorted(x.name for x in ls if _own(x, is_focal, focal_species)), "other_focal": sum(1 for x in ls if is_focal(x.name) and not _own(x, is_focal, focal_species)),
                    "species": sorted(labels.get(x, x) for x in sp) if len(sp) <= 3 else None}
        if not n.children:
            sym = symbols.get(n.name, "")
            return {"k": "leaf", "name": n.name, "species": n.sp, "group": group_of(n.sp) if n.sp else "other", "own": _own(n, is_focal, focal_species), "symbol": sym}
        kids = [emit(c) for c in n.children]
        if len(kids) == 1:
            return kids[0]
        s = sup(n)
        return {"k": "node", "support": s[0] if s else None, "dup": s[1] if s else None, "kids": kids, "split": bool(an.split and n is an.split.node)}
    return emit(g)


def _count_nodes(t) -> int:
    return 1 + sum(_count_nodes(c) for c in t.get("kids", []))


def build_record(g: Node, *, method: str, family: str, clade: str, focal_species: str, is_focal, clades, panel: list, species_labels: dict, support=None, symbols=None,
                 placements=None, calls=None, tree_mode: str = "place", root_source: str = "", overlap_t: float = 0.5, participates=None, highlight=(), gar: str = GAR,
                 reedfish: str = REEDFISH) -> dict:
    """The record of the clade whose root in a rooted tree is `g`. `panel` is the list of reference species, in drawing order. `placements` is {gene: {"lwr", "status", "n_placements"}},
    `calls` is {gene: {column: value}} from the gene table, `symbols` is {protein: symbol} from the raw dump."""
    symbols, placements, calls = symbols or {}, placements or {}, calls or {}
    def order_key(c):                                   # a lineage is named by its zebrafish gene (irx1a before irx1b), then by size: the letters mean the same in every method
        zs = sorted((symbols.get(x.name) or x.name) for x in leaves(c) if x.sp == "danio_rerio" and not is_focal(x.name))
        return (zs[0] if zs else "~", -len(leaves(c)), min(x.name for x in leaves(c)))
    an = analyse(g, is_focal, clades, support=support, overlap_t=overlap_t, participates=participates, order_key=order_key)
    group_of = lambda sp: group_label(sp, clades)
    lv = leaves(g)
    own = [x for x in lv if _own(x, is_focal, focal_species)]
    lab = _labels(g, an, is_focal, focal_species, group_of)
    tree = _collapse(g, an, lab, support, is_focal, focal_species, group_of, symbols, species_labels)
    k = len(an.letters)
    zeros = [0] * (k + 2)
    strip = {"species": list(panel), "labels": [species_labels.get(s, s) for s in panel], "groups": [group_of(s) for s in panel], "counts": [an.counts.get(s, zeros) for s in panel]}
    focal_sp = sorted(k[6:] for k in an.counts if k.startswith("focal:"))
    focal_strip = {"species": focal_sp, "labels": [species_labels.get(s, s) for s in focal_sp], "own": [s == focal_species for s in focal_sp],
                   "counts": [an.counts["focal:" + s] for s in focal_sp]}
    # the genes of the focal species: the lineage of their position and the weight of their placement
    genes = []
    for x in sorted(own, key=lambda x: x.name):
        m = an.member[id(x)]
        pl = placements.get(x.name, {})
        genes.append({"gene": x.name, "species": x.sp, "lineage": an.letters[m] if isinstance(m, int) else ("stray" if m == "stray" else "outside"),
                      "lwr": pl.get("lwr"), "status": pl.get("status"), "n_placements": pl.get("n_placements"), "calls": calls.get(x.name, {})})
    # anchors: the zebrafish and human genes of each lineage, which name a lineage the same way in every method
    anchors = {}
    for sp_key, sp_name in (("zebrafish", "danio_rerio"), ("human", "homo_sapiens")):
        for x in lv:
            if x.sp == sp_name and not is_focal(x.name):
                m = an.member[id(x)]
                key = an.letters[m] if isinstance(m, int) else "outside"
                anchors.setdefault(key, {}).setdefault(sp_key, []).append({"gene": x.name, "symbol": symbols.get(x.name, "")})
    # supports of the nodes that matter, for the second panel
    items, sp_support = [], None
    if an.split is not None and support:
        s = support.get(id(an.split.node))
        sp_support = s
        items.append({"label": "the split" + (" (Ensembl: duplication)" if s and s[1] == "Y" else ""), "value": s[0] if s else None, "kind": "split"})
    for l in an.lineages:
        if an.split is not None:
            items.append({"label": f"lineage {l['id']}: crown of {l['species']} species", "value": l["crown_support"], "kind": "crown"})
    for o in an.outgroup_inside:
        items.append({"label": f"{species_labels.get(o['species'], o['species'])} attaches inside {o['lineage']}", "value": o["attach_support"], "kind": "outgroup"})
    stems = [l["stem"] for l in an.lineages if l.get("stem")]
    crowns = [i["value"] for i in items if i["kind"] == "crown" and i["value"] is not None]
    triage = {"overlap": round(an.split.overlap, 3) if an.split else None,
              "shared_species": (len(an.S[id(an.split.lineages[0])] & an.S[id(an.split.lineages[1])]) if an.split else None),
              "split_support": sp_support[0] if sp_support else None, "weakest_crown": min(crowns) if crowns else None, "outgroup_inside": len(an.outgroup_inside),
              "stem_ratio": round(max(stems) / min(stems), 2) if len(stems) >= 2 and min(stems) > 0 else None}
    deep = deep_order(g, group_of, is_focal, gar, reedfish)
    triage["deep_order_score"] = deep["tolerant"]["score"]
    triage["deep_order_score_strict"] = deep["score"]
    flags = []
    if deep["tolerant"]["tested"] and deep["tolerant"]["passed"] < deep["tolerant"]["tested"]:
        flags.append("deep_order_inconsistent")
    if an.mode != "split":
        flags.append("no_split" if an.mode == "single" else "no_teleost_genes")
    if an.outgroup_inside:
        flags.append("outgroup_inside_lineage")
    if triage["weakest_crown"] is not None and triage["weakest_crown"] < 50:
        flags.append("weak_crown_support")
    if an.split and sum(c[-1] for c in an.counts.values()) > 0:
        flags.append("stray_teleost_genes")
    if not support:
        flags.append("no_support_values")
    return {"schema": SCHEMA, "method": method, "family": family, "clade": clade, "tree_mode": tree_mode, "root_source": root_source, "focal_species": focal_species,
            "n_leaves": len(lv), "n_own": len(own), "mode": an.mode,
            "split": ({"kind": an.split.kind, "level": an.split.level, "overlap": round(an.split.overlap, 3), "n_lineages": len(an.letters)} if an.split else None),
            "lineages": [{k_: v for k_, v in l.items() if k_ != "crown_node"} for l in an.lineages], "letters": an.letters,
            "outgroup": {"inside": an.outgroup_inside, "outside": an.outgroup_outside}, "strip": strip, "focal_strip": focal_strip, "genes": genes, "anchors": anchors,
            "deep_order": deep, "tree": tree, "support_items": items, "support_available": bool(support), "triage": triage, "flags": flags}
