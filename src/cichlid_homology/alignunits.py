"""The alignment unit of a gene or a pair of genes: the smallest clade that holds a human gene (0.24.0).

For analysing the protein alignment at the level of a gene, a pair of duplicates or any set of our genes, the sequences to look at are the ones in the smallest clade of the family tree that holds the gene
(or both genes of the pair) and a human gene. Rules, in order:

1. **human**: the smallest clade of the rooted family tree that holds the gene(s) and at least one human gene;
2. **tetrapod**: when the family has no human gene, the smallest clade that holds the gene(s) and a tetrapod gene (a sarcopterygian that is not the coelacanth);
3. **family**: when the family has neither, the whole family.

The tree is the family tree of the method with our genes in it (`tree.nwk` of the family folder: the reference tree with our genes placed on it, or a de novo tree), rooted like the clade split roots it: like the
root of the reference tree when there is one, otherwise by reconciliation with the master species tree. The clade is a node of that tree and is not cut at ancient duplications (a clade of the clade table can
be larger or smaller): for a gene whose own clade of the clade table has no human gene because the human gene is on the other side of an ancient duplication, the unit climbs past the duplication to the
smallest clade that has one. Every unit lists all the leaves below its node: our genes of every project species and the reference proteins, with a flag for the reference leaves of the shadowed species
(the Ensembl copy of a genome that we annotate too). The leaf names are the names in `combined.aln` of the family, so a unit is a list of alignment rows. The residue step (which column is which position)
is project 4's; this module only says which rows.

The units are shared: a family that has 40 of our genes in the same part of its tree gives one unit. A unit id is `<family>:<H|T|F><k>`: H, T or F for the rule, k the rank of the unit by size within the family
and rule (1 is the smallest)."""
from __future__ import annotations

import collections
import csv
import gzip
import json
import statistics
from pathlib import Path

LETTER = {"human": "H", "tetrapod": "T", "family": "F"}
UNIT_COLS = ["unit_id", "family", "basis", "whole_family", "tree_mode", "rooting", "n_genes", "n_species", "n_family_genes", "share_of_family", "n_human_genes", "n_tetrapod_genes", "n_focal_genes",
             "n_project_genes", "n_reference_genes", "n_shadowed_reference", "genes_requested", "pairs_requested", "focal_genes", "human_genes", "alignment_file", "tree_file"]
MEMBER_COLS = ["unit_id", "gene_id", "species", "kind", "ladder_level", "is_human", "is_tetrapod", "shadowed_reference"]
GENE_COLS = ["gene_id", "family", "clade_id", "unit_id", "basis", "unit_genes", "whole_family", "status"]
PAIR_COLS = ["gene_a", "gene_b", "group", "family_a", "family_b", "unit_id", "basis", "unit_genes", "whole_family", "unit_a", "unit_b", "status"]


class UnitTree:
    """A rooted tree with the leaves in order, the interval of leaves below each node, and prefix counts of the human and the tetrapod leaves."""

    def __init__(self, root, species_of, is_human, is_tetrapod):
        self.root = root
        pre, stack = [], [root]
        while stack:
            n = stack.pop()
            pre.append(n)
            stack.extend(reversed(n.children))
        self.leaves = []
        for n in pre:
            if not n.children:
                n.lo = len(self.leaves)
                n.hi = n.lo + 1
                self.leaves.append(n)
        self.par = {}
        for n in pre:
            for c in n.children:
                self.par[id(c)] = n
        for n in reversed(pre):
            if n.children:
                n.lo, n.hi = min(c.lo for c in n.children), max(c.hi for c in n.children)
        self.by_name = {x.name: x for x in self.leaves}
        self.sp = [species_of(x.name) for x in self.leaves]
        self.hpre, self.tpre = [0], [0]
        for s in self.sp:
            self.hpre.append(self.hpre[-1] + (1 if s is not None and is_human(s) else 0))
            self.tpre.append(self.tpre[-1] + (1 if s is not None and is_tetrapod(s) else 0))

    def lca(self, leaves):
        lo, hi = min(x.lo for x in leaves), max(x.hi for x in leaves)
        n = leaves[0]
        while not (n.lo <= lo and n.hi >= hi):
            n = self.par[id(n)]
        return n

    def unit(self, names):
        """(node, basis) of the unit of the genes `names`, or (None, reason)."""
        names = list(dict.fromkeys(names))
        if not names:
            return None, "no gene"
        miss = [x for x in names if x not in self.by_name]
        if miss:
            return None, "gene not in the family tree"
        node = self.lca([self.by_name[x] for x in names])
        for basis, pre in (("human", self.hpre), ("tetrapod", self.tpre)):
            n = node
            while n is not None and pre[n.hi] - pre[n.lo] == 0:
                n = self.par.get(id(n))
            if n is not None:
                return n, basis
        return self.root, "family"


def _kind(name, focal_prefix, is_project):
    if name.startswith(focal_prefix + "_"):
        return "focal"
    return "project" if is_project(name) else "reference"


def family_units(fid, root, how, tree_mode, species_of, is_human, is_tetrapod, is_project, shadowed, level_of, focal_prefix, genes, pairs, files=("", "")):
    """(unit rows, member rows, {gene: (unit key or None, status)}, {pair index: (unit key or None, status)}) of one family. `genes`: [gene ids]; `pairs`: [(index, gene_a, gene_b)]."""
    T = UnitTree(root, species_of, is_human, is_tetrapod)
    units = {}                                                       # (lo, hi, basis) -> {"node", "genes", "pairs"}
    gmap, pmap = {}, {}
    for g in genes:
        n, b = T.unit([g])
        if n is None:
            gmap[g] = (None, b)
            continue
        k = (n.lo, n.hi, b)
        units.setdefault(k, {"node": n, "genes": 0, "pairs": 0})["genes"] += 1
        gmap[g] = (k, "ok")
    for i, a, b_ in pairs:
        n, b = T.unit([a, b_])
        if n is None:
            pmap[i] = (None, b)
            continue
        k = (n.lo, n.hi, b)
        units.setdefault(k, {"node": n, "genes": 0, "pairs": 0})["pairs"] += 1
        pmap[i] = (k, "ok")
    ids, count = {}, collections.Counter()
    for k in sorted(units, key=lambda k: (k[2] != "human", k[2] != "tetrapod", k[1] - k[0], k[0])):
        count[k[2]] += 1
        ids[k] = f"{fid}:{LETTER[k[2]]}{count[k[2]]}"
    nfam = len(T.leaves)
    urows, mrows = [], []
    for k, u in units.items():
        lo, hi, basis = k
        lv = T.leaves[lo:hi]
        kinds = [_kind(x.name, focal_prefix, is_project) for x in lv]
        sps = T.sp[lo:hi]
        humans = [x.name for x, s in zip(lv, sps) if s is not None and is_human(s)]
        n_tet = sum(1 for s in sps if s is not None and is_tetrapod(s))
        shad = [s is not None and kd == "reference" and s in shadowed for s, kd in zip(sps, kinds)]
        urows.append({"unit_id": ids[k], "family": fid, "basis": basis, "whole_family": "yes" if (lo, hi) == (0, nfam) else "no", "tree_mode": tree_mode, "rooting": how,
                      "n_genes": len(lv), "n_species": len({s for s in sps if s is not None}), "n_family_genes": nfam, "share_of_family": round(len(lv) / nfam, 4), "n_human_genes": len(humans),
                      "n_tetrapod_genes": n_tet, "n_focal_genes": kinds.count("focal"), "n_project_genes": kinds.count("project") + kinds.count("focal"), "n_reference_genes": kinds.count("reference"),
                      "n_shadowed_reference": sum(shad), "genes_requested": u["genes"], "pairs_requested": u["pairs"], "focal_genes": ";".join(x.name for x, kd in zip(lv, kinds) if kd == "focal"),
                      "human_genes": ";".join(humans), "alignment_file": files[0], "tree_file": files[1]})
        for x, s, kd, sh in zip(lv, sps, kinds, shad):
            mrows.append({"unit_id": ids[k], "gene_id": x.name, "species": s or "", "kind": kd, "ladder_level": (level_of(s) if s else "") or "", "is_human": "yes" if (s and is_human(s)) else "no",
                          "is_tetrapod": "yes" if (s and is_tetrapod(s)) else "no", "shadowed_reference": "yes" if sh else "no"})
    gres = {g: ((ids[v[0]], v[0][2], len(T.leaves[v[0][0]:v[0][1]]), v[0][:2] == (0, nfam)) if v[0] else None, v[1]) for g, v in gmap.items()}
    pres = {i: ((ids[v[0]], v[0][2], len(T.leaves[v[0][0]:v[0][1]]), v[0][:2] == (0, nfam)) if v[0] else None, v[1]) for i, v in pmap.items()}
    return urows, mrows, gres, pres


def _q(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))] if v else None


def run_units(*, rows: dict, pairs: list, froot, ref_files: dict, prefixes: set, focal_prefix: str, load_family, load_ls, rooter, clades_info, out_dir, progress=print, limit: int = 0) -> dict:
    """`rows`: the gene table by gene_id; `pairs`: [{gene_a, gene_b, group}]; `froot`: the folder of the family folders; `ref_files`: {family: {"tree": path}} of the reference trees;
    `load_family(folder, prefixes, family_id)`; `load_ls(names)` returns the leaf-to-species table; `rooter(family, family_obj, ref_text, ls)` returns (rooted root, how);
    `clades_info`: (is_human, is_tetrapod, level_of) callables on a species."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    is_human, is_tetrapod, level_of = clades_info
    summary = {"genes_in_table": len(rows), "genes_with_a_family": 0, "families": 0, "skipped": collections.Counter(), "errors": []}
    fam_of = {g: r.get("family") for g, r in rows.items()}
    want_g, want_p = collections.defaultdict(list), collections.defaultdict(list)
    pair_rows = []
    for i, p in enumerate(pairs):
        a, b = p["gene_a"], p["gene_b"]
        fa, fb = fam_of.get(a), fam_of.get(b)
        pair_rows.append({"gene_a": a, "gene_b": b, "group": p.get("group", ""), "family_a": fa or "", "family_b": fb or "", "unit_id": "", "basis": "", "unit_genes": "", "whole_family": "", "unit_a": "", "unit_b": "", "status": ""})
        if not fa or not fb or "NoHit" in (fa, fb):
            pair_rows[-1]["status"] = "a gene has no family"
        elif fa != fb:
            pair_rows[-1]["status"] = "the genes are in different families"
        else:
            want_p[fa].append((i, a, b))
    for g, f in fam_of.items():
        if f and f != "NoHit":
            want_g[f].append(g)
            summary["genes_with_a_family"] += 1
    fams_all = sorted(set(want_g) | set(want_p))
    if limit:
        fams_all = fams_all[:limit]
    fams, names = {}, set()
    for fam in fams_all:
        folder = Path(froot) / fam
        if not (folder / "done.json").exists():
            summary["skipped"]["no finished family folder"] += len(want_g[fam])
            continue
        try:
            fams[fam] = load_family(folder, prefixes, fam)
        except Exception:
            summary["skipped"]["family not readable"] += len(want_g[fam])
            continue
        names |= set(fams[fam].leaves)
    ls = load_ls(names)
    shadowed = set(ls.shadowed())
    gene_rows, unit_rows, mem_rows = [], [], []
    gene_out, pair_out = {}, {}
    done = 0
    mfh = gzip.open(out_dir / "alignment_unit_members.tsv.gz", "wt", newline="")
    mw = csv.DictWriter(mfh, fieldnames=MEMBER_COLS, delimiter="\t")
    mw.writeheader()
    try:
        for fam in sorted(fams):
            f = fams[fam]
            try:
                ref = ref_files.get(fam, {}).get("tree")
                root, how = rooter(fam, f, ref.read_text(errors="replace") if ref else None, ls)
                files = (str(Path(froot) / fam / "combined.aln"), str(Path(froot) / fam / "tree.nwk"))
                u, m, g, p = family_units(fam, root, how, f.tree_mode, ls.species, is_human, is_tetrapod, ls.is_project_gene, shadowed, level_of, focal_prefix, want_g[fam],
                                          [(i, a, b) for i, a, b in want_p[fam]], files)
            except Exception as e:
                summary["errors"].append({"family": fam, "error": f"{type(e).__name__}: {e}"[:200]})
                continue
            unit_rows.extend(u)
            mw.writerows(m)
            gene_out.update({k: (fam, v) for k, v in g.items()})
            pair_out.update(p)
            done += 1
            if done % 2000 == 0:
                progress(f"  {done} families")
    finally:
        mfh.close()
    summary["families"] = done
    # the gene and pair tables
    for g, r in sorted(rows.items()):
        f = r.get("family")
        row = {"gene_id": g, "family": f or "", "clade_id": r.get("clade_id", ""), "unit_id": "", "basis": "", "unit_genes": "", "whole_family": "", "status": ""}
        if not f or f == "NoHit":
            row["status"] = "no family"
        elif g in gene_out:
            res, st_ = gene_out[g][1]
            if res:
                row.update(unit_id=res[0], basis=res[1], unit_genes=res[2], whole_family="yes" if res[3] else "no", status="ok")
            else:
                row["status"] = st_
        else:
            row["status"] = "family not read"
        gene_rows.append(row)
    gby = {r["gene_id"]: r for r in gene_rows}
    for i, pr in enumerate(pair_rows):
        if i in pair_out:
            res, st_ = pair_out[i]
            if res:
                pr.update(unit_id=res[0], basis=res[1], unit_genes=res[2], whole_family="yes" if res[3] else "no", status="ok")
            else:
                pr["status"] = st_
        elif not pr["status"]:
            pr["status"] = "family not read"
        pr["unit_a"], pr["unit_b"] = gby.get(pr["gene_a"], {}).get("unit_id", ""), gby.get(pr["gene_b"], {}).get("unit_id", "")
    for name, cols, rs in (("alignment_units", UNIT_COLS, sorted(unit_rows, key=lambda r: r["unit_id"])), ("gene_to_unit", GENE_COLS, gene_rows), ("pair_to_unit", PAIR_COLS, pair_rows)):
        with open(out_dir / f"{name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rs)
    summary["units"] = len(unit_rows)
    by_basis = collections.defaultdict(list)
    for r in unit_rows:
        by_basis[r["basis"]].append(r["n_genes"])
    summary["unit_size_by_basis"] = {b: {"units": len(v), "min": min(v), "median": statistics.median(v), "p90": _q(v, 0.9), "max": max(v)} for b, v in by_basis.items()}
    summary["genes_by_basis"] = dict(collections.Counter(r["basis"] or r["status"] for r in gene_rows))
    summary["pairs_by_basis"] = dict(collections.Counter(r["basis"] or r["status"] for r in pair_rows))
    summary["genes_whose_unit_is_the_whole_family"] = sum(1 for r in gene_rows if r["whole_family"] == "yes")
    summary["units_over_1000_genes"] = sum(1 for r in unit_rows if r["n_genes"] > 1000)
    summary["skipped"] = dict(summary["skipped"])
    (out_dir / "alignment_units.json").write_text(json.dumps(summary, indent=1, default=str))
    (out_dir / "alignment_units.md").write_text(report_text(summary, len(pairs)))
    return summary


def report_text(s: dict, n_pairs: int) -> str:
    L = ["# Alignment units: the smallest clade that holds a human gene", "",
         f"{s['genes_with_a_family']} of {s['genes_in_table']} genes have a family; {s['families']} family trees were read; {s['units']} units. Rules: the smallest clade with a human gene; else with a tetrapod gene; else the whole family.", "",
         "| rule | units | smallest | median | 90th percentile | largest |", "|---|---|---|---|---|---|"]
    for b in ("human", "tetrapod", "family"):
        if b in s["unit_size_by_basis"]:
            d = s["unit_size_by_basis"][b]
            L.append(f"| {b} | {d['units']} | {d['min']} | {d['median']} | {d['p90']} | {d['max']} |")
    L += ["", f"Genes by rule: {s['genes_by_basis']}", f"Pairs by rule ({n_pairs} pairs): {s['pairs_by_basis']}",
          f"Genes whose unit is the whole family: {s['genes_whose_unit_is_the_whole_family']}. Units with more than 1000 genes: {s['units_over_1000_genes']}.", f"Skipped: {s['skipped']}. Errors: {len(s['errors'])}.", ""]
    return "\n".join(L)
