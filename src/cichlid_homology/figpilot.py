"""The figure pilot on the Compara data (0.21.0): one summary record and one figure per clade, for a stratified sample of genes (or for all of them).

For each family of the sample: the reference tree of the family is split into clades exactly as `clade-split` does (so the clade ids of the gene table apply); the placed tree of the family (the reference
tree with our genes placed on it) is rooted where the reference tree is rooted; the clade is found in it by its reference leaves; Ensembl's support values are read from the raw dump when it is given; the
record of the clade is built (`figrecord`), and the figure is drawn from the record (`figrender`). Records are the product (`clade_records.jsonl`), figures are for looking. TreeFam is not used."""
from __future__ import annotations

import collections
import csv
import html
import json
import random
import re
import time
import traceback
from pathlib import Path

from . import nhx
from .cladesplit import split_family
from .famsurvey import reference_files
from .figrecord import build_record
from .figrender import render_many
from .parsimony import ParsimonyContext
from .rooting import root_tree
from .suppmap import attach_support, reroot_like
from .tree import Node, leaves, parse_newick, postorder

PAIR_STATUSES = ("paired", "paired (resolved by synteny)", "pair retained in few teleosts (synteny-supported)", "possible pair (the other teleosts lost a copy)")
CALL_COLS = ("symbol", "b1_category", "b1_confidence", "ohnolog_status", "ohnolog_symbol", "human_symbol", "zebrafish_symbol", "clade_type", "teleost_trend", "rel_human", "rel_zebrafish")

STRATA = {
    "published_pair": "a high-confidence pair that matches a published pair",
    "ab_not_ohnolog": "named a/b, but both genes are called single copy",
    "pair_other": "a pair by synteny, a retained or a possible pair",
    "human_1to1_clean": "1:1 with human, clean human side",
    "human_1to2": "1:2 with human, not a pair",
    "out_paralog": "an out-paralog clade",
    "later_duplicate": "a later duplicate among our genes",
    "expansion": "an expansion",
    "ambiguous_partner": "a partner that is ambiguous or not found",
    "no_human_clade": "a clade with no human gene",
    "low_placement": "a placement weight under 0.5",
    "denovo": "a gene in a de novo tree",
}


def _is_real(sym: str) -> bool:
    return bool(sym) and not sym.upper().startswith("LOC")


def ab_pairs(rows: dict) -> list:
    """Pairs of genes whose symbols are a stem plus a and the same stem plus b, one gene each, both called single copy."""
    by = collections.defaultdict(list)
    for g, r in rows.items():
        s = (r.get("symbol") or "").lower()
        if _is_real(s) and len(s) > 3 and s[-1] in "ab":
            by[s[:-1]].append((s[-1], g))
    out = []
    for stem, lst in by.items():
        d = collections.defaultdict(list)
        for suf, g in lst:
            d[suf].append(g)
        if len(d.get("a", [])) == 1 and len(d.get("b", [])) == 1:
            a, b = d["a"][0], d["b"][0]
            if rows[a].get("ohnolog_status") == rows[b].get("ohnolog_status") == "single copy (duplicate lost)":
                out.append((a, b))
    return out


def select_genes(rows: dict, per_stratum: int = 10, seed: int = 1, placement: dict | None = None, explicit=None) -> list:
    """[{gene, stratum}] from the gene table (`rows`: gene_id -> row). Every stratum is sampled without replacement; a gene keeps the first stratum that drew it."""
    placement = placement or {}
    cands = {k: [] for k in STRATA}
    abp = ab_pairs(rows)
    cands["ab_not_ohnolog"] = [g for p in abp for g in p]
    for g, r in rows.items():
        st = r.get("ohnolog_status", "")
        if st == "paired" and r.get("ohnolog_confidence") == "high" and (r.get("published_pair") or "").startswith("the two"):
            cands["published_pair"].append(g)
        elif st in PAIR_STATUSES[1:]:
            cands["pair_other"].append(g)
        if r.get("rel_human") == "1:1" and r.get("human_clean_tier") == "clean":
            cands["human_1to1_clean"].append(g)
        if r.get("rel_human") == "1:2" and st not in PAIR_STATUSES:
            cands["human_1to2"].append(g)
        if r.get("out_paralog") == "yes":
            cands["out_paralog"].append(g)
        if st in ("later duplicate", "later duplicates (three or more focal genes)"):
            cands["later_duplicate"].append(g)
        if st == "expansion":
            cands["expansion"].append(g)
        if st in ("partner ambiguous (three or more focal genes)", "partner not found on a partner chromosome"):
            cands["ambiguous_partner"].append(g)
        if (r.get("clade_type") or "").startswith(("fish-only clade", "tetrapod clade, no human")):
            cands["no_human_clade"].append(g)
        p = placement.get(g)
        if p:
            try:
                if p.get("lwr") not in ("", None) and float(p["lwr"]) < 0.5:
                    cands["low_placement"].append(g)
            except ValueError:
                pass
            if p.get("mode") == "denovo":
                cands["denovo"].append(g)
    chosen, seen = [], set()
    rng = random.Random(seed)
    for name in STRATA:
        pool = [g for g in sorted(set(cands[name])) if g not in seen and rows[g].get("family") not in ("", "NoHit", None)]
        take = rng.sample(pool, min(per_stratum, len(pool))) if name != "ab_not_ohnolog" else _take_pairs(pool, abp, per_stratum, rng, seen)
        for g in sorted(take):
            chosen.append({"gene": g, "stratum": name})
            seen.add(g)
    for g in explicit or []:
        if g in rows and g not in seen:
            chosen.append({"gene": g, "stratum": "requested"})
            seen.add(g)
    return chosen


def _take_pairs(pool, abp, n, rng, seen):
    pairs = [p for p in abp if p[0] in pool and p[1] in pool and p[0] not in seen and p[1] not in seen]
    return [g for p in rng.sample(pairs, min(n, len(pairs))) for g in p]


def panel_order(master_newick: str, keep: set) -> list:
    """The reference species in species-tree order, the smaller branch first at every node (so the tetrapods follow the cartilaginous fish and the teleosts come last)."""
    root = parse_newick(master_newick)
    size = {}
    for n in postorder(root):
        size[id(n)] = 1 if not n.children else sum(size[id(c)] for c in n.children)
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        if n.children:
            stack.extend(sorted(n.children, key=lambda c: size[id(c)], reverse=True))        # the stack pops the smallest first
        elif n.name in keep:
            out.append(n.name)
    return out


def find_clade_node(root: Node, names: frozenset, is_focal=lambda name: False):
    """The smallest node that holds every one of `names` that is in the tree. Returns (node, exact): exact when every leaf below the node is one of `names` or one of our own genes."""
    present = {x.name for x in leaves(root)} & names
    if not present:
        return None, False
    best, sets = None, {}
    for n in postorder(root):
        s = frozenset([n.name]) & present if not n.children else frozenset().union(*(sets[id(c)] for c in n.children))
        sets[id(n)] = s
        if s == present and best is None:
            best = n
    exact = best is not None and all((x.name in names) or is_focal(x.name) for x in leaves(best))
    return best, exact


def _read_tsv(path) -> list:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run_figure_pilot(*, genes: list, rows: dict, prepared, froot, placement: dict, prefixes: set, load_family, load_ls, st, clades, master_newick: str, out_dir, dump=None,
                     workers: int = 1, dpi: int = 110, fmt: str = "png", method: str = "Compara", progress=print) -> dict:
    """`genes`: [{gene, stratum}]; `rows`: the gene table by gene_id; `prepared`: the folder with the reference trees; `froot`: the folder of the family folders of the method;
    `load_family(folder, prefixes, family_id)`; `load_ls(names)` returns the leaf-to-species table for those leaf names."""
    t0 = time.time()
    out_dir = Path(out_dir)
    (out_dir / "figures").mkdir(parents=True, exist_ok=True)
    summary = {"method": method, "genes_requested": len(genes), "skipped": collections.Counter(), "errors": []}
    ref_files, _ = reference_files(Path(prepared) / "trees")
    want = collections.defaultdict(list)                                    # family -> [(gene, stratum)]
    for g in genes:
        r = rows.get(g["gene"], {})
        fam = r.get("family")
        if not fam or fam == "NoHit" or not (Path(froot) / fam / "done.json").exists():
            summary["skipped"]["no finished family"] += 1
            continue
        want[fam].append((g["gene"], g["stratum"]))
    fams = {}
    for fam in sorted(want):
        try:
            fams[fam] = load_family(Path(froot) / fam, prefixes, fam)
        except Exception:
            summary["skipped"]["family not readable"] += len(want[fam])
            continue
    refs = {}
    for fam in fams:
        if fam in ref_files and "tree" in ref_files[fam]:
            refs[fam] = ref_files[fam]["tree"].read_text(errors="replace")
    names = set()
    for fam, f in fams.items():
        names |= set(f.leaves)
    for text in refs.values():
        names |= {x.name for x in leaves(parse_newick(text))}
    ls = load_ls(names)
    is_focal = ls.is_project_gene
    ctx = ParsimonyContext(st, clades)
    support_all = {}
    if dump:
        wanted = {fam: {x.name for x in leaves(parse_newick(t))} for fam, t in refs.items()}
        progress(f"reading support values for {len(wanted)} families from {dump} ...")
        support_all = nhx.collect_support(dump, wanted, progress=lambda n, k: progress(f"  {n} trees read, {k} of {len(wanted)} families found"))
    summary["families"] = len(fams)
    summary["support_families_found"] = len(support_all)
    # ---- one rooted tree per family, and the clade nodes
    todo = []                                                                # (family, clade id, node, support, symbols, root_source, flags)
    for fam, f in sorted(fams.items()):
        try:
            root = parse_newick(f.tree)
            for x in leaves(root):
                x.sp = ls.species(x.name)
            ref_root = None
            gl = []
            if fam in refs and f.tree_mode == "place":
                ref_text = refs[fam]
                try:
                    split_family(fam, ref_text, ls, st, clades, ctx, group_leaves=gl)
                except Exception:
                    gl = []
                ref_root = parse_newick(ref_text)
            rooted, how = None, ""
            if ref_root is not None:
                rooted = reroot_like(root, ref_root)
                how = "like the reference root" if rooted else ""
            if rooted is None:
                rr = root_tree(root, st, f.rooted, k=1)
                rooted, how = rr.roots[0], rr.how
            sup = symbols = None
            if fam in support_all:
                universe = {x.name for x in leaves(ref_root)} if ref_root is not None else {x.name for x in leaves(rooted) if not is_focal(x.name)}
                sup = attach_support(rooted, support_all[fam]["support"], universe)
                symbols = support_all[fam]["symbols"]
            need = {}
            for gene, stratum in want[fam]:
                cid = rows[gene].get("clade_id") or ""
                key = cid if (cid and gl) else f"{fam}:whole"
                need.setdefault(key, []).append(gene)
            gmap = {f"{fam}:{name}": s for name, s in gl}
            for key, gl_genes in need.items():
                flags = []
                if key in gmap:
                    node, exact = find_clade_node(rooted, gmap[key], is_focal)
                    if node is None:
                        summary["skipped"]["clade not found in the placed tree"] += len(gl_genes)
                        continue
                    if not exact:
                        flags.append("clade_not_exact_in_placed_tree")
                else:
                    node = rooted
                    flags.append("whole_tree" if fam in refs and gl else "no_reference_clades")
                todo.append((fam, key, node, sup, symbols, how, flags, gl_genes, f.tree_mode))
        except Exception:
            summary["errors"].append({"family": fam, "error": traceback.format_exc()[-500:]})
    # ---- the species axis, the records
    ref_species = set()
    for fam, key, node, *_ in todo:
        ref_species |= {x.sp for x in leaves(node) if x.sp and not is_focal(x.name)}
    panel = panel_order(master_newick, ref_species)
    labels = {s: s.replace("_", " ").capitalize() for s in ref_species | set(panel)}
    records, jobs, index = [], [], []
    stratum_of = {g["gene"]: g["stratum"] for g in genes}
    for fam, key, node, sup, symbols, how, flags, gl_genes, tmode in todo:
        focal_species = ls.species(gl_genes[0])
        labels[focal_species] = gl_genes[0].split("_")[0]
        for x in leaves(node):
            if x.sp:
                labels.setdefault(x.sp, x.sp.replace("_", " ").capitalize())
        in_clade = {x.name for x in leaves(node)}
        place = {g: placement.get(g, {}) for g in in_clade if is_focal(g)}
        calls = {g: {c: rows[g].get(c, "") for c in CALL_COLS} for g in in_clade if g in rows}
        rec = build_record(node, method=method, family=fam, clade=key, focal_species=focal_species, is_focal=is_focal, clades=clades, panel=panel, species_labels=labels, support=sup,
                           symbols=symbols, placements=place, calls=calls, tree_mode=tmode, root_source=how, participates=(None if tmode == "place" else (lambda x: True)), highlight=gl_genes)
        rec["flags"] = rec["flags"] + flags
        for g in gl_genes:
            if g not in in_clade:                                          # placed on the branch that leads to the clade, outside its node
                rec["genes"].append({"gene": g, "species": focal_species, "lineage": "outside the clade", "lwr": placement.get(g, {}).get("lwr"), "status": placement.get(g, {}).get("status"),
                                     "n_placements": placement.get(g, {}).get("n_placements"), "calls": {c: rows[g].get(c, "") for c in CALL_COLS}})
                rec["flags"].append("gene_outside_clade_node")
        records.append(rec)
        fname = re.sub(r"[^A-Za-z0-9_.-]", "_", key) + "." + fmt
        jobs.append((rec, out_dir / "figures" / fname))
        for g in gl_genes:
            index.append({"gene": g, "stratum": stratum_of.get(g, ""), "clade": key, "family": fam, "figure": f"figures/{fname}", "record": rec})
    with open(out_dir / "clade_records.jsonl", "w") as fh:
        for rec in records:
            fh.write(json.dumps(rec, separators=(",", ":")) + "\n")
    progress(f"rendering {len(jobs)} figures with {workers} worker(s) ...")
    results = render_many(jobs, workers, dpi, fmt)
    bad = [r for r in results if r[1]]
    summary["figure_errors"] = bad[:20]
    _write_index(out_dir, index, rows, summary)
    flags = collections.Counter(f for r in records for f in r["flags"])
    summary.update(clades=len(records), figures=len(results) - len(bad), flags=dict(flags), modes=dict(collections.Counter(r["mode"] for r in records)),
                   splits=dict(collections.Counter((r["split"] or {}).get("kind", "none") for r in records)), seconds=round(time.time() - t0, 1))
    summary["skipped"] = dict(summary["skipped"])
    (out_dir / "figure_pilot.json").write_text(json.dumps(summary, indent=1, default=str))
    return summary


def _write_index(out_dir: Path, index: list, rows: dict, summary: dict):
    with open(out_dir / "figure_index.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "symbol", "stratum", "family", "clade", "mode", "split_kind", "overlap", "split_support", "weakest_crown", "outgroup_inside", "stem_ratio", "flags", "figure"])
        for x in sorted(index, key=lambda x: (x["stratum"], x["gene"])):
            r, t = x["record"], x["record"]["triage"]
            w.writerow([x["gene"], rows.get(x["gene"], {}).get("symbol", ""), x["stratum"], x["family"], x["clade"], r["mode"], (r["split"] or {}).get("kind", ""), t["overlap"], t["split_support"],
                        t["weakest_crown"], t["outgroup_inside"], t["stem_ratio"], ",".join(r["flags"]), x["figure"]])
    parts = ["<html><head><meta charset='utf-8'><title>Figure pilot</title><style>body{font:14px sans-serif;margin:20px}table{border-collapse:collapse}td{border:1px solid #ccc;padding:6px;vertical-align:top}"
             "img{width:560px}h2{margin-top:28px}.f{color:#b00020}.m{color:#666}</style></head><body><h1>Figure pilot (Compara)</h1>"]
    by = collections.defaultdict(list)
    for x in index:
        by[x["stratum"]].append(x)
    for name in list(STRATA) + ["requested"]:
        lst = sorted(by.get(name, []), key=lambda x: x["gene"])
        if not lst:
            continue
        parts.append(f"<h2>{html.escape(name)}: {html.escape(STRATA.get(name, 'requested genes'))} ({len(lst)})</h2><table>")
        for x in lst:
            r, t, g = x["record"], x["record"]["triage"], rows.get(x["gene"], {})
            lines = [f"<b>{html.escape(g.get('symbol') or x['gene'])}</b> {html.escape(x['gene'])}", html.escape(g.get("b1_category", "")), html.escape(g.get("ohnolog_status", "")),
                     f"clade {html.escape(x['clade'])}", f"mode: {r['mode']}" + (f", {r['split']['kind']} split" if r["split"] else ""),
                     f"overlap {t['overlap']}, split support {t['split_support']}, weakest crown {t['weakest_crown']}, outgroup inside {t['outgroup_inside']}, stem ratio {t['stem_ratio']}"]
            parts.append(f"<tr><td style='width:340px'>{'<br>'.join(l for l in lines if l)}<br><span class='f'>{html.escape(', '.join(r['flags']))}</span></td>"
                         f"<td><a href='{x['figure']}'><img loading='lazy' src='{x['figure']}'></a></td></tr>")
        parts.append("</table>")
    parts.append("</body></html>")
    (out_dir / "figure_index.html").write_text("\n".join(parts))
