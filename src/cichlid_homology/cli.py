"""The cichlid-homology command line."""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

from . import __version__
from .paths import default_out
from .read import load_family, read_tsv, sample_families


def _med(v):
    return sorted(v)[len(v) // 2] if v else None


def _survey(args) -> int:
    from cichlid_data import Catalog, get_layout

    layout = get_layout()
    cat = Catalog.load(layout)
    prefixes = {s.prefix for s in cat.species(release="v1")}
    report = {"cichlid_homology": __version__}
    for part in ("TreeFam", "Compara"):
        root, chosen, rows = sample_families(layout, part, args.species, args.n)
        by_fam = collections.defaultdict(list)
        for r in rows:
            by_fam[r["best_family"]].append(r["gene_id"])
        modes, rooted, kids = collections.Counter(), collections.Counter(), collections.Counter()
        keys, probs = collections.Counter(), collections.Counter()
        nleaf, nfocal, refs = [], [], []
        for f in chosen:
            fam = load_family(root / f, prefixes, f)
            keys.update(fam.done.keys())
            modes[fam.tree_mode] += 1
            rooted[fam.rooted] += 1
            kids[fam.root_children] += 1
            probs.update(fam.problems)
            nleaf.append(len(fam.leaves))
            nfocal.append(len(fam.focal))
            if len(refs) < 10:
                refs += fam.reference[:2]
            leafset = set(fam.leaves)
            if any(g not in leafset for g in by_fam[f]):
                probs["assigned_gene_not_a_leaf"] += 1
        info = {
            "assign_columns": list(rows[0].keys()) if rows else [],
            "assign_rows": len(rows),
            "flags": dict(collections.Counter(
                x for r in rows for x in (r["flags"] or "").split(",") if x)),
            "n_sampled": len(chosen),
            "done_json_keys": dict(keys),
            "tree_mode": dict(modes),
            "rooted": {str(k): v for k, v in rooted.items()},
            "root_children": {str(k): v for k, v in kids.items()},
            "leaves_min_med_max": [min(nleaf), _med(nleaf), max(nleaf)] if nleaf else None,
            "focal_min_med_max": [min(nfocal), _med(nfocal), max(nfocal)] if nfocal else None,
            "reference_leaf_examples": refs,
            "problems": dict(probs),
            "native_summary": {},
        }
        summ = Path(layout.results("v1", part).native().local) / "summary"
        for name in ("gene_placement.tsv", "gene_status.tsv", "skipped_families.tsv"):
            p = summ / name
            if p.exists():
                t = read_tsv(p)
                d = {"rows": len(t), "columns": list(t[0].keys()) if t else []}
                if t and "status" in t[0]:
                    d["status_counts"] = dict(collections.Counter(r["status"] for r in t))
                info["native_summary"][name] = d
            else:
                info["native_summary"][name] = "MISSING"
        report[part] = info
    text = json.dumps(report, indent=1)
    out = Path(args.out).expanduser() if args.out else default_out("survey", ".json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(text)
    print(f"\nwritten to {out}")
    return 0


def _survey2(args) -> int:
    from cichlid_data import Catalog, get_layout
    from cichlid_data.clades import Clades
    from .survey2 import survey_part

    layout = get_layout()
    cat = Catalog.load(layout)
    clades = Clades.load()
    report = {"cichlid_homology": __version__, "species": args.species}
    for part in ("TreeFam", "Compara"):
        try:
            report[part] = survey_part(layout, cat, clades, part, args.species, args.n, args.n_denovo)
        except Exception:
            import traceback
            report[part] = {"error": traceback.format_exc()}
    text = json.dumps(report, indent=1, default=str)
    out = Path(args.out).expanduser() if args.out else default_out("survey2", ".json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(text)
    print(f"\nwritten to {out}")
    return 0


def _family_survey(args) -> int:
    """Count human, teleost and other genes in every finished family of a method, or in the coordinator's prepared reference trees."""
    import csv

    from cichlid_data import Catalog, get_layout
    from cichlid_data.clades import Clades

    from . import famsurvey as F
    from .species import load_leaf_species

    layout = get_layout()
    cat = Catalog.load(layout)
    prefixes = {s.prefix for s in cat.species(release="v1")}
    res = layout.results("v1", args.part)
    root = res.family("X").local.parent
    out = Path(args.out).expanduser() if args.out else default_out("family_survey")
    out.mkdir(parents=True, exist_ok=True)
    clades = Clades.load()
    if args.source == "built":
        fams = F.finished_families(root)
        if args.limit:
            fams = fams[: args.limit]
        print(f"{args.part}: {len(fams)} finished families under {root}; loading the leaf-to-species table ...", flush=True)
        ls = load_leaf_species(layout, args.part)
        rows = F.run_survey(root, fams, ls, clades, prefixes, args.species, progress=lambda i, n: print(f"  {i}/{n}", flush=True))
        cols, summ, text = F.COLS, F.summarize(rows), None
    else:
        if args.prepared:
            prepared = Path(args.prepared).expanduser()
        else:
            fr = layout.family_resource("Compara116" if args.part == "Compara" else "TreeFam9")
            prepared = Path(getattr(fr, "local", fr)) / "prepared"
        print(f"reference trees under {prepared}; loading the leaf-to-species table ...", flush=True)
        built = F.finished_families(root)
        hit, focal = set(), set()
        for sp in cat.species(release="v1"):
            try:
                fam = {r["best_family"] for r in read_tsv(res.assign(sp.prefix).local) if r["best_family"] not in ("", "NoHit")}
            except Exception:
                continue
            hit |= fam
            if sp.prefix == args.species:
                focal = fam
        ls = load_leaf_species(layout, args.part)
        rows, diag = F.run_reference_survey(prepared, ls, clades, built, hit, focal, progress=lambda i, n: print(f"  {i}/{n}", flush=True),
                                            limit=args.limit)
        cols, summ = F.REF_COLS, F.summarize_reference(rows, diag)
        text = F.report_reference_text(summ, args.part)
    tag = f"{args.part}.{args.source}"
    with open(out / f"family_survey.{tag}.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    summ["cichlid_homology"], summ["part"], summ["source"] = __version__, args.part, args.source
    (out / f"family_survey.{tag}.json").write_text(json.dumps(summ, indent=1))
    text = text or F.report_text(summ, args.part)
    (out / f"family_survey.{tag}.md").write_text(text)
    print(text)
    print(f"\nwritten to {out}")
    return 0


def _clade_split(args) -> int:
    """Split every reference family with a human gene at its ancient duplications and score each clade (read-only)."""
    import csv

    from cichlid_data import get_layout
    from cichlid_data.clades import Clades

    from . import cladesplit as C
    from .master import load_species_tree
    from .parsimony import ParsimonyContext
    from .species import load_leaf_species

    layout = get_layout()
    clades, st = Clades.load(), load_species_tree()
    if args.prepared:
        prepared = Path(args.prepared).expanduser()
    else:
        fr = layout.family_resource("Compara116" if args.part == "Compara" else "TreeFam9")
        prepared = Path(getattr(fr, "local", fr)) / "prepared"
    print(f"reference trees under {prepared}; loading the leaf-to-species table ...", flush=True)
    ls = load_leaf_species(layout, args.part)
    fams, groups, errors, skipped = C.run_split(prepared, ls, st, clades, ParsimonyContext(st, clades), pair_overlap=args.pair_overlap,
                                                only_human=not args.all_families, limit=args.limit, already_rooted=(args.rooting == "auto"),
                                                progress=lambda i, n: print(f"  {i} families", flush=True))
    out = Path(args.out).expanduser() if args.out else default_out("clade_split")
    out.mkdir(parents=True, exist_ok=True)
    for name, cols, rows in (("families", C.FAMILY_COLS, fams), ("groups", C.GROUP_COLS, groups)):
        with open(out / f"clade_split.{args.part}.{name}.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    summ = C.summarize_split(fams, groups, errors, skipped)
    summ.update(cichlid_homology=__version__, part=args.part, pair_overlap=args.pair_overlap, rooting=args.rooting)
    (out / f"clade_split.{args.part}.json").write_text(json.dumps(summ, indent=1, default=str))
    text = C.report_split_text(summ, args.part)
    (out / f"clade_split.{args.part}.md").write_text(text)
    print(text)
    print(f"\nwritten to {out}")
    return 0


def _pilot(args) -> int:
    from cichlid_data import Catalog, get_layout
    from cichlid_data.clades import Clades
    from .dating import Params
    from .pilot import run_pilot

    layout = get_layout()
    base = Path(args.out).expanduser() if args.out else default_out("pilot")
    prm = Params(overlap_t=args.overlap, min_levels_each_side=args.min_levels,
                 mostly_one_to_one=args.mostly, k_roots=args.k_roots,
                 kept_frac=args.kept, lost_frac=args.lost, tandem_max_gap=args.tandem_gap,
                 parsimony_grid=args.parsimony_grid, gain_cost=args.gain_cost, loss_cost=args.loss_cost,
                 ohnologs_dir=args.ohnologs or "", pair_overlap=args.pair_overlap)
    species = [x.strip() for x in args.species.split(",") if x.strip()]
    alt = [] if args.all else [float(x) for x in args.alt_overlap.split(",") if x]
    if args.all:
        prm.k_roots = 1          # the rooting-stability check belongs to the old relationship columns
    cat, clades = Catalog.load(layout), Clades.load()
    for sp in species:
        out = base / sp if len(species) > 1 else base
        run_pilot(layout, cat, clades, sp, args.n, prm, alt, out,
                  n_denovo=args.n_denovo, all_genes=args.all)
        if len(species) == 1:
            print((out / "pilot_report.md").read_text())
        print(f"\n{sp}: written to {out}")
        print(f"report: {out / 'pilot_report.md'}   full summary: {out / 'pilot_summary.json'}")
    return 0


def _explain(args) -> int:
    """Print the evidence behind one gene's call in one method."""
    from cichlid_data import Catalog, get_layout
    from cichlid_data.clades import Clades
    from .calls import REFS
    from .dating import Params
    from .evidence import explain_gene
    from .master import load_species_tree
    from .pilot import prepare
    from .species import load_leaf_species

    layout = get_layout()
    prefixes = {s.prefix for s in Catalog.load(layout).species(release="v1")}
    res = layout.results("v1", args.part)
    row = next((r for r in read_tsv(res.assign(args.gene.split("_")[0]).local)
                if r["gene_id"] == args.gene), None)
    if row is None or row["best_family"] in ("", "NoHit"):
        print("this gene has no family in", args.part)
        return 1
    fam = load_family(res.family(row["best_family"]).local, prefixes, row["best_family"])
    ls = load_leaf_species(layout, args.part, only=set(fam.leaves))
    prm = Params(overlap_t=args.overlap, k_roots=1)
    prep = prepare(fam, ls, set(ls.shadowed()), load_species_tree(), prm)
    rooting, rts, dropped, _ = prep
    print(f"{args.gene} in {fam.family_id} ({args.part}, {fam.tree_mode}); rooting {rooting.how}: "
          f"{rooting.cost} duplications, {rooting.losses} losses, {rooting.n_tied} tied roots; "
          f"{len(dropped)} shadowed leaves dropped; {len(rts[0].leaves)} leaves")
    g = next(x for x in rts[0].leaves if x.name == args.gene)
    for key in ([args.ref] if args.ref else list(REFS)):
        print()
        print("\n".join(explain_gene(rts[0], g, REFS[key], g.sp, args.overlap, Clades.load())))
    return 0


def _probe(args) -> int:
    """Print what the coordinator objects look like, so the adapters can be fixed if they are wrong."""
    import cichlid_data
    from cichlid_data import get_layout

    print("cichlid_data", cichlid_data.__version__)
    print("public names:", sorted(n for n in dir(cichlid_data) if not n.startswith("_")))
    try:
        from cichlid_data.master_tree import load_master_tree
        mt = load_master_tree()
        print("master tree:", type(mt).__name__, [n for n in dir(mt) if not n.startswith("_")][:40])
        from .master import newick_of
        nw = newick_of(mt)
        print("master newick ok, starts:", nw[:80], "... length", len(nw))
    except Exception as e:
        print("master tree problem:", repr(e))
    try:
        from cichlid_data import LeafSpecies
        ls = LeafSpecies.load(get_layout(), "Compara")
        print("LeafSpecies:", [n for n in dir(ls) if not n.startswith("_")])
        print("shadowed:", ls.shadowed())
    except Exception as e:
        print("LeafSpecies problem:", repr(e))
    try:
        from cichlid_data import ReferenceNames
        t = ReferenceNames.load(get_layout(), "TreeFam")
        print("ReferenceNames:", [n for n in dir(t) if not n.startswith("_")])
        sample = next(iter(t.proteins("homo_sapiens")))
        rec = t.get(sample)
        print("  sample protein", sample, "->", type(rec).__name__, [n for n in dir(rec) if not n.startswith("_")])
        from .names import info_of
        print("  as read by cichlid-homology:", info_of(rec))
        print("  status counts (human):", t.status_counts("homo_sapiens"))
    except Exception as e:
        print("ReferenceNames problem:", repr(e))
    try:
        from cichlid_data.clades import Clades
        c = Clades.load()
        print("Clades:", [n for n in dir(c) if not n.startswith("_")])
        print("common_level(Mzeb,Drer):", c.common_level(["Mzeb", "Drer"]))
        print("lineage(lepisosteus_oculatus):", c.lineage("lepisosteus_oculatus"))
    except Exception as e:
        print("Clades problem:", repr(e))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="cichlid-homology")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("survey", help="read-only survey of sampled families (tasks 1-2)")
    s.add_argument("--species", default="Mzeb")
    s.add_argument("--n", type=int, default=300)
    s.add_argument("--out", default=None)
    s.set_defaults(func=_survey)
    s2 = sub.add_parser("survey2", help="species of leaves, shadowed leaves, coverage, support, fusions")
    s2.add_argument("--species", default="Mzeb")
    s2.add_argument("--n", type=int, default=300)
    s2.add_argument("--n-denovo", type=int, default=100)
    s2.add_argument("--out", default=None)
    s2.set_defaults(func=_survey2)
    cs = sub.add_parser("clade-split", help="split the reference families at their ancient duplications and score each clade (read-only)")
    cs.add_argument("--part", default="Compara", choices=["Compara", "TreeFam"])
    cs.add_argument("--prepared", default=None, help="the prepared/ folder (default: from the layout)")
    cs.add_argument("--rooting", default="auto", choices=["auto", "reconcile"],
                    help="auto: keep a tree's own root when it has two children; reconcile: always root by reconciliation")
    cs.add_argument("--pair-overlap", type=float, default=0.5, help="share of species two sides must have in common to be a duplication")
    cs.add_argument("--all-families", action="store_true", help="also the families without a human gene (default: only those with one)")
    cs.add_argument("--limit", type=int, default=0, help="only the first N families (a quick test)")
    cs.add_argument("--out", default=None, help="default: ~/share/clade_split_<version>_<timestamp>")
    cs.set_defaults(func=_clade_split)
    fs = sub.add_parser("family-survey", help="human, teleost and other genes in every finished family of a method (read-only)")
    fs.add_argument("--part", default="Compara", choices=["TreeFam", "Compara", "OrthoFinder"])
    fs.add_argument("--species", default="Mzeb", help="prefix of the focal species (marks the families that hold one of its genes)")
    fs.add_argument("--source", default="built", choices=["built", "reference"],
                    help="built: the finished family folders under Results; reference: the coordinator's prepared reference trees")
    fs.add_argument("--prepared", default=None, help="with --source reference: the prepared/ folder (default: from the layout)")
    fs.add_argument("--limit", type=int, default=0, help="only the first N families (a quick test)")
    fs.add_argument("--out", default=None, help="default: ~/share/family_survey_<version>_<timestamp>")
    fs.set_defaults(func=_family_survey)
    p = sub.add_parser("pilot", help="rooting, relationships, dating, TreeFam vs Compara (writes scratch files)")
    p.add_argument("--species", default="Mzeb", help="one prefix, or several separated by commas (Mzeb,Acal)")
    p.add_argument("--n", type=int, default=600, help="genes sampled (the same in both methods)")
    p.add_argument("--n-denovo", type=int, default=100, help="extra genes from de novo families, per method")
    p.add_argument("--overlap", type=float, default=0.0, help="species-overlap threshold (0 = any shared species)")
    p.add_argument("--alt-overlap", default="0.1,0.25", help="thresholds for the sensitivity table")
    p.add_argument("--min-levels", type=int, default=2)
    p.add_argument("--mostly", type=float, default=0.8)
    p.add_argument("--k-roots", type=int, default=3)
    p.add_argument("--all", action="store_true", help="every gene of the species (not a sample)")
    p.add_argument("--tandem-gap", type=int, default=10, help="genes between two copies at or below which they are tandem")
    p.add_argument("--kept", type=float, default=0.5, help="share of species with 2+ copies counted as kept")
    p.add_argument("--lost", type=float, default=0.3, help="share below which the extra copy counts as absent")
    p.add_argument("--parsimony-grid", action="store_true",
                   help="also label every gene by parsimony under a grid of gain/loss costs (for calibration; slower)")
    p.add_argument("--pair-overlap", type=float, default=0.5,
                   help="tree dating: share of species two lineages must have in common to count as a pair")
    p.add_argument("--gain-cost", type=int, default=2, help="parsimony: cost of a gained copy (a loss costs --loss-cost)")
    p.add_argument("--loss-cost", type=int, default=1)
    p.add_argument("--ohnologs", default=None,
                   help="folder of OHNOLOGS v2 3R pair files (<species>_Pairs_<Strict|Relaxed>_3R.txt); independent evidence for the tiers")
    p.add_argument("--out", default=None, help="default: ~/share/pilot_<version>_<timestamp>")
    p.set_defaults(func=_pilot)
    ex = sub.add_parser("explain", help="show the evidence behind one gene's call")
    ex.add_argument("gene")
    ex.add_argument("--part", default="TreeFam", choices=["TreeFam", "Compara"])
    ex.add_argument("--ref", default=None, choices=["human", "zebrafish"])
    ex.add_argument("--overlap", type=float, default=0.0)
    ex.set_defaults(func=_explain)
    pr = sub.add_parser("probe", help="print the shape of the coordinator objects this code relies on")
    pr.set_defaults(func=_probe)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
