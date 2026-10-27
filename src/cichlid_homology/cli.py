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


def _pilot(args) -> int:
    from cichlid_data import Catalog, get_layout
    from cichlid_data.clades import Clades
    from .dating import Params
    from .pilot import run_pilot

    layout = get_layout()
    out = Path(args.out).expanduser() if args.out else default_out("pilot")
    prm = Params(overlap_t=args.overlap, min_levels_each_side=args.min_levels,
                 mostly_one_to_one=args.mostly, k_roots=args.k_roots,
                 kept_frac=args.kept, lost_frac=args.lost)
    s = run_pilot(layout, Catalog.load(layout), Clades.load(), args.species, args.n, prm,
                  [float(x) for x in args.alt_overlap.split(",") if x], out, n_denovo=args.n_denovo)
    print((out / "pilot_report.md").read_text())
    print(f"\nwritten to {out}")
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
    p = sub.add_parser("pilot", help="rooting, relationships, dating, TreeFam vs Compara (writes scratch files)")
    p.add_argument("--species", default="Mzeb")
    p.add_argument("--n", type=int, default=600, help="genes sampled (the same in both methods)")
    p.add_argument("--n-denovo", type=int, default=100, help="extra genes from de novo families, per method")
    p.add_argument("--overlap", type=float, default=0.0, help="species-overlap threshold (0 = any shared species)")
    p.add_argument("--alt-overlap", default="0.1,0.25", help="thresholds for the sensitivity table")
    p.add_argument("--min-levels", type=int, default=2)
    p.add_argument("--mostly", type=float, default=0.8)
    p.add_argument("--k-roots", type=int, default=3)
    p.add_argument("--kept", type=float, default=0.5, help="share of species with 2+ copies counted as kept")
    p.add_argument("--lost", type=float, default=0.3, help="share below which the extra copy counts as absent")
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
