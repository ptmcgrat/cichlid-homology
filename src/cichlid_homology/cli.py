"""The cichlid-homology command line."""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

from . import __version__
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
    if args.out:
        Path(args.out).expanduser().write_text(text)
    print(text)
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
    if args.out:
        Path(args.out).expanduser().write_text(text)
    print(text)
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
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
