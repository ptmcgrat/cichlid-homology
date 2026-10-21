#!/usr/bin/env python
"""Read-only survey of TreeFam / Compara families for one focal species.

Usage:  python survey_families.py [--species Mzeb] [--n 300] [--out survey.json]

Writes nothing except the JSON report you name. Uses only the coordinator's public API.
Leaf names are pulled with a regex (no tree library), so quoted names with commas would be miscounted.
"""
import argparse, csv, json, re, sys, random, collections
from pathlib import Path

csv.field_size_limit(sys.maxsize)
from cichlid_data import get_layout, Catalog, schemas
import cichlid_data

LEAF = re.compile(r"[(,]\s*([^(),:;\s]+)(?=\s*[:),])")


def leaves(nwk):
    return LEAF.findall(nwk)


def fasta_names(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def root_children(nwk):
    """Number of children at the root (top-level commas + 1), or None if unparseable."""
    s = nwk.strip().rstrip(";")
    if not s.startswith("("):
        return None
    depth = 0
    n = 1
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


def survey_part(layout, cat, part, species, n, rng):
    res = layout.results("v1", part)
    out = {"part": part}
    prefixes = set(p.prefix for p in cat.species(release="v1"))
    out["catalog_prefixes_n"] = len(prefixes)

    # --- run record
    try:
        out["run_json"] = json.load(open(res.run_record().local))
    except Exception as e:
        out["run_json"] = f"unreadable: {e}"

    # --- assignment table
    a = read_tsv(res.assign(species).local)
    out["assign_columns"] = list(a[0].keys()) if a else []
    out["assign_rows"] = len(a)
    out["best_family_counts_top"] = collections.Counter(r["best_family"] for r in a).most_common(5)
    out["flags"] = collections.Counter(f for r in a for f in (r["flags"] or "").split(",") if f)
    out["no_hit"] = sum(1 for r in a if r["best_family"] in ("", "NoHit"))
    out["fusion_partner_filled"] = sum(1 for r in a if r.get("fusion_partner"))
    out["fusion_partner_example"] = [
        {k: r[k] for k in r if k.startswith("fusion") or k in ("gene_id", "flags", "best_family")}
        for r in a if r.get("fusion_partner")][:3]
    fam_of_gene = {r["gene_id"]: r["best_family"] for r in a}

    # --- native per-gene tables
    summ = Path(res.native().local) / "summary"
    tables = {}
    for name in ("gene_placement.tsv", "gene_status.tsv", "skipped_families.tsv"):
        p = summ / name
        if p.exists():
            rows = read_tsv(p)
            info = {"rows": len(rows), "columns": list(rows[0].keys()) if rows else []}
            if rows and "status" in rows[0]:
                info["status_counts"] = collections.Counter(r["status"] for r in rows)
            if rows and "mode" in rows[0]:
                info["mode_counts"] = collections.Counter(r["mode"] for r in rows)
            tables[name] = info
        else:
            tables[name] = "MISSING"
    out["native_summary_tables"] = tables
    out["native_summary_listing"] = sorted(x.name for x in summ.glob("*")) if summ.exists() else "no summary dir"

    # --- families on disk
    fam_root = res.family("X").local.parent
    on_disk = sorted(d.name for d in fam_root.iterdir() if d.is_dir()) if fam_root.exists() else []
    out["families_on_disk"] = len(on_disk)
    ours = sorted(set(f for f in fam_of_gene.values() if f and f != "NoHit"))
    out["families_with_focal_genes"] = len(ours)
    out["families_with_focal_genes_and_folder"] = len([f for f in ours if f in set(on_disk)])
    out["families_with_focal_genes_no_folder"] = [f for f in ours if f not in set(on_disk)][:10]

    # --- sample
    pool = [f for f in ours if (fam_root / f / "done.json").exists()]
    sample = rng.sample(pool, min(n, len(pool)))
    modes, rooted, rooted_how = collections.Counter(), collections.Counter(), collections.Counter()
    done_keys = collections.Counter()
    sizes, nfocal, nref = [], [], []
    leaf_prefix = collections.Counter()
    ref_examples, bad = [], collections.Counter()
    root_kids = collections.Counter()
    placed_branch = []
    extra_files = collections.Counter()
    for f in sample:
        d = fam_root / f
        try:
            dj = json.load(open(d / "done.json"))
        except Exception:
            bad["done_json_unreadable"] += 1
            continue
        done_keys.update(dj.keys())
        modes[dj.get("tree_mode")] += 1
        try:
            r = schemas.family_rooting(d)
            rooted[r[0]] += 1
            rooted_how[str(r[1])] += 1
        except Exception as e:
            bad[f"family_rooting_error:{type(e).__name__}"] += 1
        extra_files.update(x.name for x in d.iterdir())
        nwk = (d / "tree.nwk").read_text()
        lv = leaves(nwk)
        root_kids[root_children(nwk)] += 1
        sizes.append(len(lv))
        if len(set(lv)) != len(lv):
            bad["duplicate_leaf"] += 1
        focal = [x for x in lv if x.split("_", 1)[0] in prefixes and "_" in x]
        ref = [x for x in lv if x not in set(focal)]
        nfocal.append(len(focal))
        nref.append(len(ref))
        for x in focal:
            leaf_prefix[x.split("_", 1)[0]] += 1
        if len(ref_examples) < 12:
            ref_examples += ref[:2]
        mine = [g for g, fam in fam_of_gene.items() if fam == f]
        if any(g not in set(lv) for g in mine):
            bad["species_gene_assigned_but_not_a_leaf"] += 1
        aln = fasta_names(d / "combined.aln")
        if set(aln) - set(lv):
            bad["aln_row_without_leaf"] += 1
        if set(lv) - set(aln):
            bad["leaf_without_aln_row"] += 1
        if (d / "epa_result.jplace").exists():
            pass
    out["sample"] = {
        "n_requested": n, "n_available_with_done_json": len(pool), "n_read": len(sample),
        "done_json_keys": dict(done_keys), "tree_mode": dict(modes),
        "rooted": {str(k): v for k, v in rooted.items()}, "rooting_text": dict(rooted_how),
        "root_children": {str(k): v for k, v in root_kids.items()},
        "leaves_per_tree": [min(sizes), sorted(sizes)[len(sizes) // 2], max(sizes)] if sizes else None,
        "focal_leaves_per_tree_min_med_max": [min(nfocal), sorted(nfocal)[len(nfocal) // 2], max(nfocal)] if nfocal else None,
        "focal_leaf_prefixes": dict(leaf_prefix),
        "reference_leaf_examples": ref_examples,
        "files_in_family_folders": dict(extra_files),
        "problems": dict(bad),
    }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="Mzeb")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--out", default="survey.json")
    args = ap.parse_args()
    layout = get_layout()
    cat = Catalog.load(layout)
    rel = cat.releases["v1"]
    report = {"coordinator_version": cichlid_data.__version__, "deliver": rel.get("deliver"),
              "parts": rel.get("parts", {}).keys() and list(rel["parts"].keys())}
    rng = random.Random(1)
    for part in ("TreeFam", "Compara"):
        try:
            report[part] = survey_part(layout, cat, part, args.species, args.n, rng)
        except Exception as e:
            import traceback
            report[part] = {"error": traceback.format_exc()}
    # master tree
    try:
        from cichlid_data.master_tree import load_master_tree, canonical
        mt = load_master_tree()
        report["master_tree"] = {"loaded": True, "type": type(mt).__name__,
                                 "canonical_examples": {n: canonical(n) for n in
                                     ["canis_familiaris", "Homo sapiens", "danio_rerio", "metriaclima_zebra"]}}
    except Exception as e:
        report["master_tree"] = {"loaded": False, "error": repr(e)}
    try:
        from cichlid_data.clades import Clades
        cl = Clades.load()
        report["clades"] = {"mzeb": cl.level("Mzeb"), "common_cichlid_zebrafish":
                            cl.common_level(["Mzeb", "Drer"]), "common_cichlid_human":
                            cl.common_level(["Mzeb", "Hsap"])}
    except Exception as e:
        report["clades"] = {"error": repr(e)}
    json.dump(report, open(args.out, "w"), indent=1, default=lambda o: dict(o) if hasattr(o, "items") else str(o))
    print(json.dumps(report, indent=1, default=lambda o: dict(o) if hasattr(o, "items") else str(o)))


if __name__ == "__main__":
    main()
