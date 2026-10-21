"""Second read-only survey: species of leaves, shadowed leaves, coverage of human / zebrafish /
gar / reedfish, placement support, fusion columns, and a targeted sample of de novo families."""
from __future__ import annotations

import collections
import random
from pathlib import Path

from .read import load_family, read_tsv, sample_families
from .species import drop_shadowed, leaf_species_map, load_leaf_species

TELEOST_LEVELS = {"teleost", "percomorph", "ovalentaria", "cichlid"}


def _mmm(v):
    v = sorted(v)
    return [v[0], v[len(v) // 2], v[-1]] if v else None


def _frac(v, k):
    return round(k / v, 4) if v else None


def placement_support(rows, species):
    """rows of gene_placement.tsv; only genes of `species`."""
    mine = [r for r in rows if r["gene_id"].startswith(species + "_")]
    out = {"genes": len(mine),
           "status": dict(collections.Counter(r["status"] for r in mine)),
           "mode": dict(collections.Counter(r["mode"] for r in mine))}
    lw = []
    for r in mine:
        if r["status"] == "placed":
            try:
                lw.append(float(r["lwr"]))
            except (TypeError, ValueError):
                out.setdefault("unparsable_lwr", 0)
                out["unparsable_lwr"] += 1
    if lw:
        out["placed_n"] = len(lw)
        out["lwr_lt_0.5"] = _frac(len(lw), sum(1 for x in lw if x < 0.5))
        out["lwr_lt_0.9"] = _frac(len(lw), sum(1 for x in lw if x < 0.9))
        out["lwr_quartiles"] = [sorted(lw)[int(q * (len(lw) - 1))] for q in (0.25, 0.5, 0.75)]
    npl = []
    for r in mine:
        try:
            npl.append(int(r["n_placements"]))
        except (TypeError, ValueError):
            pass
    out["n_placements_min_med_max"] = _mmm(npl)
    return out


def fusion_columns(rows):
    flagged = [r for r in rows if "possible_fusion" in (r["flags"] or "").split(",")]
    filled = [r for r in rows if r.get("fusion_partner")]
    cols = ["fusion_partner", "fusion_partner_score", "fusion_partner_region",
            "fusion_partner_family_coverage"]
    return {
        "flagged": len(flagged),
        "partner_filled": len(filled),
        "flagged_and_filled": sum(1 for r in flagged if r.get("fusion_partner")),
        "filled_not_flagged": sum(1 for r in filled if r not in flagged),
        "flagged_not_filled": sum(1 for r in flagged if not r.get("fusion_partner")),
        "flagged_missing_each_column": {c: sum(1 for r in flagged if not r.get(c)) for c in cols},
        "examples": [{k: r[k] for k in ["gene_id", "best_family", "flags"] + cols}
                     for r in flagged[:3]],
    }


def survey_part(layout, cat, clades, part, species, n, n_denovo, seed=1):
    prefixes = {s.prefix for s in cat.species(release="v1")}
    res = layout.results("v1", part)
    summ = Path(res.native().local) / "summary"
    gp = read_tsv(summ / "gene_placement.tsv")
    out = {"fusion": None, "placement": placement_support(gp, species)}

    root, chosen, rows = sample_families(layout, part, species, n, seed)
    out["fusion"] = fusion_columns(rows)

    # targeted sample: families built de novo that hold a gene of this species
    dn = sorted({r["family"] for r in gp
                 if r["gene_id"].startswith(species + "_") and r["mode"] == "denovo"})
    dn = [f for f in dn if (root / f / "done.json").exists()]
    out["denovo_families_with_species_gene"] = len(dn)
    dn_pick = sorted(random.Random(seed).sample(dn, min(n_denovo, len(dn))))
    fams = {f: load_family(root / f, prefixes, f) for f in sorted(set(chosen) | set(dn_pick))}

    ls = load_leaf_species(layout, part, only={x for f in fams.values() for x in f.leaves})
    shadowed = set(ls.shadowed())
    out["shadowed"] = sorted(shadowed)

    level_cache, lineage_cache = {}, {}
    no_level = set()

    def level(sp):
        if sp not in level_cache:
            try:
                level_cache[sp] = clades.level(sp)
            except Exception:
                level_cache[sp] = None
                no_level.add(sp)
        return level_cache[sp]

    def lineage(sp):
        if sp not in lineage_cache:
            try:
                lineage_cache[sp] = clades.lineage(sp)
            except Exception:
                lineage_cache[sp] = None
        return lineage_cache[sp]

    def stats(ids):
        n_sp, n_tel, n_drop, unknown = [], [], [], collections.Counter()
        has = collections.Counter()
        dropped_fams = 0
        for f in ids:
            fam = fams[f]
            sp_of = leaf_species_map(ls, fam.leaves)
            kept, dropped = drop_shadowed(sp_of, shadowed, ls.is_project_gene)
            n_drop.append(len(dropped))
            dropped_fams += bool(dropped)
            for leaf, sp in kept.items():
                if sp is None:
                    unknown[leaf] += 1
            sps = {sp for sp in kept.values() if sp}
            n_sp.append(len(sps))
            lv = {sp: level(sp) for sp in sps}
            n_tel.append(sum(1 for v in lv.values() if v in TELEOST_LEVELS))
            hum = "homo_sapiens" in sps
            zeb = "danio_rerio" in sps
            has["human"] += hum
            has["zebrafish"] += zeb
            has["human_and_zebrafish"] += hum and zeb
            has["non_teleost_actinopterygian"] += any(
                lineage(sp) == "non-teleost actinopterygian" for sp in sps)
            has["cyclostome_leaf"] += any(v == "cyclostome" for v in lv.values())
            has["amphioxus_leaf"] += "branchiostoma_floridae" in sps
        k = len(ids)
        return {
            "families": k,
            "species_per_family_min_med_max": _mmm(n_sp),
            "teleost_species_per_family_min_med_max": _mmm(n_tel),
            "shadowed_leaves_dropped_total": sum(n_drop),
            "families_with_shadowed_leaves": dropped_fams,
            "unknown_species_leaves": sum(unknown.values()),
            "unknown_species_examples": list(unknown)[:10],
            "fraction_with": {k2: _frac(k, v) for k2, v in has.items()},
        }

    out["random_sample"] = stats(chosen)
    out["denovo_sample"] = stats(dn_pick) if dn_pick else "no de novo family with a gene of this species"
    out["species_without_clade_level"] = sorted(no_level)
    out["denovo_root_children"] = dict(collections.Counter(
        fams[f].root_children for f in dn_pick))
    out["denovo_rooted"] = dict(collections.Counter(str(fams[f].rooted) for f in dn_pick))
    out["problems"] = dict(collections.Counter(p for f in fams.values() for p in f.problems))
    return out
