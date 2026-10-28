"""Split a reference family at its ancient duplications and score each clade (0.16.0).

A family with several human genes is usually several ancient paralogs (the vertebrate genome duplications), sometimes with
recent duplicates on top. The ohnolog question is simplest inside a clade where the tetrapods mostly have one gene and the
teleosts mostly have at least one. So: root the family tree by reconciliation with the master species tree, find the
duplications older than the tetrapod-fish split, and cut there; each clade left is scored on its own.

An ancient duplication is a node where (1) the species-tree node of its species equals that of one of its children (a
reconciliation duplication), (2) its two largest children share at least `pair_overlap` (0.5) of their species (the rule
that separates a duplication from a gene tree that disagrees with the species tree), and (3) species of both tetrapods
(sarcopterygians without the coelacanth) lie below it, and either reference teleosts lie below it too or the duplication is
older than the amniote split (its species-tree node is the amniote ancestor or older: a paralog lost in every fish is still
cut). Duplications inside the tetrapods younger than that (human-specific duplicates, for one) and inside the teleosts only
(the teleost genome duplication) are not cut: they stay in the clade and are counted. No other cutoff is used: the scores are
reported, not thresholded.

Rooting: with `already_rooted` the tree's own root is kept when it has exactly two children, otherwise (and always without it)
the tree is rooted by reconciliation with the master species tree, which can put the fish outside a tetrapod duplication when a
paralog is missing in fish; `rooting` in the family row says which was used.

Teleost trend (0.17): the median number of genes per teleost species present in the clade, the shares of species with one, two and three or more
genes, the same median over the early teleosts and over the percomorphs and below, the number of zebrafish genes, and a label: 1:1 (median 1), 1:2
(median 2), 1:3+ (3 or more), or mixed (a half-integer median: the species are split). The median needs no cutoff.

Each clade row lists the protein IDs of its human and zebrafish genes (`human_proteins`, `zebrafish_proteins`, `;`-joined), so a gene that
has a closest human or zebrafish gene in the family joins to its clade through that protein.

Scores of a clade: human genes; tetrapod species and genes, the mean genes per tetrapod species, the share of tetrapod species
with exactly one gene; the same for teleosts, and the share of the family's teleost (and tetrapod) species that are in the clade;
supported duplications inside it by kind; and from the branch lengths of the tree, root-to-tip distances (from the clade's
root) of the tetrapods: median, spread (coefficient of variation) and the human distance relative to the other tetrapods."""
from __future__ import annotations

import collections
import csv
import json
import statistics
from pathlib import Path

from .famsurvey import NAMES, _median, reference_files
from .orthology import RTree
from .rooting import root_tree
from .tree import leaves, parse_newick, postorder, prune
from .treedate import _bits, _Tree

GROUP_COLS = ["family", "group", "n_human", "human_proteins", "zebrafish_proteins", "n_genes", "n_species", "n_tetrapod_genes", "n_tetrapod_species", "tetrapod_mean_genes",
              "tetrapod_single_copy_share", "n_teleost_genes", "n_teleost_species", "teleost_mean_genes", "teleost_median_genes", "teleost_share_one", "teleost_share_two",
              "teleost_share_three_plus", "early_teleost_median_genes", "percomorph_median_genes", "zebrafish_genes", "teleost_trend",
              "teleost_species_share_of_family",
              "tetrapod_species_share_of_family", "n_gar_reedfish", "tetrapod_clade_dups", "tetrapod_species_dups", "teleost_clade_dups",
              "mixed_dups", "tetrapod_rtt_median", "tetrapod_rtt_cv", "human_rtt", "human_rel_rate", "lengths_missing"]
FAMILY_COLS = ["family", "n_genes", "n_human", "n_groups", "n_groups_with_human", "n_ancient_dups", "n_family_tetrapod_species",
               "n_family_teleost_species", "root_cost", "rooting"]


def trend_label(median) -> str:
    """The teleost trend of a clade from the median number of genes per teleost species present (no cutoff: the median is an integer or
    a half): 1 is 1:1, 2 is 1:2, 3 or more is 1:3+, and in between the species are split."""
    if median is None or median == "":
        return ""
    m = float(median)
    if m <= 1.0:
        return "1:1"
    if m < 2.0:
        return "mixed 1-2"
    if m == 2.0:
        return "1:2"
    if m < 3.0:
        return "mixed 2-3"
    return "1:3+"


def _med(v):
    return statistics.median(v) if v else ""


def _overlap(T, v) -> float:
    inf = [c for c in v.children if T.M[id(c)] is not None]
    if len(inf) < 2:
        return 0.0
    big = sorted(inf, key=lambda c: -T.n_species(c))[:2]
    union = _bits(T.smask[id(big[0])] | T.smask[id(big[1])])
    return _bits(T.smask[id(big[0])] & T.smask[id(big[1])]) / union if union else 0.0


def _rtt(g, T):
    """Root-to-tip distances (from the clade's root, not counting the branch above it) to the informative leaves: {leaf: distance}; and
    whether any branch length was missing."""
    out, missing = {}, False
    if not g.children:
        return ({g: 0.0} if T.M[id(g)] is not None else {}), False
    stack = [(c, 0.0) for c in g.children]
    while stack:
        n, d = stack.pop()
        if n.length is None:
            missing = True
        d += n.length or 0.0
        if n.children:
            stack.extend((c, d) for c in n.children)
        elif T.M[id(n)] is not None:
            out[n] = d
    return out, missing


def _stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return None, None
    med = statistics.median(vals)
    cv = statistics.pstdev(vals) / (sum(vals) / len(vals)) if len(vals) > 1 and sum(vals) > 0 else (0.0 if len(vals) > 1 else None)
    return med, cv


def split_family(fid, newick, ls, st, clades, ctx, names=NAMES, pair_overlap=0.5, already_rooted=True):
    """(family row, [group rows]) for one reference family tree, or None when nothing is left of it."""
    root = parse_newick(newick)
    sp_of = {x.name: ls.species(x.name) for x in leaves(root)}
    unknown = {n for n, s in sp_of.items() if s is None}
    if unknown:
        root = prune(root, unknown)
    if root is None:
        return None
    for x in leaves(root):
        x.sp = sp_of[x.name]
    nleaf = len(leaves(root))
    if nleaf < 2:                                                      # a one-member family: nothing to root or split
        rooting, rt = None, None
    else:
        rooting = root_tree(root, st, already_rooted, 1)
        rt = RTree(rooting.roots[0])
    if rt is None:
        sp = root.sp
        row = {c: 0 for c in GROUP_COLS}
        row.update(family=fid, group="G1", n_genes=1, n_species=1, n_human=int(sp == names["human"]), lengths_missing=0,
                   human_proteins=root.name if sp == names["human"] else "", zebrafish_proteins=root.name if sp == names["zebrafish"] else "")
        return ({"family": fid, "n_genes": 1, "n_human": row["n_human"], "n_groups": 1, "n_groups_with_human": row["n_human"],
                 "n_ancient_dups": 0, "n_family_tetrapod_species": 0, "n_family_teleost_species": 0, "root_cost": "", "rooting": "single"}, [row])
    T = _Tree(rt, ctx, "denovo", None)
    tet = {sp for sp in T.rank if clades.lineage(sp) == "sarcopterygian" and sp != names.get("coelacanth")}
    tel = {sp for sp, r in T.rank.items() if r >= 1}
    tet_mask = sum(1 << rt.bit[sp] for sp in tet if sp in rt.bit)
    tel_mask = sum(1 << rt.bit[sp] for sp in tel if sp in rt.bit)
    gar = names["gar_reed"]
    amn = None
    if names.get("chicken") in st.leaf_id and names["human"] in st.leaf_id:
        amn = st.lca(st.id_of(names["human"]), st.id_of(names["chicken"]))            # the amniote ancestor

    def deep(v):                                                      # its species-tree node is the amniote ancestor or older
        return amn is not None and T.M[id(v)] is not None and st.lca(T.M[id(v)], amn) == T.M[id(v)]

    def supported(v):
        return bool(v.children) and T.M[id(v)] is not None and T.dup(v) and _overlap(T, v) >= pair_overlap

    ancient, has = {}, {}
    n_anc = 0
    for v in postorder(rt.root):
        a = supported(v) and bool(T.smask[id(v)] & tet_mask) and (bool(T.smask[id(v)] & tel_mask) or deep(v))
        n_anc += a
        ancient[id(v)] = a
        has[id(v)] = a or any(has[id(c)] for c in v.children)
    groups, stack = [], [rt.root]
    while stack:
        v = stack.pop()
        if not has[id(v)]:
            groups.append(v)
            continue
        stack.extend(reversed(v.children))
    fam_tet = len({x.sp for x in rt.leaves if x.sp in tet and T.M[id(x)] is not None})
    fam_tel = len({x.sp for x in rt.leaves if x.sp in tel and T.M[id(x)] is not None})
    rows = []
    for g in groups:
        lv = [x for x in rt.leaves_in(g) if T.M[id(x)] is not None]
        if not lv:
            continue                                                  # a clade with no jawed vertebrate (an outgroup)
        per = collections.Counter(x.sp for x in lv)
        tet_sp = [s for s in per if s in tet]
        tel_sp = [s for s in per if s in tel]
        n_tet, n_tel = sum(per[s] for s in tet_sp), sum(per[s] for s in tel_sp)
        tcounts = [per[s] for s in tel_sp]
        tmed = _med(tcounts)
        early_c = [per[s] for s in tel_sp if T.rank.get(s) == 1]
        perc_c = [per[s] for s in tel_sp if T.rank.get(s, 0) >= 2]
        dists, missing = _rtt(g, T)
        tet_d = [d for x, d in dists.items() if x.sp in tet]
        hum_d = [d for x, d in dists.items() if x.sp == names["human"]]
        oth_d = [d for x, d in dists.items() if x.sp in tet and x.sp != names["human"]]
        med, cv = _stats(tet_d)
        omed = statistics.median(oth_d) if len(oth_d) >= 3 else None
        hd = sum(hum_d) / len(hum_d) if hum_d else None
        kinds = collections.Counter()
        if g.children:
            for v in postorder(g):
                if v.children and supported(v):
                    m = T.smask[id(v)]
                    if m & ~tet_mask == 0:
                        kinds["tet_species" if _bits(m) == 1 else "tet_clade"] += 1
                    elif m & ~tel_mask == 0:
                        kinds["tel_clade" if _bits(m) > 1 else "tel_species"] += 1
                    else:
                        kinds["mixed"] += 1
        rows.append({"family": fid, "group": f"G{len(rows) + 1}", "n_human": per.get(names["human"], 0),
                     "human_proteins": ";".join(sorted(x.name for x in lv if x.sp == names["human"])),
                     "zebrafish_proteins": ";".join(sorted(x.name for x in lv if x.sp == names["zebrafish"])),
                     "n_genes": len(lv), "n_species": len(per),
                     "n_tetrapod_genes": n_tet, "n_tetrapod_species": len(tet_sp),
                     "tetrapod_mean_genes": round(n_tet / len(tet_sp), 3) if tet_sp else "",
                     "tetrapod_single_copy_share": round(sum(1 for s in tet_sp if per[s] == 1) / len(tet_sp), 3) if tet_sp else "",
                     "n_teleost_genes": n_tel, "n_teleost_species": len(tel_sp),
                     "teleost_mean_genes": round(n_tel / len(tel_sp), 3) if tel_sp else "",
                     "teleost_median_genes": tmed,
                     "teleost_share_one": round(sum(1 for c in tcounts if c == 1) / len(tcounts), 3) if tcounts else "",
                     "teleost_share_two": round(sum(1 for c in tcounts if c == 2) / len(tcounts), 3) if tcounts else "",
                     "teleost_share_three_plus": round(sum(1 for c in tcounts if c >= 3) / len(tcounts), 3) if tcounts else "",
                     "early_teleost_median_genes": _med(early_c), "percomorph_median_genes": _med(perc_c),
                     "zebrafish_genes": per.get(names["zebrafish"], 0), "teleost_trend": trend_label(tmed),
                     "teleost_species_share_of_family": round(len(tel_sp) / fam_tel, 3) if fam_tel else "",
                     "tetrapod_species_share_of_family": round(len(tet_sp) / fam_tet, 3) if fam_tet else "",
                     "n_gar_reedfish": sum(per[s] for s in per if s in gar), "tetrapod_clade_dups": kinds["tet_clade"],
                     "tetrapod_species_dups": kinds["tet_species"], "teleost_clade_dups": kinds["tel_clade"], "mixed_dups": kinds["mixed"],
                     "tetrapod_rtt_median": round(med, 4) if med is not None else "", "tetrapod_rtt_cv": round(cv, 3) if cv is not None else "",
                     "human_rtt": round(hd, 4) if hd is not None else "",
                     "human_rel_rate": round(hd / omed, 3) if (hd is not None and omed) else "", "lengths_missing": int(missing)})
    fam = {"family": fid, "n_genes": sum(1 for x in rt.leaves if T.M[id(x)] is not None),
           "n_human": sum(1 for x in rt.leaves if x.sp == names["human"]), "n_groups": len(rows),
           "n_groups_with_human": sum(1 for r in rows if r["n_human"] > 0), "n_ancient_dups": n_anc,
           "n_family_tetrapod_species": fam_tet, "n_family_teleost_species": fam_tel, "root_cost": rooting.cost, "rooting": rooting.how}
    return fam, rows


def run_split(prepared, ls, st, clades, ctx, names=NAMES, pair_overlap=0.5, only_human=True, limit=0, progress=None, already_rooted=True):
    prepared = Path(prepared)
    files, _ = reference_files(prepared / "trees")
    ids = sorted(f for f, d in files.items() if "tree" in d)
    fams, groups, errors, skipped = [], [], [], 0
    done = 0
    for fid in ids:
        text = files[fid]["tree"].read_text(errors="replace")
        if only_human:
            from .read import leaf_names
            import re
            lv = leaf_names(re.sub(r"\[[^\]]*\]", "", text)) or [text.strip().split(":")[0].split("[")[0].rstrip(";")]
            if not any(ls.species(x) == names["human"] for x in lv):
                skipped += 1
                continue
        try:
            r = split_family(fid, text, ls, st, clades, ctx, names, pair_overlap, already_rooted)
        except Exception as e:
            errors.append((fid, f"{type(e).__name__}: {e}"[:160]))
            continue
        if r:
            fams.append(r[0])
            groups.extend(r[1])
        done += 1
        if progress and done % 1000 == 0:
            progress(done, len(ids))
        if limit and done >= limit:
            break
    return fams, groups, errors, skipped


TET_BINS = [("1", None), ("1 to 1.25", 1.25), ("1.25 to 1.5", 1.5), ("1.5 to 2", 2.0), ("over 2", float("inf"))]
SHARE_BINS = [("1 (every tetrapod species has one gene)", 1.0), ("0.9 to 1", 0.9), ("0.7 to 0.9", 0.7), ("0.5 to 0.7", 0.5), ("under 0.5", -1)]


def _bin_mean(x):
    if x <= 1.0 + 1e-9:
        return TET_BINS[0][0]
    for n, hi in TET_BINS[1:]:
        if x <= hi:
            return n


def _bin_share(x):
    if x >= 1.0 - 1e-9:
        return SHARE_BINS[0][0]
    if x >= 0.9:
        return SHARE_BINS[1][0]
    if x >= 0.7:
        return SHARE_BINS[2][0]
    if x >= 0.5:
        return SHARE_BINS[3][0]
    return SHARE_BINS[4][0]


def _q(v):
    v = sorted(x for x in v if x != "")
    if not v:
        return None
    pick = lambda f: round(v[min(len(v) - 1, int(f * (len(v) - 1) + 0.5))], 3)
    return {"p10": pick(0.1), "p25": pick(0.25), "median": pick(0.5), "p75": pick(0.75), "p90": pick(0.9)}


def _profile(gs) -> dict:
    tm = [g["tetrapod_mean_genes"] for g in gs if g["tetrapod_mean_genes"] != ""]
    ts = [g["tetrapod_single_copy_share"] for g in gs if g["tetrapod_single_copy_share"] != ""]
    teo = [g for g in gs if g["n_teleost_species"] > 0]
    return {"groups": len(gs), "with_tetrapods": len(tm),
            "tetrapod_mean_genes_bins": {n: sum(1 for x in tm if _bin_mean(x) == n) for n, _ in TET_BINS},
            "tetrapod_single_copy_share_bins": {n: sum(1 for x in ts if _bin_share(x) == n) for n, _ in SHARE_BINS},
            "with_teleosts": len(teo), "without_teleosts": len(gs) - len(teo),
            "teleost_mean_genes_quantiles": _q([g["teleost_mean_genes"] for g in teo]),
            "teleost_species_share_of_family_quantiles": _q([g["teleost_species_share_of_family"] for g in teo]),
            "tetrapod_species_median": _median([g["n_tetrapod_species"] for g in gs if g["n_tetrapod_species"] > 0]),
            "with_a_tetrapod_clade_duplication": sum(1 for g in gs if g["tetrapod_clade_dups"] > 0),
            "human_rel_rate_quantiles": _q([g["human_rel_rate"] for g in gs]),
            "tetrapod_rtt_cv_quantiles": _q([g["tetrapod_rtt_cv"] for g in gs])}


TRENDS = ("1:1", "mixed 1-2", "1:2", "mixed 2-3", "1:3+")


def _trend_table(groups) -> dict:
    """For clades with exactly one (and two) human genes and teleosts: the teleost trend (median genes per teleost species) against the
    number of zebrafish genes in the clade, and the early-teleost and percomorph medians for the 1:1 and 1:2 clades."""
    out = {}
    for k in (1, 2):
        sub = [g for g in groups if g["n_human"] == k and g["n_teleost_species"] > 0]
        d = {"clades_with_teleosts": len(sub), "trends": {}}
        for tr in TRENDS:
            gs = [g for g in sub if g["teleost_trend"] == tr]
            d["trends"][tr] = {"clades": len(gs), "zebrafish_genes": {str(z) if z < 3 else "3+": sum(1 for g in gs if min(g["zebrafish_genes"], 3) == z) for z in (0, 1, 2, 3)},
                               "median_teleost_species": _median([g["n_teleost_species"] for g in gs]),
                               "median_share_one": _q([g["teleost_share_one"] for g in gs]), "median_share_two": _q([g["teleost_share_two"] for g in gs]),
                               "early_vs_percomorph_median": dict(collections.Counter(f"{g['early_teleost_median_genes']} / {g['percomorph_median_genes']}" for g in gs).most_common(6))}
        out[str(k)] = d
    return out


def _kbin(k):
    return str(k) if k <= 4 else ("5-9" if k <= 9 else "10+")


def summarize_split(fams, groups, errors, skipped) -> dict:
    out = {"families": len(fams), "skipped_without_human": skipped, "errors": len(errors), "error_examples": errors[:5]}
    by_k = collections.defaultdict(collections.Counter)
    for f in fams:
        by_k[_kbin(f["n_human"])][min(f["n_groups"], 5)] += 1
    out["groups_per_family_by_human_genes"] = {k: {("5+" if g == 5 else str(g)): n for g, n in sorted(c.items())}
                                               for k, c in sorted(by_k.items(), key=lambda kv: (len(kv[0]), kv[0]))}
    out["rooting_used"] = dict(collections.Counter(f["rooting"] for f in fams))
    out["groups"] = len(groups)
    out["groups_by_human_genes"] = dict(sorted(collections.Counter(min(g["n_human"], 3) for g in groups).items()))
    fam_k = {f["family"]: f["n_human"] for f in fams}
    one = [g for g in groups if g["n_human"] == 1]
    out["groups_with_one_human_gene"] = {
        "all": _profile(one),
        "from_families_with_one_human_gene": _profile([g for g in one if fam_k[g["family"]] == 1]),
        "from_families_with_several_human_genes": _profile([g for g in one if fam_k[g["family"]] > 1])}
    out["groups_with_two_human_genes"] = _profile([g for g in groups if g["n_human"] == 2])
    out["teleost_trend"] = _trend_table(groups)
    four = [f["family"] for f in fams if f["n_human"] == 4][:8]
    out["examples_families_with_four_human_genes"] = {fid: [{k: g[k] for k in ("group", "n_human", "n_tetrapod_species", "tetrapod_mean_genes",
                                                                        "tetrapod_single_copy_share", "n_teleost_species", "teleost_mean_genes",
                                                                        "teleost_species_share_of_family", "tetrapod_clade_dups", "human_rel_rate")}
                                                          for g in groups if g["family"] == fid] for fid in four}
    return out


def _profile_text(name, p) -> list:
    L = [f"### {name}: {p['groups']} clades", ""]
    if not p["groups"]:
        return L + ["(none)", ""]
    m, s = p["tetrapod_mean_genes_bins"], p["tetrapod_single_copy_share_bins"]
    L += ["| tetrapods | clades | share of tetrapod species with exactly one gene | clades |", "|---|---|---|---|"]
    for (a, x), (b, y) in zip(m.items(), s.items()):
        L.append(f"| mean {a} | {x} | {b} | {y} |")
    L += ["", f"With teleosts: {p['with_teleosts']}; without: {p['without_teleosts']}. Median tetrapod species: {p['tetrapod_species_median']}. "
              f"Teleost mean genes per species {p['teleost_mean_genes_quantiles']}; share of the family's teleost species in the clade "
              f"{p['teleost_species_share_of_family_quantiles']}. Clades with a duplication inside the tetrapods: {p['with_a_tetrapod_clade_duplication']}. "
              f"Human root-to-tip distance relative to other tetrapods {p['human_rel_rate_quantiles']}; spread (CV) of tetrapod distances "
              f"{p['tetrapod_rtt_cv_quantiles']}.", ""]
    return L


def report_split_text(s: dict, part: str) -> str:
    L = [f"# Clades of the reference families, split at their ancient duplications: {part}", "",
         f"{s['families']} families with a human gene split ({s['skipped_without_human']} without human skipped; {s['errors']} errors "
         f"{s['error_examples']}). Rooting used: {s['rooting_used']}. {s['groups']} clades; by number of human genes (3 means 3 or more): {s['groups_by_human_genes']}.", "",
         "## Clades per family, by number of human genes in the family", "", "| human genes | clades per family |", "|---|---|"]
    L += [f"| {k} | {v} |" for k, v in s["groups_per_family_by_human_genes"].items()]
    L += ["", "## Clades with exactly one human gene", ""]
    one = s["groups_with_one_human_gene"]
    L += _profile_text("all", one["all"]) + _profile_text("from families with one human gene", one["from_families_with_one_human_gene"])
    L += _profile_text("from families with several human genes (the split ones)", one["from_families_with_several_human_genes"])
    L += ["## Clades with two human genes", ""] + _profile_text("all", s["groups_with_two_human_genes"])
    L += ["## Teleost trend: the median number of genes per teleost species in the clade", "",
          "1:1 = median 1, 1:2 = median 2, 1:3+ = median 3 or more; the mixed classes have a median of 1.5 or 2.5 (the species are split). "
          "Early teleosts and percomorphs are the two halves of the teleost species.", ""]
    for k, d in s["teleost_trend"].items():
        L += [f"### Clades with {k} human gene{'s' if k != '1' else ''} and teleosts: {d['clades_with_teleosts']}", "",
              "| trend | clades | zebrafish 0 | 1 | 2 | 3+ | median teleost species | early / percomorph medians (most common) |", "|---|---|---|---|---|---|---|---|"]
        for tr, v in d["trends"].items():
            z = v["zebrafish_genes"]
            L.append(f"| {tr} | {v['clades']} | {z['0']} | {z['1']} | {z['2']} | {z['3+']} | {v['median_teleost_species']} | {v['early_vs_percomorph_median']} |")
        L += [""]
    L += ["## Examples: families with four human genes", ""]
    for fid, gs in s["examples_families_with_four_human_genes"].items():
        L += [f"{fid}: " + " | ".join(
            f"{g['group']}: {g['n_human']} human, tetrapod {g['n_tetrapod_species']} sp mean {g['tetrapod_mean_genes']} single {g['tetrapod_single_copy_share']}, "
            f"teleost {g['n_teleost_species']} sp mean {g['teleost_mean_genes']} ({g['teleost_species_share_of_family']} of family)" for g in gs), ""]
    return "\n".join(L)
