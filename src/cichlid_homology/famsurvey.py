"""Read-only survey of every finished family of a method (0.14.0): how many genes of human, mouse, zebrafish, gar and
reedfish, how many teleost genes and species, how many of our own placed genes, family by family.

Why: the ohnolog analysis is simplest where the human protein is one gene with clear mammalian orthologs, and a
method's families are not the same thing as the families that hold an M. zebra gene. Compara splits a family into a
fish tree and a tetrapod tree, so many families hold human genes and no teleosts, and none of them is reached from our
gene table. This counts them. Nothing is written into an input folder and nothing is built: trees are read as text."""
from __future__ import annotations

import collections
import json
from pathlib import Path

from .read import leaf_names

TELEOST_LEVELS = {"teleost", "percomorph", "ovalentaria", "cichlid"}
NAMES = {"human": "homo_sapiens", "mouse": "mus_musculus", "zebrafish": "danio_rerio",       # canonical master-tree names
         "gar_reed": {"lepisosteus_oculatus", "erpetoichthys_calabaricus"}}
COLS = ["family", "tree_mode", "n_leaves", "n_species", "n_human", "n_mouse", "n_zebrafish", "n_gar_reedfish",
        "n_teleost_ref_leaves", "n_teleost_ref_species", "n_nonteleost_gnathostome_leaves", "n_sarcopterygian_leaves",
        "n_other_leaves", "n_unknown_leaves", "n_project_leaves", "has_focal_gene"]
HIST_BINS = [("0", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3", 3, 3), ("4", 4, 4), ("5-9", 5, 9), ("10-19", 10, 19), ("20+", 20, 10 ** 9)]


def survey_family(folder, ls, clades, prefixes, focal_prefix="Mzeb", names=NAMES) -> dict:
    folder = Path(folder)
    try:
        mode = json.loads((folder / "done.json").read_text()).get("tree_mode", "")
    except Exception:
        mode = ""
    leaves = leaf_names((folder / "tree.nwk").read_text())
    row = {c: 0 for c in COLS}
    row.update(family=folder.name, tree_mode=mode, n_leaves=len(leaves), has_focal_gene=0)
    species, tel = set(), set()
    cache = {}

    def info(sp):
        if sp not in cache:
            try:
                lv = clades.level(sp)
            except Exception:
                lv = None
            try:
                lin = clades.lineage(sp)
            except Exception:
                lin = None
            cache[sp] = (lv, lin)
        return cache[sp]

    for leaf in leaves:
        if leaf.split("_")[0] in prefixes:
            row["n_project_leaves"] += 1
            if leaf.startswith(focal_prefix + "_"):
                row["has_focal_gene"] = 1
            continue
        sp = ls.species(leaf)
        if sp is None:
            row["n_unknown_leaves"] += 1
            continue
        species.add(sp)
        lv, lin = info(sp)
        if sp == names["human"]:
            row["n_human"] += 1
        elif sp == names["mouse"]:
            row["n_mouse"] += 1
        elif sp == names["zebrafish"]:
            row["n_zebrafish"] += 1
        if sp in names["gar_reed"]:
            row["n_gar_reedfish"] += 1
        if lv in TELEOST_LEVELS:
            row["n_teleost_ref_leaves"] += 1
            tel.add(sp)
        elif lv == "gnathostome":
            row["n_nonteleost_gnathostome_leaves"] += 1
        else:
            row["n_other_leaves"] += 1
        if lin == "sarcopterygian":
            row["n_sarcopterygian_leaves"] += 1
    row["n_species"], row["n_teleost_ref_species"] = len(species), len(tel)
    return row


def _hist(values) -> dict:
    out = {}
    for name, lo, hi in HIST_BINS:
        out[name] = sum(1 for v in values if lo <= v <= hi)
    return out


def _median(v):
    v = sorted(v)
    return v[len(v) // 2] if v else None


def summarize(rows) -> dict:
    n = len(rows)
    out = {"families": n, "tree_mode": dict(collections.Counter(r["tree_mode"] for r in rows))}
    human = [r for r in rows if r["n_human"] > 0]
    out["with_human"] = len(human)
    out["human_genes_per_family_histogram"] = _hist([r["n_human"] for r in rows])
    out["human_genes_per_family_histogram_families_with_human"] = _hist([r["n_human"] for r in human])
    out["exactly_one_human_gene"] = sum(1 for r in human if r["n_human"] == 1)
    out["median_human_genes_where_present"] = _median([r["n_human"] for r in human])
    no_tel = [r for r in human if r["n_teleost_ref_leaves"] == 0]
    out["with_human_and_no_teleost"] = len(no_tel)
    out["with_human_and_teleost"] = len(human) - len(no_tel)
    out["with_human_and_one_or_two_teleost_species"] = sum(1 for r in human if 1 <= r["n_teleost_ref_species"] <= 2)
    classes = [("1", 1, 1), ("2", 2, 2), ("3-4", 3, 4), ("5+", 5, 10 ** 9)]
    out["human_genes_by_teleost_presence"] = {
        nm: {"no_teleost": sum(1 for r in no_tel if lo <= r["n_human"] <= hi),
             "with_teleost": sum(1 for r in human if r["n_teleost_ref_leaves"] > 0 and lo <= r["n_human"] <= hi)}
        for nm, lo, hi in classes}
    tel = [r for r in rows if r["n_teleost_ref_leaves"] > 0]
    out["with_teleost"] = len(tel)
    out["with_teleost_and_no_human"] = sum(1 for r in tel if r["n_human"] == 0)
    out["no_human_and_no_teleost"] = sum(1 for r in rows if r["n_human"] == 0 and r["n_teleost_ref_leaves"] == 0)
    out["human_no_teleost_families"] = {
        "median_leaves": _median([r["n_leaves"] for r in no_tel]), "median_species": _median([r["n_species"] for r in no_tel]),
        "with_mouse": sum(1 for r in no_tel if r["n_mouse"] > 0),
        "with_gar_or_reedfish": sum(1 for r in no_tel if r["n_gar_reedfish"] > 0),
        "with_our_placed_genes": sum(1 for r in no_tel if r["n_project_leaves"] > 0),
        "with_a_focal_gene": sum(1 for r in no_tel if r["has_focal_gene"]),
        "median_sarcopterygian_leaves": _median([r["n_sarcopterygian_leaves"] for r in no_tel])}
    foc = [r for r in rows if r["has_focal_gene"]]
    fh = [r for r in foc if r["n_human"] > 0]
    out["families_with_a_focal_gene"] = {
        "families": len(foc), "with_human": len(fh), "without_human": len(foc) - len(fh),
        "human_genes_per_family_histogram": _hist([r["n_human"] for r in foc]),
        "exactly_one_human_gene": sum(1 for r in fh if r["n_human"] == 1)}
    return out


def report_text(s: dict, part: str) -> str:
    L = [f"# Family survey: {part}", "",
         f"{s['families']} finished families. Tree modes: {s['tree_mode']}.", "",
         "## Human genes per family", ""]
    h, hh = s["human_genes_per_family_histogram"], s["human_genes_per_family_histogram_families_with_human"]
    mx = max(h.values()) or 1
    L += ["| human genes | families | |", "|---|---|---|"]
    L += [f"| {k} | {v} | {'#' * max(0, int(40 * v / mx))} |" for k, v in h.items()]
    L += ["", f"{s['with_human']} families hold at least one human gene; exactly one in {s['exactly_one_human_gene']} of them; "
              f"median {s['median_human_genes_where_present']} where present. Among families with human: {hh}.", "",
          "## Families with human genes and no teleosts", "",
          f"{s['with_human_and_no_teleost']} of the {s['with_human']} families with human genes have no reference teleost leaf; "
          f"{s['with_human_and_teleost']} have at least one ({s['with_human_and_one_or_two_teleost_species']} of them with only "
          f"one or two teleost species). Families with a teleost and no human: {s['with_teleost_and_no_human']} of "
          f"{s['with_teleost']}; with neither: {s['no_human_and_no_teleost']}.", "",
          "| human genes | no teleost | with teleost |", "|---|---|---|"]
    L += [f"| {k} | {v['no_teleost']} | {v['with_teleost']} |" for k, v in s["human_genes_by_teleost_presence"].items()]
    L += ["", f"The families with human and no teleost: {s['human_no_teleost_families']}.", "",
          "## Families that hold a focal-species gene", "", f"{s['families_with_a_focal_gene']}", ""]
    return "\n".join(L)


def run_survey(root, family_ids, ls, clades, prefixes, focal_prefix="Mzeb", progress=None, names=NAMES) -> list:
    rows = []
    for i, fid in enumerate(family_ids, 1):
        try:
            rows.append(survey_family(Path(root) / fid, ls, clades, prefixes, focal_prefix, names))
        except Exception as e:                      # a family that cannot be read is reported, not fatal
            rows.append({**{c: 0 for c in COLS}, "family": fid, "tree_mode": f"unreadable:{type(e).__name__}"})
        if progress and i % 2000 == 0:
            progress(i, len(family_ids))
    return rows


def finished_families(root) -> list:
    root = Path(root)
    return sorted(p.name for p in root.iterdir() if p.is_dir() and (p / "done.json").exists() and (p / "tree.nwk").exists())
