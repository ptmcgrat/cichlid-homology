"""Read-only survey of every finished family of a method (0.14.0): how many genes of human, mouse, zebrafish, gar and
reedfish, how many teleost genes and species, how many of our own placed genes, family by family.

Why: the ohnolog analysis is simplest where the human protein is one gene with clear mammalian orthologs, and a
method's families are not the same thing as the families that hold an M. zebra gene. Compara splits a family into a
fish tree and a tetrapod tree, so many families hold human genes and no teleosts, and none of them is reached from our
gene table. This counts them. Nothing is written into an input folder and nothing is built: trees are read as text."""
from __future__ import annotations

import collections
import csv
import json
import re
from pathlib import Path

from .read import leaf_names

TELEOST_LEVELS = {"teleost", "percomorph", "ovalentaria", "cichlid"}
NAMES = {"human": "homo_sapiens", "mouse": "mus_musculus", "zebrafish": "danio_rerio",       # canonical master-tree names
         "gar_reed": {"lepisosteus_oculatus", "erpetoichthys_calabaricus"}, "coelacanth": "latimeria_chalumnae"}
COLS = ["family", "tree_mode", "n_leaves", "n_species", "n_human", "n_mouse", "n_zebrafish", "n_gar_reedfish",
        "n_teleost_ref_leaves", "n_teleost_ref_species", "n_nonteleost_gnathostome_leaves", "n_sarcopterygian_leaves",
        "n_other_leaves", "n_unknown_leaves", "n_project_leaves", "has_focal_gene", "n_tetrapod_leaves", "n_tetrapod_species"]
HIST_BINS = [("0", 0, 0), ("1", 1, 1), ("2", 2, 2), ("3", 3, 3), ("4", 4, 4), ("5-9", 5, 9), ("10-19", 10, 19), ("20+", 20, 10 ** 9)]


def survey_family(folder, ls, clades, prefixes, focal_prefix="Mzeb", names=NAMES) -> dict:
    folder = Path(folder)
    try:
        mode = json.loads((folder / "done.json").read_text()).get("tree_mode", "")
    except Exception:
        mode = ""
    return survey_tree(folder.name, (folder / "tree.nwk").read_text(), mode, ls, clades, prefixes, focal_prefix, names)


def survey_tree(name, nwk, mode, ls, clades, prefixes, focal_prefix="Mzeb", names=NAMES) -> dict:
    text = re.sub(r"\[[^\]]*\]", "", nwk)                          # NHX or other bracketed comments are dropped
    leaves = leaf_names(text)
    if not leaves:                                                  # a one-member family is a bare name: "Hsap_1;" or "Hsap_1:0.0;"
        m = _BARE.match(text)
        leaves = [m.group(1)] if m else []
    return survey_leaves(name, leaves, mode, ls, clades, prefixes, focal_prefix, names)


def survey_leaves(name, leaves, mode, ls, clades, prefixes, focal_prefix="Mzeb", names=NAMES) -> dict:
    row = {c: 0 for c in COLS}
    row.update(family=name, tree_mode=mode, n_leaves=len(leaves), has_focal_gene=0)
    species, tel, tet = set(), set(), set()
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
            if sp != names.get("coelacanth"):                     # tetrapods: the sarcopterygians without the coelacanth
                row["n_tetrapod_leaves"] += 1
                tet.add(sp)
    row["n_species"], row["n_teleost_ref_species"], row["n_tetrapod_species"] = len(species), len(tel), len(tet)
    return row


def _hist(values) -> dict:
    out = {}
    for name, lo, hi in HIST_BINS:
        out[name] = sum(1 for v in values if lo <= v <= hi)
    return out


def _median(v):
    v = sorted(v)
    return v[len(v) // 2] if v else None


MEAN_BINS = [("1 (one gene per species)", None), ("1 to 1.25", 1.25), ("1.25 to 1.5", 1.5), ("1.5 to 2", 2.0),
             ("2 to 3", 3.0), ("3 to 5", 5.0), ("over 5", float("inf"))]


def _mean_bin(x: float) -> str:
    if x <= 1.0 + 1e-9:
        return MEAN_BINS[0][0]
    for name, hi in MEAN_BINS[1:]:
        if x <= hi:
            return name
    return MEAN_BINS[-1][0]


def _quantiles(v):
    v = sorted(v)
    if not v:
        return None
    q = lambda f: round(v[min(len(v) - 1, int(f * (len(v) - 1) + 0.5))], 3)
    return {"min": round(v[0], 3), "p10": q(0.1), "p25": q(0.25), "median": q(0.5), "p75": q(0.75), "p90": q(0.9), "max": round(v[-1], 3)}


def copies_per_species(rows, human_genes=(1, 2)) -> dict:
    """For families with exactly k human genes: the distribution of the mean number of genes per species among the tetrapod species
    present, and among the reference teleost species present (genes / species per family)."""
    out = {}
    for k in human_genes:
        sub = [r for r in rows if r["n_human"] == k]
        d = {"families": len(sub)}
        for grp, lv, sp in (("tetrapod", "n_tetrapod_leaves", "n_tetrapod_species"), ("teleost", "n_teleost_ref_leaves", "n_teleost_ref_species")):
            vals = [r[lv] / r[sp] for r in sub if r[sp] > 0]
            bins = collections.Counter(_mean_bin(x) for x in vals)
            d[grp] = {"families_with_species": len(vals), "families_with_none": len(sub) - len(vals),
                      "mean_of_the_means": round(sum(vals) / len(vals), 3) if vals else None, "quantiles": _quantiles(vals),
                      "bins": {name: bins.get(name, 0) for name, _ in MEAN_BINS},
                      "median_species_present": _median([r[sp] for r in sub if r[sp] > 0])}
        out[str(k)] = d
    return out


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
    out["mean_genes_per_species_by_human_genes"] = copies_per_species(rows)
    out["families_with_a_focal_gene"] = {
        "families": len(foc), "with_human": len(fh), "without_human": len(foc) - len(fh),
        "human_genes_per_family_histogram": _hist([r["n_human"] for r in foc]),
        "exactly_one_human_gene": sum(1 for r in fh if r["n_human"] == 1)}
    return out



def _mean_section(s) -> list:
    L = ["", "## Mean number of genes per species, by number of human genes in the family", "",
         "For each family: genes divided by species, over the tetrapod species present (the coelacanth is left out) and over the reference "
         "teleost species present. 1 means one gene per species.", ""]
    for k, d in s["mean_genes_per_species_by_human_genes"].items():
        t, e = d["tetrapod"], d["teleost"]
        L += [f"### Families with {k} human gene{'s' if k != '1' else ''}: {d['families']} families", "",
              "| mean genes per species | tetrapods | teleosts |", "|---|---|---|"]
        for name in t["bins"]:
            a, b = t["bins"][name], e["bins"][name]
            L.append(f"| {name} | {a} ({100 * a / max(1, t['families_with_species']):.0f}%) | {b} ({100 * b / max(1, e['families_with_species']):.0f}%) |")
        L += [f"| no species of the group in the family | {t['families_with_none']} | {e['families_with_none']} |", "",
              f"Quartiles (min, p10, p25, median, p75, p90, max): tetrapods {t['quantiles']}; teleosts {e['quantiles']}. Mean of the means: "
              f"tetrapods {t['mean_of_the_means']}, teleosts {e['mean_of_the_means']}. Median number of species present: tetrapods "
              f"{t['median_species_present']}, teleosts {e['median_species_present']}.", ""]
    return L


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
    return "\n".join(L + _mean_section(s))


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


REF_EXTRA = ["audit_status", "built", "hit_by_project", "leaf_source"]
REF_COLS = COLS + REF_EXTRA
TREE_EXT = {".nhx", ".nwk", ".newick", ".tree", ".tre"}
ALN_EXT = {".aln", ".fa", ".fasta", ".afa", ".sto", ".stk"}


def reference_files(trees_dir):
    """({family id: {"tree": path, "aln": path}}, Counter of every file extension). The id is the name up to the first dot."""
    by, ext = {}, collections.Counter()
    for p in sorted(Path(trees_dir).iterdir()):
        if not p.is_file():
            continue
        e = p.suffix.lower()
        ext[e or "(none)"] += 1
        d = by.setdefault(p.name.split(".")[0], {})
        if e in TREE_EXT:
            d["tree"] = p
        elif e in ALN_EXT:
            d["aln"] = p
    return by, ext


_RANGE = re.compile(r"/\d+-\d+$")
_BARE = re.compile(r"^\s*([^(),:;\s]+)\s*(?::[-+0-9.eE]+)?\s*;?\s*$")


def alignment_names(path) -> list:
    """Sequence names of an alignment: FASTA headers, or the first column of Stockholm lines; a trailing /start-end is removed."""
    text = Path(path).read_text(errors="replace")
    out, seen = [], set()
    if text.lstrip().startswith("# STOCKHOLM"):
        for line in text.splitlines():
            if line and not line.startswith(("#", "//")):
                n = _RANGE.sub("", line.split()[0])
                if n not in seen:
                    seen.add(n)
                    out.append(n)
        return out
    for line in text.splitlines():
        if line.startswith(">"):
            out.append(_RANGE.sub("", line[1:].split()[0]) if line[1:].strip() else "")
    return [n for n in out if n]


def read_audit(path):
    """({family: status}, header) from family_audit.tsv: the columns named family and status, else the first two."""
    path = Path(path)
    if not path.is_file():
        return {}, []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        head = list(rd.fieldnames or [])
        fam = "family" if "family" in head else (head[0] if head else None)
        st = "status" if "status" in head else (head[1] if len(head) > 1 else None)
        return {r[fam]: (r[st] if st else "") for r in rd if fam and r.get(fam)}, head


def run_reference_survey(prepared, ls, clades, built=(), hit=(), focal_hit=(), audit=None, progress=None,
                         names=NAMES, limit=0) -> tuple:
    """Survey the coordinator's prepared reference families: the tree when there is one, else the names in the alignment.
    `built`: the families with a result folder; `hit`: the families any project species was assigned to; `focal_hit`: those the
    focal species was assigned to. Returns (rows, diagnostics)."""
    prepared = Path(prepared)
    files, ext = reference_files(prepared / "trees")
    ids = sorted(files)[:limit] if limit else sorted(files)
    audit = audit if audit is not None else read_audit(prepared / "family_audit.tsv")[0]
    built, hit, focal_hit = set(built), set(hit), set(focal_hit)
    rows = []
    for i, fid in enumerate(ids, 1):
        f = files[fid]
        try:
            if "tree" in f:
                r = survey_tree(fid, f["tree"].read_text(errors="replace"), "reference", ls, clades, set(), "", names)
                src = "tree"
                if r["n_leaves"] == 0 and "aln" in f:               # an empty or unparsable tree: the alignment names the members
                    r = survey_leaves(fid, alignment_names(f["aln"]), "reference", ls, clades, set(), "", names)
                    src = "alignment_because_tree_empty"
            elif "aln" in f:
                r = survey_leaves(fid, alignment_names(f["aln"]), "reference", ls, clades, set(), "", names)
                src = "alignment"
            else:
                raise ValueError("no tree or alignment file")
        except Exception as e:
            r, src = {**{c: 0 for c in COLS}, "family": fid, "tree_mode": f"unreadable:{type(e).__name__}"}, "none"
        r.update(audit_status=audit.get(fid, "not_in_audit"), built=int(fid in built), hit_by_project=int(fid in hit),
                 has_focal_gene=int(fid in focal_hit), leaf_source=src)
        rows.append(r)
        if progress and i % 5000 == 0:
            progress(i, len(ids))
    diag = {"families_with_files": len(files), "files_by_extension": dict(ext),
            "families_with_a_tree_file": sum(1 for f in files.values() if "tree" in f),
            "families_with_an_alignment_only": sum(1 for f in files.values() if "tree" not in f and "aln" in f),
            "families_read_from_the_alignment_because_the_tree_was_empty": sum(1 for r in rows if r["leaf_source"] == "alignment_because_tree_empty"),
            "first_files": [str(p.name) for f in list(files.values())[:3] for p in f.values()],
            "audit_rows": len(audit), "audit_statuses": dict(collections.Counter(audit.values())),
            "audited_families_without_a_file": sum(1 for f in audit if f not in files),
            "unknown_leaves_total": sum(r["n_unknown_leaves"] for r in rows),
            "leaves_total": sum(r["n_leaves"] for r in rows),
            "unreadable": sum(1 for r in rows if str(r["tree_mode"]).startswith("unreadable")),
            "families_with_no_leaves": sum(1 for r in rows if r["n_leaves"] == 0 and not str(r["tree_mode"]).startswith("unreadable"))}
    return rows, diag


def summarize_reference(rows, diag) -> dict:
    out = summarize(rows)
    out["diagnostics"] = diag
    human = [r for r in rows if r["n_human"] > 0]
    ht = [r for r in human if r["n_teleost_ref_leaves"] > 0]
    nt = [r for r in human if r["n_teleost_ref_leaves"] == 0]
    nh = [r for r in ht if not r["hit_by_project"]]
    out["reference"] = {
        "human_and_teleost_not_hit_details": {
            "families": len(nh), "audit_statuses": dict(collections.Counter(r["audit_status"] for r in nh)),
            "median_leaves": _median([r["n_leaves"] for r in nh]), "median_teleost_species": _median([r["n_teleost_ref_species"] for r in nh]),
            "with_one_or_two_teleost_species": sum(1 for r in nh if r["n_teleost_ref_species"] <= 2),
            "with_exactly_one_human_gene": sum(1 for r in nh if r["n_human"] == 1), "with_gar_or_reedfish": sum(1 for r in nh if r["n_gar_reedfish"] > 0)},
        "human_and_teleost": {"families": len(ht), "built": sum(r["built"] for r in ht), "hit_by_a_project_species": sum(r["hit_by_project"] for r in ht),
                              "hit_by_the_focal_species": sum(r["has_focal_gene"] for r in ht),
                              "not_hit_by_any_project_species": sum(1 for r in ht if not r["hit_by_project"])},
        "human_and_no_teleost": {"families": len(nt), "built": sum(r["built"] for r in nt), "hit_by_a_project_species": sum(r["hit_by_project"] for r in nt),
                                 "audit_statuses": dict(collections.Counter(r["audit_status"] for r in nt)),
                                 "median_leaves": _median([r["n_leaves"] for r in nt]), "median_species": _median([r["n_species"] for r in nt])},
        "teleost_and_no_human": {"families": sum(1 for r in rows if r["n_teleost_ref_leaves"] > 0 and r["n_human"] == 0),
                                 "built": sum(r["built"] for r in rows if r["n_teleost_ref_leaves"] > 0 and r["n_human"] == 0)},
        "by_audit_status": {st: {"families": sum(1 for r in rows if r["audit_status"] == st),
                                 "with_human": sum(1 for r in human if r["audit_status"] == st)}
                            for st in sorted({r["audit_status"] for r in rows})},
        "human_genes_histogram_by_hit": {"hit_by_a_project_species": _hist([r["n_human"] for r in rows if r["hit_by_project"]]),
                                         "not_hit": _hist([r["n_human"] for r in rows if not r["hit_by_project"]])}}
    return out


def report_reference_text(s: dict, part: str) -> str:
    d, ref = s["diagnostics"], s["reference"]
    L = [f"# Reference-tree survey: {part} (the coordinator's prepared input)", "",
         f"{s['families']} reference families read; files by extension {d['files_by_extension']}; {d['families_with_a_tree_file']} families have a tree "
         f"file and {d['families_with_an_alignment_only']} an alignment only (first files {d['first_files']}). family_audit.tsv: {d['audit_rows']} rows, "
         f"statuses {d['audit_statuses']}; {d['audited_families_without_a_file']} audited families have no file. Leaves read: {d['leaves_total']}, "
         f"with no species: {d['unknown_leaves_total']}; read from the alignment because the tree was empty: "
         f"{d['families_read_from_the_alignment_because_the_tree_was_empty']}; families with no leaves: {d['families_with_no_leaves']}; unreadable: {d['unreadable']}.", "",
         "## Human genes per reference family", ""]
    h = s["human_genes_per_family_histogram"]
    mx = max(h.values()) or 1
    L += ["| human genes | families | |", "|---|---|---|"] + [f"| {k} | {v} | {'#' * max(0, int(40 * v / mx))} |" for k, v in h.items()]
    L += ["", f"{s['with_human']} reference families hold at least one human gene; exactly one in {s['exactly_one_human_gene']}; median "
              f"{s['median_human_genes_where_present']} where present.", "",
          "## Human genes and teleosts, and which of them the project built", "",
          f"Human and at least one reference teleost: {ref['human_and_teleost']}.", "",
          f"Of those, the ones no project species was assigned to: {ref['human_and_teleost_not_hit_details']}.", "",
          f"Human and no reference teleost: {ref['human_and_no_teleost']}.", "",
          f"Teleost and no human: {ref['teleost_and_no_human']}.", "",
          "| human genes | no teleost | with teleost |", "|---|---|---|"]
    L += [f"| {k} | {v['no_teleost']} | {v['with_teleost']} |" for k, v in s["human_genes_by_teleost_presence"].items()]
    L += ["", f"By audit status: {ref['by_audit_status']}.", "",
          f"Human genes per family, families hit by a project species against not hit: {ref['human_genes_histogram_by_hit']}.", ""]
    return "\n".join(L + _mean_section(s))
