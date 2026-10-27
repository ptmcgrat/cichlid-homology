"""A readable Markdown report of the 0.7.0 summary (overall numbers and agreement)."""
from __future__ import annotations

from pathlib import Path

from .clade_counts import CATS


def _table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    out += ["| " + " | ".join("" if c is None else str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _cat_rows(A, key):
    t, c = A["overall"]["TreeFam"][key]["categories"], A["overall"]["Compara"][key]["categories"]
    return [[cat, f"{t[cat]['n']} ({t[cat]['pct']}%)", f"{c[cat]['n']} ({c[cat]['pct']}%)"] for cat in CATS]


def _size_rows(A, key):
    out = []
    for m in ("TreeFam", "Compara"):
        for lab, d in A["overall"][m][key]["by_family_size"].items():
            tot = sum(d.values()) or 1
            out.append([m, lab, sum(d.values())] + [f"{100 * d.get(cat, 0) / tot:.0f}%" for cat in CATS])
    return out


def write_report(summary: dict, path: Path) -> None:
    A = summary["A"]
    names = {"human": "Human (A1)", "zebrafish": "Zebrafish (A2)"}
    L = [f"# Pilot report: cichlid-homology {summary['cichlid_homology']}, focal species {summary['species']}", "",
         "Categories count reference genes : focal genes in the smallest clade holding the gene and a reference "
         "gene. Focal = the focal species only. `None` = no gene of that reference species in this method's tree.", "",
         "## Overall numbers", "",
         _table(["", "TreeFam", "Compara"],
                [["genes called", A["genes_called"]["TreeFam"], A["genes_called"]["Compara"]],
                 ["genes called in both methods", A["genes_in_both_methods"], A["genes_in_both_methods"]]]), ""]
    for key in ("human", "zebrafish"):
        L += [f"### {names[key]}: categories", "", _table(["Category", "TreeFam", "Compara"], _cat_rows(A, key)), "",
              f"By number of focal genes in the family (share of genes in each category):", "",
              _table(["Method", "Genes in family", "n"] + list(CATS), _size_rows(A, key)), "",
              "Median clade size (leaves) by category: TreeFam "
              + str(A["overall"]["TreeFam"][key]["median_clade_leaves"]) + "; Compara "
              + str(A["overall"]["Compara"][key]["median_clade_leaves"]), "",
              f"Best gene tied on path length: TreeFam {A['overall']['TreeFam'][key]['best_gene_tied_on_path_length']} of "
              f"{A['overall']['TreeFam'][key]['genes_with_a_reference_gene']}; Compara "
              f"{A['overall']['Compara'][key]['best_gene_tied_on_path_length']} of "
              f"{A['overall']['Compara'][key]['genes_with_a_reference_gene']}. "
              f"Low name confidence: TreeFam {A['overall']['TreeFam'][key]['best_gene_name_confidence_low']}, "
              f"Compara {A['overall']['Compara'][key]['best_gene_name_confidence_low']}.", ""]
    L += ["## Agreement between the methods", ""]
    for key in ("human", "zebrafish"):
        g = A["agreement"][key]
        L += [f"### {names[key]}", "",
              _table(["", ""], [
                  ["genes in both methods", g["genes_in_both_methods"]],
                  ["both gave a category", g["both_gave_a_category"]],
                  ["same category", f"{g['same_category']} ({g['same_category_pct']}%)"],
                  ["same reference-gene count class (1 vs 2+)", f"{g['reference_side_same_(1_vs_2plus)']}%"],
                  ["same focal-gene count class (1, 2, 3+)", f"{g['focal_side_same_(1_vs_2_vs_3plus)']}%"],
                  ["only one method has a reference gene", g["only_one_method_has_a_reference_gene"]],
                  ["neither has one", g["neither_has_one"]],
                  ["labels", g["labels"]]]), "",
              "Most common category pairs (TreeFam, Compara):", "",
              _table(["TreeFam", "Compara", "genes"], [[p[0][0], p[0][1], p[1]] for p in g["category_pairs_TreeFam_Compara"][:12]]), ""]
        b = A["best_gene"][key]
        L += [f"Best reference gene (official ID, reliable names only): "
              + ", ".join(f"{k} {v}" for k, v in b.items() if k not in ("adjudication_by_own_ncbi_symbol",)), "",
              f"Where the best genes differ, the gene's own NCBI symbol supports: {b['adjudication_by_own_ncbi_symbol']}", ""]
    L += ["## NCBI symbol check (best gene against the gene's own symbol)", ""]
    for m in ("TreeFam", "Compara"):
        s = A["symbol_check"][m]
        L += [f"{m}: {s['genes_with_a_real_own_symbol']} of {s['genes']} genes have a real own symbol.", ""]
        rows = []
        for key in ("human", "zebrafish"):
            for cat in CATS[:-1]:
                d = s[key][cat]
                if d["checkable"]:
                    rows.append([key, cat, d["checkable"], f"{d['exact_pct']}%", f"{d['exact_or_stem_pct']}%"])
        L += [_table(["Reference", "Category", "checkable", "exact", "exact or stem"], rows), ""]
    L += ["## Goal B: teleost copy numbers on the human clade", ""]
    for m in ("TreeFam", "Compara"):
        gb = A["goal_B"][m]
        L += [f"### {m}", "", _table(["Human category", "genes", "B1 labels"],
                                     [[c, gb[f"human_{c}"]["genes"], gb[f"human_{c}"]["labels"]] for c in CATS]), ""]
        for cat in ("OneToOne", "OneToTwo"):
            d = gb[f"distribution_human_{cat}"]
            L += [f"{cat} ({d['genes_scored']} genes scored). Share of non-cichlid teleost species with 2 or more copies, "
                  f"10 bins from 0 to 1: {d['histogram_share_of_non_cichlid_teleost_species_with_2plus_copies']}. "
                  f"Cichlid species: {d['histogram_share_of_cichlid_species_with_2plus_copies']}. "
                  f"Patterns: {d['patterns_(teleost,percomorph,ovalentaria,cichlid)']}.", ""]
    ga = A["goal_B_agreement"]
    L += [f"B1 label agreement for genes OneToOne or OneToTwo in both methods: {ga['genes_OneToOne_or_OneToTwo_in_both_methods']} genes, "
          f"{ga['labels']}, agree {ga['agree_pct_when_both_labelled']}% when both gave a label.", "",
          f"Category splits listed in `pilot_category_splits.tsv`: {A['category_splits_written']} rows.", ""]
    Path(path).write_text("\n".join(L))
