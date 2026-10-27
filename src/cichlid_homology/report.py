"""A readable Markdown report of the summary: coverage, overall numbers, agreement between the methods,
the symbol check, the teleost distribution, where the extra copies sit, and the cichlid-specific candidates."""
from __future__ import annotations

from pathlib import Path

from .clade_counts import GROUPS
from .summary_a import ORIGIN_CLASSES

NAMES = {"human": "Human (A1)", "zebrafish": "Zebrafish (A2)"}


def _table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    out += ["| " + " | ".join("" if c is None else str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _cat_rows(A, key):
    t, c = A["overall"]["TreeFam"][key]["categories"], A["overall"]["Compara"][key]["categories"]
    return [[cat, f"{t[cat]['n']} ({t[cat]['pct']}%)", f"{c[cat]['n']} ({c[cat]['pct']}%)"] for cat in GROUPS[key]]


def _size_rows(A, key):
    out = []
    for m in ("TreeFam", "Compara"):
        for lab, d in A["overall"][m][key]["by_family_size"].items():
            tot = sum(d.values()) or 1
            out.append([m, lab, sum(d.values())] + [f"{100 * d.get(cat, 0) / tot:.0f}%" for cat in GROUPS[key]])
    return out


def _origin_rows(d):
    return [[k, v["genes"]] + [f"{v[c]['n']} ({v[c]['pct']}%)" for c in ORIGIN_CLASSES[:4]]
            for k, v in [("all genes with 2+ copies", d["all"])] + list(d["by_teleost_profile_label"].items())]


def write_report(summary: dict, path: Path) -> None:
    A = summary["A"]
    run = summary.get("run", {})
    L = [f"# Pilot report: cichlid-homology {summary['cichlid_homology']}, focal species {summary['species']}", "",
         f"Mode: **{run.get('mode', 'sample')}**. The focal species has {run.get('species_genes_total', '?')} genes in "
         f"the assignment table; {run.get('finished_family_genes', {})} have a finished family in each method. "
         "Categories count reference genes : focal genes in the smallest clade holding the gene and a reference gene "
         "(focal = the focal species only). `None` = no gene of that reference species in this method's tree. "
         "Human: Other = ManyToOne + ManyToMany. Zebrafish: OneToMany = 1:2 and 1:many, Other = many:1 and many:many.", "",
         "## Overall numbers", "",
         _table(["", "TreeFam", "Compara"],
                [["genes called", A["genes_called"]["TreeFam"], A["genes_called"]["Compara"]],
                 ["genes called in both methods", A["genes_in_both_methods"], A["genes_in_both_methods"]]]), ""]
    for key in ("human", "zebrafish"):
        L += [f"### {NAMES[key]}: categories", "", _table(["Category", "TreeFam", "Compara"], _cat_rows(A, key)), "",
              "By number of focal genes in the family (share of genes in each category):", "",
              _table(["Method", "Genes in family", "n"] + list(GROUPS[key]), _size_rows(A, key)), "",
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
        L += [f"### {NAMES[key]}", "",
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
        L += ["Best reference gene (official ID, reliable names only): "
              + ", ".join(f"{k} {v}" for k, v in b.items() if k != "adjudication_by_own_ncbi_symbol"), "",
              f"Where the best genes differ, the gene's own NCBI symbol supports: {b['adjudication_by_own_ncbi_symbol']}", ""]
    L += ["## NCBI symbol check (best gene against the gene's own symbol)", ""]
    for m in ("TreeFam", "Compara"):
        s = A["symbol_check"][m]
        L += [f"{m}: {s['genes_with_a_real_own_symbol']} of {s['genes']} genes have a real own symbol.", ""]
        rows = []
        for key in ("human", "zebrafish"):
            for cat in GROUPS[key][:-1]:
                d = s[key][cat]
                if d["checkable"]:
                    rows.append([key, cat, d["checkable"], f"{d['exact_pct']}%", f"{d['exact_or_stem_pct']}%"])
        L += [_table(["Reference", "Category", "checkable", "exact", "exact or stem"], rows), ""]
    L += ["## Goal B: teleost copy numbers on the human clade", ""]
    for m in ("TreeFam", "Compara"):
        gb = A["goal_B"][m]
        L += [f"### {m}", "", _table(["Human category", "genes", "B1 labels"],
                                     [[c, gb[f"human_{c}"]["genes"], gb[f"human_{c}"]["labels"]] for c in GROUPS["human"]]), ""]
        for cat in ("OneToOne", "OneToTwo"):
            d = gb[f"distribution_human_{cat}"]
            L += [f"{cat} ({d['genes_scored']} genes scored). Share of non-cichlid teleost species with 2 or more copies, "
                  f"10 bins from 0 to 1: {d['histogram_share_of_non_cichlid_teleost_species_with_2plus_copies']}. "
                  f"Cichlid species: {d['histogram_share_of_cichlid_species_with_2plus_copies']}. "
                  f"Patterns: {d['patterns_(teleost,percomorph,ovalentaria,cichlid)']}.", ""]
    ga = A["goal_B_agreement"]
    L += [f"B1 label agreement for genes OneToOne or OneToTwo in both methods: {ga['genes_OneToOne_or_OneToTwo_in_both_methods']} genes, "
          f"{ga['labels']}, agree {ga['agree_pct_when_both_labelled']}% when both gave a label.", ""]
    L += ["## Lost duplicate or retained ohnolog by parsimony (0.9.0 pilot)", "",
          "Copy numbers of the jawed-vertebrate species in the anchor clade are explained by gains and losses on "
          "the species tree; the teleost ancestor having 1, 2 or 3 copies is compared by the fewest events. No share "
          "thresholds. `ambiguous` = a tie, or the call changes when a loss costs 2 gains or the reverse.", ""]
    for m in ("TreeFam", "Compara"):
        p = A["parsimony"][m]
        if not p.get("scored"):
            L += [f"### {m}: nothing scored ({p.get('status')})", ""]
            continue
        lab = ", ".join("%s %d (%s%%)" % (k, v["n"], v["pct"]) for k, v in p["labels"].items())
        L += [f"### {m}: {p['scored']} of {p['genes']} genes scored ({p['status']})", "",
              f"Labels: {lab}. Margin (events separating the best and second-best history), 0 to 6+: "
              f"{p['margin_0_1_2_3_4_5_6plus']}. Robust to the relative cost of a loss: {p['robust_to_cost_changes_pct']}%. "
              f"Zero counts read as losses or as missing data give the same label for "
              f"{p['zero_as_loss_vs_missing_agree_pct']}% of {p['zero_modes_compared']} genes.", "",
              "Against the threshold-based B1 label (rows: parsimony label):", "",
              _table(["parsimony label", "B1 labels"], [[k, v] for k, v in p["vs_threshold_B1_labels"].items()]), "",
              "By human group: " + str(p["by_human_group"]), "",
              "Check 1, zebrafish a/b names (share of best zebrafish gene names ending in a or b): "
              + str(p["check_zebrafish_ab_names_pct"]), "",
              "Check 2, position of the nearest other copy: " + str(p["check_copy_positions"]), "",
              "Later gains below the teleost ancestor (the third B1 outcome): "
              + str(p["later_duplicates_after_the_teleost_ancestor"]), "",
              f"Gain events by ladder level (all scored genes): {p['gain_events_by_level_all_scored_genes']}. "
              f"Of {p['retained_ohnologs_total']} retained ohnologs, {p['retained_ohnologs_with_a_species_level_loss']}% "
              f"have a species-level loss; most often: {p['species_most_often_lost_in_retained_ohnologs']}.", ""]
    pa = A["parsimony_agreement"]
    L += [f"TreeFam and Compara parsimony labels (ambiguous counted as no call): {pa['labels']}; agree "
          f"{pa['agree_pct_when_both_called']}% when both called. Pairs: {pa['pairs_TreeFam_Compara']}", ""]
    L += ["## Where the extra copies sit", "",
          "For genes with two or more focal-species copies in the reference clade: the nearest other copy is "
          "adjacent (no gene between), near (at most the set number of genes between), distant on the same "
          "sequence, or on a different sequence. Teleost-wide duplicates are expected on different sequences; local "
          "duplications next to each other.", ""]
    for m in ("TreeFam", "Compara"):
        o = A["origin"][m]
        L += [f"### {m}: {o['genes_with_2plus_copies_in_the_reference_clade']} genes with 2+ copies", "",
              _table(["Teleost profile label", "genes", "adjacent", "near", "same seq., distant", "different seq."],
                     _origin_rows(o)), ""]
    L += ["## Cichlid-specific candidates", "",
          "Genes labelled `cichlid_specific_duplication`: the extra copies are in the cichlids and not in the other "
          "teleosts. Risk flags are signs of an annotation or placement artifact.", ""]
    for m in ("TreeFam", "Compara"):
        c = A["candidates"][m]
        near = ", ".join("%s %d" % (k, c["origin"][k]["n"]) for k in ORIGIN_CLASSES[:4])
        L += [f"### {m}: {c['genes']} genes", "",
              f"Share of cichlid species carrying the extra copy: {c['carriers']}. "
              f"Nearest other copy: {near}. "
              f"No risk flag: {c['with_no_risk_flag']}; two or more: {c['with_2plus_risk_flags']}. "
              f"Flags: {c['risk_flag_counts']}.", ""]
    L += [f"Candidates in both methods: {A['candidates_in_both_methods']} "
          f"({A['candidates_in_both_methods_with_no_risk_flag_in_either']} with no risk flag in either).", "",
          f"Category splits listed in `pilot_category_splits.tsv`: {A['category_splits_written']} rows. "
          "Candidate tables: `pilot_candidates.TreeFam.tsv`, `pilot_candidates.Compara.tsv`.", ""]
    Path(path).write_text("\n".join(L))
