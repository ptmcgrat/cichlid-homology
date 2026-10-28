from cichlid_homology.genetable import COLUMNS, COLUMN_SPEC, build_gene_table, columns_markdown, human_tier, summarize_table


def clade(fam="F1", grp="G1", human="H1", zfish="Z1;Z2", trend="1:2", med=2, early=2, perc=2, nh=None, nz=None, tel=30, tet=11, share=1.0, rate=0.95, tshare=1.0):
    hs = [p for p in human.split(";") if p]
    zs = [p for p in zfish.split(";") if p]
    return {"family": fam, "group": grp, "n_human": nh if nh is not None else len(hs), "human_proteins": human, "zebrafish_proteins": zfish,
            "zebrafish_genes": nz if nz is not None else len(zs), "n_teleost_species": tel, "teleost_trend": trend, "teleost_median_genes": med,
            "early_teleost_median_genes": early, "percomorph_median_genes": perc, "n_tetrapod_species": tet, "tetrapod_single_copy_share": share,
            "human_rel_rate": rate, "teleost_species_share_of_family": tshare}


def gene(gid, fam="F1", z="Z1", h="H1", lineage="A", near="", cls="", cat="Both copies kept (best guess)", tier="core", **kw):
    r = {"gene_id": gid, "own_symbol": f"sym_{gid}", "gt_description": "d", "family": fam, "a_human_best": h, "a_human_best_symbol": "HSYM", "a_human_best_id": "HGNC:1",
         "a_human_best_name_reliable": "True", "a_human_group": "OneToTwo", "a_zebrafish_best": z, "a_zebrafish_best_symbol": f"zs_{z}", "a_zebrafish_best_id": "ZDB",
         "a_zebrafish_group": "OneToTwo", "d_gene_lineage": lineage, "orig_nearest_copy": near, "orig_class": cls, "t_category": cat, "t_tier": tier, "t_basis": "b",
         "t_vote_thr": "pair", "t_vote_par": "pair", "t_vote_tree": "pair", "t_call": "retained_pair", "t_when": "", "d_call": "pair", "expanded_family": "False",
         "risk_flags": "", "a_human_symbol_match": "exact", "a_zebrafish_symbol_match": "none", "t_pairs_relaxed_n": 0}
    r.update(kw)
    return r


FAM = {"F1": {"family": "F1", "n_human": 1}, "F2": {"family": "F2", "n_human": 0}, "F3": {"family": "F3", "n_human": 2}}


def table(rows, groups, **kw):
    return {r["gene_id"]: r for r in build_gene_table(rows, kw.pop("treefam", []), groups, FAM, version="9.9", run_date="d", **kw)}


def test_the_table_has_66_documented_columns_in_order():
    assert len(COLUMNS) == len(COLUMN_SPEC) == 66 and COLUMNS[0] == "gene_id" and COLUMNS[-1] == "parameters"
    assert "## Ohnolog" in columns_markdown()


def test_a_pair_in_a_1_to_2_clade_with_all_three_signals_is_high_confidence():
    g = [gene("g1", z="Z1", lineage="A", near="g2", cls="different_seqid"), gene("g2", z="Z2", lineage="B", near="g1", cls="different_seqid")]
    t = table(g, [clade()])
    a, b = t["g1"], t["g2"]
    assert (a["ohnolog_status"], a["ohnolog_gene"], b["ohnolog_gene"], a["pair_id"] == b["pair_id"] == "F1:G1") == ("paired", "g2", "g1", True)
    assert (a["ev_lineages"], a["ev_position"], a["ev_zebrafish_split"], a["ohnolog_confidence"]) == ("different", "different chromosomes", "yes", "high")
    assert (a["ohnolog_zebrafish_gene"], a["ohnolog_lineage"], b["ohnolog_lineage"], a["ohnolog_symbol"]) == ("Z2", "A", "B", "sym_g2")
    assert (a["rel_human"], a["rel_zebrafish"], a["teleost_trend"], a["copy_vs_trend"], a["early_vs_percomorph"]) == ("1:2", "2:2", "1:2", "matches trend", "2/2")
    assert (a["clade_assignment"], a["clade_id"], a["clade_n_focal"], a["clade_type"], a["out_paralog"]) == ("zebrafish", "F1:G1", 2, "human clade (1 human gene)", "no")
    assert a["human_clean_tier"] == "clean" and a["pipeline_version"] == "9.9" and a["human_agree_treefam"] == "not in TreeFam"


def test_the_confidence_falls_with_the_signals_and_tandem_copies_count_against():
    same_zeb = [gene("g1", z="Z1", lineage="A", near="g2", cls="different_seqid"), gene("g2", z="Z1", lineage="B", near="g1", cls="different_seqid")]
    assert table(same_zeb, [clade()])["g1"]["ohnolog_confidence"] == "medium"                       # two signals, none against
    tandem = [gene("g1", z="Z1", lineage="A", near="g2", cls="tandem_adjacent"), gene("g2", z="Z2", lineage="B", near="g1", cls="tandem_adjacent")]
    assert (table(tandem, [clade()])["g1"]["ev_position"], table(tandem, [clade()])["g1"]["ohnolog_confidence"]) == ("tandem", "low")
    same_lin = [gene("g1", z="Z1", lineage="A", near="g2", cls="different_seqid"), gene("g2", z="Z2", lineage="A", near="g1", cls="different_seqid")]
    assert (table(same_lin, [clade()])["g1"]["ev_lineages"], table(same_lin, [clade()])["g1"]["ohnolog_confidence"]) == ("same", "low")
    far = [gene("g1", near="x", cls="different_seqid"), gene("g2", z="Z2", lineage="B", near="y", cls="different_seqid")]
    assert table(far, [clade()])["g1"]["ev_position"] == "unknown"                                  # the partner is not the nearest copy and there are no positions


def test_positions_decide_the_chromosome_signal_when_given():
    g = [gene("g1", z="Z1", lineage="A"), gene("g2", z="Z2", lineage="B")]
    pos = {"g1": ("chr1", 100, 5), "g2": ("chr1", 900, 40)}
    assert table(g, [clade()], positions=pos)["g1"]["ev_position"] == "same chromosome, distant"
    assert table(g, [clade()], positions={"g1": ("chr1", 100, 5), "g2": ("chr1", 200, 8)})["g1"]["ev_position"] == "tandem"
    t = table(g, [clade()], positions={"g1": ("chr1", 100, 5), "g2": ("chr7", 200, 8)})
    assert (t["g1"]["ev_position"], t["g1"]["seqid"], t["g1"]["gene_index"]) == ("different chromosomes", "chr1", 5)


def test_one_focal_gene_in_a_1_to_2_clade_means_the_partner_was_lost_and_names_the_zebrafish_gene():
    t = table([gene("g3", z="Z1")], [clade()])
    assert (t["g3"]["ohnolog_status"], t["g3"]["copy_vs_trend"], t["g3"]["partner_elsewhere"], t["g3"]["ohnolog_gene"]) == ("partner lost in the focal species", "fewer copies", "Z2", "")


def test_three_or_more_focal_genes_in_a_1_to_2_clade_are_ambiguous():
    t = table([gene(f"g{i}", z="Z1") for i in range(3)], [clade()])
    assert (t["g0"]["ohnolog_status"], t["g0"]["copy_vs_trend"], t["g0"]["ohnolog_gene"]) == ("partner ambiguous (three or more focal genes)", "extra copies", "")


def test_the_1_to_1_clades_single_copy_later_duplicate_and_possible_pair():
    c = clade(trend="1:1", med=1, early=1, perc=1, zfish="Z1")
    single = table([gene("a", cat="Lost duplicate")], [c])["a"]
    assert (single["ohnolog_status"], single["rel_human"], single["rel_zebrafish"], single["copy_vs_trend"]) == ("single copy (duplicate lost)", "1:1", "1:1", "matches trend")
    tan = table([gene("a", near="b", cls="tandem_near"), gene("b", near="a", cls="tandem_near")], [c])
    assert tan["a"]["ohnolog_status"] == "later duplicate" and tan["a"]["copy_vs_trend"] == "extra copies" and tan["a"]["ohnolog_gene"] == ""
    dis = table([gene("a", near="b", cls="different_seqid"), gene("b", near="a", cls="different_seqid")], [c])
    assert dis["a"]["ohnolog_status"] == "possible pair (the other teleosts lost a copy)" and dis["a"]["ohnolog_gene"] == "b"
    many = table([gene(x) for x in "abc"], [c])
    assert many["a"]["ohnolog_status"] == "later duplicates (three or more focal genes)"


def test_expansion_mixed_and_out_paralog_clades():
    assert table([gene("a")], [clade(trend="1:3+", med=3)])["a"]["ohnolog_status"] == "expansion"
    m = table([gene("a")], [clade(trend="mixed 1-2", med=1.5)])["a"]
    assert (m["ohnolog_status"], m["copy_vs_trend"], m["teleost_median_genes"]) == ("mixed trend", "mixed trend", "1.5")
    fish = clade(fam="F3", human="", zfish="Z1", trend="1:1", med=1, tet=0, nh=0)
    fish.update(n_tetrapod_species=0)
    o = table([gene("a", fam="F3", h="")], [fish])["a"]
    assert (o["out_paralog"], o["clade_type"], o["rel_human"], o["human_clean_tier"]) == ("yes", "fish-only clade, family has a human gene (out-paralog)", "no human gene; no tetrapods", "not applicable")
    fam2 = table([gene("a", fam="F2", h="")], [dict(fish, family="F2")])["a"]
    assert (fam2["out_paralog"], fam2["clade_type"]) == ("no", "fish-only clade, family has no human gene")
    tet = table([gene("a", h="")], [clade(human="", zfish="Z1", trend="1:1", med=1, nh=0, tet=8)])["a"]
    assert (tet["clade_type"], tet["rel_human"]) == ("tetrapod clade, no human gene", "no human gene; tetrapods present")


def test_placement_through_the_human_protein_needs_teleosts_and_otherwise_the_gene_is_unplaced():
    c = clade(zfish="", trend="1:1", med=1)
    t = table([gene("a", z="")], [c])["a"]
    assert (t["clade_assignment"], t["clade_id"]) == ("human", "F1:G1")
    none = table([gene("a", z="")], [dict(c, n_teleost_species=0)])["a"]
    assert (none["clade_assignment"], none["clade_id"], none["clade_type"], none["ohnolog_status"]) == ("unplaced", "", "unplaced", "not applicable")
    other = table([gene("a", z="Z1")], [dict(clade(zfish="Z1", human="H9"), group="G2"), clade(grp="G3", human="H1", zfish="Z5")])["a"]
    assert other["clade_assignment"] == "zebrafish (closest human gene is in another clade)"


def test_the_human_clean_tier_and_treefam_agreement():
    assert human_tier(clade()) == "clean" and human_tier(clade(nh=2)) == "not applicable"
    assert human_tier(clade(share=0.75, tet=6)) == "acceptable" and human_tier(clade(share=0.5)) == "messy" and human_tier(clade(rate=1.5)) == "messy"
    assert human_tier(clade(tshare=0.85)) == "acceptable" and human_tier(clade(tet=3)) == "messy"
    tf = [{"gene_id": "a", "a_human_best_id": "HGNC:1", "a_human_group": "OneToOne", "a_zebrafish_best_id": "ZDB", "a_zebrafish_group": "OneToOne"},
          {"gene_id": "b", "a_human_best_id": "HGNC:2", "a_human_group": "OneToOne", "a_zebrafish_best_id": "ZDB2", "a_zebrafish_group": "None"}]
    t = table([gene("a"), gene("b"), gene("c")], [clade(trend="1:1", med=1, zfish="Z1")], treefam=tf)
    assert (t["a"]["human_agree_treefam"], t["b"]["human_agree_treefam"], t["c"]["human_agree_treefam"]) == ("same", "different", "not in TreeFam")
    assert (t["a"]["zebrafish_agree_treefam"], t["b"]["zebrafish_agree_treefam"]) == ("same", "TreeFam has none")


def test_the_published_pairs_are_an_outside_check_and_flags_pass_through():
    g = [gene("g1", z="Z1", lineage="A", near="g2", cls="different_seqid", t_pairs_relaxed_n=2), gene("g2", z="Z2", lineage="B", near="g1", cls="different_seqid", t_pairs_relaxed_n=2)]
    pub = {"zebrafish": {"zs_z1": {"zs_z2"}, "zs_z2": {"zs_z1"}}}
    assert table(g, [clade()], partners=pub)["g1"]["published_pair"] == "the two zebrafish names are a published pair"
    assert table(g, [clade()])["g1"]["published_pair"] == "gene in a published pair in 2 of 4 species"
    u = table([gene("u", cat="Unclassified", tier="unclassified", t_call="many_copies", expanded_family="True", risk_flags="partial_gene_model")], [clade()])["u"]
    assert (u["unclassified_reason"], u["expanded_family"], u["risk_flags"], u["symbol_match_human"]) == ("many_copies", "True", "partial_gene_model", "exact")


def test_the_summary_counts_statuses_and_confidence():
    g = [gene("g1", z="Z1", lineage="A", near="g2", cls="different_seqid"), gene("g2", z="Z2", lineage="B", near="g1", cls="different_seqid"), gene("x", fam="F9", z="", h="")]
    rows = build_gene_table(g, [], [clade()], FAM, version="9")
    s = summarize_table(rows)
    assert s["genes"] == 3 and s["ohnolog_status"] == {"paired": 2, "not applicable": 1} and s["pair_confidence"] == {"high": 2} and s["clade_assignment"]["unplaced"] == 1


def test_both_genes_of_a_pair_get_the_same_evidence_even_when_only_one_points_to_the_other():
    g = [gene("g1", z="Z1", lineage="A", near="g2", cls="different_seqid"), gene("g2", z="Z2", lineage="B", near="x", cls="tandem_near")]   # g2's nearest copy is another gene
    t = table(g, [clade()])
    assert t["g1"]["ev_position"] == t["g2"]["ev_position"] == "different chromosomes"
    assert t["g1"]["ohnolog_confidence"] == t["g2"]["ohnolog_confidence"] == "high"
    assert (t["g1"]["ohnolog_zebrafish_gene"], t["g2"]["ohnolog_zebrafish_gene"]) == ("Z2", "Z1")


def test_a_possible_pair_in_a_1_to_1_clade_has_the_same_status_for_both_genes():
    c = clade(trend="1:1", med=1, zfish="Z1")
    one_sided = table([gene("a", near="b", cls="different_seqid"), gene("b", near="x", cls="tandem_near")], [c])
    assert one_sided["a"]["ohnolog_status"] == one_sided["b"]["ohnolog_status"] == "possible pair (the other teleosts lost a copy)"
    assert one_sided["a"]["ohnolog_gene"] == "b" and one_sided["b"]["ohnolog_gene"] == "a"
    tandem = table([gene("a", near="b", cls="tandem_near"), gene("b", near="x", cls="different_seqid")], [c])
    assert tandem["a"]["ohnolog_status"] == tandem["b"]["ohnolog_status"] == "later duplicate"
