import csv

from cichlid_homology.compare import compare_clades, compare_gene_tables, read, report_text, write_outputs

COLS = ["gene_id", "human_gene", "zebrafish_gene", "b1_category", "b1_confidence", "tree_call", "ohnolog_status", "ohnolog_gene", "ohnolog_confidence",
        "teleost_trend", "rel_human", "rel_zebrafish", "clade_id", "clade_type", "out_paralog", "later_dup_level"]


def row(gid, **kw):
    r = {c: "" for c in COLS}
    r.update(gene_id=gid, human_gene="H1", zebrafish_gene="Z1", b1_category="Lost duplicate", tree_call="lost_duplicate", ohnolog_status="single copy (duplicate lost)",
           teleost_trend="1:1", rel_human="1:1", rel_zebrafish="1:1", clade_id="F:G1", clade_type="human clade (1 human gene)", out_paralog="no")
    r.update(kw)
    return r


def test_the_headline_counts_the_genes_that_change_a_call_and_the_closest_genes_are_apart():
    old = [row("a"), row("b"), row("c"), row("d", b1_category="Disputed")]
    new = [row("a"), row("b", b1_category="Both copies kept (best guess)", tree_call="pair"), row("c", human_gene="H2"), row("d", b1_category="Lost duplicate")]
    r = compare_gene_tables(old, new)
    assert (r["genes_compared"], r["only_in_old"], r["only_in_new"]) == (4, 0, 0)
    assert r["genes_with_a_changed_call"] == 2 and r["genes_with_a_changed_call_share"] == 0.5                  # b and d; c changed only a closest gene
    assert r["genes_changed_by_group"]["closest genes"] == 1 and r["genes_changed_by_group"]["B1 call"] == 2
    assert r["columns"]["b1_category"]["changed"] == 2 and ["Disputed", "Lost duplicate", 1] in r["columns"]["b1_category"]["transitions"]
    assert r["changed_genes"]["b"]["tree_call"] == ["lost_duplicate", "pair"]


def test_genes_in_only_one_table_are_counted_apart():
    r = compare_gene_tables([row("a"), row("x")], [row("a"), row("y")])
    assert (r["genes_compared"], r["only_in_old"], r["only_in_new"], r["genes_with_a_changed_call"]) == (1, 1, 1, 0)


def test_pairs_are_kept_lost_gained_or_re_paired():
    def pr(g, p, st="paired"):
        return row(g, ohnolog_status=st, ohnolog_gene=p)
    old = [pr("a", "b"), pr("b", "a"), pr("c", "d"), pr("d", "c"), pr("e", "f"), pr("f", "e"), row("g"), row("h")]
    new = [pr("a", "b", "paired (resolved by synteny)"), pr("b", "a", "paired (resolved by synteny)"), row("c"), row("d"), pr("e", "g"), pr("g", "e"), pr("f", "h"), pr("h", "f")]
    p = compare_gene_tables(old, new)["pairs"]
    assert (p["old"], p["new"], p["kept"], p["lost"], p["gained"], p["kept_with_a_different_status"]) == (3, 3, 1, 2, 2, 1)
    assert p["genes_with_a_different_partner"] == 2                                                         # e: f to g; f: e to h


def test_clades_and_the_files(tmp_path):
    og = [{"family": "F", "group": "G1", "teleost_trend": "1:1", "n_human": "1"}, {"family": "F", "group": "G2", "teleost_trend": "1:2", "n_human": "1"}, {"family": "H", "group": "G1", "teleost_trend": "1:1", "n_human": "1"}]
    ng = [{"family": "F", "group": "G1", "teleost_trend": "1:2", "n_human": "1"}, {"family": "H", "group": "G1", "teleost_trend": "1:1", "n_human": "1"}]
    c = compare_clades(og, ng)
    assert (c["clades_old"], c["clades_new"], c["clades_in_both"], c["families_with_a_different_number_of_clades"], c["changed_teleost_trend"]) == (3, 2, 2, 1, 1)
    old, new = [row("a"), row("b")], [row("a", b1_category="Disputed"), row("b")]
    res = compare_gene_tables(old, new)
    write_outputs(tmp_path, res, c)
    assert "1 genes (50.0%) change at least one call" in (tmp_path / "compare.md").read_text()
    with open(tmp_path / "changed_genes.tsv") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    assert rows[0] == ["gene_id", "b1_category_old", "b1_category_new"] and rows[1] == ["a", "Lost duplicate", "Disputed"]
    assert read(tmp_path / "changed_genes.tsv")[0]["gene_id"] == "a" and "changed_genes" in res                      # the result is left whole


def test_a_gene_that_lost_its_clade_is_counted_apart_and_the_clean_headline_leaves_it_out():
    old = [row("a"), row("b"), row("c"), row("d", b1_category="Disputed")]
    new = [row("a"), row("b", clade_id="", clade_type="unplaced", teleost_trend="", rel_human="", rel_zebrafish="", out_paralog="", ohnolog_status="not applicable"),
           row("c", clade_id="", clade_type="unplaced", teleost_trend="", rel_human="", rel_zebrafish="", out_paralog="", ohnolog_status="not applicable", b1_category="Disputed"),
           row("d", b1_category="Lost duplicate")]
    r = compare_gene_tables(old, new)
    assert r["genes_with_a_changed_call"] == 3 and r["genes_that_lost_their_clade"] == 2 and r["genes_that_gained_a_clade"] == 0
    h = r["clean_headline"]
    assert (h["pilot_calls_changed"], h["clade_calls_changed_same_coverage"], h["any"], h["share"]) == (2, 0, 2, 0.5)        # c and d change the B1 call; b only lost its clade
    assert r["warnings"] and "differ in coverage" in r["warnings"][0]
    assert "Clean headline: 2 genes (50.0%)" in report_text({k: v for k, v in r.items() if k != "changed_genes"})


def test_pairs_lost_only_because_a_clade_is_missing_and_the_coverage_check():
    def pr(g, p):
        return row(g, ohnolog_status="paired", ohnolog_gene=p)
    old = [pr("a", "b"), pr("b", "a"), pr("c", "d"), pr("d", "c")]
    new = [pr("a", "b"), pr("b", "a"), row("c", clade_id="", ohnolog_status="not applicable"), row("d", clade_id="", ohnolog_status="not applicable")]
    p = compare_gene_tables(old, new)["pairs"]
    assert (p["old"], p["new"], p["kept"], p["lost"], p["lost_only_because_a_clade_is_missing"]) == (2, 1, 1, 1, 1)
    from cichlid_homology.compare import clade_coverage
    clades = [{"family": "F", "group": "G1"}]
    assert clade_coverage([row("a"), row("b", clade_id="F:G2")], clades) == {"clade_ids_in_the_gene_table": 2, "not_in_the_clade_table": 1}


def test_a_b_named_pairs_are_classified_before_and_after():
    from cichlid_homology.compare import compare_name_pairs
    def g(gid, sym, fam, clade, partner="", status="single copy (duplicate lost)"):
        return row(gid, symbol=sym, family=fam, clade_id=clade, ohnolog_gene=partner, ohnolog_status=status)
    old = [g("a1", "gnb1a", "F1", "F1:G1"), g("b1", "gnb1b", "F1", "F1:G2"), g("a2", "ascl1a", "F2", "F2:G1"), g("b2", "ascl1b", "F3", "F3:G1"), g("a3", "xyz1a", "F4", "F4:G1"), g("a3x", "xyz1b", "F4", "F4:G1"),
           g("a4", "locb1a", "F5", "F5:G1"), g("a5", "LOC1a", "F6", "F6:G1")]
    new = [g("a1", "gnb1a", "F1", "F1:G1", "b1", "paired"), g("b1", "gnb1b", "F1", "F1:G1", "a1", "paired"), g("a2", "ascl1a", "F2", "F2:G1"), g("b2", "ascl1b", "F3", "F3:G1"),
           g("a3", "xyz1a", "F4", "F4:G1"), g("a3x", "xyz1b", "F4", "F4:G1"), g("a4", "locb1a", "F5", "F5:G1"), g("a5", "LOC1a", "F6", "F6:G1")]
    r = compare_name_pairs(old, new)
    assert r["pairs"] == 3 and r["old"] == {"same family, different clades": 1, "different families": 1, "same clade, not partners": 1}
    assert [["same family, different clades", "partners (a pair of ours)", 1]] == [t for t in r["transitions"] if t[1].startswith("partners")]
    res = compare_gene_tables(old, new)
    assert res["name_pairs"]["pairs"] == 3 and "Genes named a and b" in report_text({k: v for k, v in res.items() if k != "changed_genes"})
