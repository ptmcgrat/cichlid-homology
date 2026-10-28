from helpers import PREFIX2SP, FakeClades, st
from cichlid_homology.cladesplit import run_split, split_family, summarize_split, report_split_text, GROUP_COLS, FAMILY_COLS
from cichlid_homology.parsimony import ParsimonyContext

NAMES = {"human": "homo_sapiens", "mouse": "mmus", "zebrafish": "danio_rerio", "gar_reed": {"locu", "ecal"}, "coelacanth": "ecal"}


class Clades(FakeClades):
    def lineage(self, s):
        return "sarcopterygian" if s in ("homo_sapiens", "mmus", "ecal") else None


class LS:
    def species(self, leaf):
        return PREFIX2SP.get(leaf.split("_")[0])


def run(newick, fid="F"):
    sp, cl = st(), Clades()
    return split_family(fid, newick, LS(), sp, cl, ParsimonyContext(sp, cl), NAMES)


def test_two_paralogs_that_both_span_tetrapods_and_fish_are_cut_into_two_clades():
    fam, groups = run("(((Hsap_1,Mmus_1),(Aang_1,Drer_1)),((Hsap_2,Mmus_2),(Aang_2,Drer_2)));")
    assert set(fam) == set(FAMILY_COLS) and set(groups[0]) == set(GROUP_COLS)
    assert (fam["n_human"], fam["n_groups"], fam["n_ancient_dups"], fam["n_groups_with_human"]) == (2, 2, 1, 2)
    for g in groups:
        assert (g["n_human"], g["n_tetrapod_species"], g["tetrapod_mean_genes"], g["tetrapod_single_copy_share"]) == (1, 2, 1.0, 1.0)
        assert (g["n_teleost_species"], g["teleost_mean_genes"], g["teleost_species_share_of_family"]) == (2, 1.0, 1.0)


def test_a_clean_single_gene_family_is_one_clade_without_a_cut():
    fam, groups = run("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,Drer_1)));")
    assert (fam["n_groups"], fam["n_ancient_dups"]) == (1, 0) and groups[0]["n_gar_reedfish"] == 1 and groups[0]["n_genes"] == 5


def test_a_human_only_duplicate_is_counted_inside_the_clade_not_cut():
    fam, groups = run("(((Hsap_1,Hsap_2),Mmus_1),(Locu_1,(Aang_1,Drer_1)));")
    g = groups[0]
    assert (fam["n_groups"], g["n_human"], g["n_tetrapod_genes"], g["n_tetrapod_species"], g["tetrapod_mean_genes"]) == (1, 2, 3, 2, 1.5)
    assert (g["tetrapod_single_copy_share"], g["tetrapod_species_dups"], g["tetrapod_clade_dups"]) == (0.5, 1, 0)


def test_a_kept_teleost_duplication_is_counted_inside_the_clade_not_cut():
    fam, groups = run("((Hsap_1,Mmus_1),(Locu_1,((Aang_1,Drer_1),(Aang_2,Drer_2))));")
    g = groups[0]
    assert (fam["n_groups"], fam["n_ancient_dups"], g["teleost_mean_genes"], g["teleost_clade_dups"]) == (1, 0, 2.0, 1)


def test_a_paralog_lost_in_fish_is_still_cut_and_gives_a_clade_without_teleosts():
    fam, groups = run("(((Hsap_1,Mmus_1),(Aang_1,Drer_1)),(Hsap_2,Mmus_2));")
    assert fam["n_groups"] == 2
    with_tel = [g for g in groups if g["n_teleost_species"] > 0][0]
    no_tel = [g for g in groups if g["n_teleost_species"] == 0][0]
    assert with_tel["teleost_species_share_of_family"] == 1.0 and no_tel["teleost_mean_genes"] == "" and no_tel["n_human"] == 1


def test_branch_lengths_give_root_to_tip_distances_from_the_clade_root():
    fam, groups = run("((Hsap_1:0.4,Mmus_1:0.1):0.1,(Locu_1:0.1,(Aang_1:0.1,Drer_1:0.1):0.1):0.1);")
    g = groups[0]
    assert 0.45 <= g["human_rtt"] <= 0.55 and abs(g["tetrapod_rtt_median"] - 0.35) < 0.06 and g["lengths_missing"] == 0
    assert g["human_rel_rate"] == ""                                  # fewer than three other tetrapods
    assert run("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,Drer_1)));")[1][0]["lengths_missing"] == 1


def test_a_one_member_family_and_an_unknown_leaf():
    fam, groups = run("Hsap_1;")
    assert (fam["n_groups"], groups[0]["n_human"]) == (1, 1)
    fam2, groups2 = run("((Hsap_1,Mmus_1),(Zzzz_1,(Aang_1,Drer_1)));")
    assert groups2[0]["n_genes"] == 4


def test_the_run_and_the_summary(tmp_path):
    prep = tmp_path / "prepared"
    (prep / "trees").mkdir(parents=True)
    trees = {"P1": "(((Hsap_1,Mmus_1),(Aang_1,Drer_1)),((Hsap_2,Mmus_2),(Aang_2,Drer_2)));", "P2": "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,Drer_1)));",
             "P3": "(Aang_1,Drer_1);", "P4": "(((Hsap_1,Hsap_2),Mmus_1),(Locu_1,(Aang_1,Drer_1)));"}
    for f, t in trees.items():
        (prep / "trees" / f"{f}.nhx").write_text(t)
    sp, cl = st(), Clades()
    fams, groups, errors, skipped = run_split(prep, LS(), sp, cl, ParsimonyContext(sp, cl), NAMES)
    assert [f["family"] for f in fams] == ["P1", "P2", "P4"] and skipped == 1 and errors == []
    s = summarize_split(fams, groups, errors, skipped)
    assert s["groups_per_family_by_human_genes"] == {"1": {"1": 1}, "2": {"1": 1, "2": 1}}
    assert s["groups_by_human_genes"] == {1: 3, 2: 1}
    one = s["groups_with_one_human_gene"]
    assert one["all"]["groups"] == 3 and one["from_families_with_several_human_genes"]["groups"] == 2
    assert one["all"]["tetrapod_mean_genes_bins"]["1"] == 3 and one["all"]["tetrapod_single_copy_share_bins"]["1 (every tetrapod species has one gene)"] == 3
    assert "Clades with exactly one human gene" in report_split_text(s, "Compara")


def test_a_tetrapod_only_duplication_older_than_the_amniote_split_is_cut_even_with_the_fish_outside():
    from cichlid_homology.speciestree import SpeciesTree
    sp = SpeciesTree("(bflo,(cmil,(((homo_sapiens,mmus),ggal),((locu,ecal),(aang,(danio_rerio,((trub,gacu),(olat,(maylandia_zebra,acal)))))))));")

    class C2(Clades):
        def level(self, s):
            return "gnathostome" if s == "ggal" else super().level(s)

        def lineage(self, s):
            return "sarcopterygian" if s == "ggal" else super().lineage(s)

    class LS2(LS):
        def species(self, leaf):
            return "ggal" if leaf.startswith("Ggal_") else super().species(leaf)

    cl = C2()
    names = {**NAMES, "chicken": "ggal"}
    # the fish are outside; the two tetrapod paralogs each hold mammals and a bird, so the duplication predates the amniote split
    nwk = "((Aang_1,Drer_1),(((Hsap_2,Mmus_2),Ggal_2),((Hsap_1,Mmus_1),Ggal_1)));"
    fam, groups = split_family("F", nwk, LS2(), sp, cl, ParsimonyContext(sp, cl), names, 0.5, True)
    assert (fam["n_groups"], fam["n_ancient_dups"]) == (3, 1)
    assert sorted(g["n_human"] for g in groups) == [0, 1, 1] and sorted(g["n_teleost_species"] for g in groups) == [0, 0, 2]
    # a duplication inside the mammals only (no bird in both sides) is not older than the amniote split: not cut
    fam2, g2 = split_family("F", "((Aang_1,Drer_1),(((Hsap_2,Mmus_2),(Hsap_1,Mmus_1)),Ggal_1));", LS2(), sp, cl, ParsimonyContext(sp, cl), names, 0.5, True)
    assert fam2["n_groups"] == 1 and g2[0]["tetrapod_clade_dups"] == 1


def test_the_rooting_option_keeps_the_trees_own_root_or_reconciles():
    fam_rec, _ = run("(((Hsap_1,Mmus_1),(Aang_1,Drer_1)),(Hsap_2,Mmus_2));")
    assert fam_rec["rooting"] == "recorded" and fam_rec["n_groups"] == 2
    sp, cl = st(), Clades()
    fam, groups = split_family("F", "(((Hsap_1,Mmus_1),(Aang_1,Drer_1)),(Hsap_2,Mmus_2));", LS(), sp, cl, ParsimonyContext(sp, cl), NAMES, 0.5, False)
    assert fam["rooting"] == "reconciliation"


def test_each_clade_lists_its_human_and_zebrafish_protein_ids_so_genes_can_be_joined_to_it():
    fam, groups = run("(((Hsap_1,Mmus_1),(Aang_1,Drer_1)),((Hsap_2,Mmus_2),(Aang_2,Drer_2)));")
    by = {g["human_proteins"]: g for g in groups}
    assert set(by) == {"Hsap_1", "Hsap_2"} and by["Hsap_1"]["zebrafish_proteins"] == "Drer_1" and by["Hsap_2"]["zebrafish_proteins"] == "Drer_2"
    one = run("((Hsap_1,Hsap_2),(Mmus_1,(Aang_1,Drer_1)));")[1][0]
    assert one["human_proteins"] == "Hsap_1;Hsap_2" and one["zebrafish_proteins"] == "Drer_1"
    solo = run("Hsap_1;")[1][0]
    assert (solo["human_proteins"], solo["zebrafish_proteins"]) == ("Hsap_1", "")


def test_the_teleost_median_the_shares_and_the_trend_label():
    from cichlid_homology.cladesplit import trend_label
    assert [trend_label(x) for x in (None, "", 1, 1.0, 1.5, 2, 2.5, 3, 4)] == ["", "", "1:1", "1:1", "mixed 1-2", "1:2", "mixed 2-3", "1:3+", "1:3+"]
    # aang 1 gene, drer 2, trub 3: the median is 2; early teleosts (aang, drer) 1.5; percomorphs (trub) 3
    fam, groups = run("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,Drer_2,(Trub_1,Trub_2,Trub_3)))));")
    g = groups[0]
    assert (g["teleost_median_genes"], g["teleost_trend"], g["zebrafish_genes"]) == (2, "1:2", 2)
    assert (g["teleost_share_one"], g["teleost_share_two"], g["teleost_share_three_plus"]) == (0.333, 0.333, 0.333)
    assert (g["early_teleost_median_genes"], g["percomorph_median_genes"]) == (1.5, 3)
    clean = run("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,Drer_1)));")[1][0]
    assert (clean["teleost_median_genes"], clean["teleost_trend"], clean["teleost_share_one"]) == (1, "1:1", 1.0)
    none = run("(((Hsap_1,Mmus_1),(Aang_1,Drer_1)),(Hsap_2,Mmus_2));")[1]
    assert [x["teleost_trend"] for x in none if x["n_teleost_species"] == 0] == [""]


def test_the_summary_tabulates_the_trend_against_the_zebrafish_genes(tmp_path):
    prep = tmp_path / "prepared"
    (prep / "trees").mkdir(parents=True)
    trees = {"Q1": "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,Drer_1)));",                                  # 1:1, zebrafish 1
             "Q2": "((Hsap_1,Mmus_1),(Locu_1,((Aang_1,Drer_1),(Aang_2,Drer_2))));",                # 1:2, zebrafish 2
             "Q3": "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,Drer_2))));"}                         # aang 1, drer 2: median 1.5, mixed
    for f, tr in trees.items():
        (prep / "trees" / f"{f}.nhx").write_text(tr)
    sp, cl = st(), Clades()
    fams, groups, errors, skipped = run_split(prep, LS(), sp, cl, ParsimonyContext(sp, cl), NAMES)
    s = summarize_split(fams, groups, errors, skipped)
    one = s["teleost_trend"]["1"]
    assert one["clades_with_teleosts"] == 3
    assert one["trends"]["1:1"]["zebrafish_genes"] == {"0": 0, "1": 1, "2": 0, "3+": 0}
    assert one["trends"]["1:2"]["zebrafish_genes"]["2"] == 1 and one["trends"]["mixed 1-2"]["clades"] == 1
    assert "Teleost trend" in report_split_text(s, "Compara")
