from helpers import FakeClades, PREFIXES, leaf, rtree, tree
from cichlid_homology.calls import call_gene
from cichlid_homology.clade_counts import CATS, category, clade_counts
from cichlid_homology.dating import Params
from cichlid_homology.orthology import RTree
from cichlid_homology.profile import clade_profile_v2

OWN = "maylandia_zebra"
HUM, ZEB = "homo_sapiens", "danio_rerio"


def counts(newick, ref, gene="Mzeb_1"):
    rt = rtree(newick)
    return clade_counts(rt, leaf(rt, gene), ref, OWN)


FISH1 = "(Aang_1,(Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))))"
FISH2 = "(Aang_2,(Drer_2,(Gacu_2,(Olat_2,(Mzeb_2,Acal_2)))))"


def test_category_table():
    assert [category(r, f) for r, f in [(0, 1), (1, 1), (1, 2), (1, 3), (1, 9), (2, 1), (5, 1), (2, 2), (3, 4)]] == [
        "None", "OneToOne", "OneToTwo", "OneToMany", "OneToMany", "ManyToOne", "ManyToOne", "ManyToMany", "ManyToMany"]
    assert set(CATS) == {"OneToOne", "OneToTwo", "OneToMany", "ManyToOne", "ManyToMany", "None"}


def test_one_to_one_for_human_and_zebrafish():
    t = "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Mzeb_1,Acal_1)))));"
    h, z = counts(t, HUM), counts(t, ZEB)
    assert (h["cat"], h["ref_n"], h["focal_n"], h["best"]) == ("OneToOne", 1, 1, "Hsap_1")
    assert (z["cat"], z["clade_leaves"]) == ("OneToOne", 3)


def test_retained_pair_is_one_to_two_for_human_and_one_to_one_for_zebrafish():
    t = f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},{FISH2})));"
    assert counts(t, HUM)["cat"] == "OneToTwo" and counts(t, ZEB)["cat"] == "OneToOne"
    assert counts(t, ZEB, "Mzeb_2")["best"] == "Drer_2"


def test_one_to_many_counts_only_the_focal_species():
    t = ("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Gacu_1,(Olat_1,"
         "((Mzeb_1,Acal_1),((Mzeb_2,Acal_2),(Mzeb_3,Acal_3)))))))));")
    h = counts(t, HUM)
    assert h["cat"] == "OneToMany" and h["focal_n"] == 3          # the three Acal genes are not counted


def test_many_to_one_and_many_to_many():
    assert counts("((Hsap_1,Hsap_2),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));", HUM)["cat"] == "ManyToOne"
    assert counts("((Hsap_1,Hsap_2),(Locu_1,(Drer_1,((Mzeb_1,Acal_1),(Mzeb_2,Acal_2)))));", HUM)["cat"] == "ManyToMany"


def test_none_when_the_tree_has_no_gene_of_the_reference_species():
    h = counts("(Locu_1,(Drer_1,(Mzeb_1,Acal_1)));", HUM)
    assert h["cat"] == "None" and h["best"] == "" and h["ref_n"] == 0


def test_best_gene_is_the_shortest_path_and_ties_are_flagged():
    t = ("((Hsap_1:0.5,Hsap_2:0.1):0.2,(Locu_1:0.1,(Drer_1:0.1,(Mzeb_1:0.1,Acal_1:0.1):0.1):0.1):0.2);")
    h = counts(t, HUM)
    assert h["best"] == "Hsap_2" and h["tie"] is False
    t2 = counts("((Hsap_1,Hsap_2),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));", HUM)
    assert t2["best"] == "Hsap_1" and t2["tie"] is True            # no lengths: tie, broken by name


def test_an_old_paralog_is_just_a_larger_clade_now():
    t = "(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));"
    h = counts(t, HUM)
    assert h["cat"] == "OneToTwo" and h["clade_leaves"] == 10 and h["best"] == "Hsap_2"


def test_one_stray_species_no_longer_demotes_an_ortholog():
    # Aang has a gene on both sides of the node that joins Mzeb to zebrafish. The overlap rule called
    # zebrafish a paralog; the count does not care.
    t = "((Hsap_1,Mmus_1),(Locu_1,((Drer_1,Aang_1),((Mzeb_1,Acal_1),Aang_2))));"
    rt = rtree(t)
    row = call_gene(rt, leaf(rt, "Mzeb_1"), OWN, Params(), FakeClades(), "denovo", PREFIXES)
    assert row["zebrafish_relationship"] == "paralog"                 # the old rule
    assert row["a_zebrafish_cat"] == "OneToOne" and row["a_zebrafish_best"] == "Drer_1"


def test_goal_b_uses_the_human_clade_without_an_orthology_test():
    # the human gene is a (distant) paralog, so the old profile skipped it; the human clade is used now
    t = "(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));"
    rt = rtree(t)
    g = leaf(rt, "Mzeb_1")
    node = clade_counts(rt, g, HUM, OWN)["node"]
    r = clade_profile_v2(rt, g, FakeClades(), Params(), node)
    assert r["b_anchor"] == "human" and r["b_clade_leaves"] == 10


def test_goal_b_falls_back_to_another_non_teleost_and_to_nothing():
    rt = rtree(f"(Locu_1,{FISH1});")
    assert clade_profile_v2(rt, leaf(rt, "Mzeb_1"), FakeClades(), Params(), None)["b_anchor"] == "other_gnathostome"
    rt2 = rtree("(Drer_1,(Mzeb_1,Acal_1));")
    r = clade_profile_v2(rt2, leaf(rt2, "Mzeb_1"), FakeClades(), Params(), None)
    assert r["b_anchor"] == "none" and r["b_label"] == "no_anchor"


def test_goal_b_labels_on_the_human_clade():
    t = f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},{FISH2})));"
    rt = rtree(t)
    g = leaf(rt, "Mzeb_1")
    r = clade_profile_v2(rt, g, FakeClades(), Params(), clade_counts(rt, g, HUM, OWN)["node"])
    assert r["b_label"] == "retained_pair" and r["b_pattern"] == "++++"


def test_groups_follow_the_owners_scheme():
    from cichlid_homology.clade_counts import GROUPS, group
    assert [group("human", c) for c in CATS] == ["OneToOne", "OneToTwo", "OneToMany", "Other", "Other", "None"]
    assert [group("zebrafish", c) for c in CATS] == ["OneToOne", "OneToMany", "OneToMany", "Other", "Other", "None"]
    assert GROUPS["zebrafish"] == ("OneToOne", "OneToMany", "Other", "None")


def test_clade_counts_keeps_the_focal_leaves():
    t = ("((Hsap_1,Mmus_1),(Locu_1,(Drer_1,((Mzeb_1,Acal_1),(Mzeb_2,Acal_2)))));")
    h = counts(t, HUM)
    assert h["focal_leaves"] == ["Mzeb_1", "Mzeb_2"]
