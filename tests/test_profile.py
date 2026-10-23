from helpers import FakeClades, PREFIXES, leaf, rtree
from cichlid_homology.profile import clade_profile, label

C = FakeClades()


def prof(newick, gene="Mzeb_1"):
    rt = rtree(newick)
    return clade_profile(rt, leaf(rt, gene), C, 0.0)


FISH1 = "(Aang_1,(Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))))"
FISH2 = "(Aang_2,(Drer_2,(Gacu_2,(Olat_2,(Mzeb_2,Acal_2)))))"


def test_dominant_single_copy():
    r = prof(f"((Hsap_1,Mmus_1),(Locu_1,{FISH1}));")
    assert r["b_anchor"] == "human" and r["b_label"] == "single"
    assert r["b_teleost_1"] == 2 and r["b_cichlid_1"] == 2 and r["b_teleost_2"] == 0


def test_retained_pair_counts_two_in_every_species():
    r = prof(f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},{FISH2})));")
    assert r["b_label"] == "retained_pair"
    assert r["b_teleost_2"] == 2 and r["b_percomorph_2"] == 1 and r["b_ovalentaria_2"] == 1 and r["b_cichlid_2"] == 2


def test_the_clade_rooted_at_zebrafish_would_miss_the_second_copy_but_human_does_not():
    # the smallest clade with g and the zebrafish gene holds only one copy; the human-rooted clade holds both
    r = prof(f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},{FISH2})));")
    assert r["b_clade_leaves"] > 12


def test_cichlid_specific_duplication():
    t = "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Gacu_1,(Olat_1,((Mzeb_1,Acal_1),(Mzeb_2,Acal_2))))))));"
    r = prof(t)
    assert r["b_label"] == "cichlid_specific_duplication"


def test_cichlid_lost_one_copy():
    t = f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},(Aang_2,(Drer_2,(Gacu_2,Olat_2))))));"
    r = prof(t)
    assert r["b_label"] == "cichlid_lost_one_copy"


def test_falls_back_to_another_non_teleost_when_there_is_no_human():
    r = prof(f"(Locu_1,{FISH1});")
    assert r["b_anchor"] == "other_gnathostome"


def test_paralog_anchor_is_not_used():
    # the only human genes are paralogs of g: human is skipped, the gar gene (an ortholog) is used
    t = "(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));"
    r = prof(t)
    assert r["b_anchor"] == "other_gnathostome"


def test_no_anchor_at_all():
    r = prof("(Drer_1,(Mzeb_1,Acal_1));")
    assert r["b_label"] == "no_orthologous_anchor" and r["b_anchor"] == "none"


def test_label_needs_enough_species():
    row = {f"b_{lv}_{k}": 0 for lv in ("teleost", "percomorph", "ovalentaria", "cichlid") for k in ("n", "1", "2", "3p")}
    row["b_cichlid_n"], row["b_cichlid_1"] = 2, 2
    assert label(row) == "insufficient_species"
