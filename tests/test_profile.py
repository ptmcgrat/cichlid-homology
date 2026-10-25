from helpers import FakeClades, leaf, rtree
from cichlid_homology.dating import Params
from cichlid_homology.profile import clade_profile, label, state

C = FakeClades()
P = Params()


def prof(newick, gene="Mzeb_1"):
    rt = rtree(newick)
    return clade_profile(rt, leaf(rt, gene), C, 0.0, P)


FISH1 = "(Aang_1,(Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))))"
FISH2 = "(Aang_2,(Drer_2,(Gacu_2,(Olat_2,(Mzeb_2,Acal_2)))))"


def test_dominant_single_copy():
    r = prof(f"((Hsap_1,Mmus_1),(Locu_1,{FISH1}));")
    assert r["b_anchor"] == "human" and r["b_label"] == "single"
    assert r["b_pattern"] == "...."
    assert r["b_other_share2p"] == 0 and r["b_cichlid_share2p"] == 0


def test_retained_pair_counts_two_in_every_species():
    r = prof(f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},{FISH2})));")
    assert r["b_label"] == "retained_pair" and r["b_pattern"] == "++++"
    assert r["b_teleost_2"] == 2 and r["b_percomorph_2"] == 1 and r["b_cichlid_2"] == 2


def test_the_human_rooted_clade_holds_both_copies():
    r = prof(f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},{FISH2})));")
    assert r["b_clade_leaves"] > 12


def test_cichlid_specific_duplication_and_the_b2_carriers():
    t = "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Gacu_1,(Olat_1,((Mzeb_1,Acal_1),(Mzeb_2,Acal_2))))))));"
    r = prof(t)
    assert r["b_label"] == "cichlid_specific_duplication" and r["b_pattern"] == "...+"
    assert r["b_cichlid_with2p"] == 2 and r["b_cichlid_n"] == 2
    assert r["b_cichlid_extra_species"] == "acal;maylandia_zebra"


def test_three_or_more_cichlid_copies_count_as_carrying_the_duplicate():
    t = ("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Gacu_1,(Olat_1,"
         "((Mzeb_1,Acal_1),((Mzeb_2,Acal_2),(Mzeb_3,Acal_3)))))))));")
    r = prof(t)
    assert r["b_cichlid_with3p"] == 2 and r["b_label"] == "cichlid_specific_duplication"


def test_cichlid_lost_one_copy():
    t = f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},(Aang_2,(Drer_2,(Gacu_2,Olat_2))))));"
    assert prof(t)["b_label"] == "cichlid_lost_one_copy"


def test_falls_back_to_another_non_teleost_when_there_is_no_human():
    assert prof(f"(Locu_1,{FISH1});")["b_anchor"] == "other_gnathostome"


def test_paralog_anchor_is_not_used():
    t = "(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));"
    assert prof(t)["b_anchor"] == "other_gnathostome"


def test_no_anchor_at_all():
    r = prof("(Drer_1,(Mzeb_1,Acal_1));")
    assert r["b_label"] == "no_orthologous_anchor" and r["b_anchor"] == "none"


def test_states_follow_the_thresholds():
    assert state(0.0, P) == "lost" and state(0.29, P) == "lost"
    assert state(0.3, P) == "partial" and state(0.49, P) == "partial"
    assert state(0.5, P) == "kept" and state(1.0, P) == "kept" and state("", P) == "none"


def row(o_share, c_share, o3=0.0, n_o=6, n_c=4):
    return {"b_other_n": n_o, "b_cichlid_n": n_c, "b_other_share2p": o_share, "b_other_share3p": o3,
            "b_cichlid_share2p": c_share, "b_other_state": state(o_share, P), "b_cichlid_state": state(c_share, P)}


def test_labels_for_every_combination():
    assert label(row(0.0, 0.0), P) == "single"
    assert label(row(0.0, 0.6), P) == "cichlid_specific_duplication"
    assert label(row(0.0, 0.4), P) == "cichlid_specific_duplication"          # partial: B2 counts the carriers
    assert label(row(0.9, 0.9), P) == "retained_pair"
    assert label(row(0.9, 0.4), P) == "retained_pair"
    assert label(row(0.9, 0.0), P) == "cichlid_lost_one_copy"
    assert label(row(0.4, 0.9), P) == "partial_retention"
    assert label(row(0.9, 0.9, o3=0.5), P) == "multi_copy"


def test_label_needs_enough_species():
    assert label(row(0.0, 0.0, n_o=2), P) == "insufficient_species"
    assert label(row(0.0, 0.0, n_c=0), P) == "insufficient_species"
