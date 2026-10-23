from helpers import FakeClades, PREFIXES, leaf, rtree
from cichlid_homology.calls import call_gene
from cichlid_homology.dating import Params
from cichlid_homology.evidence import explain_gene

C = FakeClades()
P = Params()
OWN = "maylandia_zebra"


def test_explain_prints_the_ladder_level_and_the_species_that_sets_it():
    t = "(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));"
    rt = rtree(t)
    text = "\n".join(explain_gene(rt, leaf(rt, "Mzeb_1"), "homo_sapiens", OWN, 0.0, C))
    assert "level gnathostome" in text and "set by" in text and "DUPLICATION" in text


def test_paralog_overlap_level_is_reported_beside_the_all_species_level():
    t = "(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));"
    rt = rtree(t)
    r = call_gene(rt, leaf(rt, "Mzeb_1"), OWN, P, C, "denovo", PREFIXES)
    assert r["human_paralog_level"] == "gnathostome"
    assert r["human_paralog_level_overlap"] == "gnathostome"      # locu is shared by both sides


def test_overlap_level_can_be_younger_than_the_all_species_level():
    # an old outgroup (Bflo) sits inside one side only: the all-species level is the outermost
    # level present, the overlap level only counts species that both sides hold
    t = ("(((Hsap_2,Bflo_2),(Drer_2,(Mzeb_2,Acal_2))),(Drer_1,(Mzeb_1,Acal_1)));")
    rt = rtree(t)
    r = call_gene(rt, leaf(rt, "Mzeb_1"), OWN, P, C, "denovo", PREFIXES)
    assert r["human_paralog_level"] == "outside_vertebrates"
    assert r["human_paralog_level_overlap"] == "teleost"


def test_reference_only_levels_ignore_our_placed_genes():
    t = ("((Hsap_1,Mmus_1),(Locu_1,((Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))),"
         "(Drer_2,(Gacu_2,(Olat_2,(Mzeb_2,Acal_2)))))));")
    rt = rtree(t)
    r = call_gene(rt, leaf(rt, "Mzeb_1"), OWN, P, C, "denovo", PREFIXES)
    a, b = map(int, r["human_dup_levels_sides"].split("/"))
    ra, rb = map(int, r["human_dup_ref_levels_sides"].split("/"))
    assert (a, b) == (4, 4)
    assert ra <= a and rb <= b        # reference-only counts can only be equal or smaller
