from helpers import FakeClades, PREFIXES, leaf, rtree, st, tree
from cichlid_homology.calls import call_gene
from cichlid_homology.dating import Params
from cichlid_homology.evidence import explain_gene
from cichlid_homology.pilot import choose_genes
from cichlid_homology.rooting import root_tree
from cichlid_homology.tree import leaves

C = FakeClades()
P = Params()


def test_losses_break_ties_between_rootings_with_equal_duplications():
    r = root_tree(tree("(Mzeb_1,Mzeb_2,Drer_1);"), st(), False, k=3)
    assert r.cost == 1 and r.n_tied_dup == 3 and r.n_tied == 1
    top = {frozenset(x.name for x in leaves(c)) for c in r.roots[0].children}
    assert frozenset({"Mzeb_1", "Mzeb_2"}) in top      # the duplication stays inside the cichlids


def test_outgroup_is_counted_only_in_the_nearest_clade_outside_the_duplication():
    # the ohnolog pair has one gar copy next to it; a second gar copy belongs to an unrelated older
    # paralog far away and must not make the outgroup "multi"
    t = ("(((Hsap_1,Mmus_1),(Locu_1,((Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))),"
         "(Drer_2,(Gacu_2,(Olat_2,(Mzeb_2,Acal_2))))))),(Hsap_9,(Locu_9,Locu_8)));")
    rt = rtree(t)
    r = call_gene(rt, leaf(rt, "Mzeb_1"), "maylandia_zebra", P, C, "denovo", PREFIXES)
    assert r["human_dup_class"] == "teleost_ohnolog" and r["human_dup_outgroup_found"] is True


def test_teleost_pattern_has_a_fraction_of_all_species_and_a_coverage():
    rt = rtree("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Mzeb_1,Acal_1)))));")
    r = call_gene(rt, leaf(rt, "Mzeb_1"), "maylandia_zebra", P, C, "denovo", PREFIXES)
    assert r["teleost_tel_frac_one_of_all"] == 1.0 and r["teleost_tel_coverage"] == 1.0


def test_explain_names_speciation_and_duplication():
    rt = rtree("((Hsap_1,Mmus_1),(Locu_1,((Drer_1,(Mzeb_1,Acal_1)),(Drer_2,(Mzeb_2,Acal_2)))));")
    text = "\n".join(explain_gene(rt, leaf(rt, "Mzeb_1"), "homo_sapiens", "maylandia_zebra", 0.0))
    assert "DUPLICATION" in text and "SPECIATION" in text and "relationship 1:many" in text


def test_choose_genes_is_deterministic_and_bounded():
    ids = [f"Mzeb_{i}" for i in range(100)]
    a, b = choose_genes(ids, 10, seed=1), choose_genes(ids, 10, seed=1)
    assert a == b and len(a) == 10 and len(choose_genes(ids, 500)) == 100
