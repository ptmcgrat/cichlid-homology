from helpers import FakeClades, PREFIXES, leaf, rtree, st, tree
from cichlid_homology.calls import call_gene
from cichlid_homology.consensus import label_agreement
from cichlid_homology.dating import Params
from cichlid_homology.orthology import RTree, relationship
from cichlid_homology.rooting import root_tree
from cichlid_homology.species import drop_shadowed
from cichlid_homology.tree import leaves, parse_newick, prune, to_newick

C = FakeClades()
P = Params()


def call(newick, gene="Mzeb_1", mode="denovo"):
    rt = rtree(newick)
    return call_gene(rt, leaf(rt, gene), "maylandia_zebra", P, C, mode, PREFIXES)


def test_parser_roundtrip_and_prune():
    t = parse_newick("((a:1,b:2)0.9:3,c,(d,e));")
    assert [x.name for x in leaves(t)] == ["a", "b", "c", "d", "e"]
    t2 = prune(t, {"a", "d"})
    assert [x.name for x in leaves(t2)] == ["b", "c", "e"]
    assert to_newick(parse_newick("((a,b),c);")) == "((a,b),c);"


def test_one_to_one_ortholog():
    r = call("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Mzeb_1,Acal_1)))));")
    assert r["human_relationship"] == "1:1" and r["human_closest"] == "Hsap_1"
    assert r["zebrafish_relationship"] == "1:1"
    assert r["teleost_tel_frac_one_of_present"] == 1.0


def test_retained_teleost_ohnolog_is_1_to_2_for_human_and_1_to_1_for_zebrafish():
    t = ("((Hsap_1,Mmus_1),(Locu_1,((Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))),"
         "(Drer_2,(Gacu_2,(Olat_2,(Mzeb_2,Acal_2)))))));")
    r = call(t)
    assert r["human_relationship"] == "1:many" and r["human_n_focal"] == 2
    assert r["human_dup_level"] == "teleost" and r["human_dup_class"] == "teleost_ohnolog"
    assert r["human_dup_resolved"] is True
    assert r["zebrafish_relationship"] == "1:1" and r["zebrafish_closest"] == "Drer_1"


def test_recent_duplication_is_dated_to_cichlid():
    t = ("((Hsap_1,Mmus_1),(Locu_1,(Drer_1,(Gacu_1,(Olat_1,((Mzeb_1,Acal_1),(Mzeb_2,Acal_2)))))));")
    r = call(t)
    assert r["human_relationship"] == "1:many"
    assert r["human_dup_level"] == "cichlid" and r["human_dup_class"] == "lineage_specific:cichlid"


def test_paralog_and_its_ladder_level():
    t = ("(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));")
    r = call(t, "Mzeb_1")
    assert r["human_relationship"] == "paralog" and r["human_closest_is_ortholog"] is False
    assert r["human_paralog_level"] == "gnathostome"
    r2 = call(t, "Mzeb_2")
    assert r2["human_relationship"] == "1:1"


def test_none_means_no_gene_in_this_family():
    r = call("(Locu_1,(Drer_1,(Mzeb_1,Acal_1)));")
    assert r["human_relationship"] == "none" and r["human_n_ref_in_family"] == 0


def test_unrooted_tree_is_rooted_by_reconciliation_not_midpoint():
    # trifurcating root: the only rooting without a duplication is the edge to the human/mouse clade
    root = tree("((Hsap_1,Mmus_1),Locu_1,(Drer_1,(Mzeb_1,Acal_1)));")
    rooting = root_tree(root, st(), already_rooted=False, k=3)
    assert rooting.how == "reconciliation" and rooting.cost == 0 and rooting.n_tied == 1
    top = {frozenset(x.name for x in leaves(c)) for c in rooting.roots[0].children}
    assert frozenset({"Hsap_1", "Mmus_1"}) in top


def test_recorded_root_is_kept():
    root = tree("((Hsap_1,Mmus_1),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));")
    r = root_tree(root, st(), already_rooted=True)
    assert r.how == "recorded" and r.roots[0] is root


def test_rooting_a_wrongly_rooted_tree_recovers_the_cost_zero_root():
    # same tree as the 1:1 test but drawn rooted inside the fish clade: reroot must find cost 0
    root = tree("(Mzeb_1,(Acal_1,(Drer_1,(Locu_1,(Hsap_1,Mmus_1)))));")
    r = root_tree(root, st(), already_rooted=False)
    assert r.cost == 0
    rt_cost0 = {frozenset(x.name for x in leaves(c)) for c in r.roots[0].children}
    assert frozenset({"Hsap_1", "Mmus_1"}) in rt_cost0


def test_shadowed_reference_leaf_makes_a_false_duplication_unless_dropped():
    newick = ("((Hsap_1,Mmus_1),(Locu_1,(Drer_1,((Gacu_1,GacuREF_1),(Olat_1,(Mzeb_1,Acal_1))))));")
    root = parse_newick(newick)
    sp = {x.name: ("gacu" if x.name.startswith("Gacu") else
                   {"Hsap": "homo_sapiens", "Drer": "danio_rerio", "Mzeb": "maylandia_zebra"}.get(
                       x.name.split("_")[0], x.name.split("_")[0].lower())) for x in leaves(root)}
    for x in leaves(root):
        x.sp = sp[x.name]
    assert root_tree(root, st(), already_rooted=True).cost is None     # recorded root: no scoring
    unrooted = parse_newick(newick.replace("((Hsap_1,Mmus_1),(", "((Hsap_1,Mmus_1),(", 1))
    for x in leaves(unrooted):
        x.sp = sp[x.name]
    # force scoring: pass as not rooted
    cost_with_ref = root_tree(unrooted, st(), already_rooted=False).cost
    kept, dropped = drop_shadowed(sp, {"gacu"}, lambda n: n == "Gacu_1")
    assert dropped == ["GacuREF_1"]
    pruned = prune(parse_newick(newick), set(dropped))
    for x in leaves(pruned):
        x.sp = sp[x.name]
    cost_dropped = root_tree(pruned, st(), already_rooted=False).cost
    assert cost_with_ref == 1 and cost_dropped == 0
    # and the pair is a duplication only while the reference leaf is present
    rt_bad = RTree(unrooted)
    g = leaf(rt_bad, "Gacu_1")
    assert "GacuREF_1" in [x.name for (_, _, c, spec) in rt_bad.walk(g, 0.0) if not spec
                           for x in rt_bad.leaves_in(c)]


def test_place_mode_cannot_resolve_focal_copies_on_one_branch():
    t = "((Hsap_1,Mmus_1),(Locu_1,(Drer_1,(Gacu_1,(Olat_1,Mzeb_1,Mzeb_2)))));"
    r = call(t, mode="place")
    assert r["human_relationship"] == "1:many"
    assert r["human_dup_resolved"] is False and r["human_dup_class"] == "unresolved_place_mode"
    r2 = call(t, mode="denovo")
    assert r2["human_dup_resolved"] is True


def test_cyclostome_level_is_never_used_to_date():
    class Cyc(FakeClades):
        def common_level(self, names):
            return "cyclostome"
    rt = rtree("(((Hsap_2,Mmus_2),(Locu_2,(Drer_2,(Mzeb_2,Acal_2)))),(Locu_1,(Drer_1,(Mzeb_1,Acal_1))));")
    r = call_gene(rt, leaf(rt, "Mzeb_1"), "maylandia_zebra", P, Cyc(), "denovo", PREFIXES)
    assert r["human_paralog_level"] == "cyclostome_disputed"


def test_agreement_labels():
    assert label_agreement({"TreeFam": "1:1", "Compara": "1:1"})[0] == "agree"
    assert label_agreement({"TreeFam": "1:1", "Compara": "1:many"})[0] == "split"
    assert label_agreement({"TreeFam": "1:1", "Compara": None})[0] == "single"
    assert label_agreement({"a": "1:1", "b": "1:1", "c": "paralog"})[0] == "majority"
    assert label_agreement({"a": None, "b": ""})[0] == "none"
