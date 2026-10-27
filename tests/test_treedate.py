from helpers import FakeClades, PREFIXES, leaf, rtree, st
from cichlid_homology.parsimony import ParsimonyContext
from cichlid_homology.profile import anchor_for
from cichlid_homology.treedate import D_COLS, treedate_profile, treedate_summary

OUT = "(Hsap_1,Locu_1)"


def early(i):
    return f"(Aang_{i},(Drer_{i},"


def fish(i, cich=None):
    c = cich or f"(Mzeb_{i},Acal_{i})"
    return f"(Aang_{i},(Drer_{i},((Trub_{i},Gacu_{i}),(Olat_{i},{c}))))"


def run(newick, gene, mode="denovo", prefixes=None):
    rt = rtree(newick)
    g = leaf(rt, gene)
    ctx = ParsimonyContext(st(), FakeClades(), 6)
    return treedate_profile(rt, g, anchor_for(rt, g, FakeClades(), None), ctx, mode, prefixes)


def test_a_lost_duplicate_is_a_speciation_at_the_teleost_root():
    r = run(f"({OUT},{fish(1)});", "Mzeb_1")
    assert r["d_status"] == "ok" and r["d_root"] == "single" and r["d_call"] == "lost_duplicate"
    assert r["d_later_levels"] == "" and set(r) >= set(D_COLS)


def test_a_pair_at_the_teleost_root_has_two_lineages_each_with_early_teleosts():
    r = run(f"({OUT},({fish(1)},{fish(2)}));", "Mzeb_1")
    assert (r["d_root"], r["d_root_level"], r["d_support"], r["d_call"]) == ("pair", "teleost", "both_early", "pair")
    assert (r["d_gene_lineage"], r["d_focal_lineages"], r["d_lineages"]) == ("A", 2, 2)
    assert run(f"({OUT},({fish(1)},{fish(2)}));", "Mzeb_2")["d_gene_lineage"] == "B"


def test_a_pair_kept_by_percomorphs_only_in_one_lineage_has_weaker_support():
    perc = "((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))"
    r = run(f"({OUT},({fish(1)},{perc}));", "Mzeb_2")
    assert (r["d_root"], r["d_dup_level"], r["d_support"], r["d_call"]) == ("pair", "teleost", "early_within", "pair")
    assert 0.5 < r["d_overlap"] < 1.0                      # all the early-teleost genes are inside one lineage: the second copy was lost there


def test_a_new_duplicate_in_the_percomorphs_is_a_later_duplicate_not_a_pair():
    p1, p2 = "((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1)))", "((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))"
    r = run(f"({OUT},(Aang_1,(Drer_1,({p1},{p2}))));", "Mzeb_1")
    assert (r["d_root"], r["d_call"], r["d_first_later_level"]) == ("single", "later_duplicate_before_cichlids", "percomorph")
    assert r["d_later_levels"].startswith("percomorph:")


def test_lineages_that_share_no_species_are_disjoint_not_a_pair():
    """Reciprocal losses and a gene tree that disagrees with the species tree look the same: no shared species."""
    a, b = "(Aang_1,(Trub_1,Mzeb_1))", "(Drer_1,(Olat_1,Acal_1))"
    r = run(f"({OUT},({a},{b}));", "Mzeb_1")
    assert (r["d_root"], r["d_support"], r["d_overlap"], r["d_call"]) == ("disjoint", "both_early", 0.0, "lost_duplicate")
    assert r["d_focal_lineages"] == 1                         # M. zebra itself holds one of the two lineages


def test_a_cichlid_duplicate_lists_the_species_that_carry_both_copies():
    r = run(f"({OUT},{fish(1, '((Mzeb_1,Acal_1),(Mzeb_2,Acal_2))')});", "Mzeb_1")
    assert (r["d_call"], r["d_first_later_level"], r["d_carriers"], r["d_carrier_n"]) == ("later_duplicate_in_cichlids", "cichlid", "acal;maylandia_zebra", 2)
    one = run(f"({OUT},{fish(1, '((Mzeb_1,Acal_1),Mzeb_2)')});", "Mzeb_1")
    assert one["d_carriers"] == "maylandia_zebra"             # Acanthopterygii partner lost: only M. zebra carries both


def test_a_duplicate_inside_one_species_is_species_level_and_not_a_later_duplicate():
    r = run(f"({OUT},{fish(1, '((Mzeb_1,Mzeb_2),Acal_1)')});", "Mzeb_1")
    assert r["d_species_dups"] == 1 and r["d_first_later_level"] == "" and r["d_call"] == "lost_duplicate"


def test_other_paralog_lineages_with_their_tetrapod_genes_do_not_get_in_the_way():
    """A clade holding two 2R paralog lineages, each with a tetrapod gene: g's own teleost clade is what is dated."""
    rt = rtree(f"(({OUT},{fish(1)}),((Hsap_2,Locu_2),{fish(2)}));")
    g = leaf(rt, "Mzeb_1")
    ctx = ParsimonyContext(st(), FakeClades(), 6)
    r = treedate_profile(rt, g, ("human", rt.root), ctx, "denovo", None)
    assert r["d_status"] == "ok" and r["d_root"] == "single" and r["d_call"] == "lost_duplicate"


def test_a_gene_with_no_teleost_neighbour_is_not_dated():
    r = run(f"({OUT},(Locu_2,Mzeb_1));", "Mzeb_1")
    assert r["d_status"] == "one_teleost_gene" and r["d_call"] == ""


def test_in_a_placement_tree_project_genes_are_assigned_to_a_lineage_but_do_not_create_duplications():
    pl = f"({OUT},((Aang_1,(Drer_1,Mzeb_1)),(Aang_2,(Drer_2,Mzeb_2))));"
    r = run(pl, "Mzeb_1", "place", PREFIXES)
    assert (r["d_call"], r["d_support"], r["d_gene_lineage"], r["d_focal_lineages"]) == ("pair", "both_early", "A", 2)
    assert run(pl, "Mzeb_2", "place", PREFIXES)["d_gene_lineage"] == "B"
    two = f"({OUT},(Aang_1,(Drer_1,(Mzeb_1,Mzeb_2))));"          # two placed genes on one branch: nothing to date
    r2 = run(two, "Mzeb_1", "place", PREFIXES)
    assert r2["d_call"] == "lost_duplicate" and r2["d_later_levels"] == ""


def test_summary_counts_calls_and_support():
    rows = [dict(run(f"({OUT},({fish(1)},{fish(2)}));", "Mzeb_1"), t_category="Both copies kept (best guess)"),
            dict(run(f"({OUT},{fish(1)});", "Mzeb_1"), t_category="Lost duplicate"),
            dict(run(f"({OUT},(Locu_2,Mzeb_1));", "Mzeb_1"))]
    s = treedate_summary(rows)
    assert s["scored"] == 2 and s["calls"] == {"pair": 1, "lost_duplicate": 1} and s["pair_support"] == {"both_early": 1}
    assert s["pairs_with_focal_species_in_two_or_more_lineages"] == 1
    assert s["by_category"]["Lost duplicate"] == {"lost_duplicate": 1}


def test_the_overlap_separates_a_retained_pair_from_reciprocal_losses_and_from_a_stray_gene():
    pair = run(f"({OUT},({fish(1)},{fish(2)}));", "Mzeb_1")
    assert pair["d_overlap"] == 1.0
    recip = run(f"({OUT},((Aang_1,(Trub_1,Mzeb_1)),(Drer_1,(Olat_1,Acal_1))));", "Mzeb_1")
    assert recip["d_overlap"] == 0.0
    # a percomorph-only duplicate with one stray early-teleost gene (Aang_2) inside a percomorph clade: a false pair at the root
    p1, p2 = "((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1),Aang_2))", "((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))"
    stray = run(f"({OUT},(Aang_1,(Drer_1,({p1},{p2}))));", "Mzeb_1")
    assert stray["d_root"] == "thin"                                          # the root has almost nothing to support it
    assert (stray["d_call"], stray["d_dup_level"]) == ("later_duplicate_before_cichlids", "percomorph")   # the stray does not make it a teleost pair
    # a stray early-teleost gene of another species (Drer_2) does not make a pair: it only lowers the level of the later duplicate
    p1b = "((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1),Drer_2))"
    weak = run(f"({OUT},(Aang_1,(Drer_1,({p1b},{p2}))));", "Mzeb_1")
    assert weak["d_root"] == "single" and weak["d_first_later_level"] == "percomorph"


def test_the_overlap_cutoff_is_a_parameter():
    one = run(f"({OUT},({fish(1)},((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))));", "Mzeb_2")
    assert one["d_root"] == "pair"                                         # overlap 5 / 7
    rt = rtree(f"({OUT},({fish(1)},((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))));")
    g = leaf(rt, "Mzeb_2")
    ctx = ParsimonyContext(st(), FakeClades(), 6)
    strict = treedate_profile(rt, g, anchor_for(rt, g, FakeClades(), None), ctx, "denovo", None, 0.9)
    assert strict["d_root"] == "thin" and strict["d_call"] == "thin_pair"


def test_a_node_where_the_gene_tree_disagrees_with_the_species_tree_is_not_a_later_duplication():
    """(Trub,Mzeb) against (Gacu,Acal): reconciliation says duplication, but no species has both: nothing was duplicated."""
    scrambled = "(((Trub_1,Mzeb_1),(Gacu_1,Acal_1)),Olat_1)"
    r = run(f"({OUT},(Aang_1,(Drer_1,{scrambled})));", "Mzeb_1")
    assert r["d_later_levels"] == "" and r["d_call"] == "lost_duplicate"
    real = "(((Trub_1,Mzeb_1),(Trub_2,Mzeb_2)),Olat_1)"                 # the same shape with species in both halves: a duplication
    r2 = run(f"({OUT},(Aang_1,(Drer_1,{real})));", "Mzeb_1")
    assert r2["d_call"] == "later_duplicate_before_cichlids" and r2["d_first_later_level"] == "percomorph"


def test_a_real_pair_below_a_top_that_disagrees_with_the_species_tree_is_found():
    """Eel (Aang) is single and sits above the pair of (zebrafish + percomorphs) lineages: the pair node is below the root."""
    p1, p2 = "((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1)))", "((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))"
    r = run(f"({OUT},(Aang_1,((Drer_1,{p1}),(Drer_2,{p2}))));", "Mzeb_1")
    assert (r["d_root"], r["d_call"], r["d_dup_level"], r["d_dup_at_root"]) == ("single", "pair", "teleost", "no")
    assert (r["d_support"], r["d_overlap"], r["d_focal_lineages"]) == ("both_early", 1.0, 2)


def test_a_duplicate_of_percomorphs_alone_is_unresolved_when_no_early_teleost_is_in_the_clade():
    p1, p2 = "((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1)))", "((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))"
    r = run(f"({OUT},({p1},{p2}));", "Mzeb_1")
    assert (r["d_call"], r["d_dup_level"], r["d_support"], r["d_dup_at_root"]) == ("pair_unresolved_depth", "percomorph", "none_early", "yes")
