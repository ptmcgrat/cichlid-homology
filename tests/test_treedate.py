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
    assert (r["d_root"], r["d_root_level"], r["d_support"], r["d_call"]) == ("pair", "teleost", "one_early", "pair")
    assert 0.5 < r["d_overlap"] < 1.0


def test_a_new_duplicate_in_the_percomorphs_is_a_later_duplicate_not_a_pair():
    p1, p2 = "((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1)))", "((Trub_2,Gacu_2),(Olat_2,(Mzeb_2,Acal_2)))"
    r = run(f"({OUT},(Aang_1,(Drer_1,({p1},{p2}))));", "Mzeb_1")
    assert (r["d_root"], r["d_call"], r["d_first_later_level"]) == ("single", "later_duplicate_before_cichlids", "percomorph")
    assert r["d_later_levels"].startswith("percomorph:")


def test_reciprocal_losses_make_a_pair_although_every_species_has_one_copy():
    a, b = "(Aang_1,(Trub_1,Mzeb_1))", "(Drer_1,(Olat_1,Acal_1))"
    r = run(f"({OUT},({a},{b}));", "Mzeb_1")
    assert (r["d_root"], r["d_support"], r["d_call"]) == ("pair", "both_early", "pair")
    assert r["d_focal_lineages"] == 1                         # M. zebra itself holds one of the two lineages


def test_a_cichlid_duplicate_lists_the_species_that_carry_both_copies():
    r = run(f"({OUT},{fish(1, '((Mzeb_1,Acal_1),(Mzeb_2,Acal_2))')});", "Mzeb_1")
    assert (r["d_call"], r["d_first_later_level"], r["d_carriers"], r["d_carrier_n"]) == ("later_duplicate_in_cichlids", "cichlid", "acal;maylandia_zebra", 2)
    one = run(f"({OUT},{fish(1, '((Mzeb_1,Acal_1),Mzeb_2)')});", "Mzeb_1")
    assert one["d_carriers"] == "maylandia_zebra"             # Acanthopterygii partner lost: only M. zebra carries both


def test_a_duplicate_inside_one_species_is_species_level_and_not_a_later_duplicate():
    r = run(f"({OUT},{fish(1, '((Mzeb_1,Mzeb_2),Acal_1)')});", "Mzeb_1")
    assert r["d_species_dups"] == 1 and r["d_first_later_level"] == "" and r["d_call"] == "lost_duplicate"


def test_teleosts_that_are_not_one_clade_are_reported_not_dated():
    r = run("((Hsap_1,(Aang_1,Locu_1)),(Drer_1,((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1)))));", "Mzeb_1")
    assert r["d_status"] == "teleosts_not_monophyletic" and r["d_call"] == ""


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
            dict(run("((Hsap_1,(Aang_1,Locu_1)),(Drer_1,(Mzeb_1,Acal_1)));", "Mzeb_1"))]
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
    assert stray["d_root"] == "pair" and stray["d_overlap"] < 0.2             # the pair is there, with almost nothing to support it
    # a stray early-teleost gene of another species (Drer_2) does not make a pair: it only lowers the level of the later duplicate
    p1b = "((Trub_1,Gacu_1),(Olat_1,(Mzeb_1,Acal_1),Drer_2))"
    weak = run(f"({OUT},(Aang_1,(Drer_1,({p1b},{p2}))));", "Mzeb_1")
    assert weak["d_root"] == "single" and weak["d_first_later_level"] == "teleost"
