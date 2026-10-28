from cichlid_homology.ohnologs import load_pair_sets, support
from cichlid_homology.tiers import assign, tiers_summary

HEAD = "Ohno1\tOhno2\tSymbol1\tSymbol2\tOg weighted global mean\tSelf weighted global mean\tOutgroup Support\tDuplication time\tGene type\n"


def write_pairs(d):
    rows = {"drerio": [("phf12a", "phf12b", "FishWGD", "protein_coding"), ("mir1", "mir2", "FishWGD", "miRNA"),
                       ("x1", "x2", "Vertebrata", "protein_coding"), ("clstn2", "CLSTN2 (1 to many)", "FishWGD", "protein_coding")],
            "olatipes": [("capns1a", "capns1b", "FishWGD", "protein_coding")]}
    for sp, lst in rows.items():
        for lvl, keep in (("Strict", lst[:1]), ("Relaxed", lst)):
            (d / f"{sp}_Pairs_{lvl}_3R.txt").write_text(HEAD + "".join(
                f"ENS1\tENS2\t{a}\t{b}\t1e-9\t1e-5\t7\t{t}\t{g}\n" for a, b, t, g in keep))


def test_the_pair_files_keep_only_protein_coding_fish_genome_duplication_pairs_with_real_symbols(tmp_path):
    write_pairs(tmp_path)
    s = load_pair_sets(tmp_path)
    assert set(s) == {"zebrafish", "medaka"}
    assert s["zebrafish"]["strict"] == {"phf12a", "phf12b"} and s["zebrafish"]["relaxed"] == {"phf12a", "phf12b", "clstn2"}


def test_support_counts_the_species_whose_published_pairs_hold_the_gene(tmp_path):
    write_pairs(tmp_path)
    s = load_pair_sets(tmp_path)
    row = {"a_zebrafish_best_symbol": "PHF12A", "a_zebrafish_best_name_reliable": "True", "own_symbol": "capns1b"}
    assert support(row, s) == {"t_pairs_relaxed_n": 2, "t_pairs_strict_n": 2, "t_pairs_species": "zebrafish;medaka"}
    unreliable = dict(row, a_zebrafish_best_name_reliable="False", own_symbol="LOC123")
    assert support(unreliable, s)["t_pairs_relaxed_n"] == 0
    assert support(row, {}) == {"t_pairs_relaxed_n": "", "t_pairs_strict_n": "", "t_pairs_species": ""}


def base(**kw):
    """A gene all three calls call lost: thresholds single, parsimony lost, the tree a lost duplicate."""
    r = {"p_status": "ok", "expanded_family": "False", "b_label_raw": "single", "p_label": "lost_duplicate", "d_status": "ok",
         "d_call": "lost_duplicate", "p_gain_levels": "", "p_loss_species": "", "orig_class": "single_copy", "orig_copies_n": 1,
         "t_pairs_relaxed_n": ""}
    r.update(kw)
    return r


def pair(**kw):
    return base(b_label_raw="retained_pair", p_label="retained_ohnolog", d_call="pair", orig_class="different_seqid", orig_copies_n=2, **kw)


def test_unclassified():
    assert assign(base(p_status="no_anchor"))["t_call"] == "no_anchor"
    assert assign(base(expanded_family="True"))["t_call"] == "expanded_family"
    assert assign(base(p_label="multi_copy"))["t_category"] == "Unclassified"
    assert assign(base(p_status="no_anchor"))["t_tier"] == "unclassified"


def test_three_votes_make_a_core_call():
    a = assign(base())
    assert (a["t_category"], a["t_tier"], a["t_call"], a["t_row"]) == ("Lost duplicate", "core", "lost_duplicate", "R1")
    assert (a["t_vote_thr"], a["t_vote_par"], a["t_vote_tree"]) == ("lost", "lost", "lost")
    b = assign(pair())
    assert (b["t_category"], b["t_tier"], b["t_call"]) == ("Both copies kept (best guess)", "core", "retained_pair")
    assert assign(pair(p_loss_species="danio_rerio"))["t_category"] == "Both copies kept (best guess)"     # zebrafish lost is no category
    assert assign(dict(pair(), b_label_raw="cichlid_lost_one_copy"))["t_tier"] == "core"                    # the teleost ancestor kept both


def test_two_of_three_is_probable_whether_the_third_abstains_or_disagrees():
    assert assign(dict(pair(), d_call="thin_pair"))["t_tier"] == "probable"                    # the tree abstains
    a = assign(dict(pair(), d_call="lost_duplicate"))
    assert (a["t_category"], a["t_tier"]) == ("Both copies kept (best guess)", "probable") and "tree says lost" in a["t_basis"]
    assert assign(dict(pair(), p_label="lost_duplicate"))["t_tier"] == "probable"
    assert assign(dict(pair(), b_label_raw="single"))["t_tier"] == "probable"
    lost = assign(dict(base(), d_call="thin_pair"))
    assert (lost["t_category"], lost["t_tier"]) == ("Lost duplicate", "probable")
    assert "thresholds says a pair" in assign(dict(base(), d_call="pair"))["t_basis"] or assign(dict(base(), d_call="pair"))["t_category"] == "Lost duplicate"
    assert assign(dict(base(), d_call="pair"))["t_tier"] == "probable"                         # the tree alone is outvoted


def test_no_majority_is_disputed():
    a = assign(dict(base(), p_label="retained_ohnolog", d_call="thin_pair"))                   # one pair, one lost, one abstains
    assert (a["t_category"], a["t_tier"], a["t_row"]) == ("Disputed", "disputed", "other") and "no majority" in a["t_basis"]
    b = assign(dict(base(), b_label_raw="partial_retention", p_label="ambiguous", d_call="thin_pair"))
    assert b["t_category"] == "Disputed"


def test_an_undated_gene_can_only_reach_probable():
    a = assign(dict(pair(), d_status="no_teleost", d_call=""))
    assert (a["t_category"], a["t_tier"]) == ("Both copies kept (best guess)", "probable")


def test_the_published_pairs_are_not_part_of_the_calls():
    one = assign(dict(pair(), d_call="lost_duplicate", p_label="lost_duplicate"))             # one pair vote against two lost
    assert one["t_category"] == "Lost duplicate"
    assert assign(dict(pair(), d_call="lost_duplicate", p_label="lost_duplicate", t_pairs_relaxed_n=4))["t_category"] == "Lost duplicate"


def test_a_later_duplicate_in_the_cichlids_needs_two_calls():
    c = base(b_label_raw="cichlid_specific_duplication", p_gain_levels="cichlid:1")
    a = assign(c)
    assert (a["t_category"], a["t_tier"], a["t_row"]) == ("Later duplicate, in cichlids", "core", "R5")
    t = assign(dict(c, p_label="ambiguous", d_call="later_duplicate_in_cichlids"))
    assert (t["t_category"], t["t_tier"]) == ("Later duplicate, in cichlids", "probable")
    assert assign(dict(c, p_gain_levels="cichlid(species):1"))["t_category"] == "Lost duplicate"     # one call is not enough


def test_a_later_duplicate_before_the_cichlids_needs_a_clade_level_duplication_tandem_copies_and_no_published_pair():
    one = base(p_gain_levels="percomorph:1", orig_class="tandem_adjacent", orig_copies_n=2)
    a = assign(one)
    assert (a["t_call"], a["t_category"], a["t_tier"], a["t_when"]) == ("later_duplicate_before_cichlids", "Later duplicate, before cichlids", "possible", "percomorph")
    both = assign(dict(one, d_call="later_duplicate_before_cichlids", d_dup_level="percomorph"))
    assert both["t_tier"] == "probable" and both["t_when"] == "percomorph"
    tree = assign(dict(base(), d_call="later_duplicate_before_cichlids", d_dup_level="ovalentaria", orig_class="tandem_near", orig_copies_n=3))
    assert (tree["t_category"], tree["t_when"]) == ("Later duplicate, before cichlids", "ovalentaria")
    assert assign(dict(one, orig_class="different_seqid"))["t_category"] == "Lost duplicate"       # dispersed copies look like a pair
    assert assign(dict(one, t_pairs_relaxed_n=1))["t_category"] == "Lost duplicate"                # a published pair: an ohnolog
    assert assign(dict(one, orig_copies_n=1))["t_category"] == "Lost duplicate"
    assert assign(dict(pair(), p_gain_levels="percomorph:1", orig_class="tandem_adjacent"))["t_category"] == "Both copies kept (best guess)"   # a pair majority wins


def test_a_many_copy_gene_with_a_pair_at_the_teleost_root_is_a_pair_but_a_thin_pair_is_not():
    many = base(p_label="multi_copy", b_label_raw="multi_copy")
    assert assign(many)["t_category"] == "Unclassified" and assign(many)["t_call"] == "many_copies"
    a = assign(dict(many, d_call="pair"))
    assert (a["t_category"], a["t_tier"]) == ("Both copies kept (best guess)", "possible") and "three or more copies" in a["t_basis"]
    assert assign(dict(many, d_call="thin_pair"))["t_category"] == "Unclassified"
    assert assign(dict(many, d_call="pair", expanded_family="True"))["t_category"] == "Unclassified"   # expanded families stay out


def test_summary_counts_categories_tiers_and_votes():
    rows = []
    for g in (base(), pair(), dict(pair(), d_call="lost_duplicate"), base(p_status="no_anchor"), base(expanded_family="True"),
              dict(base(), p_label="retained_ohnolog", d_call="thin_pair")):
        r = dict(g)
        r.update(assign(r))
        rows.append(r)
    t = tiers_summary(rows)
    assert t["placed"] == 4 and t["core"] == 2 and t["probable"] == 1 and t["disputed"] == 1
    assert t["by_category"]["Lost duplicate"] == {"core": 1} and t["by_category"]["Both copies kept (best guess)"] == {"core": 1, "probable": 1}
    assert t["unclassified"] == {"no_anchor": 1, "expanded_family": 1} and t["votes_for_a_pair"] == {0: 1, 1: 1, 2: 1, 3: 1}


def test_the_pair_files_may_use_dots_or_underscores_and_a_wrong_folder_is_reported(tmp_path):
    from cichlid_homology.ohnologs import find_files
    d = tmp_path / "dots"
    d.mkdir()
    for lvl in ("Strict", "Relaxed"):
        (d / f"drerio.Pairs.{lvl}.3R.txt").write_text(HEAD + "E1\tE2\tphf12a\tphf12b\t1e-9\t1e-5\t7\tFishWGD\tprotein_coding\n")
    (d / "notes.txt").write_text("not a pair file")
    assert set(find_files(d)) == {("drerio", "Strict"), ("drerio", "Relaxed")}
    assert load_pair_sets(d)["zebrafish"]["strict"] == {"phf12a", "phf12b"}
    other = tmp_path / "empty"
    other.mkdir()
    (other / "readme.txt").write_text("x")
    assert load_pair_sets(other) == {}
    try:
        find_files(tmp_path / "missing")
        raise AssertionError("a missing folder must raise")
    except FileNotFoundError:
        pass
