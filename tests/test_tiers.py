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
    r = {"p_status": "ok", "expanded_family": "False", "b_label_raw": "single", "p_label": "lost_duplicate",
         "p_gain_levels": "", "p_loss_species": "", "orig_class": "single_copy", "t_pairs_relaxed_n": ""}
    r.update(kw)
    return r


def test_unclassified():
    assert assign(base(p_status="no_anchor"))["t_call"] == "no_anchor"
    assert assign(base(expanded_family="True"))["t_call"] == "expanded_family"
    assert assign(base(p_label="multi_copy"))["t_tier"] == "unclassified"


def test_row_1_core_and_the_disputes_around_it():
    a = assign(base())
    assert (a["t_row"], a["t_tier"], a["t_call"]) == ("R1", "core", "lost_duplicate")
    assert assign(base(p_gain_levels="percomorph:1"))["t_tier"] == "disputed"        # a later duplicate in a clade
    assert assign(base(p_label="ambiguous"))["t_tier"] == "disputed"


def test_pairs_are_core_when_both_agree_and_row_3_when_zebrafish_lacks_the_copy():
    pair = base(b_label_raw="retained_pair", p_label="retained_ohnolog", orig_class="different_seqid")
    assert (assign(pair)["t_row"], assign(pair)["t_tier"]) == ("R2", "core")
    z = assign(dict(pair, p_loss_species="callorhinchus_milii;danio_rerio"))
    assert (z["t_row"], z["t_tier"]) == ("R3", "core")
    assert assign(base(b_label_raw="cichlid_lost_one_copy", p_label="retained_ohnolog"))["t_tier"] == "core"


def test_one_call_for_a_pair_needs_independent_evidence_to_be_probable():
    one = base(b_label_raw="retained_pair", p_label="lost_duplicate")
    assert assign(dict(one, orig_class="different_seqid"))["t_tier"] == "possible"       # positions alone are weak
    assert assign(dict(one, t_pairs_relaxed_n=2))["t_tier"] == "probable"                # a published pair is not
    tandem = assign(dict(one, orig_class="tandem_adjacent"))
    assert tandem["t_tier"] == "disputed" and "tandem" in tandem["t_basis"]
    assert assign(dict(one, orig_class="tandem_adjacent", t_pairs_relaxed_n=1))["t_tier"] == "probable"
    assert assign(base(b_label_raw="partial_retention", p_label="retained_ohnolog", orig_class="different_seqid"))["t_tier"] == "possible"
    assert assign(base(b_label_raw="retained_pair", p_label="lost_duplicate"))["t_tier"] == "disputed"


def test_row_5_core_needs_both_calls_and_a_cichlid_level_gain():
    c = base(b_label_raw="cichlid_specific_duplication", p_label="lost_duplicate", p_gain_levels="cichlid:1")
    assert (assign(c)["t_row"], assign(c)["t_tier"]) == ("R5", "core")
    assert assign(dict(c, p_gain_levels="cichlid(species):1"))["t_tier"] == "disputed"


def test_summary_counts_tiers_by_row():
    rows = [dict(base(), **assign(base())), dict(base(b_label_raw="retained_pair", p_label="retained_ohnolog"),
                                                    **assign(base(b_label_raw="retained_pair", p_label="retained_ohnolog"))),
            dict(base(p_status="no_anchor"), **assign(base(p_status="no_anchor")))]
    t = tiers_summary(rows)
    assert t["placed"] == 2 and t["core"] == 2 and t["unclassified"] == {"no_anchor": 1}
    assert t["by_row_and_tier"]["R1"] == {"core": 1} and t["published_pair_evidence_loaded"] is False


def test_a_later_duplicate_before_the_cichlids_needs_a_clade_gain_tandem_copies_and_no_published_pair():
    one = base(p_gain_levels="percomorph:1", orig_class="tandem_adjacent", orig_copies_n=2)
    a = assign(one)
    assert (a["t_call"], a["t_category"], a["t_tier"], a["t_when"]) == ("later_duplicate_before_cichlids", "Later duplicate, before cichlids", "possible", "percomorph")
    assert assign(dict(one, p_gain_levels="teleost:1;ovalentaria:1"))["t_when"] == "teleost;ovalentaria"
    assert assign(dict(one, orig_class="different_seqid"))["t_category"] == "Disputed"      # dispersed copies look like a pair
    assert assign(dict(one, t_pairs_relaxed_n=1))["t_category"] == "Disputed"               # a published pair
    assert assign(dict(one, orig_copies_n=1))["t_category"] == "Disputed"                   # one copy: nothing to see
    assert assign(dict(one, p_gain_levels="percomorph(species):1"))["t_category"] == "Lost duplicate"
    pair = assign(dict(one, b_label_raw="retained_pair"))                                    # a pair call overruled by tandem copies
    assert pair["t_category"] == "Later duplicate, before cichlids"
    assert assign(dict(one, p_label="ambiguous", b_label_raw="partial_retention"))["t_category"] == "Disputed"


def test_the_simplified_categories():
    pair = base(b_label_raw="retained_pair", p_label="retained_ohnolog")
    assert assign(pair)["t_category"] == "Both copies kept (best guess)"
    assert assign(dict(pair, p_loss_species="danio_rerio"))["t_category"] == "Both copies kept (best guess)"      # zebrafish lost is not its own category
    assert assign(dict(base(b_label_raw="retained_pair", p_label="lost_duplicate"), t_pairs_relaxed_n=2))["t_category"] == "Both copies kept (best guess)"
    assert assign(dict(base(b_label_raw="retained_pair", p_label="lost_duplicate"), orig_class="different_seqid"))["t_category"] == "Both copies kept (best guess)"
    assert assign(base())["t_category"] == "Lost duplicate"
    assert assign(base(b_label_raw="cichlid_specific_duplication", p_gain_levels="cichlid:1"))["t_category"] == "Later duplicate, in cichlids"
    for kw in (dict(p_status="no_anchor"), dict(expanded_family="True"), dict(p_label="multi_copy")):
        assert assign(base(**kw))["t_category"] == "Unclassified"


def test_summary_counts_the_simplified_categories():
    rows = []
    for kw in (dict(), dict(b_label_raw="retained_pair", p_label="retained_ohnolog"), dict(p_status="no_anchor"), dict(expanded_family="True")):
        r = base(**kw)
        r.update(assign(r))
        rows.append(r)
    t = tiers_summary(rows)["by_category"]
    assert t["Lost duplicate"] == {"core": 1} and t["Both copies kept (best guess)"] == {"core": 1} and t["Unclassified"] == {"unclassified": 2}


def test_a_pair_in_the_gene_tree_is_weak_evidence_for_the_possible_tier():
    one = base(b_label_raw="retained_pair", p_label="lost_duplicate", orig_class="tandem_adjacent")
    assert assign(one)["t_tier"] == "disputed"
    a = assign(dict(one, d_call="pair"))
    assert (a["t_tier"], a["t_category"]) == ("possible", "Both copies kept (best guess)") and "gene tree shows a pair" in a["t_basis"]
    assert assign(dict(one, d_call="pair", t_pairs_relaxed_n=1))["t_tier"] == "probable"     # a published pair still wins
    assert assign(dict(one, d_call="thin_pair"))["t_tier"] == "disputed"


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
