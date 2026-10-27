from helpers import FakeClades, st
from cichlid_homology.parsimony import ParsimonyContext, classify, trans_matrix, INF

OUT = ["homo_sapiens", "mmus", "locu", "ecal", "cmil"]
TEL = ["aang", "danio_rerio", "trub", "gacu", "olat", "maylandia_zebra", "acal"]
SCOPE = frozenset(OUT + TEL)


def ctx():
    return ParsimonyContext(st(), FakeClades(), kcap=6)


def counts(**kw):
    d = {s: 1 for s in SCOPE}
    for k, v in kw.items():
        d[{"drer": "danio_rerio", "mzeb": "maylandia_zebra", "human": "homo_sapiens"}.get(k, k)] = v
    return d


def test_transition_costs_forbid_birth_from_nothing():
    t = trans_matrix(4, 1, 2)
    assert t[0][1] == INF and t[0][0] == 0
    assert t[1][3] == 2 and t[3][1] == 4 and t[2][2] == 0


def test_all_single_copy_is_a_lost_duplicate():
    r = ctx().score(counts(), SCOPE)
    assert r["p_status"] == "ok" and r["p_label"] == "lost_duplicate" and r["p_cost1"] == 0
    assert r["p_margin"] > 0 and r["p_robust"] is True and r["p_gains_n"] == 0


def test_every_teleost_with_two_copies_is_a_retained_ohnolog_with_one_event():
    c = {s: (2 if s in TEL else 1) for s in SCOPE}
    r = ctx().score(c, SCOPE)
    assert r["p_label"] == "retained_ohnolog" and r["p_cost2"] == 1 and r["p_cost1"] == 2 and r["p_margin"] == 1
    assert r["p_gain_levels"] == "teleost_ancestor:1" and r["p_losses_n"] == 0


def test_with_thin_early_teleost_sampling_one_loss_ties_with_separate_gains():
    # only aang and danio_rerio branch before the percomorphs: one loss in zebrafish (2 events) costs the same as
    # a gain in aang plus a gain in the percomorph ancestor (2 events). The method says so instead of guessing.
    c = {s: (2 if s in TEL else 1) for s in SCOPE}
    c["danio_rerio"] = 1
    r = ctx().score(c, SCOPE)
    assert r["p_cost1"] == r["p_cost2"] == 2 and r["p_label"] == "ambiguous" and r["p_margin"] == 0


def test_with_denser_sampling_a_zebrafish_loss_is_found():
    from helpers import LEVEL
    from cichlid_homology.speciestree import SpeciesTree

    class Clades(FakeClades):
        def level(self, s):
            return {"pkin": "teleost", "amex": "teleost", "char": "teleost"}.get(s) or LEVEL[s]

    newick = ("(bflo,(cmil,((homo_sapiens,mmus),((locu,ecal),(pkin,(aang,(amex,(danio_rerio,(char,"
              "((trub,gacu),(olat,(maylandia_zebra,acal))))))))))));")
    tel = TEL + ["pkin", "amex", "char"]
    scope = frozenset(OUT + tel)
    c = {s: (2 if s in tel else 1) for s in scope}
    c["danio_rerio"] = 1
    r = ParsimonyContext(SpeciesTree(newick), Clades(), 6).score(c, scope)
    assert r["p_label"] == "retained_ohnolog" and r["p_robust"] is True
    assert r["p_loss_species"] == "danio_rerio" and r["p_losses_below_T"] == 1 and r["p_cost2"] == 2 and r["p_cost1"] == 3


def test_a_cichlid_only_duplicate_is_one_event_on_the_cichlid_branch_and_the_teleost_call_stays_single():
    c = {s: 1 for s in SCOPE}
    c["maylandia_zebra"] = c["acal"] = 2
    r = ctx().score(c, SCOPE)
    assert r["p_label"] == "lost_duplicate" and r["p_cost1"] == 1
    assert r["p_gain_levels"] == "cichlid:1" and r["p_gains_below_T"] == 1


def test_a_single_species_duplicate_is_one_species_level_gain():
    c = {s: 1 for s in SCOPE}
    c["maylandia_zebra"] = 2
    r = ctx().score(c, SCOPE)
    assert r["p_label"] == "lost_duplicate" and r["p_gain_levels"] == "cichlid(species):1"


def test_many_copies_everywhere_is_multi_copy():
    c = {s: (3 if s in TEL else 1) for s in SCOPE}
    assert ctx().score(c, SCOPE)["p_label"] == "multi_copy"


def test_half_the_teleosts_with_two_copies_is_ambiguous_or_cost_dependent():
    c = {s: 1 for s in SCOPE}
    for s in ("aang", "danio_rerio", "trub", "gacu"):
        c[s] = 2
    r = ctx().score(c, SCOPE)
    assert r["p_label"] == "ambiguous" or r["p_margin"] <= 2


def test_zero_counts_as_loss_or_as_missing_can_differ():
    # a retained pair seen only in two teleosts: with zeros as losses it needs many losses
    c = {s: 0 for s in SCOPE}
    c.update({"homo_sapiens": 1, "mmus": 1, "locu": 1, "ecal": 1, "cmil": 1, "danio_rerio": 2, "aang": 2})
    r = ctx().score(c, SCOPE)
    assert r["p_label_missing"] == "retained_ohnolog"
    assert r["p_zero_modes_agree"] in (True, False, "")


def test_no_outgroup_and_too_few_teleosts():
    only_tel = frozenset(TEL)
    assert ctx().score({s: 1 for s in only_tel}, only_tel)["p_status"] == "no_outgroup"
    few = frozenset(OUT + ["danio_rerio"])
    assert ctx().score({s: 1 for s in few}, few)["p_status"] == "insufficient_teleosts"


def test_classify_ties_are_ambiguous():
    assert classify({1: 3, 2: 3, 3: 9}) == ("ambiguous", 0)
    assert classify({1: 0, 2: 4, 3: 9}) == ("lost_duplicate", 4)
