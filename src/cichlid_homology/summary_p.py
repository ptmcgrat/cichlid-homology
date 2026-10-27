"""Summary of the parsimony labels (lost duplicate / retained ohnolog): label counts, confidence, how they
compare with the threshold-based B1 labels, two checks that use nothing from the labels (zebrafish a/b
names, positions of the copies), the events of the best histories, and agreement between the methods."""
from __future__ import annotations

import collections
import re

from .consensus import label_agreement

LABELS = ("lost_duplicate", "retained_ohnolog", "multi_copy", "ambiguous")
B1 = ("single", "retained_pair", "cichlid_specific_duplication", "cichlid_lost_one_copy", "partial_retention",
      "multi_copy", "insufficient_species", "no_anchor")
_AB = re.compile(r"[a-z0-9\-]{3,}[ab]")


def _pct(k, n):
    return round(100 * k / n, 1) if n else None


def _i(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def _ab(sym):
    s = sym or ""
    return bool(_AB.fullmatch(s)) and not s.startswith(("si:", "zgc"))


def parsimony_summary(rows):
    C = collections.Counter
    ok = [r for r in rows if r.get("p_status") == "ok"]
    out = {"genes": len(rows), "status": dict(C(r.get("p_status", "") for r in rows)), "scored": len(ok)}
    if not ok:
        return out
    out["labels"] = {k: {"n": sum(1 for r in ok if r["p_label"] == k), "pct": _pct(sum(1 for r in ok if r["p_label"] == k), len(ok))}
                     for k in LABELS}
    m = C(min(_i(r["p_margin"]), 6) for r in ok)
    out["margin_0_1_2_3_4_5_6plus"] = [m.get(i, 0) for i in range(7)]
    out["robust_to_cost_changes_pct"] = _pct(sum(1 for r in ok if r["p_robust"] in (True, "True")), len(ok))
    both = [r for r in ok if r["p_zero_modes_agree"] in (True, False, "True", "False")]
    out["zero_as_loss_vs_missing_agree_pct"] = _pct(sum(1 for r in both if r["p_zero_modes_agree"] in (True, "True")), len(both))
    out["zero_modes_compared"] = len(both)
    out["vs_threshold_B1_labels"] = {k: dict(C(r.get("b_label_raw") or r.get("b_label") for r in ok if r["p_label"] == k).most_common())
                                     for k in LABELS}
    out["by_human_group"] = {g: dict(C(r["p_label"] for r in ok if r["a_human_group"] == g))
                             for g in ("OneToOne", "OneToTwo", "OneToMany", "Other", "None")}
    # check 1: the zebrafish a/b names (ZFIN marks duplicated pairs this way), zebrafish OneToOne, reliable names, small clades
    base = [r for r in ok if r["a_zebrafish_group"] == "OneToOne" and r.get("a_zebrafish_best_symbol")
            and str(r.get("a_zebrafish_best_name_reliable")) == "True" and _i(r["a_zebrafish_clade_leaves"]) <= 100]
    out["check_zebrafish_ab_names_pct"] = {k: {"genes": sum(1 for r in base if r["p_label"] == k),
                                               "pct_a_or_b": _pct(sum(1 for r in base if r["p_label"] == k and _ab(r["a_zebrafish_best_symbol"])),
                                                                  sum(1 for r in base if r["p_label"] == k))} for k in LABELS[:3]}
    # check 2: where the nearest other copy of the gene sits
    multi = [r for r in ok if _i(r.get("orig_copies_n")) >= 2]
    out["check_copy_positions"] = {}
    for k in LABELS[:3]:
        sub = [r for r in multi if r["p_label"] == k]
        c = C(r.get("orig_class") for r in sub)
        out["check_copy_positions"][k] = {"genes": len(sub),
                                         "tandem_adjacent_or_near_pct": _pct(c["tandem_adjacent"] + c["tandem_near"], len(sub)),
                                         "different_sequence_pct": _pct(c["different_seqid"], len(sub))}
    # events of the best histories
    gains = C()
    for r in ok:
        for tok in (r.get("p_gain_levels") or "").split(";"):
            if tok:
                lvl, _, n = tok.rpartition(":")
                gains[lvl] += _i(n)
    out["gain_events_by_level_all_scored_genes"] = dict(gains.most_common())
    ret = [r for r in ok if r["p_label"] == "retained_ohnolog"]
    loss = C(s for r in ret for s in (r.get("p_loss_species") or "").split(";") if s)
    out["retained_ohnologs_with_a_species_level_loss"] = _pct(sum(1 for r in ret if r.get("p_loss_species")), len(ret))
    out["species_most_often_lost_in_retained_ohnologs"] = dict(loss.most_common(12))
    out["retained_ohnologs_total"] = len(ret)
    # the three outcomes of B1: a later duplicate shows as gains below the teleost ancestor
    later = {}
    for k in LABELS[:2]:
        sub = [r for r in ok if r["p_label"] == k]
        with_gain = [r for r in sub if _i(r.get("p_gains_below_T")) > 0]
        lv = C()
        for r in with_gain:
            for tok in (r.get("p_gain_levels") or "").split(";"):
                if tok and not tok.startswith("teleost_ancestor"):
                    lv[tok.rpartition(":")[0]] += 1
        later[k] = {"genes": len(sub), "with_a_later_gain": len(with_gain), "pct": _pct(len(with_gain), len(sub)),
                    "levels_of_the_later_gains": dict(lv.most_common(8))}
    out["later_duplicates_after_the_teleost_ancestor"] = later
    return out


def parsimony_agreement(by, genes):
    A, B = "TreeFam", "Compara"
    labs, pairs = collections.Counter(), collections.Counter()
    for g in genes:
        va = by[A][g].get("p_label") if by[A][g].get("p_status") == "ok" else None
        vb = by[B][g].get("p_label") if by[B][g].get("p_status") == "ok" else None
        va = None if va == "ambiguous" else va
        vb = None if vb == "ambiguous" else vb
        lab, _, _ = label_agreement({A: va, B: vb})
        labs[lab] += 1
        pairs[(va, vb)] += 1
    both = labs["agree"] + labs["split"]
    return {"labels": dict(labs), "agree_pct_when_both_called": _pct(labs["agree"], both),
            "pairs_TreeFam_Compara": [[list(map(str, k)), v] for k, v in pairs.most_common(12)]}
