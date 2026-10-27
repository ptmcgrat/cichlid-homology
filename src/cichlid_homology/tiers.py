"""Confidence tiers for goal B1: for each gene, which of the histories it fits and how sure we are.

Two independent-looking calls are combined, both read from the same clade copy numbers: the share-based call
(`b_label`, thresholds 0.5 and 0.3) and the parsimony call (`p_label`, gain 2 loss 1; calibrated against published
pairs in 0.9.1: the thresholds gave 84% recall at 0.8% false positives, parsimony 82% at 5.8%, and requiring both
gave 79% at 0.7%). Rows follow the seven histories of docs/NOTES.md:

    R1  a duplicate lost early                    R2  both copies kept (zebrafish has both)
    R3  both kept, zebrafish lost one             R2/3/4  a pair whose depth the counts cannot fix (rows 3 and 4 give
    R5  a later duplicate in cichlids                      identical counts; rows 6 and 7 cannot be told from them)

Tiers: `core` (both calls agree); `probable` (exactly one calls a pair and the gene is in a published pair, which is
independent evidence); `possible` (exactly one calls a pair and either the copies are on different sequences, which is weak
evidence: among genes that are not pairs, 45 to 58% of those with two copies are also on different sequences, or the
gene tree shows a pair, `d_call` = pair: two lineages at the teleost level sharing half their species, which on the genes
where the counts disagree was in a published pair 2.5 times as often as for genes with no evidence);
`disputed` (the calls conflict and nothing settles it, or the copies are tandem); `unclassified` (not scored, expanded
family, or many copies). No share cutoffs are added here.

The simplified category (`t_category`) groups these: Lost duplicate; Both copies kept (best guess) (R2, R3, probable and
possible, with the confidence in `t_tier`); Later duplicate, in cichlids (R5); Later duplicate, before cichlids (a
duplication after the teleost root but before the cichlids, not an ohnolog); Disputed; Unclassified (one group, the
reason in `t_call`). A later duplicate before the cichlids needs parsimony to place a clade-level gain at the
teleost, percomorph or ovalentarian level, the gene to have two or more copies that are tandem, and no published
pair: counts alone cannot tell such a duplicate from a teleost-duplication pair kept in the percomorphs, so the
rule is strict and the group is small. `t_when` gives the level of the gain."""
from __future__ import annotations

from collections import Counter

TIER_COLS = ["t_category", "t_row", "t_call", "t_tier", "t_basis", "t_positions", "t_when"]
CATEGORIES = ("Lost duplicate", "Both copies kept (best guess)", "Later duplicate, in cichlids",
              "Later duplicate, before cichlids", "Disputed", "Unclassified")
CATEGORY_OF = {"lost_duplicate": CATEGORIES[0], "retained_pair": CATEGORIES[1], "cichlid_specific_duplicate": CATEGORIES[2],
               "later_duplicate_before_cichlids": CATEGORIES[3], "unresolved": CATEGORIES[4]}
LADDER = ("teleost", "percomorph", "ovalentaria")
CLADE_LEVELS = ("teleost", "percomorph", "ovalentaria", "cichlid")
TIERS = ("core", "probable", "possible", "disputed", "unclassified")
PAIR_THR = ("retained_pair", "cichlid_lost_one_copy")      # the teleost ancestor kept both copies


def _later_levels(row) -> set:
    out = set()
    for t in str(row.get("p_gain_levels") or "").split(";"):
        if t:
            lv = t.rpartition(":")[0]
            if lv in CLADE_LEVELS:
                out.add(lv)
    return out


def _positions(row) -> str:
    c = row.get("orig_class") or ""
    if c == "different_seqid":
        return "dispersed"
    if c in ("tandem_adjacent", "tandem_near"):
        return "tandem"
    return ""


def _n(v) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def _before_cichlids(row, cl, p, pos, pub):
    older = [x for x in LADDER if x in cl]
    ok = p == "lost_duplicate" and older and _n(row.get("orig_copies_n")) >= 2 and pos == "tandem" and not pub
    return ";".join(older) if ok else ""


def assign(row: dict) -> dict:
    out = _assign(row)
    out["t_category"] = "Unclassified" if out["t_tier"] == "unclassified" else CATEGORY_OF.get(out["t_call"], "Disputed")
    out.setdefault("t_when", "")
    return out


def _assign(row: dict) -> dict:
    thr = row.get("b_label_raw") or row.get("b_label") or ""
    p = row.get("p_label") or ""
    pos = _positions(row)
    out = {"t_positions": pos}
    if str(row.get("p_status")) != "ok":
        return {**out, "t_row": "unclassified", "t_call": "no_anchor", "t_tier": "unclassified", "t_basis": "no outgroup anchor in the tree"}
    if str(row.get("expanded_family")) == "True":
        return {**out, "t_row": "unclassified", "t_call": "expanded_family", "t_tier": "unclassified", "t_basis": "family with many focal genes"}
    if thr == "multi_copy" or p == "multi_copy":
        return {**out, "t_row": "unclassified", "t_call": "many_copies", "t_tier": "unclassified", "t_basis": "three or more copies in many teleosts"}
    pub = _n(row.get("t_pairs_relaxed_n"))
    cl = _later_levels(row)
    pair_t, pair_p = thr in PAIR_THR, p == "retained_ohnolog"
    zloss = "danio_rerio" in str(row.get("p_loss_species") or "").split(";")
    when = _before_cichlids(row, cl, p, pos, pub)

    def res(r, call, tier, basis):
        return {**out, "t_row": r, "t_call": call, "t_tier": tier, "t_basis": basis}

    def before():
        return {**out, "t_row": "R4", "t_call": "later_duplicate_before_cichlids", "t_tier": "possible", "t_when": when,
                "t_basis": f"parsimony places a clade-level duplication at {when.replace(';', ' and ')}; copies are tandem and in no published pair"}

    if pair_t and pair_p:
        return res("R3" if zloss else "R2", "retained_pair", "core", "both calls: a pair" + ("; zebrafish lacks the second copy" if zloss else ""))
    if pair_t or pair_p:
        who = "thresholds" if pair_t else "parsimony"
        other = f"parsimony says {p or 'nothing'}" if pair_t else f"thresholds say {thr or 'nothing'}"
        if pub:
            return res("R2/3/4", "retained_pair", "probable", f"{who} say a pair, {other}; in a published pair in {pub} of 4 species")
        tree = str(row.get("d_call")) == "pair"
        if pos == "dispersed" or tree:
            ev = " and ".join(x for x, y in (("copies on different sequences", pos == "dispersed"), ("the gene tree shows a pair", tree)) if y)
            return res("R2/3/4", "retained_pair", "possible", f"{who} say a pair, {other}; {ev} (weak evidence)")
        if when:
            return before()
        why = "copies are tandem and in no published pair" if pos == "tandem" else "no independent evidence"
        return res("R2/3/4", "unresolved", "disputed", f"{who} say a pair, {other}; {why}")
    if thr == "single":
        if p == "lost_duplicate" and not cl:
            return res("R1", "lost_duplicate", "core", "both calls: lost early, no later duplicate in a clade")
        if when:
            return before()
        why = "parsimony sees a later duplicate in a clade" if p == "lost_duplicate" else f"parsimony says {p or 'nothing'}"
        return res("R1", "unresolved", "disputed", f"thresholds say single; {why}")
    if thr == "cichlid_specific_duplication":
        if p == "lost_duplicate" and "cichlid" in cl:
            return res("R5", "cichlid_specific_duplicate", "core", "both calls: a later duplicate in the cichlids")
        if when:
            return before()
        return res("R5", "unresolved", "disputed", f"thresholds say a cichlid duplicate; parsimony says {p or 'nothing'}")
    if when:
        return before()
    if p == "lost_duplicate" and "cichlid" in cl:
        return res("R5", "unresolved", "disputed", f"parsimony sees a cichlid duplicate; thresholds say {thr or 'nothing'}")
    return res("other", "unresolved", "disputed", f"thresholds say {thr or 'nothing'}; parsimony says {p or 'nothing'}")


def tiers_summary(rows) -> dict:
    C = Counter
    scored = [r for r in rows if r.get("t_tier")]
    by = C((r["t_row"], r["t_tier"]) for r in scored)
    rows_in = sorted({k[0] for k in by})
    out = {"genes": len(scored), "by_row_and_tier": {r: {t: by[(r, t)] for t in TIERS if by[(r, t)]} for r in rows_in}}
    placed = [r for r in scored if r["t_tier"] != "unclassified"]
    out["placed"] = len(placed)
    out["core"] = sum(1 for r in placed if r["t_tier"] == "core")
    out["probable"] = sum(1 for r in placed if r["t_tier"] == "probable")
    out["possible"] = sum(1 for r in placed if r["t_tier"] == "possible")
    out["disputed"] = sum(1 for r in placed if r["t_tier"] == "disputed")
    out["unclassified"] = {k: v for k, v in C(r["t_call"] for r in scored if r["t_tier"] == "unclassified").items()}
    out["by_category"] = {c: {t: n for t in TIERS if (n := sum(1 for r in scored if r["t_category"] == c and r["t_tier"] == t))}
                          for c in CATEGORIES}
    pair_like = [r for r in placed if r["t_row"] in ("R2", "R3", "R2/3/4")]
    have = [r for r in pair_like if str(r.get("t_pairs_relaxed_n")) not in ("", "None")]
    out["published_pair_evidence_loaded"] = bool(have)
    if have:
        out["pair_like_genes_in_a_published_pair_pct"] = round(100 * sum(1 for r in have if _n(r.get("t_pairs_relaxed_n")) >= 1) / len(have), 1)
    return out
