"""Confidence tiers for goal B1: for each gene, which of the histories it fits and how sure we are.

Two independent-looking calls are combined, both read from the same clade copy numbers: the share-based call
(`b_label`, thresholds 0.5 and 0.3) and the parsimony call (`p_label`, gain 2 loss 1; calibrated against published
pairs in 0.9.1: the thresholds gave 84% recall at 0.8% false positives, parsimony 82% at 5.8%, and requiring both
gave 79% at 0.7%). Rows follow the seven histories of docs/NOTES.md:

    R1  a duplicate lost early                    R2  both copies kept (zebrafish has both)
    R3  both kept, zebrafish lost one             R2/3/4  a pair whose depth the counts cannot fix (rows 3 and 4 give
    R5  a later duplicate in cichlids                      identical counts; rows 6 and 7 cannot be told from them)

Tiers: `core` (both calls agree); `probable` (exactly one calls a pair and the gene is in a published pair, which is
independent evidence); `possible` (exactly one calls a pair and the copies are on different sequences, which is weak
evidence: among genes that are not pairs, 45 to 58% of those with two copies are also on different sequences);
`disputed` (the calls conflict and nothing settles it, or the copies are tandem); `outside` (not scored, expanded
family, or many copies). No share cutoffs are added here."""
from __future__ import annotations

from collections import Counter

TIER_COLS = ["t_row", "t_call", "t_tier", "t_basis", "t_positions"]
CLADE_LEVELS = ("teleost", "percomorph", "ovalentaria", "cichlid")
TIERS = ("core", "probable", "possible", "disputed", "outside")
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


def assign(row: dict) -> dict:
    thr = row.get("b_label_raw") or row.get("b_label") or ""
    p = row.get("p_label") or ""
    pos = _positions(row)
    out = {"t_positions": pos}
    if str(row.get("p_status")) != "ok":
        return {**out, "t_row": "outside", "t_call": "no_anchor", "t_tier": "outside", "t_basis": "no outgroup anchor in the tree"}
    if str(row.get("expanded_family")) == "True":
        return {**out, "t_row": "outside", "t_call": "expanded_family", "t_tier": "outside", "t_basis": "family with many focal genes"}
    if thr == "multi_copy" or p == "multi_copy":
        return {**out, "t_row": "outside", "t_call": "many_copies", "t_tier": "outside", "t_basis": "three or more copies in many teleosts"}
    pub = _n(row.get("t_pairs_relaxed_n"))
    cl = _later_levels(row)
    pair_t, pair_p = thr in PAIR_THR, p == "retained_ohnolog"
    zloss = "danio_rerio" in str(row.get("p_loss_species") or "").split(";")

    def res(r, call, tier, basis):
        return {**out, "t_row": r, "t_call": call, "t_tier": tier, "t_basis": basis}

    if pair_t and pair_p:
        return res("R3" if zloss else "R2", "retained_pair", "core", "both calls: a pair" + ("; zebrafish lacks the second copy" if zloss else ""))
    if pair_t or pair_p:
        who = "thresholds" if pair_t else "parsimony"
        other = f"parsimony says {p or 'nothing'}" if pair_t else f"thresholds say {thr or 'nothing'}"
        if pub:
            return res("R2/3/4", "retained_pair", "probable", f"{who} say a pair, {other}; in a published pair in {pub} of 4 species")
        if pos == "dispersed":
            return res("R2/3/4", "retained_pair", "possible", f"{who} say a pair, {other}; copies on different sequences (weak evidence)")
        why = "copies are tandem and in no published pair" if pos == "tandem" else "no independent evidence"
        return res("R2/3/4", "unresolved", "disputed", f"{who} say a pair, {other}; {why}")
    if thr == "single":
        if p == "lost_duplicate" and not cl:
            return res("R1", "lost_duplicate", "core", "both calls: lost early, no later duplicate in a clade")
        why = "parsimony sees a later duplicate in a clade" if p == "lost_duplicate" else f"parsimony says {p or 'nothing'}"
        return res("R1", "unresolved", "disputed", f"thresholds say single; {why}")
    if thr == "cichlid_specific_duplication":
        if p == "lost_duplicate" and "cichlid" in cl:
            return res("R5", "cichlid_specific_duplicate", "core", "both calls: a later duplicate in the cichlids")
        return res("R5", "unresolved", "disputed", f"thresholds say a cichlid duplicate; parsimony says {p or 'nothing'}")
    if p == "lost_duplicate" and "cichlid" in cl:
        return res("R5", "unresolved", "disputed", f"parsimony sees a cichlid duplicate; thresholds say {thr or 'nothing'}")
    return res("other", "unresolved", "disputed", f"thresholds say {thr or 'nothing'}; parsimony says {p or 'nothing'}")


def tiers_summary(rows) -> dict:
    C = Counter
    scored = [r for r in rows if r.get("t_tier")]
    by = C((r["t_row"], r["t_tier"]) for r in scored)
    rows_in = sorted({k[0] for k in by})
    out = {"genes": len(scored), "by_row_and_tier": {r: {t: by[(r, t)] for t in TIERS if by[(r, t)]} for r in rows_in}}
    placed = [r for r in scored if r["t_tier"] != "outside"]
    out["placed"] = len(placed)
    out["core"] = sum(1 for r in placed if r["t_tier"] == "core")
    out["probable"] = sum(1 for r in placed if r["t_tier"] == "probable")
    out["possible"] = sum(1 for r in placed if r["t_tier"] == "possible")
    out["disputed"] = sum(1 for r in placed if r["t_tier"] == "disputed")
    out["outside"] = {k: v for k, v in C(r["t_call"] for r in scored if r["t_tier"] == "outside").items()}
    pair_like = [r for r in placed if r["t_row"] in ("R2", "R3", "R2/3/4")]
    have = [r for r in pair_like if str(r.get("t_pairs_relaxed_n")) not in ("", "None")]
    out["published_pair_evidence_loaded"] = bool(have)
    if have:
        out["pair_like_genes_in_a_published_pair_pct"] = round(100 * sum(1 for r in have if _n(r.get("t_pairs_relaxed_n")) >= 1) / len(have), 1)
    return out
