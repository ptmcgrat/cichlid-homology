"""Categories and confidence for goal B1, from the votes of three calls (0.13.0).

Three calls read the same clade of the same family in different ways:
  thresholds  the share of non-cichlid teleosts with two or more copies (`b_label`; cutoffs 0.5 and 0.3)
  parsimony   copy-number changes on the species tree (`p_label`; a gain costs 2, a loss 1)
  tree        the gene tree reconciled with the species tree: the deepest duplication with half-shared species (`d_call`)
Each votes `pair` (the teleost ancestor kept both copies of the teleost-genome duplication), `lost` (one copy) or abstains.
On Compara, scored against 2,637 genes in published teleost-duplication pairs and 8,192 genes with no evidence of one:
3 of 3 votes for a pair held 1,877 positives and 34 negatives, 2 of 3 held 365 and 27, 1 of 3 held 203 and 427, none 192
and 7,704. The published pairs are no longer part of the calls: they stay in the table as an outside check
(`t_pairs_*`).

Categories (`t_category`): Lost duplicate; Both copies kept (best guess) (R2, R3 and R4 of docs/NOTES.md cannot be told
apart by counts); Later duplicate, in cichlids; Later duplicate, before cichlids (after the teleost root, before the
cichlids: not an ohnolog); Disputed (no majority); Unclassified (one group, the reason in `t_call`).
Confidence (`t_tier`): `core` = all three calls agree (for later duplicates in cichlids, the two count calls agree);
`probable` = two of the three agree; `possible` = a later duplicate before the cichlids on one line of evidence, or a pair in
a many-copy family on the tree alone; `disputed` = no majority; `unclassified`.
A later duplicate before the cichlids needs parsimony or the tree to place a clade-level duplication at the teleost,
percomorph or ovalentarian level, two or more M. zebra copies that are tandem, and no published pair (an ohnolog): counts
alone cannot tell such a duplicate from a teleost pair kept in the percomorphs, so the rule is strict and the group is small."""
from __future__ import annotations

from collections import Counter

TIER_COLS = ["t_category", "t_row", "t_call", "t_tier", "t_basis", "t_vote_thr", "t_vote_par", "t_vote_tree",
             "t_positions", "t_when"]
CATEGORIES = ("Lost duplicate", "Both copies kept (best guess)", "Later duplicate, in cichlids",
              "Later duplicate, before cichlids", "Disputed", "Unclassified")
CATEGORY_OF = {"lost_duplicate": CATEGORIES[0], "retained_pair": CATEGORIES[1], "cichlid_specific_duplicate": CATEGORIES[2],
               "later_duplicate_before_cichlids": CATEGORIES[3], "unresolved": CATEGORIES[4]}
LADDER = ("teleost", "percomorph", "ovalentaria")
CLADE_LEVELS = ("teleost", "percomorph", "ovalentaria", "cichlid")
TIERS = ("core", "probable", "possible", "disputed", "unclassified")
PAIR_THR = ("retained_pair", "cichlid_lost_one_copy")      # the teleost ancestor kept both copies
LOST_THR = ("single", "cichlid_specific_duplication")      # one copy at the teleost ancestor


def _n(v) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


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


def votes(row) -> tuple:
    thr = row.get("b_label_raw") or row.get("b_label") or ""
    v_thr = "pair" if thr in PAIR_THR else ("lost" if thr in LOST_THR else "")
    v_par = {"retained_ohnolog": "pair", "lost_duplicate": "lost"}.get(row.get("p_label") or "", "")
    v_tree = {"pair": "pair", "lost_duplicate": "lost"}.get(str(row.get("d_call") or ""), "") if str(row.get("d_status")) == "ok" else ""
    return v_thr, v_par, v_tree


def _names(v, what):
    return ", ".join(n for n, x in zip(("thresholds", "parsimony", "tree"), v) if x == what)


def assign(row: dict) -> dict:
    out = _assign(row)
    out["t_category"] = "Unclassified" if out["t_tier"] == "unclassified" else CATEGORY_OF.get(out["t_call"], "Disputed")
    out.setdefault("t_when", "")
    return out


def _assign(row: dict) -> dict:
    thr = row.get("b_label_raw") or row.get("b_label") or ""
    p = row.get("p_label") or ""
    pos = _positions(row)
    v = votes(row)
    out = {"t_positions": pos, "t_vote_thr": v[0], "t_vote_par": v[1], "t_vote_tree": v[2]}

    def res(r, call, tier, basis, when=""):
        return {**out, "t_row": r, "t_call": call, "t_tier": tier, "t_basis": basis, "t_when": when}

    if str(row.get("p_status")) != "ok":
        return res("unclassified", "no_anchor", "unclassified", "no outgroup anchor in the tree")
    if str(row.get("expanded_family")) == "True":
        return res("unclassified", "expanded_family", "unclassified", "family with many focal genes")
    d_call = str(row.get("d_call") or "") if str(row.get("d_status")) == "ok" else ""
    if thr == "multi_copy" or p == "multi_copy":
        if d_call == "pair":           # the counts see an expansion; the gene tree still shows a pair at the teleost root
            return res("R2/3/4", "retained_pair", "possible", "three or more copies in many teleosts, but the gene tree shows a pair at the teleost root")
        return res("unclassified", "many_copies", "unclassified", "three or more copies in many teleosts")
    npair, nlost = v.count("pair"), v.count("lost")
    cl = _later_levels(row)
    pub = _n(row.get("t_pairs_relaxed_n"))
    if npair < 2:
        # a later duplicate in the cichlids: thresholds and parsimony agree on it, or the tree and one of them
        thr_c, par_c = thr == "cichlid_specific_duplication", p == "lost_duplicate" and "cichlid" in cl
        tree_c = d_call == "later_duplicate_in_cichlids"
        if thr_c + par_c + tree_c >= 2:
            tier = "core" if (thr_c and par_c) else "probable"
            who = ", ".join(n for n, x in (("thresholds", thr_c), ("parsimony", par_c), ("tree", tree_c)) if x)
            return res("R5", "cichlid_specific_duplicate", tier, f"a later duplicate in the cichlids ({who})")
        # a later duplicate before the cichlids: a clade-level duplication at the teleost, percomorph or ovalentarian level
        tree_level = str(row.get("d_dup_level") or "") if d_call == "later_duplicate_before_cichlids" else ""
        par_older = [x for x in LADDER if x in cl] if p == "lost_duplicate" else []
        older = [x for x in LADDER if x in par_older or x == tree_level]
        if older and _n(row.get("orig_copies_n")) >= 2 and pos == "tandem" and not pub:
            tier = "probable" if (par_older and tree_level) else "possible"
            src = "parsimony and the tree" if tier == "probable" else ("the tree" if tree_level else "parsimony")
            return res("R4", "later_duplicate_before_cichlids", tier,
                       f"{src} {'place' if ' and ' in src else 'places'} a clade-level duplication at {' and '.join(older)}; copies are tandem and in no published pair", ";".join(older))
    if npair >= 2:
        tier = "core" if npair == 3 else "probable"
        other = f"; {_names(v, 'lost') or 'the third call'} says lost" if nlost else ("" if npair == 3 else "; the third call abstains")
        return res("R2/3/4", "retained_pair", tier, f"{npair} of 3 calls say a pair ({_names(v, 'pair')}){other}")
    if nlost >= 2:
        tier = "core" if nlost == 3 else "probable"
        other = f"; {_names(v, 'pair')} says a pair" if npair else ("" if nlost == 3 else "; the third call abstains")
        return res("R1", "lost_duplicate", tier, f"{nlost} of 3 calls say lost ({_names(v, 'lost')}){other}")
    say = "; ".join(f"{n}: {x or 'abstains'}" for n, x in zip(("thresholds", "parsimony", "tree"), v))
    return res("other", "unresolved", "disputed", f"no majority ({say})")


def tiers_summary(rows) -> dict:
    C = Counter
    scored = [r for r in rows if r.get("t_tier")]
    by = C((r["t_row"], r["t_tier"]) for r in scored)
    rows_in = sorted({k[0] for k in by})
    out = {"genes": len(scored), "by_row_and_tier": {r: {t: by[(r, t)] for t in TIERS if by[(r, t)]} for r in rows_in}}
    placed = [r for r in scored if r["t_tier"] != "unclassified"]
    out["placed"] = len(placed)
    for t in ("core", "probable", "possible", "disputed"):
        out[t] = sum(1 for r in placed if r["t_tier"] == t)
    out["unclassified"] = {k: v for k, v in C(r["t_call"] for r in scored if r["t_tier"] == "unclassified").items()}
    out["by_category"] = {c: {t: n for t in TIERS if (n := sum(1 for r in scored if r["t_category"] == c and r["t_tier"] == t))}
                          for c in CATEGORIES}
    out["votes_for_a_pair"] = dict(sorted(C(sum(1 for x in (r.get("t_vote_thr"), r.get("t_vote_par"), r.get("t_vote_tree")) if x == "pair")
                                            for r in placed if r["t_call"] in ("retained_pair", "lost_duplicate", "unresolved")).items()))
    pair_like = [r for r in placed if r["t_call"] == "retained_pair"]
    have = [r for r in pair_like if str(r.get("t_pairs_relaxed_n")) not in ("", "None")]
    out["published_pair_evidence_loaded"] = bool(have)
    if have:
        out["pair_like_genes_in_a_published_pair_pct"] = round(100 * sum(1 for r in have if _n(r.get("t_pairs_relaxed_n")) >= 1) / len(have), 1)
    return out
