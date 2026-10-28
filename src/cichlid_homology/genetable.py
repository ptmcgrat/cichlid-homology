"""The gene table of the focal species (0.18.0): one row per gene with its closest human and zebrafish gene, the clade it sits in, the
teleost trend of that clade, the B1 call, and its ohnolog.

Inputs (all outputs of earlier steps, nothing is downloaded or computed from sequence): the pilot gene table of Compara (and of TreeFam, for
agreement), the clade table of `clade-split` (a row per clade of the reference families, with the protein IDs of its human and zebrafish genes),
the published teleost-duplication pairs (optional), and the gene positions of the focal species (optional).

A gene's clade is the clade holding its closest zebrafish protein (a teleost gene belongs with the teleost genes of its clade); when it has no
zebrafish gene, the clade holding its closest human protein if that clade has teleosts. Relationships are the clade's: the number of human (or
zebrafish) genes against the MEDIAN number of genes per teleost species in the clade (`1:1`, `1:2`, ...), not the focal species' own copy number,
which is compared with the trend in `copy_vs_trend`.

Ohnolog (the other copy of the teleost genome duplication): in a clade with a 1:2 trend and exactly two focal genes they are a pair; with one, the
partner was lost in the focal species (the zebrafish genes in the clade are listed); with three or more, the partner is ambiguous until the trees
order our own genes (OrthoFinder). In a 1:1 clade with two focal genes on different chromosomes the pair is possible (the other teleosts lost a
copy, the focal species kept it); tandem copies are a later duplicate. A pair's confidence counts three independent signals: the two genes sit in
different lineages of the duplication, on different chromosomes, and match two different zebrafish genes (high: all three; medium: two and none
against; low: otherwise). Published pairs are an outside check, not part of the call."""
from __future__ import annotations

import collections
import json

from .tiers import assign

# (name, group, meaning)
COLUMN_SPEC = [
    ("gene_id", "Identity", "gene ID of the focal species"),
    ("symbol", "Identity", "the gene's own annotated symbol"),
    ("description", "Identity", "its annotated description"),
    ("seqid", "Identity", "chromosome or scaffold (blank unless positions are given)"),
    ("start", "Identity", "start position"),
    ("gene_index", "Identity", "order of the gene along its sequence (genes between two copies are counted from it)"),
    ("family", "Family and clade", "the reference family"),
    ("clade_id", "Family and clade", "family plus clade, cut at the ancient duplications"),
    ("clade_assignment", "Family and clade", "how the gene was placed: zebrafish (closest zebrafish protein), human (closest human protein), unplaced"),
    ("clade_type", "Family and clade", "human clade (1 or 2+ human genes), tetrapod clade with no human gene, fish-only clade (out-paralog if the family has a human gene)"),
    ("clade_n_human", "Family and clade", "human genes in the clade"),
    ("clade_n_zebrafish", "Family and clade", "zebrafish genes in the clade"),
    ("clade_n_focal", "Family and clade", "focal-species genes in the clade"),
    ("human_gene", "Closest human gene", "protein ID"),
    ("human_symbol", "Closest human gene", "symbol"),
    ("human_hgnc", "Closest human gene", "HGNC ID"),
    ("human_name_reliable", "Closest human gene", "whether the name is reliable"),
    ("human_group", "Closest human gene", "one-to-one, one-to-two, one-to-many or other (focal genes against human genes in the smallest clade holding both)"),
    ("human_agree_treefam", "Closest human gene", "TreeFam's closest human gene is the same, different, or TreeFam has none or does not hold the gene"),
    ("zebrafish_gene", "Closest zebrafish gene", "protein ID"),
    ("zebrafish_symbol", "Closest zebrafish gene", "symbol"),
    ("zfin_id", "Closest zebrafish gene", "ZFIN ID"),
    ("zebrafish_group", "Closest zebrafish gene", "same grouping as for human"),
    ("zebrafish_agree_treefam", "Closest zebrafish gene", "TreeFam's closest zebrafish gene is the same, different, or TreeFam has none"),
    ("teleost_trend", "Teleost trend and relationships", "1:1, 1:2, 1:3+ or mixed, from the median genes per teleost species in the clade"),
    ("teleost_median_genes", "Teleost trend and relationships", "that median"),
    ("early_vs_percomorph", "Teleost trend and relationships", "the median in the early teleosts and in the percomorphs and below, e.g. 1/2"),
    ("rel_human", "Teleost trend and relationships", "human genes : teleost median (5a and 6), or no human gene (tetrapods present or not)"),
    ("rel_zebrafish", "Teleost trend and relationships", "zebrafish genes : teleost median (5b), or no zebrafish gene"),
    ("copy_vs_trend", "Teleost trend and relationships", "the focal copy number in the clade against the trend: matches, extra copies, fewer copies, mixed trend"),
    ("human_clean_tier", "Human-side cleanliness", "clean, acceptable or messy (clades with one human gene), else not applicable"),
    ("tetrapod_single_copy_share", "Human-side cleanliness", "share of tetrapod species with exactly one gene in the clade"),
    ("tetrapod_species", "Human-side cleanliness", "tetrapod species in the clade"),
    ("human_rel_rate", "Human-side cleanliness", "human root-to-tip distance relative to the other tetrapods"),
    ("b1_category", "B1 call", "lost duplicate, both copies kept, later duplicate in cichlids, later duplicate before cichlids, disputed, unclassified"),
    ("b1_confidence", "B1 call", "core (all three calls agree), probable (two), possible, disputed or unclassified"),
    ("b1_basis", "B1 call", "the reason in words"),
    ("vote_thresholds", "B1 call", "the share-based call's vote: pair, lost or blank"),
    ("vote_parsimony", "B1 call", "the parsimony call's vote"),
    ("vote_tree", "B1 call", "the gene-tree call's vote"),
    ("tree_call", "B1 call", "the gene-tree call"),
    ("ohnolog_status", "Ohnolog", "paired; possible pair; partner lost in the focal species; partner ambiguous; single copy (duplicate lost); later duplicate; expansion; mixed trend; not applicable"),
    ("pair_id", "Ohnolog", "shared by the two genes of a pair"),
    ("ohnolog_gene", "Ohnolog", "the partner's gene ID"),
    ("ohnolog_symbol", "Ohnolog", "the partner's symbol"),
    ("ohnolog_lineage", "Ohnolog", "this gene's lineage (A or B) in the duplication"),
    ("ohnolog_confidence", "Ohnolog", "high, medium or low, from the three signals below"),
    ("ev_lineages", "Ohnolog", "the two genes sit in different lineages of the duplication: different, same, unassigned"),
    ("ev_position", "Ohnolog", "different chromosomes, tandem, same chromosome distant, or unknown"),
    ("ev_zebrafish_split", "Ohnolog", "the two genes match two different zebrafish genes: yes, no, unknown"),
    ("ohnolog_zebrafish_gene", "Ohnolog", "the zebrafish gene matching the partner"),
    ("published_pair", "Ohnolog", "the two zebrafish names are a published pair, or the gene is in a published pair in n of 4 species (outside check)"),
    ("partner_elsewhere", "Ohnolog", "when the partner was lost in the focal species: the other zebrafish gene(s) in the clade"),
    ("later_dup_level", "Later duplicates and out-paralogs", "teleost, percomorph or ovalentarian level of a later duplicate before the cichlids"),
    ("nearest_copy_gene", "Later duplicates and out-paralogs", "the nearest other focal copy"),
    ("nearest_copy_class", "Later duplicates and out-paralogs", "tandem adjacent, tandem near, same chromosome distant, different chromosome"),
    ("out_paralog", "Later duplicates and out-paralogs", "yes when the clade is fish-only in a family that has human genes"),
    ("carriers", "Later duplicates and out-paralogs", "cichlid species carrying a later duplicate (OrthoFinder trees only)"),
    ("unclassified_reason", "Flags", "no anchor, expanded family or many copies"),
    ("expanded_family", "Flags", "many focal genes in the family"),
    ("risk_flags", "Flags", "low placement support, partial or low-quality gene model, fusion candidate"),
    ("symbol_match_human", "Flags", "the gene's own symbol against the closest human symbol: exact, stem, none, unnamed"),
    ("symbol_match_zebrafish", "Flags", "the same against the closest zebrafish symbol"),
    ("pipeline_version", "Provenance", "the code version"),
    ("run_date", "Provenance", "the date the table was made"),
    ("parameters", "Provenance", "the settings used"),
]
COLUMNS = [c for c, _, _ in COLUMN_SPEC]
DEFAULTS = {"strict": {"single": 0.9, "tetrapod_species": 10, "rate": (0.8, 1.2), "teleost_share": 0.9},
            "loose": {"single": 0.7, "tetrapod_species": 5, "rate": (0.7, 1.3), "teleost_share": 0.8}, "tandem_gap_genes": 10}
POS_OF_CLASS = {"different_seqid": "different chromosomes", "tandem_adjacent": "tandem", "tandem_near": "tandem", "same_seqid_distant": "same chromosome, distant"}


def _i(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return 0


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _fmt(m):
    m = _f(m)
    if m is None:
        return ""
    return str(int(m)) if m == int(m) else str(m)


def _key(row, side):
    if row is None:
        return ""
    i = row.get(f"a_{side}_best_id")
    s = (row.get(f"a_{side}_best_symbol") or "").lower()
    return i if i else (f"sym:{s}" if s else "")


def _agree(row, tf, side):
    if tf is None:
        return "not in TreeFam"
    if tf.get(f"a_{side}_group") == "None":
        return "TreeFam has none"
    k = _key(row, side)
    if not k or not _key(tf, side):
        return ""
    return "same" if k == _key(tf, side) else "different"


def human_tier(c, p=DEFAULTS) -> str:
    """clean, acceptable or messy for a clade with one human gene; the settings are the owner's to change (see docs/decisions.md)."""
    if c is None or _i(c["n_human"]) != 1:
        return "not applicable"

    def ok(s):
        sh, rate = _f(c.get("tetrapod_single_copy_share")), _f(c.get("human_rel_rate"))
        if sh is None or sh < s["single"] or _i(c.get("n_tetrapod_species")) < s["tetrapod_species"]:
            return False
        if rate is not None and not (s["rate"][0] <= rate <= s["rate"][1]):
            return False
        if _i(c.get("n_teleost_species")) == 0 or (_f(c.get("teleost_species_share_of_family")) or 0) < s["teleost_share"]:
            return False
        return True
    return "clean" if ok(p["strict"]) else ("acceptable" if ok(p["loose"]) else "messy")


def _pair_evidence(a, b, positions, params):
    """(lineages, position, zebrafish split, confidence) of the pair of genes a and b, the same whichever gene is asked about."""
    la, lb = a.get("d_gene_lineage"), b.get("d_gene_lineage")
    lin = "unassigned" if not (la in ("A", "B") and lb in ("A", "B")) else ("different" if la != lb else "same")
    pa, pb = (positions or {}).get(a["gene_id"]), (positions or {}).get(b["gene_id"])
    if pa and pb:
        if pa[0] != pb[0]:
            pos = "different chromosomes"
        else:
            gap = abs(int(pa[2]) - int(pb[2])) - 1 if pa[2] != "" and pb[2] != "" else None
            pos = "unknown" if gap is None else ("tandem" if gap <= params.get("tandem_gap_genes", 10) else "same chromosome, distant")
    elif a.get("orig_nearest_copy") == b["gene_id"]:
        pos = POS_OF_CLASS.get(a.get("orig_class", ""), "unknown")
    elif b.get("orig_nearest_copy") == a["gene_id"]:
        pos = POS_OF_CLASS.get(b.get("orig_class", ""), "unknown")
    else:
        pos = "unknown"
    za, zb = a.get("a_zebrafish_best"), b.get("a_zebrafish_best")
    zs = "unknown" if not (za and zb) else ("yes" if za != zb else "no")
    plus = [lin == "different", pos == "different chromosomes", zs == "yes"]
    minus = [lin == "same", pos == "tandem"]
    return lin, pos, zs, ("high" if all(plus) else ("medium" if sum(plus) >= 2 and not any(minus) else "low"))


def _clade_type(c, fam_human):
    if c is None:
        return "unplaced"
    h = _i(c["n_human"])
    if h == 1:
        return "human clade (1 human gene)"
    if h >= 2:
        return "human clade (2+ human genes)"
    if _i(c["n_tetrapod_species"]) > 0:
        return "tetrapod clade, no human gene"
    return "fish-only clade, family has a human gene (out-paralog)" if fam_human else "fish-only clade, family has no human gene"


def build_gene_table(compara, treefam, groups, families, partners=None, positions=None, version="", run_date="", params=None, focal_prefix="", reassign=False):
    """The table as a list of dicts in COLUMNS order. `compara`/`treefam`: pilot gene rows (lists of dicts); `groups`/`families`: the clade tables;
    `partners`: published pairs from ohnologs.load_pair_partners; `positions`: {gene_id: (seqid, start, index)}."""
    params = params or DEFAULTS
    tf = {r["gene_id"]: r for r in (treefam or [])}
    by_h, by_z, fam_h = {}, {}, {int(0): 0}
    for g in groups:
        for p in g["human_proteins"].split(";"):
            if p:
                by_h[(g["family"], p)] = g
        for p in g["zebrafish_proteins"].split(";"):
            if p:
                by_z[(g["family"], p)] = g
    fam_h = {f: _i(r["n_human"]) for f, r in families.items()} if isinstance(families, dict) else {r["family"]: _i(r["n_human"]) for r in families}
    rows = []
    for r in compara:
        if reassign or "t_category" not in r:                      # an older pilot table: the categories are recomputed with the current logic
            r = {**r, **assign(r)}
        zc = by_z.get((r["family"], r.get("a_zebrafish_best"))) if r.get("a_zebrafish_best") else None
        hc = by_h.get((r["family"], r.get("a_human_best"))) if r.get("a_human_best") else None
        c, how = None, "unplaced"
        if zc is not None:
            c, how = zc, "zebrafish" if (hc is None or hc is zc) else "zebrafish (closest human gene is in another clade)"
        elif hc is not None and _i(hc["n_teleost_species"]) > 0:
            c, how = hc, "human"
        rows.append({"_r": r, "_c": c, "_how": how})
    n_in = collections.Counter((x["_c"]["family"], x["_c"]["group"]) for x in rows if x["_c"] is not None)
    by_clade = collections.defaultdict(list)
    for x in rows:
        if x["_c"] is not None:
            by_clade[(x["_c"]["family"], x["_c"]["group"])].append(x)
    out = []
    partner_of = {}
    for key, xs in by_clade.items():
        trend = xs[0]["_c"]["teleost_trend"]
        if len(xs) == 2 and trend in ("1:1", "1:2"):
            partner_of[xs[0]["_r"]["gene_id"]], partner_of[xs[1]["_r"]["gene_id"]] = xs[1], xs[0]
    pub = (partners or {}).get("zebrafish", {})
    pair_cache = {}
    for x in rows:
        r, c = x["_r"], x["_c"]
        gid = r["gene_id"]
        pos = (positions or {}).get(gid)
        o = {k: "" for k in COLUMNS}
        o.update(gene_id=gid, symbol=r.get("own_symbol", ""), description=r.get("gt_description", ""), family=r.get("family", ""),
                 seqid=pos[0] if pos else "", start=pos[1] if pos else "", gene_index=pos[2] if pos else "")
        o.update(human_gene=r.get("a_human_best", ""), human_symbol=r.get("a_human_best_symbol", ""), human_hgnc=r.get("a_human_best_id", ""),
                 human_name_reliable=r.get("a_human_best_name_reliable", ""), human_group=r.get("a_human_group", ""),
                 human_agree_treefam=_agree(r, tf.get(gid), "human") if r.get("a_human_best") else "",
                 zebrafish_gene=r.get("a_zebrafish_best", ""), zebrafish_symbol=r.get("a_zebrafish_best_symbol", ""),
                 zfin_id=r.get("a_zebrafish_best_id", ""), zebrafish_group=r.get("a_zebrafish_group", ""),
                 zebrafish_agree_treefam=_agree(r, tf.get(gid), "zebrafish") if r.get("a_zebrafish_best") else "")
        o.update(b1_category=r.get("t_category", ""), b1_confidence=r.get("t_tier", ""), b1_basis=r.get("t_basis", ""),
                 vote_thresholds=r.get("t_vote_thr", ""), vote_parsimony=r.get("t_vote_par", ""), vote_tree=r.get("t_vote_tree", ""),
                 tree_call=r.get("d_call", ""), later_dup_level=r.get("t_when", ""), nearest_copy_gene=r.get("orig_nearest_copy", ""),
                 nearest_copy_class=r.get("orig_class", ""), carriers=r.get("d_carriers", ""),
                 unclassified_reason=r.get("t_call", "") if r.get("t_tier") == "unclassified" else "", expanded_family=r.get("expanded_family", ""),
                 risk_flags=r.get("risk_flags", ""), symbol_match_human=r.get("a_human_symbol_match", ""),
                 symbol_match_zebrafish=r.get("a_zebrafish_symbol_match", ""), pipeline_version=version, run_date=run_date,
                 parameters=json.dumps(params, default=list, sort_keys=True))
        o["clade_assignment"] = x["_how"]
        o["clade_type"] = _clade_type(c, fam_h.get(r.get("family"), 0) > 0)
        if c is None:
            o["out_paralog"] = ""
            o["ohnolog_status"] = "not applicable"
            out.append(o)
            continue
        nf, nh, nz = n_in[(c["family"], c["group"])], _i(c["n_human"]), _i(c["zebrafish_genes"])
        med, trend = c["teleost_median_genes"], c["teleost_trend"]
        o.update(clade_id=f"{c['family']}:{c['group']}", clade_n_human=nh, clade_n_zebrafish=nz, clade_n_focal=nf, teleost_trend=trend,
                 teleost_median_genes=_fmt(med), early_vs_percomorph=f"{_fmt(c['early_teleost_median_genes']) or '-'}/{_fmt(c['percomorph_median_genes']) or '-'}",
                 tetrapod_single_copy_share=c["tetrapod_single_copy_share"], tetrapod_species=c["n_tetrapod_species"],
                 human_rel_rate=c["human_rel_rate"], human_clean_tier=human_tier(c, params))
        o["out_paralog"] = "yes" if o["clade_type"].startswith("fish-only clade, family has a human gene") else "no"
        has_tel = _i(c["n_teleost_species"]) > 0
        m = _f(med)
        o["rel_human"] = (f"{nh}:{_fmt(med)}" if nh and has_tel else ("no human gene; tetrapods present" if _i(c["n_tetrapod_species"]) > 0 else "no human gene; no tetrapods")) if has_tel or nh == 0 else ""
        o["rel_zebrafish"] = f"{nz}:{_fmt(med)}" if nz and has_tel else "no zebrafish gene"
        if m is not None and has_tel:
            o["copy_vs_trend"] = "mixed trend" if m != int(m) else ("matches trend" if nf == m else ("extra copies" if nf > m else "fewer copies"))
        # ohnolog
        st = "not applicable"
        if trend == "1:2":
            st = "paired" if nf == 2 else ("partner lost in the focal species" if nf == 1 else "partner ambiguous (three or more focal genes)")
        elif trend == "1:1":
            if nf == 1:
                st = "single copy (duplicate lost)"
            elif nf == 2:
                pr = partner_of.get(gid)
                if pr:                                                       # decided on the pair's evidence, so both genes get the same status
                    key = frozenset((gid, pr["_r"]["gene_id"]))
                    if key not in pair_cache:
                        pair_cache[key] = _pair_evidence(r, pr["_r"], positions, params)
                    st = "possible pair (the other teleosts lost a copy)" if pair_cache[key][1] == "different chromosomes" else "later duplicate"
                else:
                    st = "later duplicate"
            else:
                st = "later duplicates (three or more focal genes)"
        elif trend == "1:3+":
            st = "expansion"
        elif trend in ("mixed 1-2", "mixed 2-3"):
            st = "mixed trend"
        o["ohnolog_status"] = st
        o["ohnolog_lineage"] = r.get("d_gene_lineage", "") if st in ("paired", "possible pair (the other teleosts lost a copy)") else ""
        pr = partner_of.get(gid)
        if pr and st in ("paired", "possible pair (the other teleosts lost a copy)"):
            pr_r = pr["_r"]
            o["pair_id"] = f"{c['family']}:{c['group']}"
            o["ohnolog_gene"], o["ohnolog_symbol"] = pr_r["gene_id"], pr_r.get("own_symbol", "")
            ev = pair_cache.get(frozenset((gid, pr_r["gene_id"])))
            if ev is None:                                                # the evidence is computed once per pair: both genes get the same
                ev = _pair_evidence(r, pr_r, positions, params)
                pair_cache[frozenset((gid, pr_r["gene_id"]))] = ev
            o["ev_lineages"], o["ev_position"], o["ev_zebrafish_split"], o["ohnolog_confidence"] = ev
            o["ohnolog_zebrafish_gene"] = pr_r.get("a_zebrafish_best") or ""
            sa, sb = (r.get("a_zebrafish_best_symbol") or "").lower(), (pr_r.get("a_zebrafish_best_symbol") or "").lower()
            if sa and sb and sb in pub.get(sa, ()):
                o["published_pair"] = "the two zebrafish names are a published pair"
            elif _i(r.get("t_pairs_relaxed_n")) >= 1:
                o["published_pair"] = f"gene in a published pair in {_i(r.get('t_pairs_relaxed_n'))} of 4 species"
        elif _i(r.get("t_pairs_relaxed_n")) >= 1:
            o["published_pair"] = f"gene in a published pair in {_i(r.get('t_pairs_relaxed_n'))} of 4 species"
        if st == "partner lost in the focal species":
            z = [p for p in c["zebrafish_proteins"].split(";") if p and p != r.get("a_zebrafish_best")]
            o["partner_elsewhere"] = ";".join(z)
        out.append(o)
    return out


def summarize_table(rows) -> dict:
    C = collections.Counter
    return {"genes": len(rows), "clade_assignment": dict(C(r["clade_assignment"] for r in rows)),
            "ohnolog_status": dict(C(r["ohnolog_status"] for r in rows).most_common()),
            "pair_confidence": dict(C(r["ohnolog_confidence"] for r in rows if r["ohnolog_status"] == "paired")),
            "rel_human": dict(C(r["rel_human"] for r in rows).most_common(8)), "rel_zebrafish": dict(C(r["rel_zebrafish"] for r in rows).most_common(8)),
            "copy_vs_trend": dict(C(r["copy_vs_trend"] for r in rows)), "human_clean_tier": dict(C(r["human_clean_tier"] for r in rows)),
            "b1_category": dict(C(r["b1_category"] for r in rows).most_common()), "out_paralog": dict(C(r["out_paralog"] for r in rows))}


def columns_markdown() -> str:
    L, grp = ["# Columns of the gene table", ""], None
    for name, g, meaning in COLUMN_SPEC:
        if g != grp:
            L += ["", f"## {g}", "", "| column | meaning |", "|---|---|"]
            grp = g
        L.append(f"| {name} | {meaning} |")
    return "\n".join(L) + "\n"
