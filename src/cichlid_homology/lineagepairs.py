"""Lineage-aware pairing of our genes in the clades that have a teleost split (0.23.0).

The gene table decides a clade's pairs from copy numbers: a 1:2 clade with two of our genes is a pair, one with three or more is "ambiguous", a 1:3+ clade is an "expansion". Those counts cannot see
that an expansion is often a pair plus an extra copy of one lineage: 78 of the published strict pairs we missed sit in clades with three or more of our genes, 95% of them on a chromosome pair that
carries about 170 of our other pairs. The records of `figure-pilot` say, for every gene of a clade that has a teleost split, which of the two lineages it sits in (from its place in the tree, not
from the counts). Here that is used.

In a clade with a teleost split: genes in lineage A and lineage B form the pair. One gene in each is a pair. With more, the pair is the cross-lineage combination whose chromosome pair carries most of
our other firm pairs (it has to carry at least `synteny_min_pairs`); a tie is broken by the zebrafish gene that anchors each lineage, then by the weight of the placement. The other genes of a lineage
are in-paralogs of the paired gene. Without synteny support the clade stays ambiguous, with its candidates listed. Genes in one lineage only are a lost partner (one gene) or in-paralogs of that
lineage (several). If gar or reedfish carry a gene in each lineage the duplication is older than the teleost duplication, and the pair is labelled that way and does not count for the synteny of the
teleost pairs. The lineage signal of the pair confidence is rebuilt from the lineages found here: different lineages with both genes placed at a weight of at least 0.5. Clades without a record,
with no teleost split, or whose gene is not in the clade node, are left as they were."""
from __future__ import annotations

import collections
import json
from pathlib import Path

FIRM = ("paired", "paired (resolved by synteny)", "pair retained in few teleosts (synteny-supported)")
POSSIBLE = "possible pair (the other teleosts lost a copy)"
RESOLVED = "paired (resolved by synteny)"
OLDER = "older duplicate pair (gar or reedfish carry both copies)"
INPARALOG = "in-paralog (extra copy of one lineage)"
ONE_LINEAGE = "later duplicates of one lineage (the other lineage is not in the focal species)"
AMBIGUOUS = "partner ambiguous (lineages known)"
CLASS_OHNOLOG = "ohnolog (teleost genome duplication)"
CLASS_OLDER = "older than the teleost duplication"
GAR, REEDFISH = "lepisosteus_oculatus", "erpetoichthys_calabaricus"
MIN_LWR = 0.5


def load_records(path) -> dict:
    """{clade: record} for the records with a teleost split."""
    out = {}
    with open(Path(path).expanduser()) as fh:
        for line in fh:
            if line.strip():
                r = json.loads(line)
                if r.get("mode") == "split" and (r.get("split") or {}).get("kind") == "teleost":
                    out[r["clade"]] = r
    return out


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _chrom_pair(a, b):
    return tuple(sorted((a["seqid"], b["seqid"]))) if a.get("seqid") and b.get("seqid") and a["seqid"] != b["seqid"] else None


def _position(a, b, params):
    if not a.get("seqid") or not b.get("seqid"):
        return None
    if a["seqid"] != b["seqid"]:
        return "different chromosomes"
    try:
        gap = abs(int(a["gene_index"]) - int(b["gene_index"])) - 1
    except (TypeError, ValueError):
        return "unknown"
    return "tandem" if gap <= params.get("tandem_gap_genes", 10) else "same chromosome, distant"


def _confidence(lin, pos, zs, syn, min_pairs):
    from .genetable import _confidence as conf
    return conf(lin, pos, zs, syn, min_pairs)


def _older_evidence(rec: dict, letters) -> list:
    ev = []
    for name, sp in (("gar", GAR), ("reedfish", REEDFISH)):
        lins = {o["lineage"] for o in rec["outgroup"]["inside"] if o["species"] == sp}
        if set(letters) <= lins:
            ev.append(name)
    return ev


def _pub_text(a, b, pub):
    sa, sb = (a.get("zebrafish_symbol") or "").lower(), (b.get("zebrafish_symbol") or "").lower()
    return "the two zebrafish names are a published pair" if sa and sb and sb in (pub or {}).get(sa, ()) else ""


def apply(rows: list, records: dict, partners=None, params=None) -> dict:
    """Changes `rows` (the gene table) in place for the clades that have a record; returns a summary with the transitions of ohnolog_status."""
    params = params or {}
    minp = params.get("synteny_min_pairs", 5)
    pub = (partners or {}).get("zebrafish", {})
    by_id = {r["gene_id"]: r for r in rows}
    by_clade = collections.defaultdict(list)
    for r in rows:
        if r.get("clade_id") in records:
            by_clade[r["clade_id"]].append(r)
    before = {r["gene_id"]: r["ohnolog_status"] for r in rows}
    # chromosome-pair counts of the firm pairs as they stand (the older duplicates are another paralogon and are not counted)
    cp = collections.Counter()
    seen = set()
    for r in rows:
        if r["ohnolog_status"] in FIRM and r.get("ohnolog_gene") in by_id:
            key = frozenset((r["gene_id"], r["ohnolog_gene"]))
            if key in seen:
                continue
            seen.add(key)
            c = _chrom_pair(r, by_id[r["ohnolog_gene"]])
            if c:
                cp[c] += 1
    summary = collections.Counter()
    new_pairs = {}
    for cid in sorted(by_clade):
        rec, genes = records[cid], by_clade[cid]
        letters = rec["letters"][:2]
        rg = {g["gene"]: g for g in rec["genes"]}
        for r in genes:
            g = rg.get(r["gene_id"])
            if g is None:
                r["lineage"], r["lineage_basis"] = "", ""
                continue
            lin = g["lineage"]
            r["lineage"] = "outside" if lin.startswith("outside") else lin
            r["lineage_basis"] = "placement"
        if len(letters) < 2:
            continue
        la = [r for r in genes if r["lineage"] == letters[0]]
        lb = [r for r in genes if r["lineage"] == letters[1]]
        anchors = {L: {a["gene"] for a in rec["anchors"].get(L, {}).get("zebrafish", [])} for L in letters}
        older = _older_evidence(rec, letters)
        if la and lb:
            if len(la) == 1 and len(lb) == 1:
                best, how = (la[0], lb[0]), "lineage"
                if best[0]["ohnolog_status"] not in FIRM + (POSSIBLE,):
                    best[0]["ohnolog_status"] = best[1]["ohnolog_status"] = "paired"
            else:
                scored = []
                for a in la:
                    for b in lb:
                        c = _chrom_pair(a, b)
                        s = cp.get(c, 0) if c else -1
                        m = (a.get("zebrafish_gene") in anchors[letters[0]]) + (b.get("zebrafish_gene") in anchors[letters[1]])
                        w = (_f(rg[a["gene_id"]].get("lwr")) or 0) + (_f(rg[b["gene_id"]].get("lwr")) or 0)
                        scored.append(((s >= minp), s, m, w, a["gene_id"], b["gene_id"]))
                scored.sort(key=lambda t: (-t[0], -t[1], -t[2], -t[3], t[4], t[5]))
                top = scored[0]
                if not top[0]:
                    for r in la + lb:
                        r["ohnolog_status"] = AMBIGUOUS
                        r["ohnolog_candidates"] = ";".join(x["gene_id"] for x in (lb if r in la else la))
                    summary["clades left ambiguous (no synteny support)"] += 1
                    continue
                best, how = (by_id[top[4]], by_id[top[5]]), "lineage and synteny" + ("" if len([t for t in scored if t[:2] == top[:2]]) == 1 else ", the synteny tied: broken by the zebrafish anchor and the placement weight")
                a, b = best
                for x, y in ((a, b), (b, a)):
                    x["ohnolog_status"] = RESOLVED
                    x["lineage_basis"] = f"placement; paired by {how}"
                for r in la + lb:
                    if r is not a and r is not b:
                        mate = a if r["lineage"] == letters[0] else b
                        r["ohnolog_status"], r["in_paralog_of"], r["ohnolog_gene"], r["ohnolog_symbol"], r["pair_id"] = INPARALOG, mate["gene_id"], "", "", ""
                        r["ohnolog_candidates"], r["ohnolog_lineage"] = "", r["lineage"]
                summary["pairs chosen among several genes"] += 1
            a, b = best
            for x, y in ((a, b), (b, a)):
                x["ohnolog_gene"], x["ohnolog_symbol"] = y["gene_id"], y.get("symbol", "")
                x["pair_id"] = cid if how == "lineage" else f"{cid}:S{a['gene_id'][-6:]}"
                x["ohnolog_lineage"] = x["lineage"]
                x["ohnolog_zebrafish_gene"] = y.get("zebrafish_gene", "")
                x["ohnolog_candidates"] = ""
                x["published_pair"] = _pub_text(x, y, pub) or x.get("published_pair", "")
                lw = [_f(rg[z["gene_id"]].get("lwr")) for z in (x, y)]
                x["ev_lineages"] = "different" if all(w is None or w >= MIN_LWR for w in lw) else "unassigned"
                pos = _position(x, y, params)
                if pos:
                    x["ev_position"] = pos
                zs = "unknown" if not (x.get("zebrafish_gene") and y.get("zebrafish_gene")) else ("yes" if x["zebrafish_gene"] != y["zebrafish_gene"] else "no")
                x["ev_zebrafish_split"] = zs
                if older:
                    x["ohnolog_status"] = OLDER
                    x["duplication_class"], x["older_evidence"] = CLASS_OLDER, " and ".join(older) + " carry a gene in each lineage"
                else:
                    x["duplication_class"], x["older_evidence"] = CLASS_OHNOLOG, ""
            new_pairs[frozenset((a["gene_id"], b["gene_id"]))] = (a, b, bool(older))
            if older:
                summary["pairs labelled older than the teleost duplication"] += 1
        elif la or lb:
            one = la or lb
            if len(one) == 1:
                if one[0]["ohnolog_status"] not in ("partner lost in the focal species",):
                    one[0]["ohnolog_status"] = "partner lost in the focal species"
                    summary["genes whose partner is lost (one lineage)"] += 1
            else:
                for r in one:
                    r["ohnolog_status"] = ONE_LINEAGE
                summary["genes in the one lineage only"] += len(one)
    # the synteny and the confidence of every pair, from the pairs as they now are
    pairs = {}
    for r in rows:
        if r.get("ohnolog_gene") in by_id and r["ohnolog_status"] in FIRM + (POSSIBLE, OLDER):
            pairs.setdefault(frozenset((r["gene_id"], r["ohnolog_gene"])), []).append(r)
    cp2, firm = collections.Counter(), set()
    for key, rs in pairs.items():
        a, b = (by_id[g] for g in tuple(key))
        r0 = rs[0]
        if r0["ohnolog_status"] == "paired" and r0["ev_position"] == "different chromosomes":
            conf = _confidence(r0["ev_lineages"], r0["ev_position"], r0["ev_zebrafish_split"], "", minp)
            if conf in ("high", "medium"):
                c = _chrom_pair(a, b)
                if c:
                    cp2[c] += 1
                    firm.add(key)
    for key, rs in pairs.items():
        a, b = (by_id[g] for g in tuple(key))
        c = _chrom_pair(a, b)
        for r in rs:
            if r["ohnolog_status"] == RESOLVED:
                syn = cp2.get(c, 0) if c else ""
                r["ev_synteny"] = syn
                r["ohnolog_confidence"] = _confidence(r["ev_lineages"], r["ev_position"], r["ev_zebrafish_split"], "", minp)
                continue
            syn = ""
            if r["ev_position"] == "different chromosomes" and c:
                syn = cp2.get(c, 0) - (1 if key in firm else 0)
            r["ev_synteny"] = syn
            r["ohnolog_confidence"] = _confidence(r["ev_lineages"], r["ev_position"], r["ev_zebrafish_split"], syn, minp)
            if r["ohnolog_status"] == POSSIBLE and syn != "" and syn >= minp:
                r["ohnolog_status"] = "pair retained in few teleosts (synteny-supported)"
    tr = collections.Counter((before[r["gene_id"]], r["ohnolog_status"]) for r in rows if before[r["gene_id"]] != r["ohnolog_status"])
    summary_out = {"clades_with_a_record": len(by_clade), "genes_in_those_clades": sum(len(v) for v in by_clade.values()), **summary,
                   "status_transitions": [[a, b, k] for (a, b), k in tr.most_common()]}
    return summary_out
