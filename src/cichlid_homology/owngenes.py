"""The focal species' own gene table (through the coordinator's layout): symbol, position and the
annotation's quality flags. Used for the symbol check, the origin of extra copies (tandem or not) and
the artifact-risk flags."""
from __future__ import annotations

import csv
import sys

csv.field_size_limit(sys.maxsize)


def truthy(v) -> bool:
    return str(v).strip().lower() in ("true", "1", "yes")


def _int(v, default=0):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return default


def own_genes(layout, prefix: str) -> dict:
    """{gene_id: {symbol, seqid, start, end, strand, partial, low_quality, internal_stops,
    protein_length, description, in_proteome}} from `<prefix>.genes.tsv`."""
    out = {}
    with open(layout.genes_table(prefix).local, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["gene_id"]] = {
                "symbol": r.get("symbol", ""), "seqid": r.get("seqid", ""),
                "start": _int(r.get("start"), -1), "end": _int(r.get("end"), -1),
                "strand": r.get("strand", ""), "partial": truthy(r.get("partial")),
                "low_quality": truthy(r.get("low_quality")), "internal_stops": _int(r.get("internal_stops")),
                "protein_length": _int(r.get("protein_length")), "description": r.get("description", ""),
                "in_proteome": bool(r.get("chosen_protein")) and not truthy(r.get("mitochondrial"))}
    return out


def gene_order(genes: dict) -> dict:
    """{gene_id: (seqid, index)}: genes of the proteome ordered by start along each sequence, so the number
    of genes between two copies can be counted."""
    by_seq: dict = {}
    for g, i in genes.items():
        if i["in_proteome"] and i["seqid"] and i["start"] >= 0:
            by_seq.setdefault(i["seqid"], []).append((i["start"], i["end"], g))
    order = {}
    for seq, lst in by_seq.items():
        for idx, (_, _, g) in enumerate(sorted(lst)):
            order[g] = (seq, idx)
    return order
