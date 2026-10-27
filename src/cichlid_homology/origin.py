"""Where do the extra copies sit? For a gene with two or more focal-species copies in its reference clade
(human clade, or the zebrafish clade when there is no human gene), classify its nearest other copy:

    tandem_adjacent      next to it on the same sequence (no gene between)
    tandem_near          same sequence, at most `max_gap` genes between
    same_seqid_distant   same sequence, further apart
    different_seqid      on another sequence
    single_copy          no other copy in the clade
    no_position          a position is missing

Teleost-wide (genome duplication) copies are expected on different sequences; local duplications next to
each other. Gene order uses the proteome's genes only (owngenes.gene_order)."""
from __future__ import annotations

RANK = {"tandem_adjacent": 0, "tandem_near": 1, "same_seqid_distant": 2, "different_seqid": 3}
CAP = 6


def relation(g: str, o: str, genes: dict, order: dict, max_gap: int):
    """(class, gap_genes, gap_bp) between two genes, or None if a position is missing."""
    a, b = genes.get(g), genes.get(o)
    if not a or not b or g not in order or o not in order:
        return None
    if order[g][0] != order[o][0]:
        return "different_seqid", "", ""
    gap = abs(order[g][1] - order[o][1]) - 1
    bp = max(0, max(a["start"], b["start"]) - min(a["end"], b["end"]))
    cls = "tandem_adjacent" if gap <= 0 else "tandem_near" if gap <= max_gap else "same_seqid_distant"
    return cls, gap, bp


def attach_origin(rows: list, genes: dict, order: dict, max_gap: int = 10) -> None:
    """Adds orig_basis, orig_copies_n, orig_class and, when there is another copy, orig_nearest_copy,
    orig_gap_genes, orig_gap_bp, orig_copies_same_seqid and orig_copy_ids."""
    for r in rows:
        if r.get("a_human_cat") not in (None, "None"):
            basis = "human"
        elif r.get("a_zebrafish_cat") not in (None, "None"):
            basis = "zebrafish"
        else:
            basis = "none"
        copies = list(r.get(f"_a_{basis}_focal_leaves") or []) if basis != "none" else []
        r["orig_basis"] = basis
        r["orig_copies_n"] = len(copies)
        others = [c for c in copies if c != r["gene_id"]]
        if not others:
            r["orig_class"] = "single_copy"
            continue
        rel = []
        for o in others:
            x = relation(r["gene_id"], o, genes, order, max_gap)
            if x is not None:
                rel.append((RANK.get(x[0], 9), x[1] if x[1] != "" else 10 ** 9, o, x))
        if not rel:
            r["orig_class"] = "no_position"
            continue
        rel.sort()
        _, _, nearest, (cls, gap, bp) = rel[0]
        r["orig_class"], r["orig_nearest_copy"], r["orig_gap_genes"], r["orig_gap_bp"] = cls, nearest, gap, bp
        r["orig_copies_same_seqid"] = sum(1 for _, _, _, x in rel if x[0] != "different_seqid")
        r["orig_copy_ids"] = ";".join(others[:CAP]) + (f";+{len(others) - CAP}" if len(others) > CAP else "")
