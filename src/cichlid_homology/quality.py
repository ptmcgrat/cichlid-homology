"""Artifact-risk flags for a gene's extra-copy call: signs that the extra copy may be an annotation or
placement artifact rather than a real duplicate. The call is kept; the flags are evidence."""
from __future__ import annotations


def _lwr(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def attach_quality(rows: list, genes: dict) -> None:
    for r in rows:
        gi = genes.get(r["gene_id"], {})
        r["gt_partial"] = gi.get("partial", "")
        r["gt_low_quality"] = gi.get("low_quality", "")
        r["gt_internal_stops"] = gi.get("internal_stops", "")
        r["gt_protein_length"] = gi.get("protein_length", "")
        r["gt_description"] = gi.get("description", "")
        flags = []
        lwr = _lwr(r.get("lwr"))
        if lwr is not None and lwr < 0.5:
            flags.append("low_placement_support")
        if "possible_fusion" in (r.get("assign_flags") or ""):
            flags.append("fusion_candidate")
        if gi.get("partial"):
            flags.append("partial_gene_model")
        if gi.get("low_quality"):
            flags.append("low_quality_model")
        if gi.get("internal_stops", 0):
            flags.append("internal_stops")
        try:
            if int(r.get("b_cichlid_with2p") or 0) <= 1 and r.get("b_label") == "cichlid_specific_duplication":
                flags.append("one_species_only")
        except (TypeError, ValueError):
            pass
        r["risk_flags"] = ";".join(flags)
        r["n_risk_flags"] = len(flags)
