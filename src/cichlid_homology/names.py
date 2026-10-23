"""Names and official IDs of the reference proteins (human, zebrafish), from the coordinator's
ReferenceNames table (TreeFam and Compara). OrthoFinder has no such table yet.

The common key across methods is the official ID (HGNC for human, ZFIN for zebrafish). A name is
`reliable` when the coordinator says so (status current, grch37 or legacy_id). A call that rests on a
replaced, legacy_symbol or unmapped protein, or on a protein with no record, gets name_confidence
"low" (the call is kept, never dropped).

The record's field names are not documented in the sheet, so `info_of` accepts attribute-style and
dict-style records and the usual column names; `cichlid-homology probe` shows what the real record has.
"""
from __future__ import annotations

RELIABLE_STATUSES = {"current", "grch37", "legacy_id"}
CAP = 10


def _field(rec, *names):
    for n in names:
        v = rec.get(n) if isinstance(rec, dict) else getattr(rec, n, None)
        if v is not None and v != "":
            return v
    return ""


def info_of(rec) -> dict:
    """Plain dict {symbol, official_id, status, reliable} from a ReferenceNames record or None."""
    if rec is None:
        return {"symbol": "", "official_id": "", "status": "absent", "reliable": False}
    status = str(_field(rec, "status"))
    rel = _field(rec, "reliable")
    if callable(rel):
        rel = rel()
    if rel == "":
        rel = status in RELIABLE_STATUSES
    return {"symbol": str(_field(rec, "symbol")), "official_id": str(_field(rec, "official_id")),
            "status": status, "reliable": bool(rel)}


def load_names(layout, part: str):
    from cichlid_data import ReferenceNames

    return ReferenceNames.load(layout, part)


def lookup(table, protein: str) -> dict:
    try:
        rec = table.get(protein)
    except Exception:
        rec = None
    return info_of(rec)


def _capped(items: list) -> str:
    return ";".join(items[:CAP]) + (f";+{len(items) - CAP}" if len(items) > CAP else "")


def attach_names(rows: list, table, refs=("human", "zebrafish")) -> None:
    """Add symbol, official ID, name status and confidence of the closest gene, and the official IDs
    and symbols of all reference orthologs, to every row. Rows must hold `<ref>_closest` and
    `_<ref>_orth_leaves` (the full list of reference orthologs)."""
    cache: dict = {}

    def info(p):
        if p not in cache:
            cache[p] = lookup(table, p)
        return cache[p]

    for r in rows:
        for key in refs:
            c = r.get(f"{key}_closest")
            if c:
                i = info(c)
                r[f"{key}_closest_symbol"] = i["symbol"]
                r[f"{key}_closest_id"] = i["official_id"]
                r[f"{key}_closest_name_status"] = i["status"]
                r[f"{key}_closest_name_reliable"] = i["reliable"]
                r[f"{key}_name_confidence"] = "ok" if i["reliable"] and i["official_id"] else "low"
            leaves = r.get(f"_{key}_orth_leaves") or []
            infos = [info(p) for p in leaves]
            ids = sorted({i["official_id"] for i in infos if i["reliable"] and i["official_id"]})
            syms = sorted({i["symbol"] for i in infos if i["symbol"]})
            r[f"_{key}_ids"] = set(ids)
            r[f"{key}_orthologs_ids"] = _capped(ids)
            r[f"{key}_orthologs_symbols"] = _capped(syms)
            r[f"{key}_orthologs_unreliable"] = sum(1 for i in infos if not i["reliable"])


def identity_agreement(ids_a: set, ids_b: set) -> str:
    """Do two methods name the same reference gene? Only reliable official IDs are compared."""
    if not ids_a or not ids_b:
        return "not_comparable"
    return "same_gene" if ids_a & ids_b else "different_gene"
