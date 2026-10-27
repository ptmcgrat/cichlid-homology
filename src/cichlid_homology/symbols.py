"""An independent check on the reference gene: does the closest human / zebrafish gene have the same
symbol as the gene's own NCBI symbol (from the proteome's genes table)?

NCBI names its genes partly by orthology and synteny, so this is evidence, not proof. Many cichlid genes
carry only a placeholder (LOC...), and those cannot be checked. Grades, strongest first:
    exact   same symbol, ignoring case
    stem    same after removing one trailing A or B (the ohnolog suffix: pax6a, pax6b, PAX6)
    none    named, but no reference symbol matches
    unnamed / no_reference_symbol   cannot be checked
`stem` also matches different paralogs with an a/b suffix (CDKN1A, CDKN1B), so report it separately.
"""
from __future__ import annotations

import csv
import sys

csv.field_size_limit(sys.maxsize)


def is_named(sym: str) -> bool:
    s = (sym or "").strip()
    return bool(s) and not s.upper().startswith("LOC") and s != "-"


def stem(sym: str) -> str:
    s = sym.upper()
    return s[:-1] if len(s) > 3 and s[-1] in "AB" else s


def match_grade(own: str, candidates: set) -> str:
    if not is_named(own):
        return "unnamed"
    cands = {c for c in candidates if c}
    if not cands:
        return "no_reference_symbol"
    o = own.upper()
    if any(o == c.upper() for c in cands):
        return "exact"
    if any(stem(own) == stem(c) for c in cands):
        return "stem"
    return "none"


def own_symbols(layout, prefix: str) -> dict:
    """{gene_id: symbol} from the species' genes table (through the coordinator's layout)."""
    out = {}
    with open(layout.genes_table(prefix).local, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["gene_id"]] = r.get("symbol", "")
    return out


def attach_symbol_check(rows: list, own: dict, refs=("human", "zebrafish")) -> None:
    """Adds own_symbol and `<ref>_symbol_match`. Candidates: the closest gene and every reference
    ortholog (needs names.attach_names to have run)."""
    for r in rows:
        sym = own.get(r["gene_id"], "")
        r["own_symbol"] = sym
        for key in refs:
            cands = set(r.get(f"_{key}_symbols", set()))
            if r.get(f"{key}_closest_symbol"):
                cands.add(r[f"{key}_closest_symbol"])
            r[f"{key}_symbol_match"] = match_grade(sym, cands)


def matched(grade: str) -> bool:
    return grade in ("exact", "stem")


def attach_symbol_check_a(rows: list, own: dict, refs=("human", "zebrafish")) -> None:
    """For the 0.7.0 counts. `a_<ref>_symbol_match` compares the gene's own symbol with the best reference
    gene only; `a_<ref>_symbol_match_any` with every reference gene in the clade."""
    for r in rows:
        sym = own.get(r["gene_id"], "")
        r.setdefault("own_symbol", sym)
        for key in refs:
            best = r.get(f"a_{key}_best_symbol", "")
            r[f"a_{key}_symbol_match"] = match_grade(sym, {best} if best else set())
            anyset = set(r.get(f"_a_{key}_symbols", set()))
            if best:
                anyset.add(best)
            r[f"a_{key}_symbol_match_any"] = match_grade(sym, anyset)
