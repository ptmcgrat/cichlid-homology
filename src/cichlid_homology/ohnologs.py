"""Published teleost-duplication pairs (OHNOLOGS v2, Singh and Isambert 2020, NAR 48:D724) as independent evidence.

The pairs come from macro-synteny against outgroups, so they do not depend on our gene trees or counts. The
files are given by the owner on the command line (a folder holding `<species>_Pairs_<Strict|Intermediate|Relaxed>_3R.txt`,
with `_` or `.` between the parts
for drerio, olatipes, gaculeatus, tnigroviridis); nothing is downloaded here. Only protein-coding pairs dated to the
fish genome duplication are used. Genes are matched by symbol: the zebrafish best gene of the family (reliable
names only) for zebrafish, and the gene's own annotated symbol for medaka, stickleback and Tetraodon. The published
lists are incomplete (about 40% of zebrafish genes named with an a or b are in them), so a gene in no pair is not
evidence against a pair; a gene in a pair is evidence for one."""
from __future__ import annotations

import csv
import re
from pathlib import Path

SPECIES = {"drerio": "zebrafish", "olatipes": "medaka", "gaculeatus": "stickleback", "tnigroviridis": "tetraodon"}
LEVELS = ("Strict", "Relaxed")
SUPPORT_COLS = ["t_pairs_relaxed_n", "t_pairs_strict_n", "t_pairs_species"]


def _symbols(path) -> set:
    out = set()
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("Gene type") != "protein_coding" or r.get("Duplication time") != "FishWGD":
                continue
            for k in ("Symbol1", "Symbol2"):
                s = (r.get(k) or "").strip().lower()
                if s and "(" not in s and not s.startswith("ens"):
                    out.add(s)
    return out


_NAME = re.compile(r"^(?P<sp>" + "|".join(SPECIES) + r")[._]Pairs[._](?P<lvl>Strict|Intermediate|Relaxed)[._]3R\.txt$", re.I)


def find_files(directory) -> dict:
    """{(species code, level): path} for the pair files in a folder; `_` or `.` may separate the parts."""
    d = Path(directory).expanduser()
    if not d.is_dir():
        raise FileNotFoundError(f"{d} is not a folder")
    out = {}
    for f in sorted(d.iterdir()):
        m = _NAME.match(f.name)
        if m:
            out[(m["sp"].lower(), m["lvl"].capitalize())] = f
    return out


def load_pair_sets(directory) -> dict:
    """{species: {"strict": set, "relaxed": set}} for the files found in the folder; {} when there are none."""
    files = find_files(directory)
    sets = {}
    for code, name in SPECIES.items():
        got = {lvl.lower(): _symbols(files[(code, lvl)]) for lvl in LEVELS if (code, lvl) in files}
        if got:
            sets[name] = {"strict": got.get("strict", set()), "relaxed": got.get("relaxed", got.get("strict", set()))}
    return sets


def support(row: dict, sets: dict) -> dict:
    """How many of the four species have this gene in a published pair (relaxed and strict) and which."""
    if not sets:
        return {c: "" for c in SUPPORT_COLS}
    syms = {}
    z = (row.get("a_zebrafish_best_symbol") or "").lower()
    if z and str(row.get("a_zebrafish_best_name_reliable")) == "True":
        syms["zebrafish"] = z
    own = (row.get("own_symbol") or "").lower()
    if own and not own.startswith("loc"):
        for sp in ("medaka", "stickleback", "tetraodon"):
            syms[sp] = own
    rel = [sp for sp, s in syms.items() if sp in sets and s in sets[sp]["relaxed"]]
    strict = [sp for sp, s in syms.items() if sp in sets and s in sets[sp]["strict"]]
    return {"t_pairs_relaxed_n": len(rel), "t_pairs_strict_n": len(strict), "t_pairs_species": ";".join(rel)}
