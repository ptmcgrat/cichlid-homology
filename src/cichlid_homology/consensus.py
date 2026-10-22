"""Agreement labels across methods (section 6 of the sheet): agree, majority, split, single, none."""
from __future__ import annotations

from collections import Counter


def label_agreement(values: dict):
    """values: {method: value or None/''}. Returns (label, value or None, n_answering)."""
    present = {m: v for m, v in values.items() if v not in (None, "")}
    n = len(present)
    if n == 0:
        return "none", None, 0
    if n == 1:
        return "single", next(iter(present.values())), 1
    top, cnt = Counter(present.values()).most_common(1)[0]
    if cnt == n:
        return "agree", top, n
    if cnt * 2 > n:
        return "majority", top, n
    return "split", None, n
