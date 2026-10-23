"""Plain-text evidence for one gene's call: what the walk up the tree saw. For hand checks."""
from __future__ import annotations

from .dating import date_mask, level_of_species
from .orthology import RTree, relationship
from .tree import Node


def _outermost(rt: RTree, mask: int, clades, level) -> list[str]:
    """The species that set the level of a node (those at the outermost level present)."""
    if clades is None or level is None:
        return []
    want = "cyclostome" if level == "cyclostome_disputed" else level
    return [s for s in rt.names_of_mask(mask) if level_of_species(clades, s) == want][:3]


def explain_gene(rt: RTree, g: Node, ref: str, own: str, t: float, clades=None,
                 extra_levels: int = 1, max_names: int = 6) -> list[str]:
    """Walk from the gene towards the root, nearest ancestor first, until the first ancestor that
    holds a `ref` leaf, plus `extra_levels` more. One block per ancestor; for each sibling subtree:
    its leaves, speciation or duplication, and the species that overlap. Each ancestor also shows the
    ladder level of the species below it and the species that set that level."""
    rel = relationship(rt, g, ref, own, t)
    lines = [f"gene {g.name} ({g.sp}); reference {ref}: {len(rt.by_species.get(ref, []))} gene(s) in this tree",
             f"relationship {rel['rel']}; closest {rel['closest'] or '-'}"
             f" (ortholog: {rel['closest_is_ortholog']}); n_ref {rel['n_ref']}, n_focal {rel['n_focal']}"]
    groups: list[tuple] = []
    for (a, c0, c, spec) in rt.walk(g, t):
        if not groups or groups[-1][0] is not a:
            groups.append((a, []))
        groups[-1][1].append((c0, c, spec))
    seen_ref, after = False, 0
    for a, items in groups:
        if seen_ref:
            after += 1
            if after > extra_levels:
                break
        lvl = date_mask(rt, a.mask, clades) if clades is not None else None
        why = _outermost(rt, a.mask, clades, lvl)
        lines.append(f"-- ancestor: {a.hi - a.lo} leaves, {bin(a.mask).count('1')} species, "
                     f"{len(a.children)} children"
                     + (f"; level {lvl}" if lvl else "")
                     + (f" (set by {', '.join(why)})" if why else ""))
        for c0, c, spec in items:
            sib = rt.leaves_in(c)
            names = [x.name for x in sib][:max_names]
            more = f" (+{len(sib) - max_names} more)" if len(sib) > max_names else ""
            has_ref = any(x.sp == ref for x in sib)
            over = rt.names_of_mask(c0.mask & c.mask)
            lines.append(f"   sibling: {', '.join(names)}{more}")
            lines.append(f"     {'SPECIATION' if spec else 'DUPLICATION'}"
                         f"{'' if spec else ' (shared species: ' + ', '.join(over[:6]) + ')'}"
                         f"{'; holds ' + ref if has_ref else ''}")
            seen_ref = seen_ref or has_ref
    return lines
