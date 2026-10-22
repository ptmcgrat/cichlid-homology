"""Plain-text evidence for one gene's call: what the walk up the tree saw. For hand checks."""
from __future__ import annotations

from .orthology import RTree, relationship
from .tree import Node


def explain_gene(rt: RTree, g: Node, ref: str, own: str, t: float, clades=None,
                 extra_levels: int = 1, max_names: int = 6) -> list[str]:
    """Walk from the gene towards the root, nearest ancestor first, until the first ancestor that
    holds a `ref` leaf, plus `extra_levels` more. One block per sibling subtree."""
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
        lines.append(f"-- ancestor with {a.hi - a.lo} leaves, {bin(a.mask).count('1')} species, "
                     f"{len(a.children)} children")
        for c0, c, spec in items:
            sib = rt.leaves_in(c)
            names = [x.name for x in sib][:max_names]
            more = f" (+{len(sib) - max_names} more)" if len(sib) > max_names else ""
            has_ref = any(x.sp == ref for x in sib)
            over = rt.names_of_mask(c0.mask & c.mask)
            lines.append(f"   sibling: {', '.join(names)}{more}")
            lines.append(f"     {'SPECIATION' if spec else 'DUPLICATION'}"
                         f"{'' if spec else ' (shared species: ' + ', '.join(over[:5]) + ')'}"
                         f"{'; holds ' + ref if has_ref else ''}")
            seen_ref = seen_ref or has_ref
    return lines
