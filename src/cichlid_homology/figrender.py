"""Drawing a clade record (0.21.0). The figure is made from the record alone: no tree, no table, no coordinator.

Panel A is the collapsed tree: branches are solid when the bootstrap is 70 or more, grey from 30, dashed below 30, and the bootstrap is written at the clade's end of the branch. A star marks
the split; a hollow star marks any other node Ensembl calls a duplication. Panel B shows the support of the nodes that matter (the split, the lineage crowns, where a non-teleost gene
attaches), or the triage numbers when the tree carries no support. Panel C is the copies per species outside the split and in each lineage, in species-tree order, with our own species after
a gap. Panel D lists our own genes of the clade with the lineage of their position, the weight of their placement and their calls."""
from __future__ import annotations

from pathlib import Path

INK, MUTED, GREY, PURPLE = "#1f1f1f", "#7a7872", "#6b6a66", "#8e44ad"
COLORS = {"A": "#2a78d6", "B": "#e0642a", "C": "#3a9d5d", "D": "#b5469b", "E": "#b8960b", "F": "#169c9c", "G": "#a04060", "H": "#607080"}
TELEOST_GROUPS = ("early teleost", "percomorph", "ovalentarian", "cichlid")
XTIP, TRI_W = 9.0, 2.2
ROW_IN = 0.40


def _color(label):
    return COLORS.get(label, MUTED)


def _flatten(t):
    rows = []

    def walk(n, d):
        n["_d"] = d
        if n["k"] == "node":
            for c in n["kids"]:
                walk(c, d + 1)
            n["_y"] = sum(c["_y"] for c in n["kids"]) / len(n["kids"])
        else:
            n["_y"] = len(rows)
            rows.append(n)
    import sys
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    walk(t, 0)
    return rows


def _tip_label(n, rec):
    anchors = rec.get("anchors", {})
    if n["k"] == "leaf":
        name = n["name"]
        if n.get("own"):
            g = next((x for x in rec["genes"] if x["gene"] == name), {})
            sym = g.get("calls", {}).get("symbol") or ""
            w = g.get("lwr")
            return f"{sym + ' ' if sym else ''}{name}" + (f"  (placement weight {float(w):.2f})" if w not in (None, "") else ""), True
        sp = (n.get("species") or "").replace("_", " ").capitalize()
        return sp + (f", {n['symbol']}" if n.get("symbol") else ""), False
    lab = n["label"]
    pl = lambda k, w: f"{k} {w}{'' if k == 1 else 's'}"
    spn = lambda k: f"{k} species"
    own = ""
    if n.get("own"):
        sym = {g["gene"]: (g.get("calls", {}).get("symbol") or g["gene"]) for g in rec["genes"]}
        own = "\n   \u2605 " + ", ".join(sym.get(o, o) for o in n["own"][:4]) + (" ..." if len(n["own"]) > 4 else "")
    if lab == "tetrapod":
        return f"{pl(n['n_genes'], 'tetrapod gene')} in {spn(n['n_species'])}", False
    if lab in ("stray", "other teleost"):
        sp = ", ".join(n["species"]) if n.get("species") else spn(n["n_species"])
        return f"{pl(n['n_genes'], 'stray teleost gene' if lab == 'stray' else 'other teleost gene')} ({sp})" + own, bool(own)
    z = ", ".join(a["symbol"] or a["gene"] for a in anchors.get(lab, {}).get("zebrafish", []))
    return f"Lineage {lab}: {n['n_species']} teleost species" + (f" (zebrafish {z})" if z else "") + own, True


def render_record(rec: dict, path, dpi: int = 110, fmt: str = "png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle

    path = Path(path)
    tree = rec["tree"]
    rows = _flatten(tree)
    R = len(rows)
    inner = [n for n in _iter(tree) if n["k"] == "node"]
    maxd = max([n["_d"] for n in inner] + [1])
    xs = min(1.05, (XTIP - TRI_W - 0.7) / maxd)
    xn = lambda n: XTIP if n["k"] != "node" else n["_d"] * xs
    xend = lambda n: (XTIP - TRI_W) if n["k"] == "block" else xn(n)
    W = 14.5
    hA = max(3.0, ROW_IN * R + 1.0)
    items = rec.get("support_items") or []
    genes = rec.get("genes") or []
    hBax = (0.52 * max(len(items), 1) + 0.6) if (rec.get("support_available") and items) else 1.5
    nl = sum(2 + (1 if (x.get("calls", {}).get("b1_category") or x.get("calls", {}).get("ohnolog_status")) else 0) + (1 if (x.get("calls", {}).get("human_symbol") or x.get("calls", {}).get("zebrafish_symbol")) else 0) for x in genes[:6]) or 1
    hD = 0.17 * nl + 0.3
    Hz = max(hA, 0.45 + hBax + 0.75 + hD, 4.0)
    Hs = 3.7
    HEAD = 1.5
    H = HEAD + Hz + Hs + 0.45
    fig = plt.figure(figsize=(W, H), dpi=dpi)
    fig.patch.set_facecolor("white")
    ax = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])
    # ---------------- header
    sp_ = rec.get("split")
    head = f"{rec['clade']}   ({rec['method']}, {rec['tree_mode']} family {rec['family']})"
    fig.text(0.015, 1 - 0.28 / H, head, fontsize=14, fontweight="bold", ha="left", va="center", color=INK)
    if sp_:
        sub = (f"{sp_['kind'].capitalize()}-level split into {sp_['n_lineages']} lineages ({', '.join(str(l['species']) for l in rec['lineages'])} teleost species); "
               f"overlap {sp_['overlap']:.2f}, {rec['triage']['shared_species']} shared species; {rec['triage']['outgroup_inside']} non-teleost gene(s) attach inside a lineage.")
    elif rec["mode"] == "single":
        sub = "No lineage split: the teleost genes form one lineage (one copy per species, or species-level duplicates only)."
    else:
        sub = "No teleost gene in this clade."
    fig.text(0.015, 1 - 0.60 / H, sub, fontsize=9.3, ha="left", va="center", color=MUTED)
    d = rec.get("deep_order") or {}
    if d.get("tested"):
        names = {"shark_outside": "shark outside the bony fishes", "coelacanth_with_tetrapods": "coelacanth with the tetrapods", "reedfish_outside_gar_teleosts": "reedfish outside gar and teleosts",
                 "gar_outside_teleosts": "gar outside the teleosts"}
        failed = [names[k] for k in names if d.get(k) is False]
        fig.text(0.015, 1 - 0.88 / H, f"Known deep order recovered: {d['passed']} of {d['tested']} tested" + (f" (not: {'; '.join(failed)})" if failed else ""), fontsize=8.6, ha="left", va="center",
                 color=("#b00020" if failed else MUTED))
    if rec["flags"]:
        fig.text(0.015, 1 - 1.14 / H, "flags: " + ", ".join(rec["flags"]), fontsize=8.2, ha="left", va="center", color="#b00020")
    # ---------------- panel A: the collapsed tree
    axA = ax(0.2, 0.45 + Hs, 9.4, Hz)
    axA.axis("off")
    axA.set_xlim(-0.3, 17.5)
    axA.set_ylim(R - 0.4 + 0.2, -0.8)
    for n in _iter(tree):
        if n["k"] != "node":
            continue
        n["_x"] = xn(n)
        for c in n["kids"]:
            axA.plot([n["_x"], n["_x"]], [n["_y"], c["_y"]], color=INK, lw=1.2, alpha=0.9, solid_capstyle="butt", zorder=2)
            b = c.get("support")
            lw, ls, al = (1.8, "-", 1) if b is None or b >= 70 else ((1.5, "-", 0.75) if b >= 30 else (1.4, (0, (3, 2)), 0.75))
            col = _color(c["label"]) if c["k"] == "block" and c["label"] in COLORS else (PURPLE if c["k"] == "leaf" and c.get("group") in ("cartilaginous fish", "coelacanth", "gar or reedfish") else (GREY if c["k"] == "block" else INK))
            if c["k"] == "leaf" and c.get("own"):
                col = "#c79200"
            axA.plot([n["_x"], xend(c)], [c["_y"], c["_y"]], color=col, lw=lw, ls=ls, alpha=al, solid_capstyle="butt", zorder=2)
            if b is not None and (c["k"] == "node" or c["k"] == "block"):
                axA.text(xend(c) + 0.1, c["_y"] - 0.2, str(b), fontsize=6.8, ha="left", va="bottom", color=("#b00020" if b < 30 else (INK if b >= 70 else MUTED)),
                         fontweight="bold" if (b < 30 or b >= 70) else None, zorder=7)
    for n in rows:
        y = n["_y"]
        txt, bold = _tip_label(n, rec)
        if n["k"] == "leaf":
            own = n.get("own")
            grp = n.get("group")
            col = "#c79200" if own else (PURPLE if grp in ("cartilaginous fish", "coelacanth", "gar or reedfish") else GREY)
            axA.plot([XTIP], [y], marker="*" if own else "o", ms=11 if own else 5.5, color=col, markeredgecolor=INK if own else col, markeredgewidth=0.6, zorder=4)
            axA.text(XTIP + 0.3, y, txt, va="center", ha="left", fontsize=7.8, color=INK, fontweight="bold" if own else None)
        else:
            h = min(0.46, 0.05 * (n["n_genes"] ** 0.5) + 0.06)
            col = _color(n["label"]) if n["label"] in COLORS else GREY
            axA.add_patch(Polygon([[XTIP - TRI_W, y], [XTIP, y - h], [XTIP, y + h]], closed=True, facecolor=col, edgecolor=col, alpha=0.85, zorder=3))
            axA.text(XTIP + 0.3, y, txt, va="center", ha="left", fontsize=7.8, color=INK, fontweight="bold" if bold else None, linespacing=1.15)
    for n in inner:
        if n.get("split"):
            axA.plot([n["_x"]], [n["_y"]], marker="*", ms=17, color="#f2c230", markeredgecolor=INK, markeredgewidth=0.8, zorder=6)
        elif n.get("dup") == "Y":
            axA.plot([n["_x"]], [n["_y"]], marker="*", ms=9, color="white", markeredgecolor=INK, markeredgewidth=0.7, zorder=6)
    fig.text(0.2 / W, (0.45 + Hs + Hz + 0.08) / H, "A   Collapsed tree", fontsize=10.5, fontweight="bold", ha="left", va="bottom", color=INK)
    # ---------------- panel B: support
    xb = 10.3
    ztop = 0.45 + Hs + Hz
    axB = ax(xb, ztop - 0.35 - hBax, W - xb - 0.35, hBax)
    axB.set_xlim(0, 100)
    n_it = max(len(items), 1)
    axB.set_ylim(n_it, -0.1)
    axB.set_yticks([])
    for s_ in ("top", "right", "left"):
        axB.spines[s_].set_visible(False)
    if rec.get("support_available") and items:
        for i, it in enumerate(items):
            c = {"split": "#d9a400", "crown": _color(it["label"].split()[1].rstrip(":")), "outgroup": PURPLE}.get(it["kind"], MUTED)
            v = it["value"]
            axB.text(0, i + 0.05, it["label"], va="top", ha="left", fontsize=8.2, color=INK)
            if v is None:
                axB.text(1, i + 0.62, "no value", va="center", fontsize=8, color=MUTED)
            else:
                axB.barh(i + 0.62, max(v, 0.6), color=c, alpha=0.9, height=0.3)
                axB.text(max(v, 0.6) + 1.5, i + 0.62, str(v), va="center", fontsize=8.8, fontweight="bold", color=INK)
        axB.axvline(70, color=INK, lw=0.9, ls=(0, (4, 3)))
        axB.set_xticks([0, 25, 50, 75, 100])
        axB.tick_params(labelsize=7.5)
        axB.set_xlabel("bootstrap of the clade (Ensembl)", fontsize=8)
    else:
        axB.axis("off")
        t = rec["triage"]
        lines = ["No node support is available for this tree." if not rec.get("support_available") else "No split: nothing to score."]
        if rec["split"]:
            lines += [f"overlap of the lineages: {t['overlap']}", f"stem ratio of the lineages: {t['stem_ratio'] if t['stem_ratio'] is not None else '-'}"]
        for i, l in enumerate(lines):
            axB.text(0, 0.35 + i * 0.42, l, fontsize=8.6, color=INK if i else MUTED, va="center")
        axB.set_ylim(3, -0.1)
    fig.text(xb / W, (ztop + 0.08) / H, "B   Support of the nodes that matter", fontsize=10.5, fontweight="bold", ha="left", va="bottom", color=INK)
    # ---------------- panel D: our genes
    yDt = ztop - 0.35 - hBax - 0.75
    yD = yDt - hD
    axD = ax(xb, yD, W - xb - 0.35, hD)
    axD.axis("off")
    axD.set_xlim(0, 1)
    axD.set_ylim(0, 1)
    gl = []
    for g in genes[:6]:
        c = g.get("calls", {})
        sym = c.get("symbol") or g["gene"]
        lw = g.get("lwr")
        gl.append((f"{sym}  ({g['gene']})", True))
        gl.append((f"lineage {g['lineage']}" + (f", placement weight {float(lw):.2f}" if lw not in (None, "") else ""), False))
        call = "; ".join(x for x in (c.get("b1_category"), (c.get("ohnolog_status") + (f" with {c['ohnolog_symbol']}" if c.get("ohnolog_symbol") else "")) if c.get("ohnolog_status") else "") if x)
        if call:
            gl.append((call, False))
        nm = ", ".join(x for x in ((f"human {c['human_symbol']}" if c.get("human_symbol") else ""), (f"zebrafish {c['zebrafish_symbol']}" if c.get("zebrafish_symbol") else "")) if x)
        if nm:
            gl.append((nm, False))
    if not gl:
        gl = [("no gene of the focal species in this clade", False)]
    top = 1.0
    for txt, bold in gl[:30]:
        axD.text(0, top, txt, fontsize=7.8, va="top", ha="left", color=INK if bold else MUTED, fontweight="bold" if bold else None)
        top -= 0.17 / hD
    fig.text(xb / W, (yDt + 0.08) / H, "D   Our genes in this clade", fontsize=10.5, fontweight="bold", ha="left", va="bottom", color=INK)
    # ---------------- panel C: copies per species
    st = rec["strip"]
    fs = rec["focal_strip"]
    cols = [(s, l, g, c, False) for s, l, g, c in zip(st["species"], st["labels"], st["groups"], st["counts"])]
    gap = 1 if fs["species"] else 0
    n_cols = len(cols) + gap + len(fs["species"])
    letters = rec["letters"]
    k = len(letters)
    has_stray = any(c[-1] for c in st["counts"]) or any(c[-1] for c in fs["counts"])
    rws = [("outside", "outside the split", GREY, 0)] + [(L, f"lineage {L}" + (f" ({', '.join(a['symbol'] or a['gene'] for a in rec['anchors'].get(L, {}).get('zebrafish', []))[:22]})" if rec['anchors'].get(L, {}).get('zebrafish') else ""), _color(L), 1 + i) for i, L in enumerate(letters)]
    if has_stray:
        rws.append(("stray", "stray teleost genes", MUTED, k + 1))
    axC = ax(2.2, 0.3, W - 2.4, Hs - 0.3)
    axC.axis("off")
    axC.set_xlim(0, n_cols + 0.4)
    axC.set_ylim(len(rws) + 3.4, -1.0)
    sm = rec["mode"] in ("split", "single")
    allcols = cols + ([None] * gap) + [(s, l, "cichlid", c, o) for s, l, c, o in zip(fs["species"], fs["labels"], fs["counts"], fs["own"])]
    for i, col in enumerate(allcols):
        if col is None:
            continue
        s, lab_, grp, cnt, own = col
        x = i + 0.5
        is_tel = grp in TELEOST_GROUPS
        for r, (key, _t, c, idx) in enumerate(rws):
            v = cnt[idx] if idx < len(cnt) else 0
            y = r * 1.05 + 0.55
            expected = (key == "outside" and not is_tel) or (key not in ("outside", "stray") and is_tel and sm)
            if v > 0:
                axC.add_patch(Rectangle((x - 0.4, y - 0.4), 0.8, 0.8, fc=c, ec=c, lw=0.8))
                if v >= 2:
                    axC.text(x, y + 0.02, str(v), ha="center", va="center", fontsize=6.5, fontweight="bold", color="white")
            elif expected:
                axC.add_patch(Rectangle((x - 0.4, y - 0.4), 0.8, 0.8, fc="white", ec="#b3b1aa", lw=0.8, ls=(0, (2, 2))))
                axC.plot([x - 0.2, x + 0.2], [y - 0.2, y + 0.2], color="#b3b1aa", lw=0.7)
                axC.plot([x - 0.2, x + 0.2], [y + 0.2, y - 0.2], color="#b3b1aa", lw=0.7)
        axC.text(x + 0.35, len(rws) * 1.05 + 0.5, lab_, rotation=62, ha="right", va="top", fontsize=6.3, color=INK, style="normal" if not own and grp == "cichlid" and col[4] is not None and col in allcols[len(cols) + gap:] else "italic",
                 fontweight="bold" if own else None, rotation_mode="anchor")
    for r, (key, t, c, idx) in enumerate(rws):
        axC.text(-0.3, r * 1.05 + 0.55, t, ha="right", va="center", fontsize=7.8, color=c if key != "outside" else GREY, fontweight="bold")
    fig.text(0.2 / W, (0.3 + Hs - 0.3 + 0.18) / H, "C   Copies per species: outside the split and in each lineage (a hollow cross is a missing copy; a number is more than one gene)", fontsize=10.5,
             fontweight="bold", ha="left", va="bottom", color=INK)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, facecolor="white", format=fmt)
    plt.close(fig)
    return path


def _iter(t):
    stack = [t]
    while stack:
        n = stack.pop()
        yield n
        if n["k"] == "node":
            stack.extend(n["kids"])


def _render_one(args):
    rec, path, dpi, fmt = args
    try:
        render_record(rec, path, dpi, fmt)
        return str(path), None
    except Exception as e:                                       # one bad clade must not stop the batch
        import traceback
        return str(path), traceback.format_exc()[-600:]


def render_many(jobs: list, workers: int = 1, dpi: int = 110, fmt: str = "png") -> list:
    """jobs: [(record, path)]. Returns [(path, error or None)]."""
    args = [(r, p, dpi, fmt) for r, p in jobs]
    if workers <= 1 or len(args) < 2:
        return [_render_one(a) for a in args]
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(_render_one, args, chunksize=4))
