import json

from cichlid_homology import lineagepairs as L
from cichlid_homology.genetable import COLUMNS

GAR, REED = L.GAR, L.REEDFISH


def row(gid, clade, status, seqid="", idx="", zf="", partner="", symbol="", **kw):
    r = {c: "" for c in COLUMNS}
    r.update(gene_id=gid, clade_id=clade, ohnolog_status=status, seqid=seqid, gene_index=idx, zebrafish_gene=zf, ohnolog_gene=partner, symbol=symbol or gid)
    r.update(kw)
    return r


def rec(clade, lineages, lwr=0.9, anchors=None, inside=(), letters=("A", "B")):
    return {"clade": clade, "mode": "split", "split": {"kind": "teleost"}, "letters": list(letters),
            "genes": [{"gene": g, "lineage": l, "lwr": str(lwr) if not isinstance(lwr, dict) else str(lwr[g])} for g, l in lineages.items()],
            "anchors": {k: {"zebrafish": [{"gene": z, "symbol": ""} for z in v]} for k, v in (anchors or {}).items()},
            "outgroup": {"inside": [{"gene": f"o{i}", "species": s, "lineage": l} for i, (s, l) in enumerate(inside)], "outside": []}}


def pair_rows(n, c1, c2, prefix="p"):
    """n clean pairs between two chromosomes: the synteny that the other clades can lean on."""
    rows = []
    for i in range(n):
        a, b = f"{prefix}{i}a", f"{prefix}{i}b"
        for x, y, sx, sy in ((a, b, c1, c2), (b, a, c2, c1)):
            rows.append(row(x, f"{prefix}{i}:G1", "paired", sx, str(i * 100), f"Z{x}", y, ev_lineages="different", ev_position="different chromosomes", ev_zebrafish_split="yes"))
    return rows


def by(rows):
    return {r["gene_id"]: r for r in rows}


def test_a_clean_pair_keeps_its_status_and_gets_a_lineage_a_class_and_a_rebuilt_lineage_signal():
    rows = [row("a", "F:G1", "paired", "c1", "10", "ZA", "b", ev_lineages="unassigned", ev_position="different chromosomes", ev_zebrafish_split="yes"),
            row("b", "F:G1", "paired", "c2", "20", "ZB", "a", ev_lineages="unassigned", ev_position="different chromosomes", ev_zebrafish_split="yes")]
    s = L.apply(rows, {"F:G1": rec("F:G1", {"a": "A", "b": "B"})}, params={"synteny_min_pairs": 5})
    r = by(rows)
    assert (r["a"]["lineage"], r["b"]["lineage"], r["a"]["ohnolog_status"], r["a"]["duplication_class"]) == ("A", "B", "paired", L.CLASS_OHNOLOG)
    assert r["a"]["ev_lineages"] == "different" and r["a"]["ohnolog_gene"] == "b" and r["a"]["lineage_basis"] == "placement"           # the old signal was "unassigned": nested outgroups
    assert s["clades_with_a_record"] == 1


def test_an_expansion_is_a_pair_plus_an_in_paralog_chosen_by_lineage_and_synteny():
    rows = pair_rows(6, "c1", "c2") + [row("g1", "F:G1", "expansion", "c1", "50", "ZA"), row("g2", "F:G1", "expansion", "c1", "51", "ZX"), row("g3", "F:G1", "expansion", "c2", "70", "ZB")]
    s = L.apply(rows, {"F:G1": rec("F:G1", {"g1": "A", "g2": "A", "g3": "B"}, anchors={"A": ["ZA"], "B": ["ZB"]})}, params={"synteny_min_pairs": 5})
    r = by(rows)
    assert (r["g1"]["ohnolog_status"], r["g3"]["ohnolog_status"], r["g1"]["ohnolog_gene"], r["g3"]["ohnolog_gene"]) == (L.RESOLVED, L.RESOLVED, "g3", "g1")
    assert (r["g2"]["ohnolog_status"], r["g2"]["in_paralog_of"], r["g2"]["ohnolog_gene"]) == (L.INPARALOG, "g1", "")                 # the tandem neighbour of g1 is an in-paralog
    assert "the synteny tied: broken by the zebrafish anchor" in r["g1"]["lineage_basis"] and r["g1"]["duplication_class"] == L.CLASS_OHNOLOG and r["g2"]["duplication_class"] == ""
    assert s["pairs chosen among several genes"] == 1 and ["expansion", L.RESOLVED, 2] in s["status_transitions"] and ["expansion", L.INPARALOG, 1] in s["status_transitions"]


def test_without_synteny_support_the_clade_stays_ambiguous_with_its_candidates():
    rows = [row("g1", "F:G1", "expansion", "c1", "1"), row("g2", "F:G1", "expansion", "c1", "2"), row("g3", "F:G1", "expansion", "c2", "3")]
    s = L.apply(rows, {"F:G1": rec("F:G1", {"g1": "A", "g2": "A", "g3": "B"})}, params={"synteny_min_pairs": 5})
    r = by(rows)
    assert r["g1"]["ohnolog_status"] == L.AMBIGUOUS and r["g1"]["ohnolog_candidates"] == "g3" and r["g3"]["ohnolog_candidates"] == "g1;g2"
    assert s["clades left ambiguous (no synteny support)"] == 1


def test_genes_in_one_lineage_are_a_lost_partner_or_in_paralogs_of_that_lineage():
    rows = [row("a", "F1:G1", "partner ambiguous (three or more focal genes)"), row("b", "F2:G1", "expansion"), row("c", "F2:G1", "expansion")]
    L.apply(rows, {"F1:G1": rec("F1:G1", {"a": "A"}), "F2:G1": rec("F2:G1", {"b": "B", "c": "B"})})
    r = by(rows)
    assert r["a"]["ohnolog_status"] == "partner lost in the focal species" and r["b"]["ohnolog_status"] == r["c"]["ohnolog_status"] == L.ONE_LINEAGE
    assert r["a"]["lineage"] == "A" and r["b"]["lineage"] == "B"


def test_a_pair_whose_lineages_both_hold_a_gar_gene_is_older_and_does_not_count_for_the_synteny_of_teleost_pairs():
    rows = pair_rows(5, "c1", "c2") + [row("o1", "F:G1", "paired", "c1", "9", "ZO1", "o2", ev_lineages="different", ev_position="different chromosomes", ev_zebrafish_split="yes"),
                                       row("o2", "F:G1", "paired", "c2", "9", "ZO2", "o1", ev_lineages="different", ev_position="different chromosomes", ev_zebrafish_split="yes")]
    s = L.apply(rows, {"F:G1": rec("F:G1", {"o1": "A", "o2": "B"}, inside=[(GAR, "A"), (GAR, "B"), (REED, "A")])}, params={"synteny_min_pairs": 5})
    r = by(rows)
    assert r["o1"]["ohnolog_status"] == L.OLDER and r["o1"]["duplication_class"] == L.CLASS_OLDER and r["o1"]["older_evidence"] == "gar carry a gene in each lineage"
    assert r["p0a"]["ev_synteny"] == 4 and s["pairs labelled older than the teleost duplication"] == 1                                 # 5 teleost pairs on c1/c2: 4 others; the older one is not counted
    assert L.apply([row("x", "F:G1", "paired")], {})["clades_with_a_record"] == 0


def test_a_low_placement_weight_makes_the_lineage_signal_unassigned_and_lowers_the_confidence():
    def pair(lwr):
        rows = pair_rows(6, "c1", "c2") + [row("a", "F:G1", "paired", "c1", "1", "ZZ", "b", ev_lineages="unassigned", ev_position="different chromosomes", ev_zebrafish_split="no"),
                                           row("b", "F:G1", "paired", "c2", "1", "ZZ", "a", ev_lineages="unassigned", ev_position="different chromosomes", ev_zebrafish_split="no")]
        L.apply(rows, {"F:G1": rec("F:G1", {"a": "A", "b": "B"}, lwr=lwr)}, params={"synteny_min_pairs": 5})
        return by(rows)["a"]
    good, weak = pair(0.9), pair(0.3)
    assert good["ev_lineages"] == "different" and good["ohnolog_confidence"] == "high"                    # lineages, chromosomes and synteny; the zebrafish split is against nothing
    assert weak["ev_lineages"] == "unassigned" and weak["ohnolog_confidence"] == "medium"                 # one signal fewer


def test_a_clade_without_a_record_or_a_gene_outside_the_clade_node_is_left_alone():
    rows = [row("a", "F:G1", "expansion"), row("b", "F:G1", "expansion"), row("c", "F:G1", "expansion")]
    L.apply(rows, {"F:G1": rec("F:G1", {"a": "A", "b": "B"})})                                  # c is not in the record: placed outside the clade node
    r = by(rows)
    assert r["c"]["ohnolog_status"] == "expansion" and r["c"]["lineage"] == "" and r["a"]["lineage"] == "A"
    rows2 = [row("z", "F2:G1", "expansion")]
    L.apply(rows2, {"F:G1": rec("F:G1", {"a": "A"})})
    assert rows2[0]["ohnolog_status"] == "expansion" and rows2[0]["lineage"] == ""


def test_only_records_with_a_teleost_split_are_loaded(tmp_path):
    p = tmp_path / "r.jsonl"
    recs = [{"clade": "A", "mode": "split", "split": {"kind": "teleost"}}, {"clade": "B", "mode": "split", "split": {"kind": "later"}}, {"clade": "C", "mode": "single", "split": None}]
    p.write_text("\n".join(json.dumps(x) for x in recs))
    assert list(L.load_records(p)) == ["A"]
