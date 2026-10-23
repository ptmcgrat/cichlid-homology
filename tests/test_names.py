import csv

from cichlid_homology.names import CAP, attach_names, identity_agreement, info_of, lookup
from cichlid_homology.pilot import _identity


class Rec:
    """Attribute-style record with `reliable` as a property."""
    def __init__(self, symbol, official_id, status):
        self.symbol, self.official_id, self.status = symbol, official_id, status

    @property
    def reliable(self):
        return self.status in ("current", "grch37", "legacy_id")


class RecMethod(Rec):
    reliable = lambda self: self.status == "current"          # noqa: E731  (a method, not a property)


class Table:
    def __init__(self, recs):
        self.recs = recs

    def get(self, protein):
        return self.recs.get(protein)


def test_info_of_accepts_attribute_dict_and_missing_records():
    assert info_of(Rec("PAX6", "HGNC:8620", "current")) == {
        "symbol": "PAX6", "official_id": "HGNC:8620", "status": "current", "reliable": True}
    d = info_of({"symbol": "pax6a", "official_id": "ZDB-GENE-1", "status": "legacy_id"})
    assert d["reliable"] is True                                # computed from the status
    assert info_of({"symbol": "x", "official_id": "y", "status": "replaced"})["reliable"] is False
    assert info_of(None)["status"] == "absent" and info_of(None)["reliable"] is False
    assert info_of(RecMethod("A", "HGNC:1", "current"))["reliable"] is True


def test_lookup_survives_a_table_that_raises():
    class Bad:
        def get(self, p):
            raise KeyError(p)
    assert lookup(Bad(), "ENSP1")["status"] == "absent"


def test_attach_names_marks_unreliable_names_low_and_keeps_the_call():
    t = Table({"H1": Rec("PAX6", "HGNC:8620", "current"), "H2": Rec("OLD", "HGNC:9", "replaced"),
               "Z1": Rec("pax6a", "ZDB-GENE-1", "legacy_id")})
    rows = [{"human_closest": "H1", "_human_orth_leaves": ["H1"],
             "zebrafish_closest": "Z1", "_zebrafish_orth_leaves": ["Z1", "Z9"]},
            {"human_closest": "H2", "_human_orth_leaves": ["H2", "H1"],
             "zebrafish_closest": "", "_zebrafish_orth_leaves": []}]
    attach_names(rows, t)
    a, b = rows
    assert a["human_closest_symbol"] == "PAX6" and a["human_name_confidence"] == "ok"
    assert a["zebrafish_name_confidence"] == "ok"                  # legacy_id counts as reliable
    assert a["zebrafish_orthologs_unreliable"] == 1               # Z9 has no record
    assert b["human_name_confidence"] == "low" and b["human_closest_name_status"] == "replaced"
    assert b["_human_ids"] == {"HGNC:8620"}                       # only reliable IDs are compared
    assert "zebrafish_name_confidence" not in b                   # no closest gene, nothing to say


def test_long_ortholog_lists_are_capped_in_the_output_but_not_in_memory():
    recs = {f"P{i}": Rec(f"G{i}", f"HGNC:{i}", "current") for i in range(25)}
    rows = [{"human_closest": "P0", "_human_orth_leaves": list(recs)}]
    attach_names(rows, Table(recs))
    assert len(rows[0]["_human_ids"]) == 25
    assert rows[0]["human_orthologs_ids"].endswith(f";+{25 - CAP}")


def test_identity_agreement():
    assert identity_agreement({"HGNC:1"}, {"HGNC:1", "HGNC:2"}) == "same_gene"
    assert identity_agreement({"HGNC:1"}, {"HGNC:2"}) == "different_gene"
    assert identity_agreement(set(), {"HGNC:2"}) == "not_comparable"


def test_identity_summary_counts_and_writes_the_differences(tmp_path):
    def row(rel, ids, fam):
        return {"human_relationship": rel, "zebrafish_relationship": "none", "family": fam,
                "_human_ids": set(ids), "_zebrafish_ids": set(),
                "human_orthologs_symbols": ",".join(ids), "human_orthologs_ids": ";".join(ids)}
    by = {"TreeFam": {"g1": row("1:1", ["A"], "TF1"), "g2": row("1:1", ["A"], "TF2"),
                      "g3": row("paralog", [], "TF3"), "g4": row("1:1", [], "TF4")},
          "Compara": {"g1": row("1:1", ["A"], "C1"), "g2": row("1:1", ["B"], "C2"),
                      "g3": row("1:1", ["A"], "C3"), "g4": row("1:1", ["A"], "C4")}}
    out = _identity(by, ["g1", "g2", "g3", "g4"], tmp_path)
    h = out["human"]
    assert h["same_gene"] == 1 and h["different_gene"] == 1 and h["not_comparable"] == 1
    assert h["not_applicable_no_ortholog_in_one_method"] == 1 and h["same_gene_rate_when_comparable"] == 0.5
    rows = list(csv.DictReader(open(tmp_path / "pilot_identity_differences.tsv"), delimiter="\t"))
    assert [r["gene_id"] for r in rows] == ["g2"]


def test_calls_feed_the_names_end_to_end():
    from helpers import FakeClades, PREFIXES, leaf, rtree
    from cichlid_homology.calls import call_gene
    from cichlid_homology.dating import Params
    rt = rtree("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Mzeb_1,Acal_1)))));")
    row = call_gene(rt, leaf(rt, "Mzeb_1"), "maylandia_zebra", Params(), FakeClades(), "denovo", PREFIXES)
    assert row["_human_orth_leaves"] == ["Hsap_1"]
    attach_names([row], Table({"Hsap_1": Rec("GENE", "HGNC:1", "current"),
                               "Drer_1": Rec("gene", "ZDB-GENE-1", "current")}))
    assert row["human_closest_symbol"] == "GENE" and row["human_closest_id"] == "HGNC:1"
    assert row["zebrafish_closest_id"] == "ZDB-GENE-1" and row["_human_ids"] == {"HGNC:1"}
