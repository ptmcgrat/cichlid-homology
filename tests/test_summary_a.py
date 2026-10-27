import csv
from pathlib import Path

from cichlid_homology.names import attach_names_a
from cichlid_homology.report import write_report
from cichlid_homology.summary_a import build, overall, symbol_summary
from cichlid_homology.symbols import attach_symbol_check_a


class Rec:
    def __init__(self, symbol, oid, status="current"):
        self.symbol, self.official_id, self.status = symbol, oid, status
        self.reliable = status in ("current", "grch37", "legacy_id")


class Table:
    def __init__(self, d):
        self.d = d

    def get(self, p):
        return self.d.get(p)


def row(gene, fam, h, z, own="", fsize=2, hbest="", zbest="", hleaves=(), zleaves=(), b="single", hn=(1, 1), zn=(1, 1)):
    return {"gene_id": gene, "family": fam, "family_mzeb_genes": fsize, "own_symbol": own,
            "a_human_cat": h, "a_human_ref_n": hn[0], "a_human_focal_n": hn[1], "a_human_clade_leaves": 9,
            "a_human_best": hbest, "a_human_tie": False, "_a_human_ref_leaves": list(hleaves),
            "a_zebrafish_cat": z, "a_zebrafish_ref_n": zn[0], "a_zebrafish_focal_n": zn[1], "a_zebrafish_clade_leaves": 5,
            "a_zebrafish_best": zbest, "a_zebrafish_tie": False, "_a_zebrafish_ref_leaves": list(zleaves),
            "b_label": b, "b_label_raw": b, "b_pattern": "....", "b_other_share2p": 0.0, "b_cichlid_share2p": 0.0,
            "b_cichlid_with2p": 0, "b_cichlid_n": 10}


TABLE = Table({"H1": Rec("PAX6", "HGNC:1"), "H2": Rec("PAX7", "HGNC:2"), "Z1": Rec("pax6a", "ZDB-1"),
               "Z2": Rec("pax6b", "ZDB-2")})


def make():
    tf = [row("g1", "TF1", "OneToOne", "OneToOne", "pax6a", hbest="H1", zbest="Z1", hleaves=["H1"], zleaves=["Z1"]),
          row("g2", "TF2", "OneToTwo", "OneToOne", "pax7", hbest="H2", zbest="Z1", hleaves=["H2"], zleaves=["Z1"], hn=(1, 2)),
          row("g3", "TF3", "None", "OneToOne", "", zbest="Z1", zleaves=["Z1"]),
          row("g4", "TF4", "OneToOne", "None", "pax6a", hbest="H1", hleaves=["H1"])]
    cp = [row("g1", "C1", "OneToOne", "OneToOne", "pax6a", hbest="H1", zbest="Z1", hleaves=["H1"], zleaves=["Z1"]),
          row("g2", "C2", "OneToOne", "OneToOne", "pax7", hbest="H2", zbest="Z2", hleaves=["H2"], zleaves=["Z2"]),
          row("g3", "C3", "None", "None", ""),
          row("g4", "C4", "OneToOne", "OneToOne", "pax6a", hbest="H2", zbest="Z2", hleaves=["H2"], zleaves=["Z2"])]
    for rows in (tf, cp):
        attach_names_a(rows, TABLE)
        attach_symbol_check_a(rows, {r["gene_id"]: r["own_symbol"] for r in rows})
    return {"TreeFam": tf, "Compara": cp}


def test_overall_numbers_and_family_size_breakdown():
    o = overall(make()["TreeFam"], "human")
    assert o["genes"] == 4 and o["categories"]["OneToOne"] == {"n": 2, "pct": 50.0}
    assert o["categories"]["None"]["n"] == 1 and o["by_family_size"]["2-3"]["OneToOne"] == 2


def test_agreement_between_methods(tmp_path):
    s = build(make(), tmp_path)
    h = s["agreement"]["human"]
    assert h["genes_in_both_methods"] == 4 and h["both_gave_a_category"] == 3          # g3 is None in both
    assert h["same_category"] == 2 and h["same_category_pct"] == 66.7                 # g2 differs
    assert h["only_one_method_has_a_reference_gene"] == 0 and h["neither_has_one"] == 1
    z = s["agreement"]["zebrafish"]
    assert z["only_one_method_has_a_reference_gene"] == 2                              # g3 and g4


def test_best_gene_agreement_and_adjudication(tmp_path):
    s = build(make(), tmp_path)
    b = s["best_gene"]["human"]
    assert b["same_best_gene"] == 2 and b["different_gene"] == 1                       # g4: H1 against H2
    assert b["same_best_gene_pct_of_comparable"] == 66.7
    assert b["adjudication_by_own_ncbi_symbol"] == {"different_gene | TreeFam_matches_own_symbol": 1}
    rows = list(csv.DictReader(open(tmp_path / "pilot_best_gene_differences.human.tsv"), delimiter="\t"))
    assert [r["gene_id"] for r in rows] == ["g4"]
    z = s["best_gene"]["zebrafish"]
    assert z["same_best_gene"] == 1 and z["different_gene"] == 1                       # g2: Z1 against Z2


def test_symbol_summary_by_category():
    t = symbol_summary(make()["TreeFam"])
    assert t["genes_with_a_real_own_symbol"] == 3
    h = t["human"]["OneToOne"]
    assert h["checkable"] == 2 and h["exact_or_stem_pct"] == 100.0                      # pax6a against PAX6, twice


def test_category_splits_file_and_report(tmp_path):
    rows = make()
    A = build(rows, tmp_path)
    assert A["category_splits_written"] == 1                                            # only g2 differs (human)
    split = list(csv.DictReader(open(tmp_path / "pilot_category_splits.tsv"), delimiter="\t"))
    assert split[0]["gene_id"] == "g2" and split[0]["TreeFam_cat"] == "OneToTwo" and split[0]["Compara_cat"] == "OneToOne"
    write_report({"cichlid_homology": "test", "species": "Mzeb", "A": A}, tmp_path / "pilot_report.md")
    text = (tmp_path / "pilot_report.md").read_text()
    for h in ("## Overall numbers", "### Human (A1): categories", "## Agreement between the methods",
              "## NCBI symbol check", "## Goal B"):
        assert h in text


def test_rows_from_the_calling_code_fill_the_columns_and_feed_the_report(tmp_path):
    from helpers import FakeClades, PREFIXES, st
    from test_v06 import FakeLS, fam
    from cichlid_homology.dating import Params
    from cichlid_homology.pilot import A_COLS, COLUMNS, process_family
    newick = ("((Hsap_1,Mmus_1),(Locu_1,((Aang_1,(Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1))))),"
              "(Aang_2,(Drer_2,(Gacu_2,(Olat_2,(Mzeb_2,Acal_2))))))));")
    table = Table({"Hsap_1": Rec("PAX6", "HGNC:1"), "Drer_1": Rec("pax6a", "ZDB-1"), "Drer_2": Rec("pax6b", "ZDB-2")})
    out = {}
    for method in ("TreeFam", "Compara"):
        rows, _ = process_family(fam(newick), FakeLS(), set(), st(), FakeClades(), Params(), "Mzeb", PREFIXES, [], None)
        for r in rows:
            r["method"] = method
        attach_names_a(rows, table)
        attach_symbol_check_a(rows, {"Mzeb_1": "pax6a", "Mzeb_2": "pax6b"})
        out[method] = rows
    r = out["TreeFam"][0]
    assert r["a_human_cat"] == "OneToTwo" and r["a_zebrafish_cat"] == "OneToOne"
    assert r["a_zebrafish_best_symbol"] == "pax6a" and r["a_zebrafish_symbol_match"] == "exact"
    assert not (set(A_COLS) - set(r)), set(A_COLS) - set(r)
    with open(tmp_path / "genes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(out["TreeFam"])
    A = build(out, tmp_path)
    assert A["agreement"]["human"]["same_category_pct"] == 100.0
    assert A["goal_B"]["TreeFam"]["distribution_human_OneToTwo"]["genes_scored"] == 2
    write_report({"cichlid_homology": "test", "species": "Mzeb", "A": A}, tmp_path / "r.md")
    assert "OneToTwo" in (tmp_path / "r.md").read_text()
