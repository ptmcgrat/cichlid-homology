import csv
from types import SimpleNamespace

from helpers import FakeClades, PREFIX2SP, PREFIXES, st
from cichlid_homology.dating import Params
from cichlid_homology.paths import default_out
from cichlid_homology.pilot import NO_B, _goal_b, _identity, _symbol_summary, process_family
from cichlid_homology.symbols import attach_symbol_check, is_named, match_grade, own_symbols
from cichlid_homology.tree import leaves, parse_newick


def test_symbol_grades():
    assert match_grade("pax6a", {"PAX6"}) == "stem"
    assert match_grade("PAX6", {"pax6"}) == "exact"
    assert match_grade("tp53", {"TP53", "TP63"}) == "exact"
    assert match_grade("cdkn1b", {"CDKN1A"}) == "stem"            # weak: different paralogs, same stem
    assert match_grade("nanos2", {"NANOS3"}) == "none"
    assert match_grade("LOC12345", {"X"}) == "unnamed" and match_grade("", {"X"}) == "unnamed"
    assert match_grade("pax6a", set()) == "no_reference_symbol"
    assert is_named("pax6a") and not is_named("LOC1") and not is_named("-")


def test_own_symbols_reads_the_genes_table(tmp_path):
    p = tmp_path / "Mzeb.genes.tsv"
    p.write_text("gene_id\tsymbol\nMzeb_1\tpax6a\nMzeb_2\tLOC2\n")
    layout = SimpleNamespace(genes_table=lambda prefix: SimpleNamespace(local=p))
    assert own_symbols(layout, "Mzeb") == {"Mzeb_1": "pax6a", "Mzeb_2": "LOC2"}


def test_attach_symbol_check_uses_closest_and_all_orthologs():
    rows = [{"gene_id": "Mzeb_1", "human_closest_symbol": "PAX6", "_human_symbols": {"PAX6"},
             "zebrafish_closest_symbol": "", "_zebrafish_symbols": set()},
            {"gene_id": "Mzeb_2", "human_closest_symbol": "ABC", "_human_symbols": {"ABC", "XYZ"},
             "zebrafish_closest_symbol": "", "_zebrafish_symbols": set()}]
    attach_symbol_check(rows, {"Mzeb_1": "pax6a", "Mzeb_2": "xyz"})
    assert rows[0]["human_symbol_match"] == "stem" and rows[0]["zebrafish_symbol_match"] == "no_reference_symbol"
    assert rows[1]["human_symbol_match"] == "exact"            # matched through the ortholog set


def test_symbol_summary_rates():
    rows = [{"own_symbol": "a", "human_relationship": "1:1", "human_symbol_match": "exact",
             "zebrafish_relationship": "1:1", "zebrafish_symbol_match": "none"},
            {"own_symbol": "b", "human_relationship": "1:1", "human_symbol_match": "stem",
             "zebrafish_relationship": "none", "zebrafish_symbol_match": "no_reference_symbol"},
            {"own_symbol": "LOC1", "human_relationship": "1:1", "human_symbol_match": "unnamed",
             "zebrafish_relationship": "none", "zebrafish_symbol_match": "unnamed"}]
    o = _symbol_summary(rows)
    assert o["own_symbol_named"] == 2
    assert o["human"]["checkable_1:1"] == 2 and o["human"]["exact_rate_1:1"] == 0.5
    assert o["human"]["exact_or_stem_rate_1:1"] == 1.0


def test_identity_differences_are_adjudicated_by_the_genes_own_symbol(tmp_path):
    def row(ids, match, fam):
        return {"human_relationship": "1:1", "zebrafish_relationship": "none", "family": fam,
                "_human_ids": set(ids), "_zebrafish_ids": set(), "own_symbol": "nanos3",
                "human_symbol_match": match, "human_orthologs_symbols": ",".join(ids),
                "human_orthologs_ids": ";".join(ids)}
    by = {"TreeFam": {"g1": row(["A"], "exact", "TF1"), "g2": row(["A"], "none", "TF2"),
                      "g3": row(["A"], "none", "TF3")},
          "Compara": {"g1": row(["B"], "none", "C1"), "g2": row(["B"], "exact", "C2"),
                      "g3": row(["B"], "none", "C3")}}
    out = _identity(by, ["g1", "g2", "g3"], tmp_path)
    adj = out["human"]["adjudication_by_own_ncbi_symbol"]
    assert adj == {"TreeFam_matches_own_symbol": 1, "Compara_matches_own_symbol": 1,
                   "neither_matches_own_symbol": 1}
    rows = list(csv.DictReader(open(tmp_path / "pilot_identity_differences.tsv"), delimiter="\t"))
    assert {r["gene_id"]: r["own_symbol"] for r in rows} == {"g1": "nanos3", "g2": "nanos3", "g3": "nanos3"}


def test_default_output_goes_to_share(tmp_path):
    p = default_out("pilot", base=str(tmp_path / "share"), now="20261023-120000")
    from cichlid_homology import __version__
    assert p.parent == tmp_path / "share" and p.name == f"pilot_{__version__}_20261023-120000"
    assert default_out("survey", ".json", now="x").name == f"survey_{__version__}_x.json"
    assert str(default_out("pilot")).startswith(str(__import__("pathlib").Path("~/share").expanduser()))


class FakeLS:
    def species(self, leaf):
        return PREFIX2SP.get(leaf.split("_")[0])

    def is_project_gene(self, leaf):
        return leaf.split("_")[0] in PREFIXES and "_" in leaf


def fam(newick):
    return SimpleNamespace(tree=newick, leaves=[x.name for x in leaves(parse_newick(newick))],
                           rooted=False, tree_mode="denovo", family_id="F1")


def test_expanded_family_replaces_ambiguous_labels_only():
    mz = ",".join(f"Mzeb_{i}" for i in range(1, 11))
    copies = lambda s: ",".join(f"{s}_{i}" for i in (1, 2, 3))
    newick = (f"((Hsap_1,Mmus_1),(Locu_1,({copies('Aang')},({copies('Drer')},"
              f"({copies('Gacu')},({copies('Olat')},({mz},Acal_1)))))));")
    rows, _ = process_family(fam(newick), FakeLS(), set(), st(), FakeClades(), Params(),
                             "Mzeb", PREFIXES, [], None)
    assert len(rows) == 10 and all(r["expanded_family"] for r in rows)
    assert {r["b_label_raw"] for r in rows} == {"multi_copy"}
    assert {r["b_label"] for r in rows} == {"expanded_family"}


def test_small_family_keeps_its_label():
    newick = "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))))));"
    rows, _ = process_family(fam(newick), FakeLS(), set(), st(), FakeClades(), Params(),
                             "Mzeb", PREFIXES, [], None)
    assert rows[0]["expanded_family"] is False and rows[0]["b_label"] == "single"


def test_goal_b_summary_counts_carriers_and_histograms():
    base = {"family_mzeb_genes": 2, "human_relationship": "1:many", "b_anchor": "human", "b_pattern": "...+",
            "b_cichlid_n": 10, "expanded_family": False}
    rows = [dict(base, b_label="cichlid_specific_duplication", b_label_raw="cichlid_specific_duplication",
                 b_other_share2p=0.0, b_cichlid_share2p=0.4, b_cichlid_with2p=4),
            dict(base, b_label="retained_pair", b_label_raw="retained_pair", b_other_share2p=0.9,
                 b_cichlid_share2p=1.0, b_cichlid_with2p=10),
            dict(base, b_label="no_orthologous_anchor", b_label_raw="no_orthologous_anchor",
                 b_other_share2p="", b_cichlid_share2p="", b_cichlid_with2p="")]
    g = _goal_b(rows)
    assert g["b2_cichlid_carriers_for_cichlid_specific_duplications"] == {"4 of 10 cichlid species": 1}
    assert sum(g["histogram_share_with_2plus_copies_non_cichlid"]) == 2     # the unanchored gene is left out
    assert "expanded_family" in NO_B
