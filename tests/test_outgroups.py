import csv
import json

from cichlid_homology import outgroups as O
from cichlid_homology.figpilot import pair_genes

GAR, REED, SHARK, COEL = "lepisosteus_oculatus", "erpetoichthys_calabaricus", "callorhinchus_milii", "latimeria_chalumnae"


def rec(clade, inside=(), outside=(), stems=(0.1, 0.3), pair=True, published=False, conf="high", kind="teleost", mode="split", tet_outside=1, tet_inside=0):
    ins = [{"gene": f"i{n}", "species": s, "lineage": l, "attach_support": 0} for n, (s, l) in enumerate(inside)]
    out = [{"gene": f"o{n}", "species": s} for n, s in enumerate(outside)]
    out += [{"gene": f"t{n}", "species": "homo_sapiens"} for n in range(tet_outside)]
    ins += [{"gene": f"ti{n}", "species": "mus_musculus", "lineage": "A", "attach_support": 5} for n in range(tet_inside)]
    calls = {"ohnolog_status": "paired" if pair else "single copy (duplicate lost)", "ohnolog_confidence": conf, "published_pair": "the two zebrafish names are a published pair" if published else ""}
    return {"clade": clade, "family": clade.split(":")[0], "mode": mode, "split": {"kind": kind, "level": "teleost", "overlap": 0.9}, "lineages": [{"id": "A", "stem": stems[0]}, {"id": "B", "stem": stems[1]}],
            "outgroup": {"inside": ins, "outside": out}, "strip": {"species": ["homo_sapiens", "mus_musculus"], "groups": ["tetrapod", "tetrapod"]}, "triage": {"stem_ratio": round(stems[1] / stems[0], 2), "deep_order_score": 1.0},
            "genes": [{"gene": "Mzeb_1", "calls": calls}, {"gene": "Mzeb_2", "calls": calls}]}


def test_each_outgroup_species_falls_in_one_of_five_cases():
    r = rec("F:G1", inside=[(GAR, "B"), (REED, "A"), (REED, "B"), (COEL, "A")], outside=[SHARK, COEL, REED][:1] + [COEL])
    assert O.species_case(r, GAR) == ("inside one lineage", ["B"])
    assert O.species_case(r, REED) == ("in both lineages", ["A", "B"])
    assert O.species_case(r, SHARK) == ("outside both lineages", None)
    assert O.species_case(r, COEL) == ("mixed (some outside, some inside)", ["A"])
    assert O.species_case(r, "homo_sapiens") == ("outside both lineages", None) and O.species_case(r, "bos_taurus") == ("absent", None)


def test_the_summary_counts_clades_and_splits_the_inside_cases_into_slower_and_faster():
    recs = [rec("F:G1", inside=[(GAR, "A")] * 30, stems=(0.1, 0.3)), rec("F2:G1", inside=[(GAR, "B")], stems=(0.1, 0.3)), rec("F3:G1", outside=[GAR, REED]), rec("F4:G1", outside=[REED]),
            rec("F5:G1", pair=False, inside=[(GAR, "A")]), rec("F6:G1", kind="later"), rec("F7:G1", mode="single")]
    s = O.summarize(recs)
    assert s["clades"] == 4 and s["pair_genes"] == 8                                                        # the clade with 30 inside genes counts once; no pair, a later split, no split: left out
    g = s["species"]["gar"]
    assert (g["clades_with_the_species"], g["outside both lineages"], g["inside one lineage"], g["absent"]) == (3, 1, 2, 1)
    assert g["inside_one_lineage_which"] == {"slower lineage": 1, "faster lineage": 1} and g["inside_one_lineage_by_letter"] == {"A": 1, "B": 1}
    assert s["species"]["reedfish"]["outside both lineages"] == 2 and s["clades_by_number_of_outgroup_species_inside"] == {0: 2, 1: 2}
    assert O.summarize(recs, pairs_only=False)["clades"] == 5


def test_tetrapods_are_one_row_and_the_quality_of_the_pair_is_compared():
    recs = [rec("F:G1", tet_inside=2, published=False, conf="low"), rec("F2:G1", published=True, conf="high"), rec("F3:G1", published=True, conf="high", inside=[(GAR, "A")])]
    s = O.summarize(recs)
    assert s["species"]["tetrapods"]["mixed (some outside, some inside)"] == 1 and s["species"]["tetrapods"]["outside both lineages"] == 2         # one gene outside, mouse genes inside
    assert s["by_pair_quality"]["published"] == {"clades": 2, "any outgroup species inside a lineage": 1, "tetrapods inside a lineage": 0}
    assert s["by_pair_quality"]["not published"] == {"clades": 1, "any outgroup species inside a lineage": 0, "tetrapods inside a lineage": 1}      # the tetrapods are counted apart from the four fish
    assert s["by_pair_quality"]["confidence low"]["clades"] == 1


def test_two_species_inside_one_lineage_are_compared_for_the_same_lineage():
    recs = [rec("F:G1", inside=[(GAR, "A"), (REED, "A")]), rec("F2:G1", inside=[(GAR, "A"), (REED, "B")])]
    assert O.summarize(recs)["clades_with_two_or_more_species_inside_one_lineage"] == {"all_in_the_same_lineage": 1, "in_different_lineages": 1}


def test_the_files_are_written_and_the_row_per_clade_says_where_each_species_is(tmp_path):
    recs = [rec("F:G1", inside=[(GAR, "A")], published=True)]
    O.write_outputs(tmp_path, recs)
    assert "Where the outgroups fall" in (tmp_path / "outgroup_summary.md").read_text() and json.loads((tmp_path / "outgroup_summary.json").read_text())["clades"] == 1
    rows = list(csv.DictReader(open(tmp_path / "outgroup_by_clade.tsv"), delimiter="\t"))
    assert rows[0]["gar"] == "inside one lineage (A)" and rows[0]["reedfish"] == "absent" and rows[0]["published"] == "True" and rows[0]["outgroup_species_inside"] == "1"


def test_every_gene_in_a_pair_is_selected_for_the_run_over_all_pairs():
    rows = {"a": {"ohnolog_status": "paired", "family": "F1"}, "b": {"ohnolog_status": "single copy (duplicate lost)", "family": "F1"}, "c": {"ohnolog_status": "possible pair (the other teleosts lost a copy)", "family": "F2"},
            "d": {"ohnolog_status": "paired", "family": "NoHit"}}
    assert [g["gene"] for g in pair_genes(rows)] == ["a", "c"]
