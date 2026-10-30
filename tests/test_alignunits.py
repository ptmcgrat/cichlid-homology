import csv
import gzip
import json
from types import SimpleNamespace

from cichlid_homology import alignunits as U
from cichlid_homology.tree import parse_newick

SP = {"Mzeb_1": "maylandia_zebra", "Mzeb_2": "maylandia_zebra", "Mzeb_3": "maylandia_zebra", "Acal_1": "astatotilapia_calliptera", "H1": "homo_sapiens", "H2": "homo_sapiens", "M1": "mus_musculus",
      "Z1": "danio_rerio", "Z2": "danio_rerio", "G1": "lepisosteus_oculatus", "S1": "gasterosteus_aculeatus", "S2": "gasterosteus_aculeatus", "C1": "latimeria_chalumnae"}
is_human = lambda s: s == "homo_sapiens"
is_tet = lambda s: s in ("homo_sapiens", "mus_musculus")
is_project = lambda n: n.startswith(("Mzeb_", "Acal_"))
LEVEL = {"maylandia_zebra": "cichlid", "astatotilapia_calliptera": "cichlid", "danio_rerio": "teleost", "homo_sapiens": "outside_vertebrates"}


def fu(newick, genes=(), pairs=(), shadowed=("gasterosteus_aculeatus",)):
    root = parse_newick(newick)
    return U.family_units("F1", root, "as given", "place", SP.get, is_human, is_tet, is_project, set(shadowed), LEVEL.get, "Mzeb", list(genes), list(pairs), ("F1/combined.aln", "F1/tree.nwk"))


def ids(rows):
    return {r["unit_id"]: sorted(m["gene_id"] for m in rows[1] if m["unit_id"] == r["unit_id"]) for r in rows[0]}


def test_the_unit_is_the_smallest_clade_with_a_human_gene():
    t = "(((Mzeb_1:1,Z1:1):1,(H1:1,Mzeb_2:1):1):1,(S1:1,G1:1):1);"
    u, m, g, p = fu(t, ["Mzeb_1", "Mzeb_2"])
    members = ids((u, m))
    assert g["Mzeb_2"][0][:2] == ("F1:H1", "human") and members["F1:H1"] == ["H1", "Mzeb_2"]                    # its sister is the human gene
    assert g["Mzeb_1"][0][:2] == ("F1:H2", "human") and members["F1:H2"] == ["H1", "Mzeb_1", "Mzeb_2", "Z1"]      # the smallest clade that has one climbs one node
    assert [r["unit_id"] for r in u if r["unit_id"] == "F1:H1"][0] == "F1:H1" and {r["unit_id"]: r["n_genes"] for r in u} == {"F1:H1": 2, "F1:H2": 4}


def test_without_a_human_gene_the_unit_has_a_tetrapod_and_without_either_it_is_the_whole_family():
    u, m, g, p = fu("(((Mzeb_1:1,M1:1):1,Z1:1):1,(S1:1,G1:1):1);", ["Mzeb_1"])
    assert g["Mzeb_1"][0][:2] == ("F1:T1", "tetrapod") and ids((u, m))["F1:T1"] == ["M1", "Mzeb_1"]
    u, m, g, p = fu("(((Mzeb_1:1,Z1:1):1,Z2:1):1,(S1:1,G1:1):1);", ["Mzeb_1"])
    assert g["Mzeb_1"][0] == ("F1:F1", "family", 5, True) and u[0]["whole_family"] == "yes" and u[0]["n_genes"] == u[0]["n_family_genes"] == 5


def test_a_pair_takes_the_smallest_clade_that_holds_both_genes_and_a_human_gene():
    t = "(((Mzeb_1:1,Z1:1):1,(H1:1,Mzeb_2:1):1):1,((Mzeb_3:1,H2:1):1,S1:1):1);"
    u, m, g, p = fu(t, ["Mzeb_1", "Mzeb_2", "Mzeb_3"], [(0, "Mzeb_1", "Mzeb_2"), (1, "Mzeb_1", "Mzeb_3"), (2, "Mzeb_2", "Mzeb_3")])
    members = ids((u, m))
    assert sorted(members) == ["F1:H1", "F1:H2", "F1:H3", "F1:H4"] and [len(members[k]) for k in sorted(members)] == [2, 2, 4, 7]               # numbered by size: the two of two, the one of four, the whole tree
    assert p[0][0][:2] == ("F1:H3", "human") and g["Mzeb_1"][0][0] == "F1:H3"                                    # both genes sit under the clade that is already Mzeb_1's own unit
    assert p[1][0] == ("F1:H4", "human", 7, True) and p[2][0] == ("F1:H4", "human", 7, True)                      # genes on the two sides of the root: the whole tree, which holds a human gene


def test_units_are_shared_numbered_by_size_and_count_what_asked_for_them():
    t = "(((Mzeb_1:1,Z1:1):1,(H1:1,Mzeb_2:1):1):1,(S1:1,G1:1):1);"
    u, m, g, p = fu(t, ["Mzeb_1", "Z1"], [(0, "Mzeb_1", "Z1")])
    assert len(u) == 1 and u[0]["unit_id"] == "F1:H1" and u[0]["genes_requested"] == 2 and u[0]["pairs_requested"] == 1                 # one unit for two genes and their pair
    assert g["Mzeb_1"][0][0] == g["Z1"][0][0] == p[0][0][0]


def test_members_say_what_each_leaf_is():
    t = "(((Mzeb_1:1,Acal_1:1):1,(H1:1,S1:1):1):1,(Z1:1,G1:1):1);"
    u, m, g, p = fu(t, ["Mzeb_1"])
    by = {r["gene_id"]: r for r in m}
    assert by["Mzeb_1"]["kind"] == "focal" and by["Acal_1"]["kind"] == "project" and by["H1"]["kind"] == "reference"
    assert by["H1"]["is_human"] == "yes" and by["H1"]["is_tetrapod"] == "yes" and by["Mzeb_1"]["is_human"] == "no"
    assert by["S1"]["shadowed_reference"] == "yes" and by["H1"]["shadowed_reference"] == "no" and by["Mzeb_1"]["ladder_level"] == "cichlid"
    r = u[0]
    assert (r["n_genes"], r["n_focal_genes"], r["n_project_genes"], r["n_reference_genes"], r["n_shadowed_reference"], r["n_human_genes"], r["human_genes"], r["focal_genes"]) == (4, 1, 2, 2, 1, 1, "H1", "Mzeb_1")
    assert r["alignment_file"] == "F1/combined.aln" and r["share_of_family"] == round(4 / 6, 4)


def test_a_gene_that_is_not_in_the_tree_has_no_unit():
    u, m, g, p = fu("((Mzeb_1:1,H1:1):1,Z1:1);", ["Mzeb_1", "Mzeb_9"], [(0, "Mzeb_1", "Mzeb_9")])
    assert g["Mzeb_9"] == (None, "gene not in the family tree") and p[0] == (None, "gene not in the family tree") and g["Mzeb_1"][1] == "ok"


def _family_fixture(tmp_path):
    froot = tmp_path / "families"
    for fam, tree in (("FAM1", "(((Mzeb_1:1,Z1:1):1,(H1:1,Mzeb_2:1):1):1,(S1:1,G1:1):1);"), ("FAM2", "((Mzeb_3:1,Z2:1):1,(S2:1,G1:1):1);")):
        (froot / fam).mkdir(parents=True)
        (froot / fam / "done.json").write_text("{}")
    trees = {"FAM1": "(((Mzeb_1:1,Z1:1):1,(H1:1,Mzeb_2:1):1):1,(S1:1,G1:1):1);", "FAM2": "((Mzeb_3:1,Z2:1):1,(S2:1,G1:1):1);"}

    def load_family(folder, prefixes, fam):
        t = trees[fam]
        return SimpleNamespace(tree=t, leaves=[x.name for x in parse_newick(t).children and __import__("cichlid_homology.tree", fromlist=["leaves"]).leaves(parse_newick(t))], tree_mode="place")
    ls = SimpleNamespace(species=SP.get, is_project_gene=is_project, shadowed=lambda: ["gasterosteus_aculeatus"])
    return froot, load_family, ls


def test_a_whole_run_writes_the_unit_gene_and_pair_tables(tmp_path):
    froot, load_family, ls = _family_fixture(tmp_path)
    rows = {"Mzeb_1": {"family": "FAM1", "clade_id": "FAM1:G1"}, "Mzeb_2": {"family": "FAM1", "clade_id": "FAM1:G1"}, "Mzeb_3": {"family": "FAM2", "clade_id": "FAM2:G1"}, "Mzeb_4": {"family": "NoHit", "clade_id": ""}}
    pairs = [{"gene_a": "Mzeb_1", "gene_b": "Mzeb_2", "group": "ohnolog pairs"}, {"gene_a": "Mzeb_1", "gene_b": "Mzeb_3", "group": "later duplicates"}, {"gene_a": "Mzeb_1", "gene_b": "Mzeb_4", "group": "x"}]
    s = U.run_units(rows=rows, pairs=pairs, froot=froot, ref_files={}, prefixes={"Mzeb"}, focal_prefix="Mzeb", load_family=load_family, load_ls=lambda names: ls,
                    rooter=lambda fam, f, ref, ls_: (parse_newick(f.tree), "as given"), clades_info=(is_human, is_tet, LEVEL.get), out_dir=tmp_path / "out", progress=lambda m: None)
    rd = lambda n: list(csv.DictReader(open(tmp_path / "out" / f"{n}.tsv"), delimiter="\t"))
    g = {r["gene_id"]: r for r in rd("gene_to_unit")}
    assert (g["Mzeb_1"]["unit_id"], g["Mzeb_1"]["basis"], g["Mzeb_2"]["unit_id"]) == ("FAM1:H2", "human", "FAM1:H1") and g["Mzeb_4"]["status"] == "no family"
    assert (g["Mzeb_3"]["basis"], g["Mzeb_3"]["whole_family"], g["Mzeb_3"]["unit_genes"]) == ("family", "yes", "4")
    p = rd("pair_to_unit")
    assert (p[0]["unit_id"], p[0]["status"], p[0]["unit_a"], p[0]["unit_b"]) == ("FAM1:H2", "ok", "FAM1:H2", "FAM1:H1")                      # the pair's unit is the bigger of the two
    assert p[1]["status"] == "the genes are in different families" and p[2]["status"] == "a gene has no family"
    units = {r["unit_id"]: r for r in rd("alignment_units")}
    assert set(units) == {"FAM1:H1", "FAM1:H2", "FAM2:F1"} and units["FAM1:H2"]["pairs_requested"] == "1" and units["FAM2:F1"]["tree_file"].endswith("FAM2/tree.nwk")
    with gzip.open(tmp_path / "out" / "alignment_unit_members.tsv.gz", "rt") as fh:
        mem = list(csv.DictReader(fh, delimiter="\t"))
    assert len(mem) == 2 + 4 + 4 and {r["gene_id"] for r in mem if r["unit_id"] == "FAM1:H1"} == {"H1", "Mzeb_2"}
    assert s["units"] == 3 and s["genes_by_basis"] == {"human": 2, "family": 1, "no family": 1} and s["genes_whose_unit_is_the_whole_family"] == 1
    assert json.loads((tmp_path / "out" / "alignment_units.json").read_text())["unit_size_by_basis"]["human"]["units"] == 2
    assert "smallest clade that holds a human gene" in (tmp_path / "out" / "alignment_units.md").read_text()


def test_a_family_that_cannot_be_read_is_counted_and_its_genes_have_a_status(tmp_path):
    froot, load_family, ls = _family_fixture(tmp_path)
    rows = {"Mzeb_1": {"family": "FAM1", "clade_id": ""}, "Mzeb_7": {"family": "FAM9", "clade_id": ""}}
    s = U.run_units(rows=rows, pairs=[], froot=froot, ref_files={}, prefixes={"Mzeb"}, focal_prefix="Mzeb", load_family=load_family, load_ls=lambda names: ls,
                    rooter=lambda fam, f, ref, ls_: (parse_newick(f.tree), "as given"), clades_info=(is_human, is_tet, LEVEL.get), out_dir=tmp_path / "out2", progress=lambda m: None)
    g = {r["gene_id"]: r for r in csv.DictReader(open(tmp_path / "out2" / "gene_to_unit.tsv"), delimiter="\t")}
    assert g["Mzeb_7"]["status"] == "family not read" and g["Mzeb_1"]["status"] == "ok" and s["skipped"] == {"no finished family folder": 1}
