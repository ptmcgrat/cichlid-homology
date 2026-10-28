import json

from helpers import PREFIX2SP, SPECIES_TREE, FakeClades, st
from cichlid_homology.cladeprofile import build_files, group_label, species_tree_json, write_files
from cichlid_homology.cladesplit import run_split, split_family
from cichlid_homology.parsimony import ParsimonyContext

NAMES = {"human": "homo_sapiens", "mouse": "mmus", "zebrafish": "danio_rerio", "gar_reed": {"locu", "ecal"}, "coelacanth": "ecal"}


class Clades(FakeClades):
    def lineage(self, s):
        return "sarcopterygian" if s in ("homo_sapiens", "mmus") else super().lineage(s)


class LS:
    def species(self, leaf):
        return PREFIX2SP.get(leaf.split("_")[0])


def profiles_of(newick):
    sp, cl, out = st(), Clades(), []
    split_family("F", newick, LS(), sp, cl, ParsimonyContext(sp, cl), NAMES, 0.5, True, out)
    return out


def test_a_pair_is_counted_outside_in_lineage_a_and_in_lineage_b():
    (p,) = profiles_of("((Hsap_1,Mmus_1),((Aang_1,Drer_1),(Aang_2,Drer_2)));")
    assert p["clade"] == "F:G1" and p["n_lineages"] == 2
    assert p["counts"] == {"homo_sapiens": [1, 0, 0], "mmus": [1, 0, 0], "aang": [0, 1, 1], "danio_rerio": [0, 1, 1]}
    assert (p["zebrafish_a"], p["zebrafish_b"], p["later"], p["later_omitted"]) == (["Drer_1"], ["Drer_2"], [], 0)


def test_a_species_that_lost_one_ohnolog_has_a_zero_in_that_lineage():
    (p,) = profiles_of("((Hsap_1,Mmus_1),((Aang_1,Drer_1,Trub_1),(Drer_2,Trub_2)));")
    assert p["counts"]["aang"] == [0, 1, 0] and p["counts"]["danio_rerio"] == [0, 1, 1] and p["counts"]["trub"] == [0, 1, 1]
    assert profiles_of("((Hsap_1,Mmus_1),((Aang_1,Drer_1),(Drer_2,Trub_1)));") == []       # the lineages share too few species to be a supported duplication


def test_a_later_duplicate_inside_one_ohnolog_is_listed_with_its_lineage_level_and_species():
    (p,) = profiles_of("((Hsap_1,Mmus_1),((Aang_1,Drer_1,((Trub_1,Gacu_1),(Trub_2,Gacu_2))),(Aang_2,Drer_2)));")
    assert p["counts"]["trub"] == [0, 2, 0] and p["counts"]["gacu"] == [0, 2, 0] and p["counts"]["aang"] == [0, 1, 1]
    assert p["later"] == [{"lineage": "A", "level": "percomorph", "species": ["gacu", "trub"]}]
    assert p["zebrafish_a"] == ["Drer_1"] and p["zebrafish_b"] == ["Drer_2"]


def test_a_clade_without_a_teleost_duplication_has_no_profile():
    assert profiles_of("((Hsap_1,Mmus_1),(Aang_1,Drer_1));") == []
    assert profiles_of("((Hsap_1,Mmus_1),((Trub_1,Gacu_1),(Trub_2,Gacu_2)));") == []        # a duplication of percomorphs only is not the teleost duplication


def test_the_files_have_the_species_in_tree_order_with_groups_and_a_pruned_tree():
    (p,) = profiles_of("((Hsap_1,Mmus_1),((Aang_1,Drer_1),(Aang_2,Drer_2)));")
    species, tree, rows = build_files([p], Clades(), SPECIES_TREE)
    assert [s["name"] for s in species] == ["homo_sapiens", "mmus", "aang", "danio_rerio"]
    assert [s["group"] for s in species] == ["tetrapod", "tetrapod", "early teleost", "early teleost"] and species[0]["label"] == "Homo sapiens"
    assert rows == [{"clade": "F:G1", "family": "F", "group": "G1", "counts": "1.0.0,1.0.0,0.1.1,0.1.1", "n_lineages": 2, "later": "[]", "later_omitted": 0,
                     "zebrafish_a": "Drer_1", "zebrafish_b": "Drer_2"}]
    assert [c.get("name", "") for c in tree["children"]] == ["", ""] and tree["children"][0]["children"][0]["name"] == "homo_sapiens"


def test_the_pruned_species_tree_collapses_unary_nodes_and_drops_numeric_labels():
    t = species_tree_json("((a:1,b:1)0.9:1,(c:1,(d:1,e:1):1)Teleostei:1);", ["a", "b", "d"])
    assert t == {"name": "", "children": [{"name": "", "children": [{"name": "a"}, {"name": "b"}]}, {"name": "d"}]}
    assert group_label("danio_rerio", Clades()) == "early teleost" and group_label("trub", Clades()) == "percomorph" and group_label("locu", Clades()) == "gar or reedfish"


def test_the_run_collects_profiles_and_the_files_are_written(tmp_path):
    prep = tmp_path / "prepared"
    (prep / "trees").mkdir(parents=True)
    (prep / "trees" / "P1.nhx").write_text("((Hsap_1,Mmus_1),((Aang_1,Drer_1),(Aang_2,Drer_2)));")
    (prep / "trees" / "P2.nhx").write_text("((Hsap_1,Mmus_1),(Aang_1,Drer_1));")
    sp, cl, profiles = st(), Clades(), []
    run_split(prep, LS(), sp, cl, ParsimonyContext(sp, cl), NAMES, profiles=profiles)
    assert [p["clade"] for p in profiles] == ["P1:G1"]
    out = tmp_path / "out"
    out.mkdir()
    assert write_files(out, "Compara", profiles, cl, SPECIES_TREE) == {"species": 4, "clades": 1}
    assert json.loads((out / "clade_species.Compara.clades.json").read_text())[0]["counts"] == "1.0.0,1.0.0,0.1.1,0.1.1"
    assert json.loads((out / "clade_species.Compara.species.json").read_text())[3]["name"] == "danio_rerio"
    assert "children" in json.loads((out / "clade_species.Compara.tree.json").read_text())
