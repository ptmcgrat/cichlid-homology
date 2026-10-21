import json
from cichlid_homology.read import leaf_names, root_children, load_family, is_focal

PREFIXES = {"Mzeb", "Acal", "Drer"}


def make(tmp_path, tree, aln_names, done):
    d = tmp_path / "TF000001"
    d.mkdir()
    (d / "tree.nwk").write_text(tree)
    (d / "combined.aln").write_text("".join(f">{n}\nAC-T\n" for n in aln_names))
    (d / "done.json").write_text(json.dumps(done))
    return d


def test_leaf_names_ignores_internal_labels_and_lengths():
    t = "((a:0.1,b_2:0.2)0.9:0.3,(Mzeb_1.2:1,ENSP0001:2):0.1,c);"
    assert leaf_names(t) == ["a", "b_2", "Mzeb_1.2", "ENSP0001", "c"]


def test_root_children():
    assert root_children("((a,b),c);") == 2
    assert root_children("(a,b,c);") == 3
    assert root_children("a;") is None


def test_is_focal_uses_text_before_first_underscore():
    assert is_focal("Mzeb_101467881", PREFIXES)
    assert not is_focal("Mzeb", PREFIXES)
    assert not is_focal("ENSDARP0001", PREFIXES)
    assert not is_focal("Drerx_1", PREFIXES)


def test_load_unrooted_place_family(tmp_path):
    names = ["Mzeb_1", "Acal_2", "ENSP1", "ENSP2"]
    d = make(tmp_path, "((Mzeb_1,Acal_2),ENSP1,ENSP2);", names,
             {"release": "v1", "rooted": False, "tree_mode": "place", "inputs": "x"})
    f = load_family(d, PREFIXES)
    assert f.tree_mode == "place" and not f.rooted and f.root_children == 3
    assert f.focal == ["Mzeb_1", "Acal_2"] and f.reference == ["ENSP1", "ENSP2"]
    assert f.problems == []


def test_rooted_true_needs_text_and_binary_root(tmp_path):
    names = ["Mzeb_1", "Acal_2", "ENSP1"]
    d = make(tmp_path, "((Mzeb_1,Acal_2),ENSP1);", names,
             {"rooted": True, "rooting": "orthofinder species tree", "tree_mode": "denovo"})
    f = load_family(d, PREFIXES)
    assert f.rooted and f.rooting == "orthofinder species tree"


def test_rooted_claim_not_usable_is_treated_unrooted(tmp_path):
    d = make(tmp_path, "(Mzeb_1,Acal_2,ENSP1);", ["Mzeb_1", "Acal_2", "ENSP1"],
             {"rooted": True, "rooting": "x", "tree_mode": "denovo"})
    f = load_family(d, PREFIXES)
    assert not f.rooted and "rooted_claim_not_usable" in f.problems


def test_alignment_leaf_mismatch_reported(tmp_path):
    d = make(tmp_path, "((Mzeb_1,Acal_2),ENSP1);", ["Mzeb_1", "ENSP1", "extra"],
             {"rooted": False, "tree_mode": "place"})
    f = load_family(d, PREFIXES)
    assert "aln_row_without_leaf" in f.problems and "leaf_without_aln_row" in f.problems
