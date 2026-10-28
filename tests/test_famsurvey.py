import json

from helpers import PREFIX2SP, PREFIXES, FakeClades
from cichlid_homology.famsurvey import COLS, finished_families, report_text, run_survey, summarize, survey_family


NAMES = {"human": "homo_sapiens", "mouse": "mmus", "zebrafish": "danio_rerio", "gar_reed": {"locu", "ecal"}}     # the helper's short names


class LS:
    def species(self, leaf):
        return PREFIX2SP.get(leaf.split("_")[0])


def make(tmp_path):
    fams = {
        "F1": "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Mzeb_1,Acal_1)))));",               # human, mouse, gar, teleosts, ours
        "F2": "(((Hsap_1,Hsap_2),Hsap_3),(Mmus_1,Ecal_1));",                                  # tetrapod only: three human genes
        "F3": "(Aang_1,(Drer_1,(Drer_2,Trub_1)));",                                           # fish only
        "F4": "((Hsap_1,Hsap_2),(Mmus_1,Mzeb_1));",                                           # a placed M. zebra gene in a tetrapod tree
        "F5": "(Bflo_1,Cmil_1);",                                                             # neither human nor teleost
    }
    for f, nwk in fams.items():
        d = tmp_path / f
        d.mkdir()
        (d / "tree.nwk").write_text(nwk)
        (d / "done.json").write_text(json.dumps({"tree_mode": "place", "rooted": False}))
    (tmp_path / "unfinished").mkdir()
    (tmp_path / "unfinished" / "tree.nwk").write_text("(a,b);")
    return fams


def test_only_finished_families_are_listed(tmp_path):
    make(tmp_path)
    assert finished_families(tmp_path) == ["F1", "F2", "F3", "F4", "F5"]


def test_one_family_is_counted_by_species_group(tmp_path):
    make(tmp_path)
    r = survey_family(tmp_path / "F1", LS(), FakeClades(), PREFIXES, "Mzeb", NAMES)
    assert set(r) == set(COLS)
    assert (r["n_human"], r["n_mouse"], r["n_gar_reedfish"], r["n_zebrafish"]) == (1, 1, 1, 1)
    assert (r["n_teleost_ref_leaves"], r["n_teleost_ref_species"], r["n_nonteleost_gnathostome_leaves"]) == (2, 2, 3)
    assert (r["n_project_leaves"], r["has_focal_gene"], r["tree_mode"]) == (2, 1, "place")
    g = survey_family(tmp_path / "F4", LS(), FakeClades(), PREFIXES, "Mzeb", NAMES)
    assert (g["n_human"], g["n_teleost_ref_leaves"], g["n_project_leaves"], g["has_focal_gene"]) == (2, 0, 1, 1)   # a placed gene is not a reference teleost


def test_the_summary_counts_human_genes_per_family_and_families_without_teleosts(tmp_path):
    make(tmp_path)
    rows = run_survey(tmp_path, finished_families(tmp_path), LS(), FakeClades(), PREFIXES, "Mzeb", names=NAMES)
    s = summarize(rows)
    assert s["families"] == 5 and s["with_human"] == 3 and s["exactly_one_human_gene"] == 1
    assert s["human_genes_per_family_histogram"] == {"0": 2, "1": 1, "2": 1, "3": 1, "4": 0, "5-9": 0, "10-19": 0, "20+": 0}
    assert (s["with_human_and_no_teleost"], s["with_human_and_teleost"]) == (2, 1)           # F2 and F4 have no reference teleost
    assert s["human_genes_by_teleost_presence"]["3-4"] == {"no_teleost": 1, "with_teleost": 0}
    assert (s["with_teleost"], s["with_teleost_and_no_human"], s["no_human_and_no_teleost"]) == (2, 1, 1)
    assert s["human_no_teleost_families"]["with_our_placed_genes"] == 1 and s["human_no_teleost_families"]["with_a_focal_gene"] == 1
    assert s["families_with_a_focal_gene"]["families"] == 2 and s["families_with_a_focal_gene"]["with_human"] == 2
    assert "Families with human genes and no teleosts" in report_text(s, "Compara")


def test_an_unreadable_family_is_reported_not_fatal(tmp_path):
    make(tmp_path)
    (tmp_path / "F3" / "tree.nwk").unlink()
    rows = run_survey(tmp_path, ["F1", "F3"], LS(), FakeClades(), PREFIXES, "Mzeb", names=NAMES)
    assert rows[1]["tree_mode"].startswith("unreadable") and rows[0]["n_human"] == 1


def make_reference(tmp_path):
    from cichlid_homology.famsurvey import run_reference_survey  # noqa: F401
    prep = tmp_path / "prepared"
    (prep / "trees").mkdir(parents=True)
    trees = {
        "C1": "((Hsap_1:0.1[&&NHX:S=Homo_sapiens],Mmus_1:0.2),(Locu_1:0.3,(Aang_1:0.1,Drer_1:0.1)));",    # human, mouse, gar, teleosts
        "C2": "((Hsap_1,Hsap_2),(Mmus_1,Ecal_1));",                                                      # human, no teleost
        "C3": "(Aang_1,(Drer_1,Trub_1));",                                                               # teleost, no human
        "C4": "((Hsap_1,Mmus_1),(Aang_1,Drer_1));",                                                      # human and teleosts, not hit by us
        "C5": "(Bflo_1,Cmil_1);",
    }
    for f, nwk in trees.items():
        (prep / "trees" / f"{f}.nwk").write_text(nwk)
    (prep / "family_audit.tsv").write_text("family\tn_members\tstatus\nC1\t7\tok\nC2\t4\tok\nC3\t3\ttiny_tree\nC4\t4\tok\nC5\t2\tno_panel_members\nC9\t0\tno_panel_members\n")
    return prep


def test_the_reference_survey_counts_human_and_teleost_genes_and_what_the_project_built(tmp_path):
    from cichlid_homology.famsurvey import REF_COLS, report_reference_text, run_reference_survey, summarize_reference
    prep = make_reference(tmp_path)
    rows, diag = run_reference_survey(prep, LS(), FakeClades(), built={"C1"}, hit={"C1", "C3"}, focal_hit={"C1"}, names=NAMES)
    assert [r["family"] for r in rows] == ["C1", "C2", "C3", "C4", "C5"] and set(rows[0]) == set(REF_COLS)
    c1 = rows[0]
    assert (c1["n_human"], c1["n_mouse"], c1["n_gar_reedfish"], c1["n_teleost_ref_leaves"], c1["n_leaves"]) == (1, 1, 1, 2, 5)   # the NHX comment is ignored
    assert (c1["audit_status"], c1["built"], c1["hit_by_project"], c1["has_focal_gene"]) == ("ok", 1, 1, 1)
    assert rows[2]["audit_status"] == "tiny_tree" and rows[3]["built"] == 0 and rows[3]["hit_by_project"] == 0
    assert diag["families_with_files"] == 5 and diag["files_by_extension"] == {".nwk": 5} and diag["families_with_a_tree_file"] == 5
    assert diag["audited_families_without_a_file"] == 1 and rows[0]["leaf_source"] == "tree" and diag["families_with_no_leaves"] == 0
    assert diag["audit_statuses"] == {"ok": 3, "tiny_tree": 1, "no_panel_members": 2}
    s = summarize_reference(rows, diag)
    assert s["with_human"] == 3 and s["with_human_and_no_teleost"] == 1 and s["with_human_and_teleost"] == 2
    ref = s["reference"]
    assert ref["human_and_teleost"] == {"families": 2, "built": 1, "hit_by_a_project_species": 1, "hit_by_the_focal_species": 1, "not_hit_by_any_project_species": 1}
    assert ref["human_and_no_teleost"]["families"] == 1 and ref["human_and_no_teleost"]["audit_statuses"] == {"ok": 1}
    assert ref["teleost_and_no_human"] == {"families": 1, "built": 0}
    assert "Reference-tree survey" in report_reference_text(s, "Compara")


def test_the_audit_columns_are_found_by_name_or_by_position(tmp_path):
    from cichlid_homology.famsurvey import read_audit
    (tmp_path / "a.tsv").write_text("family\tstatus\nX\tok\nY\ttiny_tree\n")
    (tmp_path / "b.tsv").write_text("fam_id\tresult\textra\nX\tok\t1\n")
    assert read_audit(tmp_path / "a.tsv")[0] == {"X": "ok", "Y": "tiny_tree"}
    assert read_audit(tmp_path / "b.tsv")[0] == {"X": "ok"} and read_audit(tmp_path / "b.tsv")[1][0] == "fam_id"
    assert read_audit(tmp_path / "missing.tsv") == ({}, [])


def test_a_garbled_reference_tree_is_counted_as_a_family_with_no_leaves_not_fatal(tmp_path):
    from cichlid_homology.famsurvey import run_reference_survey
    prep = make_reference(tmp_path)
    (prep / "trees" / "C3.nwk").write_bytes(b"\xff\xfe\x00")
    rows, diag = run_reference_survey(prep, LS(), FakeClades(), names=NAMES)
    assert diag["families_with_no_leaves"] == 1 and rows[2]["n_leaves"] == 0 and rows[0]["n_human"] == 1       # garbled text has no leaves


def test_a_family_with_only_an_alignment_is_read_from_its_sequence_names_and_a_tree_wins_when_both_exist(tmp_path):
    from cichlid_homology.famsurvey import alignment_names, reference_files, run_reference_survey
    prep = tmp_path / "prepared"
    (prep / "trees").mkdir(parents=True)
    (prep / "trees" / "A1.aln").write_text(">Hsap_1/1-300\nMKV--\n>Hsap_2/1-290 extra\nMKV--\n>Aang_1/5-310\nMKVAA\n>Drer_1\nMKVAA\n")
    (prep / "trees" / "A2.aln").write_text(">Hsap_1\nMK\n>Mmus_1\nMK\n")
    (prep / "trees" / "A2.nhx").write_text("(Hsap_1,(Mmus_1,Drer_1));")           # the tree has a leaf the alignment lacks: the tree wins
    (prep / "trees" / "A3.sto").write_text("# STOCKHOLM 1.0\nAang_1/1-9 MKV\nDrer_1/1-9 MKV\n#=GC RF xxx\n//\n")
    (prep / "trees" / "notes.txt").write_text("not a family file")
    by, ext = reference_files(prep / "trees")
    assert set(by) == {"A1", "A2", "A3", "notes"} and ext == {".aln": 2, ".nhx": 1, ".sto": 1, ".txt": 1}
    assert alignment_names(prep / "trees" / "A1.aln") == ["Hsap_1", "Hsap_2", "Aang_1", "Drer_1"]
    assert alignment_names(prep / "trees" / "A3.sto") == ["Aang_1", "Drer_1"]
    rows, diag = run_reference_survey(prep, LS(), FakeClades(), names=NAMES)
    r = {x["family"]: x for x in rows}
    assert (r["A1"]["n_human"], r["A1"]["n_teleost_ref_leaves"], r["A1"]["leaf_source"]) == (2, 2, "alignment")
    assert (r["A2"]["n_leaves"], r["A2"]["n_teleost_ref_leaves"], r["A2"]["leaf_source"]) == (3, 1, "tree")
    assert r["A3"]["leaf_source"] == "alignment" and r["A3"]["n_human"] == 0
    assert diag["families_with_a_tree_file"] == 1 and diag["families_with_an_alignment_only"] == 2
    assert r["notes"]["leaf_source"] == "none" and diag["unreadable"] == 1
