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
