"""The whole pilot on a fake layout: families on disk, assignment and placement tables, a genes table,
fake coordinator objects. Runs sample mode and all-genes mode and checks the outputs."""
import csv
import json
from types import SimpleNamespace
from unittest import mock

from helpers import FakeClades, PREFIX2SP, PREFIXES, st
from test_names import Rec, Table
from cichlid_homology import pilot
from cichlid_homology.dating import Params
from cichlid_homology.tree import leaves, parse_newick

TREE1 = "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Gacu_1,(Olat_1,(Mzeb_1,Acal_1)))))));"
FISH1 = "(Aang_1,(Drer_1,(Gacu_1,(Olat_1,(Mzeb_2,Acal_2)))))"
FISH2 = "(Aang_2,(Drer_2,(Gacu_2,(Olat_2,(Mzeb_3,Acal_3)))))"
TREE2 = f"((Hsap_1,Mmus_1),(Locu_1,({FISH1},{FISH2})));"
TREE3 = "((Hsap_1,Mmus_1),(Locu_1,(Aang_1,(Drer_1,(Gacu_1,(Olat_1,((Mzeb_4,Acal_4),(Mzeb_5,Acal_5))))))));"


class FakeLS:
    def species(self, leaf):
        return PREFIX2SP.get(leaf.split("_")[0])

    def is_project_gene(self, leaf):
        return leaf.split("_")[0] in PREFIXES and "_" in leaf

    def shadowed(self):
        return set()


def build_layout(tmp):
    genes = {"Mzeb_1": ("c1", 100, 200), "Mzeb_2": ("c1", 1000, 1100), "Mzeb_3": ("c9", 100, 200),
             "Mzeb_4": ("c2", 100, 200), "Mzeb_5": ("c2", 300, 400)}
    gt = tmp / "Mzeb.genes.tsv"
    with open(gt, "w") as fh:
        fh.write("gene_id\tsymbol\tseqid\tstart\tend\tstrand\tpartial\tlow_quality\tinternal_stops\tprotein_length\t"
                 "description\tchosen_protein\tmitochondrial\n")
        for g, (sq, a, b) in genes.items():
            fh.write(f"{g}\tgene{g[-1]}\t{sq}\t{a}\t{b}\t+\tFalse\tFalse\t0\t300\tdesc\tp{g}\tFalse\n")
    results = {}
    for part, ids in (("TreeFam", ("TF1", "TF2", "TF3")), ("Compara", ("C1", "C2", "C3"))):
        root = tmp / part
        (root / "families").mkdir(parents=True)
        (root / "native" / "summary").mkdir(parents=True)
        fam_genes = {ids[0]: ["Mzeb_1"], ids[1]: ["Mzeb_2", "Mzeb_3"], ids[2]: ["Mzeb_4", "Mzeb_5"]}
        for fid, tree in zip(ids, (TREE1, TREE2, TREE3)):
            d = root / "families" / fid
            d.mkdir()
            (d / "tree.nwk").write_text(tree)
            (d / "combined.aln").write_text("".join(f">{x.name}\nAC\n" for x in leaves(parse_newick(tree))))
            (d / "done.json").write_text(json.dumps({"rooted": False, "tree_mode": "place"}))
        with open(root / "assign.tsv", "w") as fh:
            fh.write("gene_id\tbest_family\tflags\tfusion_partner\n")
            for fid, gl in fam_genes.items():
                for g in gl:
                    fh.write(f"{g}\t{fid}\t\t\n")
        with open(root / "native" / "summary" / "gene_placement.tsv", "w") as fh:
            fh.write("gene_id\tfamily\tmode\tstatus\tlwr\tn_placements\n")
            for fid, gl in fam_genes.items():
                for g in gl:
                    fh.write(f"{g}\t{fid}\tplace\tplaced\t0.9\t1\n")
        results[part] = SimpleNamespace(
            assign=lambda prefix, root=root: SimpleNamespace(local=root / "assign.tsv"),
            family=lambda f, root=root: SimpleNamespace(local=root / "families" / f),
            native=lambda root=root: SimpleNamespace(local=root / "native"))
    return SimpleNamespace(results=lambda release, part: results[part],
                           genes_table=lambda prefix: SimpleNamespace(local=gt))


def run(tmp, **kw):
    layout = build_layout(tmp)
    cat = SimpleNamespace(species=lambda release=None: [SimpleNamespace(prefix=p) for p in PREFIXES])
    names = Table({"Hsap_1": Rec("GENE", "HGNC:1", "current"), "Drer_1": Rec("gene", "ZDB-1", "current"), "Drer_2": Rec("gene2", "ZDB-2", "current")})
    out = tmp / "out"
    prm = kw.pop("prm", Params())
    with mock.patch.object(pilot, "load_species_tree", lambda: st()), \
         mock.patch.object(pilot, "load_leaf_species", lambda layout, part, only=None: FakeLS()), \
         mock.patch.object(pilot, "load_names", lambda layout, part: names):
        summary = pilot.run_pilot(layout, cat, FakeClades(), "Mzeb", 100, prm, [0.25], out, n_denovo=0, **kw)
    return summary, out


def rows(path):
    return {r["gene_id"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def test_whole_pilot_in_all_genes_mode(tmp_path):
    summary, out = run(tmp_path, all_genes=True)
    assert summary["run"]["mode"] == "all genes" and summary["run"]["species_genes_total"] == 5
    t = rows(out / "pilot_genes.TreeFam.tsv")
    assert set(t) == {"Mzeb_1", "Mzeb_2", "Mzeb_3", "Mzeb_4", "Mzeb_5"}
    assert t["Mzeb_1"]["a_human_group"] == "OneToOne" and t["Mzeb_1"]["a_zebrafish_group"] == "OneToOne"
    assert t["Mzeb_2"]["a_human_group"] == "OneToTwo" and t["Mzeb_2"]["a_zebrafish_group"] == "OneToOne"
    assert t["Mzeb_4"]["a_human_group"] == "OneToMany" or t["Mzeb_4"]["a_human_cat"] == "OneToTwo"
    assert t["Mzeb_2"]["a_human_best_symbol"] == "GENE" and t["Mzeb_2"]["a_zebrafish_best_symbol"] == "gene"
    # origin: Mzeb_2 and Mzeb_3 are on different sequences (the genome-duplication signature), 4 and 5 are adjacent
    assert t["Mzeb_2"]["orig_class"] == "different_seqid" and t["Mzeb_3"]["orig_class"] == "different_seqid"
    assert t["Mzeb_4"]["orig_class"] == "tandem_adjacent" and t["Mzeb_1"]["orig_class"] == "single_copy"
    assert t["Mzeb_2"]["b_label"] == "retained_pair" and t["Mzeb_4"]["b_label"] == "cichlid_specific_duplication"
    assert "risk_flags" in t["Mzeb_4"] and t["Mzeb_4"]["gt_partial"] == "False"
    cand = rows(out / "pilot_candidates.TreeFam.tsv")
    assert set(cand) == {"Mzeb_4", "Mzeb_5"} and cand["Mzeb_4"]["orig_class"] == "tandem_adjacent"
    report = (out / "pilot_report.md").read_text()
    for h in ("Mode: **all genes**", "## Where the extra copies sit", "## Cichlid-specific candidates",
              "## Agreement between the methods"):
        assert h in report
    A = summary["A"]
    assert A["candidates_in_both_methods"] == 2 and A["agreement"]["human"]["same_category_pct"] == 100.0
    assert A["origin"]["TreeFam"]["by_teleost_profile_label"]["retained_pair"]["different_seqid"]["n"] == 2


def test_sample_mode_takes_the_whole_family_of_a_sampled_gene(tmp_path):
    summary, out = run(tmp_path, all_genes=False)
    assert summary["run"]["mode"] == "sample"
    assert set(rows(out / "pilot_genes.Compara.tsv")) == {"Mzeb_1", "Mzeb_2", "Mzeb_3", "Mzeb_4", "Mzeb_5"}


def test_parsimony_columns_labels_and_summary_in_the_whole_pilot(tmp_path):
    summary, out = run(tmp_path, all_genes=True)
    t = rows(out / "pilot_genes.TreeFam.tsv")
    assert t["Mzeb_1"]["p_label"] == "lost_duplicate" and t["Mzeb_1"]["p_robust"] == "True"
    assert t["Mzeb_2"]["p_label"] == "retained_ohnolog" and t["Mzeb_2"]["p_gain_levels"] == "teleost_ancestor:1"
    assert t["Mzeb_4"]["p_label"] == "lost_duplicate" and t["Mzeb_4"]["p_gain_levels"] == "cichlid:1"   # a later duplicate
    P = summary["A"]["parsimony"]["TreeFam"]
    assert P["scored"] == 5 and P["labels"]["retained_ohnolog"]["n"] == 2
    assert P["later_duplicates_after_the_teleost_ancestor"]["lost_duplicate"]["with_a_later_gain"] == 2
    assert P["check_copy_positions"]["retained_ohnolog"]["different_sequence_pct"] == 100.0
    assert summary["A"]["parsimony_agreement"]["agree_pct_when_both_called"] == 100.0
    assert "## Lost duplicate or retained ohnolog by parsimony" in (out / "pilot_report.md").read_text()


def test_tiers_in_the_whole_pilot_with_and_without_published_pairs(tmp_path):
    summary, out = run(tmp_path, all_genes=True)
    t = rows(out / "pilot_genes.TreeFam.tsv")
    assert (t["Mzeb_1"]["t_row"], t["Mzeb_1"]["t_tier"]) == ("R1", "core")
    assert (t["Mzeb_2"]["t_row"], t["Mzeb_2"]["t_tier"]) == ("R2", "core") and t["Mzeb_2"]["t_positions"] == "dispersed"
    assert (t["Mzeb_4"]["t_row"], t["Mzeb_4"]["t_tier"]) == ("R5", "core")
    assert [t[g]["t_category"] for g in ("Mzeb_1", "Mzeb_2", "Mzeb_4")] == ["Lost duplicate", "Both copies kept (best guess)", "Later duplicate, in cichlids"]
    T = summary["A"]["tiers"]["TreeFam"]
    assert T["placed"] == 5 and T["core"] == 5 and T["published_pair_evidence_loaded"] is False
    assert "## Confidence tiers for B1" in (out / "pilot_report.md").read_text()
    d = tmp_path / "pairs"
    d.mkdir()
    for lvl in ("Strict", "Relaxed"):
        (d / f"drerio_Pairs_{lvl}_3R.txt").write_text("Ohno1\tOhno2\tSymbol1\tSymbol2\tOg\tSelf\tOutgroup Support\tDuplication time\tGene type\n"
                                                      "E1\tE2\tgene\tgene2\t1e-9\t1e-5\t7\tFishWGD\tprotein_coding\n")
    summary2, out2 = run(tmp_path / "second", all_genes=True, prm=Params(ohnologs_dir=str(d))) if (tmp_path / "second").mkdir() is None else (None, None)
    t2 = rows(out2 / "pilot_genes.TreeFam.tsv")
    assert t2["Mzeb_2"]["t_pairs_relaxed_n"] == "1" and t2["Mzeb_2"]["t_pairs_species"] == "zebrafish"
    assert summary2["ohnologs_species_loaded"] == ["zebrafish"]
    assert summary2["A"]["tiers"]["TreeFam"]["published_pair_evidence_loaded"] is True


def test_tree_based_dating_in_the_whole_pilot(tmp_path):
    summary, out = run(tmp_path, all_genes=True)
    t = rows(out / "pilot_genes.TreeFam.tsv")
    assert t["Mzeb_1"]["d_call"] == "lost_duplicate" and t["Mzeb_1"]["d_root"] == "single"
    # the fake placement tree holds two reference fish clades with early teleosts in both: a pair at the teleost root
    assert t["Mzeb_2"]["d_call"] == "pair" and t["Mzeb_2"]["d_support"] == "both_early"
    assert {t["Mzeb_2"]["d_gene_lineage"], t["Mzeb_3"]["d_gene_lineage"]} == {"A", "B"}     # the letters are arbitrary
    assert t["Mzeb_2"]["d_focal_lineages"] == "2"
    D = summary["A"]["treedate"]["TreeFam"]
    assert D["scored"] == 5 and D["calls"]["pair"] == 2 and D["by_category"]["Both copies kept (best guess)"] == {"pair": 2}
    assert "## Tree-based dating" in (out / "pilot_report.md").read_text()
