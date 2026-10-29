import json
import re
from pathlib import Path

from helpers import PREFIX2SP, SPECIES_TREE, FakeClades, st
from cichlid_homology import figpilot
from cichlid_homology.figrecord import build_record
from cichlid_homology.read import load_family
from cichlid_homology.tree import leaves, parse_newick

REF = "((Hsap_1,Mmus_1),(Locu_1,((Aang_1,(Drer_1,Trub_1)),(Aang_2,(Drer_2,Trub_2)))));"
PLACED = "(Hsap_1,Mmus_1,(Locu_1,((Aang_1,(Drer_1,Mzeb_1,Trub_1)),(Aang_2,(Drer_2,Mzeb_2,Trub_2)))));"          # unrooted: our genes placed next to the zebrafish genes
SYM = {"Drer_1": "irx1a", "Drer_2": "irx1b"}
DUMP_TREE = re.sub(r"\)", ")[&&NHX:D=N:B=77]", REF)


class Clades(FakeClades):
    def lineage(self, s):
        return "sarcopterygian" if s in ("homo_sapiens", "mmus") else super().lineage(s)


class LS:
    def species(self, leaf):
        return PREFIX2SP.get(leaf.split("_")[0])

    def is_project_gene(self, leaf):
        return leaf.split("_")[0] in {"Mzeb", "Acal"} and "_" in leaf


def rows_of(**extra):
    base = {"family": "F1", "clade_id": "F1:G1", "symbol": "", "ohnolog_status": "single copy (duplicate lost)", "ohnolog_confidence": "", "published_pair": "", "rel_human": "", "human_clean_tier": "",
            "out_paralog": "no", "clade_type": "human clade (1 human gene)"}
    return {g: {**base, "gene_id": g, **kw} for g, kw in extra.items()}


def test_the_sample_takes_each_stratum_once_and_finds_the_a_b_named_single_copy_pairs():
    rows = rows_of(g1={"symbol": "tbx2a"}, g2={"symbol": "tbx2b"}, g3={"symbol": "foxa"}, g4={"ohnolog_status": "paired", "ohnolog_confidence": "high", "published_pair": "the two names are a published pair"},
                   g5={"ohnolog_status": "paired", "ohnolog_confidence": "high", "published_pair": "the two names are a published pair"}, g6={"rel_human": "1:1", "human_clean_tier": "clean"},
                   g7={"ohnolog_status": "expansion"}, g8={"family": ""}, g9={"symbol": "LOC1a"}, g10={"symbol": "LOC1b"})
    assert figpilot.ab_pairs(rows) == [("g1", "g2")]
    ch = figpilot.select_genes(rows, per_stratum=1, seed=1)
    by = {c["stratum"]: [x["gene"] for x in ch if x["stratum"] == c["stratum"]] for c in ch}
    assert sorted(by["ab_not_ohnolog"]) == ["g1", "g2"] and len(by["published_pair"]) == 1 and by["human_1to1_clean"] == ["g6"] and by["expansion"] == ["g7"]
    assert len({c["gene"] for c in ch}) == len(ch) and "g8" not in {c["gene"] for c in ch}                  # a gene is drawn once; a gene with no family is not drawn
    assert [c["gene"] for c in figpilot.select_genes(rows, 1, 1, explicit=["g3", "g1"]) if c["stratum"] == "requested"] == ["g3"]


def test_the_species_axis_follows_the_species_tree_smaller_branch_first():
    order = figpilot.panel_order(SPECIES_TREE, {"aang", "danio_rerio", "homo_sapiens", "mmus", "locu", "trub", "cmil"})
    assert order[:2] == ["cmil", "homo_sapiens"] or order[0] in ("homo_sapiens", "cmil")
    assert order.index("homo_sapiens") < order.index("locu") < order.index("aang") < order.index("danio_rerio") < order.index("trub")


def test_a_clade_is_found_by_its_reference_leaves_and_a_placed_gene_does_not_make_it_inexact():
    root = parse_newick(PLACED)
    node, exact = figpilot.find_clade_node(root, frozenset({"Aang_1", "Drer_1", "Trub_1"}), lambda n: n.startswith("Mzeb"))
    assert exact and {x.name for x in leaves(node)} == {"Aang_1", "Drer_1", "Mzeb_1", "Trub_1"}
    node, exact = figpilot.find_clade_node(root, frozenset({"Aang_1", "Drer_1"}), lambda n: n.startswith("Mzeb"))
    assert not exact


def build(tmp, with_dump=True, workers=1):
    prepared = tmp / "prepared"
    (prepared / "trees").mkdir(parents=True)
    (prepared / "trees" / "F1.nhx").write_text(REF)
    froot = tmp / "families"
    d = froot / "F1"
    d.mkdir(parents=True)
    (d / "tree.nwk").write_text(PLACED)
    (d / "combined.aln").write_text("".join(f">{x.name}\nAC\n" for x in leaves(parse_newick(PLACED))))
    (d / "done.json").write_text(json.dumps({"rooted": False, "tree_mode": "place"}))
    dump = None
    if with_dump:
        dump = tmp / "dump.emf"
        seqs = "".join(f"SEQ {PREFIX2SP[x.name.split('_')[0]]} {x.name} 1 1 2 1 G{x.name} {SYM.get(x.name, '')}\n" for x in leaves(parse_newick(REF)))
        dump.write_text(seqs + f"\nDATA\n{DUMP_TREE}\n//\n")
    rows = rows_of(Mzeb_1={"symbol": "geneA", "b1_category": "Both copies kept (best guess)", "ohnolog_status": "paired", "ohnolog_symbol": "geneB"},
                   Mzeb_2={"symbol": "geneB", "b1_category": "Both copies kept (best guess)", "ohnolog_status": "paired", "ohnolog_symbol": "geneA"})
    placement = {"Mzeb_1": {"lwr": "0.93", "status": "placed", "n_placements": "2", "mode": "place"}, "Mzeb_2": {"lwr": "0.41", "status": "placed", "n_placements": "3", "mode": "place"}}
    genes = [{"gene": "Mzeb_1", "stratum": "published_pair"}, {"gene": "Mzeb_2", "stratum": "published_pair"}]
    out = tmp / "out"
    summary = figpilot.run_figure_pilot(genes=genes, rows=rows, prepared=prepared, froot=froot, placement=placement, prefixes={"Mzeb", "Acal"}, load_family=load_family, load_ls=lambda names: LS(),
                                        st=st(), clades=Clades(), master_newick=SPECIES_TREE, out_dir=out, dump=dump, workers=workers, dpi=60, progress=lambda m: None)
    return summary, out


def test_the_pilot_makes_a_record_and_a_figure_per_clade_with_support_and_our_genes_in_their_lineages(tmp_path):
    summary, out = build(tmp_path)
    assert summary["genes_requested"] == 2 and summary["clades"] == 1 and summary["figures"] == 1 and summary["support_families_found"] == 1 and not summary["errors"] and not summary["figure_errors"]
    rec = json.loads((out / "clade_records.jsonl").read_text().splitlines()[0])
    assert rec["clade"] == "F1:G1" and rec["mode"] == "split" and rec["split"]["kind"] == "teleost" and rec["root_source"] == "like the reference root"
    assert {g["gene"]: g["lineage"] for g in rec["genes"]} == {"Mzeb_1": "A", "Mzeb_2": "B"}
    assert {g["gene"]: g["lwr"] for g in rec["genes"]} == {"Mzeb_1": "0.93", "Mzeb_2": "0.41"}
    assert rec["support_available"] and [l["crown_support"] for l in rec["lineages"]] == [77, 77] and rec["outgroup"]["inside"] == [] and "no_support_values" not in rec["flags"]
    assert rec["strip"]["counts"][rec["strip"]["species"].index("aang")] == [0, 1, 1, 0] and rec["focal_strip"]["species"] == ["maylandia_zebra"]
    assert rec["focal_strip"]["counts"] == [[0, 1, 1, 0]] and rec["anchors"]["A"]["zebrafish"] == [{"gene": "Drer_1", "symbol": "irx1a"}] and rec["anchors"]["B"]["zebrafish"] == [{"gene": "Drer_2", "symbol": "irx1b"}]
    assert (out / "figures" / "F1_G1.png").stat().st_size > 5000
    idx = list(__import__("csv").DictReader(open(out / "figure_index.tsv"), delimiter="\t"))
    assert [r["gene_id"] for r in idx] == ["Mzeb_1", "Mzeb_2"] and idx[0]["split_kind"] == "teleost" and idx[0]["figure"] == "figures/F1_G1.png"
    assert "published_pair" in (out / "figure_index.html").read_text()


def test_without_the_dump_there_is_no_support_and_the_record_says_so(tmp_path):
    summary, out = build(tmp_path, with_dump=False)
    rec = json.loads((out / "clade_records.jsonl").read_text().splitlines()[0])
    assert not rec["support_available"] and "no_support_values" in rec["flags"] and rec["split"]["kind"] == "teleost" and summary["figures"] == 1


def test_a_gene_without_a_finished_family_is_skipped_and_counted(tmp_path):
    prepared = tmp_path / "prepared"
    (prepared / "trees").mkdir(parents=True)
    summary = figpilot.run_figure_pilot(genes=[{"gene": "Mzeb_9", "stratum": "x"}], rows=rows_of(Mzeb_9={"family": "F9"}), prepared=prepared, froot=tmp_path / "none", placement={}, prefixes={"Mzeb"},
                                        load_family=load_family, load_ls=lambda names: LS(), st=st(), clades=Clades(), master_newick=SPECIES_TREE, out_dir=tmp_path / "o", progress=lambda m: None)
    assert summary["skipped"] == {"no finished family": 1} and summary["clades"] == 0


def test_the_record_of_a_clade_with_no_split_is_one_lineage_and_draws(tmp_path):
    root = parse_newick("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,Drer_1,Mzeb_1)));")
    for x in leaves(root):
        x.sp = PREFIX2SP[x.name.split("_")[0]]
    rec = build_record(root, method="Compara", family="F", clade="F:G1", focal_species="maylandia_zebra", is_focal=lambda n: n.startswith("Mzeb"), clades=Clades(), panel=["homo_sapiens", "locu", "aang"],
                       species_labels={}, support=None)
    assert rec["mode"] == "single" and rec["letters"] == ["A"] and "no_split" in rec["flags"] and rec["genes"][0]["lineage"] == "A"
    from cichlid_homology.figrender import render_record
    assert render_record(rec, tmp_path / "x.png", dpi=50).stat().st_size > 3000
