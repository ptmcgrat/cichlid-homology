import json
import re
from pathlib import Path

from helpers import FakeClades, PREFIX2SP
from cichlid_homology import nhx
from cichlid_homology.lineages import analyse
from cichlid_homology.suppmap import attach_support, reroot_like
from cichlid_homology.tree import leaves, parse_newick

DATA = Path(__file__).parent / "data"
FOCAL = {"Mzeb", "Acal"}
is_focal = lambda name: name.split("_")[0] in FOCAL and "_" in name


def tr(newick):
    root = parse_newick(newick)
    for x in leaves(root):
        x.sp = PREFIX2SP.get(x.name.split("_")[0])
    return root


def an(newick, **kw):
    root = tr(newick)
    return root, analyse(root, is_focal, FakeClades(), **kw)


def test_two_teleost_lineages_and_the_ancestral_copy_outside():
    root, a = an("((Hsap_1,Mmus_1),(Locu_1,((Aang_1,Drer_1,Trub_1),(Aang_2,Drer_2,Trub_2))));")
    assert a.mode == "split" and a.split.kind == "teleost" and a.letters == ["A", "B"] and a.split.overlap == 1.0
    assert a.counts["aang"] == [0, 1, 1, 0] and a.counts["homo_sapiens"] == [1, 0, 0, 0] and a.counts["locu"] == [1, 0, 0, 0]
    assert [(l["id"], l["species"], l["genes"]) for l in a.lineages] == [("A", 3, 3), ("B", 3, 3)]
    assert a.outgroup_inside == [] and {x["gene"] for x in a.outgroup_outside} == {"Hsap_1", "Mmus_1", "Locu_1"}


def test_a_gar_gene_inside_one_lineage_does_not_stop_the_split_and_stays_the_ancestral_copy():
    root, a = an("((Hsap_1,Mmus_1),((Locu_1,(Aang_1,Drer_1,Trub_1)),(Aang_2,Drer_2,Trub_2)));")
    assert a.mode == "split" and a.split.kind == "teleost" and a.split.overlap == 1.0
    assert a.counts["locu"] == [1, 0, 0, 0]                                            # counted as the ancestral copy
    assert a.outgroup_inside == [{"gene": "Locu_1", "species": "locu", "lineage": "A", "attach_support": None}]
    assert {x["gene"] for x in a.outgroup_outside} == {"Hsap_1", "Mmus_1"}


def test_no_split_is_a_single_lineage_and_a_species_level_cherry_is_not_a_split():
    root, a = an("((Hsap_1,Mmus_1),(Locu_1,(Aang_1,Drer_1,Trub_1)));")
    assert a.mode == "single" and a.letters == ["A"] and a.counts["aang"] == [0, 1, 0] and a.counts["locu"] == [1, 0, 0]
    root, a = an("((Hsap_1,Mmus_1),(Aang_1,(Drer_1,Drer_2)));")
    assert a.mode == "single"                                                          # two genes of one species: no lineages


def test_a_duplication_among_percomorphs_only_is_a_later_split():
    root, a = an("((Hsap_1,Mmus_1),(Aang_1,(Drer_1,((Trub_1,Gacu_1),(Trub_2,Gacu_2)))));")
    assert a.mode == "split" and a.split.kind == "later" and a.split.level == "percomorph" and a.counts["trub"] == [0, 1, 1, 0]


def test_a_teleost_split_wins_over_a_later_one_higher_up_or_lower_down():
    root, a = an("((Hsap_1,Mmus_1),((Aang_1,Drer_1,((Trub_1,Gacu_1),(Trub_2,Gacu_2))),(Aang_2,Drer_2)));")
    assert a.split.kind == "teleost" and a.counts["trub"] == [0, 2, 0, 0] and a.counts["aang"] == [0, 1, 1, 0]


def test_a_gene_that_shares_no_species_is_a_stray_not_a_lineage():
    root, a = an("((Hsap_1,Mmus_1),((Aang_1,Drer_1,Trub_1),(Aang_2,Drer_2,Trub_2),Gacu_1));")
    assert a.letters == ["A", "B"] and a.counts["gacu"] == [0, 0, 0, 1]


def test_our_own_genes_take_the_lineage_of_their_position_and_do_not_move_the_split():
    root, a = an("((Hsap_1,Mmus_1),((Aang_1,Drer_1,Mzeb_1),(Aang_2,Drer_2,Mzeb_2,Acal_2)));")
    assert a.letters == ["A", "B"]
    sp = {x.name: a.member[id(x)] for x in leaves(root)}
    assert sp["Mzeb_1"] == 0 and sp["Mzeb_2"] == 1 and sp["Acal_2"] == 1
    root2, b = an("((Hsap_1,Mmus_1),((Aang_1,Drer_1,Mzeb_1),(Aang_2,Drer_2,Mzeb_2,Acal_2)));", participates=lambda x: True)
    assert b.split.overlap == 0.75                                                     # a de novo tree: every leaf takes part, so Acal_2 adds a species to side B


def test_supports_are_read_for_the_crowns_and_for_where_an_outgroup_attaches():
    root = tr("((Hsap_1,Mmus_1),((Locu_1,(Aang_1,Drer_1,Trub_1)),(Aang_2,Drer_2,Trub_2)));")
    ref = {x.name for x in leaves(root)}
    sup = {frozenset({"Aang_1", "Drer_1", "Trub_1"}): (1, "N", None), frozenset({"Locu_1", "Aang_1", "Drer_1", "Trub_1"}): (25, "N", None),
           frozenset({"Aang_2", "Drer_2", "Trub_2"}): (78, "N", None)}
    m = attach_support(root, sup, ref)
    a = analyse(root, is_focal, FakeClades(), support=m)
    assert [l["crown_support"] for l in a.lineages] == [1, 78] and a.outgroup_inside[0]["attach_support"] == 25


def test_reroot_like_puts_the_root_where_the_reference_root_is():
    ref = parse_newick("((a,b),(c,(d,e)));")
    placed = parse_newick("(((a,b),F1),c,(d,(e,F2)));")                              # unrooted, with two placed genes
    rooted = reroot_like(placed, ref)
    sides = [frozenset(x.name for x in leaves(c) if x.name in set("abcde")) for c in rooted.children]
    assert sorted(map(sorted, sides)) == [["a", "b"], ["c", "d", "e"]]
    assert reroot_like(parse_newick("((a,b),(c,d));"), parse_newick("((a,c),(b,d));")) is None


def test_the_reader_and_the_real_irx1_tree_agree_with_what_we_saw_in_the_figure():
    pan = {s["name"]: s for s in json.load(open(DATA / "panel_species.json"))}
    pfx = {v: k for k, v in json.load(open(DATA / "panel_prefixes.json")).items()}

    class Clades:
        order = ["outside_vertebrates", "cyclostome", "gnathostome", "teleost", "percomorph", "ovalentaria", "cichlid"]

        def level(self, s):
            return pan[s]["level"]

        def common_level(self, names):
            return min((self.level(s) for s in names), key=self.order.index)

    root = parse_newick((DATA / "irx1_prepared.nhx").read_text().strip())
    for x in leaves(root):
        x.sp = pfx[re.match(r"^(ENS[A-Z]*?P)\d", x.name).group(1)]
    ref = {x.name for x in leaves(root)}
    blocks = list(nhx.iter_emf(DATA / "irx1_raw.emf"))
    assert len(blocks) == 1 and len(blocks[0][0]) == 253
    smap = nhx.block_supports(blocks[0][1], ref)
    sup = attach_support(root, smap, ref)
    a = analyse(root, lambda n: False, Clades(), support=sup)
    assert a.mode == "split" and a.split.kind == "teleost" and a.split.level == "teleost" and round(a.split.overlap, 2) == 0.92
    assert [(l["id"], l["species"], l["genes"]) for l in a.lineages] == [("A", 35, 35), ("B", 34, 34)]
    assert [l["crown_support"] for l in a.lineages] == [1, 78]
    assert abs(a.lineages[0]["stem"] - 0.133) < 0.001 and abs(a.lineages[1]["stem"] - 0.245) < 0.001            # the longer stem of lineage B
    assert abs(a.lineages[0]["depth"] - 0.357) < 0.001 and abs(a.lineages[1]["depth"] - 0.453) < 0.001
    assert sorted((x["species"], x["lineage"], x["attach_support"]) for x in a.outgroup_inside) == [("erpetoichthys_calabaricus", "A", 25), ("lepisosteus_oculatus", "A", 0)]
    assert a.counts["danio_rerio"] == [0, 1, 1, 0] and a.counts["ictalurus_punctatus"] == [1, 0, 0, 0] and a.counts["denticeps_clupeoides"] == [1, 0, 1, 0]      # the two genes that sit outside the split are counted outside
    assert a.counts["homo_sapiens"] == [1, 0, 0, 0] and a.counts["lepisosteus_oculatus"] == [1, 0, 0, 0]
    split_sup = sup[id(a.split.node)]
    assert split_sup[1] == "Y" and split_sup[0] == 0                                   # Ensembl calls it a duplication; its bootstrap is 0


def test_nhx_tags_are_read_on_leaves_and_inner_nodes_and_an_emf_block_is_found_by_any_of_its_proteins(tmp_path):
    t = nhx.parse_nhx("((a:0.1[&&NHX:D=N:G=g1],b:0.2)[&&NHX:D=Y:B=55]:0.3,c:0.4)[&&NHX:B=7];")
    assert t.tags == {"B": "7"} and t.children[0].tags == {"D": "Y", "B": "55"} and t.children[0].children[0].tags["G"] == "g1" and t.children[0].children[0].length == 0.1
    sup = nhx.block_supports("((a,b)[&&NHX:D=Y:B=55],(c,d)[&&NHX:D=N:B=3])[&&NHX:B=7];", {"a", "b", "c"})
    assert sup == {frozenset({"a", "b"}): (55, "Y", None), frozenset({"a", "b", "c"}): (7, None, None)}        # (c,d) holds one reference protein only
    emf = tmp_path / "x.emf"
    emf.write_text("SEQ sp1 p1 1 1 2 1 g1 s1\nSEQ sp2 p2 1 1 2 1 g2\n\nDATA\n(p1,p2)[&&NHX:B=90];\n//\nSEQ sp1 q1 1 1 2 1 h1 t1\nSEQ sp2 q2 1 1 2 1 h2 t2\n\nDATA\n(q1,q2)[&&NHX:B=10];\n//\n")
    found = nhx.collect_support(emf, {"famB": {"q2", "q1"}})
    assert list(found) == ["famB"] and found["famB"]["support"] == {frozenset({"q1", "q2"}): (10, None, None)} and found["famB"]["symbols"] == {"q1": "t1", "q2": "t2"}


def _deep(newick):
    from cichlid_homology.figrecord import deep_order
    pan = {s["name"]: s for s in json.load(open(DATA / "panel_species.json"))}

    class C:
        order = ["outside_vertebrates", "cyclostome", "gnathostome", "teleost", "percomorph", "ovalentaria", "cichlid"]

        def level(self, s):
            return pan[s]["level"]

        def lineage(self, s):
            return {"tetrapod": "sarcopterygian", "coelacanth": "sarcopterygian", "cartilaginous fish": "chondrichthyan"}.get(pan[s]["group"], "actinopterygian")
    root = parse_newick(newick)
    for x in leaves(root):
        x.sp = x.name.split("#")[0]
    from cichlid_homology.cladeprofile import group_label
    return deep_order(root, lambda sp: group_label(sp, C(), "latimeria_chalumnae"), lambda n: False)


IDEAL = ("(callorhinchus_milii#1,((latimeria_chalumnae#1,homo_sapiens#1),(erpetoichthys_calabaricus#1,(lepisosteus_oculatus#1,(danio_rerio#1,gasterosteus_aculeatus#1)))));")


def test_the_deep_order_of_an_ideal_tree_is_fully_recovered_and_each_failure_is_named():
    d = _deep(IDEAL)
    assert (d["tested"], d["passed"], d["score"]) == (4, 4, 1.0)
    d = _deep("(callorhinchus_milii#1,((latimeria_chalumnae#1,homo_sapiens#1),(erpetoichthys_calabaricus#1,(danio_rerio#1,(lepisosteus_oculatus#1,gasterosteus_aculeatus#1)))));")
    assert d["gar_outside_teleosts"] is False and d["reedfish_outside_gar_teleosts"] is True and d["score"] == 0.75          # gar nested among the teleosts
    d = _deep("(callorhinchus_milii#1,((latimeria_chalumnae#1,(danio_rerio#1,gasterosteus_aculeatus#1)),(homo_sapiens#1,(erpetoichthys_calabaricus#1,lepisosteus_oculatus#1))));")
    assert d["coelacanth_with_tetrapods"] is False
    d = _deep("((callorhinchus_milii#1,homo_sapiens#1),((latimeria_chalumnae#1,homo_sapiens#2),(erpetoichthys_calabaricus#1,(lepisosteus_oculatus#1,(danio_rerio#1,gasterosteus_aculeatus#1)))));")
    assert d["shark_outside"] is False and d["coelacanth_with_tetrapods"] is True


def test_a_relationship_is_tested_only_when_its_species_are_there_and_the_real_irx1_tree_recovers_two_of_four():
    d = _deep("((danio_rerio#1,gasterosteus_aculeatus#1),(danio_rerio#2,gasterosteus_aculeatus#2));")
    assert d["tested"] == 0 and d["score"] is None and d["shark_outside"] is None
    pan = {s["name"]: s for s in json.load(open(DATA / "panel_species.json"))}
    pfx = {v: k for k, v in json.load(open(DATA / "panel_prefixes.json")).items()}
    root = parse_newick((DATA / "irx1_prepared.nhx").read_text().strip())
    for x in leaves(root):
        x.sp = pfx[re.match(r"^(ENS[A-Z]*?P)\d", x.name).group(1)]
    from cichlid_homology.cladeprofile import group_label
    from cichlid_homology.figrecord import deep_order

    class C:
        def level(self, s):
            return pan[s]["level"]

        def lineage(self, s):
            return {"tetrapod": "sarcopterygian", "coelacanth": "sarcopterygian", "cartilaginous fish": "chondrichthyan"}.get(pan[s]["group"], "actinopterygian")
    d = deep_order(root, lambda sp: group_label(sp, C(), "latimeria_chalumnae"), lambda n: False)
    assert d["score"] == 0.5 and (d["tested"], d["passed"]) == (4, 2) and d["gar_outside_teleosts"] is False and d["reedfish_outside_gar_teleosts"] is False
    assert d["tolerant"]["score"] == 0.5 and d["tolerant"]["gar_outside_teleosts"] is False                     # a real inside-lineage case is not rescued by the tolerance


def test_one_stray_teleost_gene_breaks_the_strict_test_but_not_the_tolerant_one():
    D = ",".join(f"danio_rerio#{i}" for i in range(1, 10))
    nw = f"(callorhinchus_milii#1,((latimeria_chalumnae#1,homo_sapiens#1),((erpetoichthys_calabaricus#1,(lepisosteus_oculatus#1,({D}))),gasterosteus_aculeatus#1)));"
    d = _deep(nw)
    assert d["gar_outside_teleosts"] is False and d["reedfish_outside_gar_teleosts"] is False and d["score"] == 0.5             # strict: the stray gene puts gar inside the teleost clade
    t = d["tolerant"]
    assert t["gar_outside_teleosts"] is True and t["reedfish_outside_gar_teleosts"] is True and t["shark_outside"] is True and t["score"] == 1.0
    two = ",".join(f"danio_rerio#{i}" for i in range(1, 6)) + ",gasterosteus_aculeatus#1,gasterosteus_aculeatus#2,gasterosteus_aculeatus#3"
    d = _deep(f"(callorhinchus_milii#1,((latimeria_chalumnae#1,homo_sapiens#1),(erpetoichthys_calabaricus#1,((lepisosteus_oculatus#1,(danio_rerio#6,danio_rerio#7)),({two})))));")
    assert d["tolerant"]["gar_outside_teleosts"] is False                                                                   # gar inside half of the teleosts is not a stray gene
