import random
import time
from types import SimpleNamespace

from helpers import FakeClades, PREFIX2SP, PREFIXES, st
from cichlid_homology.dating import Params
from cichlid_homology.pilot import process_family


class FakeLS:
    def species(self, leaf):
        return PREFIX2SP.get(leaf.split("_")[0])

    def is_project_gene(self, leaf):
        return leaf.split("_")[0] in PREFIXES and "_" in leaf


def fam(newick, rooted=False, mode="place"):
    from cichlid_homology.tree import leaves, parse_newick
    return SimpleNamespace(tree=newick, leaves=[x.name for x in leaves(parse_newick(newick))],
                           rooted=rooted, tree_mode=mode, family_id="F1")


def test_process_family_end_to_end_with_shadow_and_unrooted():
    newick = ("((Hsap_1,Mmus_1),Locu_1,(Drer_1,((Gacu_1,GacuREF_1),(Olat_1,(Mzeb_1,Mzeb_2,Acal_1)))));")
    ls = FakeLS()
    ls.species = lambda leaf: ("gacu" if leaf.startswith("Gacu") else PREFIX2SP.get(leaf.split("_")[0]))
    ls.is_project_gene = lambda leaf: leaf in {"Mzeb_1", "Mzeb_2", "Acal_1", "Gacu_1", "Olat_1"}
    rows, info = process_family(fam(newick), ls, {"gacu"}, st(), FakeClades(), Params(),
                                "Mzeb", PREFIXES, [0.25], None)
    assert info["dropped"] == 1 and info["how"] == "reconciliation" and info["cost"] is not None
    assert [r["gene_id"] for r in rows] == ["Mzeb_1", "Mzeb_2"]
    assert rows[0]["human_relationship"] == "1:many"
    assert rows[0]["human_dup_class"] == "unresolved_place_mode"


def test_big_caterpillar_is_fast_and_does_not_recurse():
    rnd = random.Random(1)
    names = [f"{rnd.choice(['Hsap','Mmus','Locu','Aang','Drer','Trub','Olat','Mzeb','Acal'])}_{i}"
             for i in range(1500)]
    newick = names[0]
    for n in names[1:]:
        newick = f"({newick},{n})"
    newick += ";"
    ls = FakeLS()
    t0 = time.time()
    rows, info = process_family(fam(newick), ls, set(), st(), FakeClades(), Params(k_roots=2),
                                "Mzeb", PREFIXES, [], None)
    assert info["leaves"] == 1500 and rows
    assert time.time() - t0 < 120
    print("seconds", round(time.time() - t0, 1), "genes", len(rows))
