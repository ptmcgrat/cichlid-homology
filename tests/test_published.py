from cichlid_homology import published as PB
from cichlid_homology.genetable import COLUMNS


def row(gid, zsym, clade, fam, status="single copy (duplicate lost)", partner=""):
    r = {c: "" for c in COLUMNS}
    r.update(gene_id=gid, zebrafish_symbol=zsym, clade_id=clade, family=fam, ohnolog_status=status, ohnolog_gene=partner)
    return r


PUBS = {"Strict": {frozenset(("gnb1a", "gnb1b")), frozenset(("irx1a", "irx1b")), frozenset(("ascl1a", "ascl1b")), frozenset(("abc1a", "abc1b")), frozenset(("old1a", "old1b")), frozenset(("sol1a", "sol1b")),
                   frozenset(("xyz1a", "xyz1b")), frozenset(("same1a", "same1b"))}}


def test_recall_counts_recovered_pairs_and_says_where_the_others_are():
    rows = [row("g1", "gnb1a", "F1:G1", "F1", "paired", "g2"), row("g2", "gnb1b", "F1:G1", "F1", "paired", "g1"),                                            # recovered
            row("i1", "irx1a", "F2:G1", "F2", "expansion"), row("i2", "irx1b", "F2:G1", "F2", "expansion"),                                               # in one clade, not partners
            row("a1", "ascl1a", "F3:G1", "F3"), row("a2", "ascl1b", "F4:G1", "F4"),                                                                       # different families
            row("b1", "abc1a", "F5:G1", "F5"), row("b2", "abc1b", "F5:G2", "F5"),                                                                         # one family, different clades
            row("o1", "old1a", "F6:G1", "F6", PB.OLDER, "o2"), row("o2", "old1b", "F6:G1", "F6", PB.OLDER, "o1"),                              # labelled older
            row("s1", "sol1a", "F7:G1", "F7"),                                                                                                            # nothing for sol1b: not testable
            row("h1", "same1a", "F8:G1", "F8", "paired", "h2"), row("h2", "same1a", "F8:G1", "F8", "paired", "h1")]                                       # both closest to one zebrafish gene
    r = PB.recall(rows, PUBS["Strict"])
    assert (r["published"], r["represented"], r["recovered"], r["hidden_recoveries"]) == (8, 5, 1, 1) and r["recovered_share_of_represented"] == 0.2 and r["recovered_incl_hidden_share"] == 0.3333
    assert r["kinds"] == {"recovered": 1, "in one clade, not partners": 1, "in different families": 1, "in one family, different clades": 1, "recovered, labelled older duplicate": 1,
                          "not testable (a zebrafish gene has no gene of ours)": 2, "recovered, both genes closest to one zebrafish gene": 1}


def test_precision_separates_published_plausible_and_unrelated_pairs():
    rows = [row("g1", "gnb1a", "C:1", "F", "paired", "g2"), row("g2", "gnb1b", "C:1", "F", "paired", "g1"),
            row("p1", "foo1a", "C:2", "F", "paired", "p2"), row("p2", "foo1b", "C:2", "F", "paired", "p1"),
            row("u1", "kcna1", "C:3", "F", "paired", "u2"), row("u2", "zxcv9", "C:3", "F", "paired", "u1"),
            row("s1", "same1a", "C:4", "F", "paired", "s2"), row("s2", "same1a", "C:4", "F", "paired", "s1")]
    p = PB.precision(rows, PUBS)
    assert p["pairs"] == 4 and p["evaluable"] == 3 and p["published"] == 1 and p["published_or_a_b_stem"] == 2
    assert p["kinds"] == {"published, strict": 1, "not published; a/b names of one stem": 1, "not published; unrelated names": 1, "both closest to the same zebrafish gene": 1}
