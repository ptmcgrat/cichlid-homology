from helpers import SPECIES_TREE
from cichlid_homology.master import newick_of, node_newick
from cichlid_homology.speciestree import SpeciesTree
from cichlid_homology.tree import leaves, parse_newick


class N:
    """Stands in for the coordinator's Node: children, name, parent, length."""
    def __init__(self, name=None):
        self.name, self.children, self.parent, self.length = name, [], None, None


def build(root):
    def conv(n):
        m = N(n.name if not n.children else None)
        for c in n.children:
            k = conv(c)
            k.parent = m
            m.children.append(k)
        return m
    return conv(root)


def test_coordinator_style_node_becomes_a_species_tree():
    node = build(parse_newick(SPECIES_TREE))
    nw = newick_of(node)
    assert nw == node_newick(node)
    st = SpeciesTree(nw)
    assert set(st.leaf_id) == {x.name for x in leaves(parse_newick(SPECIES_TREE))}
    a, b = st.id_of("maylandia_zebra"), st.id_of("acal")
    assert st.lca(a, b) == st.parent[a]          # the two cichlids are sisters


def test_big_node_does_not_recurse():
    root = N()
    cur = root
    for i in range(3000):
        leaf, nxt = N(f"s{i}"), N()
        for k in (leaf, nxt):
            k.parent = cur
            cur.children.append(k)
        cur = nxt
    cur.name = "last"
    st = SpeciesTree(newick_of(root))
    assert len(st.leaf_id) == 3001
