from cichlid_homology.orthology import RTree
from cichlid_homology.speciestree import SpeciesTree
from cichlid_homology.tree import leaves, parse_newick

SPECIES_TREE = ("(bflo,(cmil,((homo_sapiens,mmus),((locu,ecal),(aang,(danio_rerio,((trub,gacu),(olat,(maylandia_zebra,acal)))))))));")
PREFIX2SP = {"Hsap": "homo_sapiens", "Mmus": "mmus", "Locu": "locu", "Ecal": "ecal", "Aang": "aang",
             "Drer": "danio_rerio", "Trub": "trub", "Gacu": "gacu", "Olat": "olat", "Mzeb": "maylandia_zebra",
             "Acal": "acal", "Bflo": "bflo", "Cmil": "cmil"}
LEVEL = {"bflo": "outside_vertebrates", "cmil": "gnathostome", "homo_sapiens": "gnathostome",
         "mmus": "gnathostome", "locu": "gnathostome", "ecal": "gnathostome", "aang": "teleost",
         "danio_rerio": "teleost", "trub": "percomorph", "gacu": "percomorph", "olat": "ovalentaria",
         "maylandia_zebra": "cichlid", "acal": "cichlid"}
ORDER = ["outside_vertebrates", "cyclostome", "gnathostome", "teleost", "percomorph",
         "ovalentaria", "cichlid"]
LINEAGE = {"locu": "non-teleost actinopterygian", "ecal": "non-teleost actinopterygian"}
PREFIXES = {"Mzeb", "Acal", "Gacu", "Olat"}


class FakeClades:
    def level(self, s):
        return LEVEL[s]

    def common_level(self, names):
        return min((LEVEL[n] for n in names), key=ORDER.index)

    def lineage(self, s):
        return LINEAGE.get(s)


def st():
    return SpeciesTree(SPECIES_TREE)


def tree(newick, t=None):
    root = parse_newick(newick)
    for x in leaves(root):
        x.sp = PREFIX2SP[x.name.split("_")[0]]
    return root


def rtree(newick):
    return RTree(tree(newick))


def leaf(rt, name):
    return next(x for x in rt.leaves if x.name == name)
