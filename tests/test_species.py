from cichlid_homology.species import drop_shadowed, leaf_species_map
from cichlid_homology.survey2 import placement_support, fusion_columns


class FakeLS:
    m = {"Mzeb_1": "maylandia_zebra", "Gacu_9": "gasterosteus_aculeatus",
         "ENSGACP1": "gasterosteus_aculeatus", "ENSDARP1": "danio_rerio",
         "ENSPNYP1": "pundamilia_nyererei"}

    def species(self, leaf):
        return self.m.get(leaf)

    def is_project_gene(self, leaf):
        return leaf.split("_", 1)[0] in {"Mzeb", "Gacu"} and "_" in leaf


def test_shadowed_reference_leaf_dropped_project_gene_kept():
    ls = FakeLS()
    sp = leaf_species_map(ls, ["Mzeb_1", "Gacu_9", "ENSGACP1", "ENSDARP1", "ENSPNYP1", "mystery"])
    kept, dropped = drop_shadowed(sp, {"gasterosteus_aculeatus"}, ls.is_project_gene)
    assert dropped == ["ENSGACP1"]
    assert "Gacu_9" in kept and "ENSDARP1" in kept
    assert "ENSPNYP1" in kept          # reference species that is not a project species stays
    assert kept["mystery"] is None      # unknown species is kept, flagged by the caller


def test_placement_support_only_counts_the_species():
    rows = [
        {"gene_id": "Mzeb_1", "family": "F", "mode": "place", "status": "placed", "lwr": "0.4", "n_placements": "3"},
        {"gene_id": "Mzeb_2", "family": "F", "mode": "place", "status": "placed", "lwr": "0.95", "n_placements": "1"},
        {"gene_id": "Mzeb_3", "family": "G", "mode": "denovo", "status": "in_denovo_tree", "lwr": "", "n_placements": ""},
        {"gene_id": "Acal_1", "family": "F", "mode": "place", "status": "placed", "lwr": "0.1", "n_placements": "9"},
    ]
    o = placement_support(rows, "Mzeb")
    assert o["genes"] == 3 and o["placed_n"] == 2
    assert o["lwr_lt_0.5"] == 0.5 and o["lwr_lt_0.9"] == 0.5
    assert o["mode"] == {"place": 2, "denovo": 1}


def test_fusion_columns_counts():
    base = {"gene_id": "Mzeb_1", "best_family": "F", "fusion_partner": "", "fusion_partner_score": "",
            "fusion_partner_region": "", "fusion_partner_family_coverage": "", "flags": ""}
    rows = [dict(base, flags="possible_fusion", fusion_partner="G"),
            dict(base, gene_id="Mzeb_2", flags="possible_fusion"),
            dict(base, gene_id="Mzeb_3", fusion_partner="H")]
    o = fusion_columns(rows)
    assert (o["flagged"], o["partner_filled"], o["flagged_and_filled"],
            o["filled_not_flagged"], o["flagged_not_filled"]) == (2, 2, 1, 1, 1)
