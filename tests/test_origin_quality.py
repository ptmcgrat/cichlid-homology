from cichlid_homology.origin import attach_origin, relation
from cichlid_homology.owngenes import gene_order, own_genes, truthy
from cichlid_homology.quality import attach_quality

GENES = {g: {"seqid": s, "start": a, "end": b, "in_proteome": True}
         for g, (s, a, b) in {"A": ("c1", 100, 200), "B": ("c1", 300, 400), "C": ("c1", 9000, 9100),
                              "D": ("c2", 100, 200), "E": ("c1", 500, 600)}.items()}
ORDER = gene_order(GENES)


def row(gene, copies, cat="OneToTwo"):
    return {"gene_id": gene, "a_human_cat": cat, "_a_human_focal_leaves": copies}


def test_gene_order_counts_genes_between_copies():
    assert ORDER["A"] == ("c1", 0) and ORDER["B"] == ("c1", 1) and ORDER["C"] == ("c1", 3) and ORDER["D"] == ("c2", 0)
    assert relation("A", "B", GENES, ORDER, 10)[0] == "tandem_adjacent"
    assert relation("A", "C", GENES, ORDER, 10)[:2] == ("tandem_near", 2)
    assert relation("A", "C", GENES, ORDER, 1)[0] == "same_seqid_distant"
    assert relation("A", "D", GENES, ORDER, 10)[0] == "different_seqid"
    assert relation("A", "Z", GENES, ORDER, 10) is None


def test_origin_classes_and_single_copy():
    rows = [row("A", ["A", "B"]), row("A", ["A", "D"]), row("A", ["A"], "OneToOne"),
            row("A", ["A", "Z"]), {"gene_id": "A", "a_human_cat": "None", "a_zebrafish_cat": "None"}]
    attach_origin(rows, GENES, ORDER, 10)
    assert [r["orig_class"] for r in rows] == ["tandem_adjacent", "different_seqid", "single_copy", "no_position", "single_copy"]
    assert rows[0]["orig_nearest_copy"] == "B" and rows[0]["orig_gap_genes"] == 0


def test_the_nearest_copy_wins():
    r = row("A", ["A", "D", "C", "B"])
    attach_origin([r], GENES, ORDER, 10)
    assert r["orig_class"] == "tandem_adjacent" and r["orig_nearest_copy"] == "B"
    assert r["orig_copies_same_seqid"] == 2 and r["orig_copies_n"] == 4          # B and C share the sequence


def test_zebrafish_clade_is_used_when_there_is_no_human_gene():
    r = {"gene_id": "A", "a_human_cat": "None", "a_zebrafish_cat": "OneToTwo", "_a_zebrafish_focal_leaves": ["A", "B"]}
    attach_origin([r], GENES, ORDER, 10)
    assert r["orig_basis"] == "zebrafish" and r["orig_class"] == "tandem_adjacent"


def test_quality_flags():
    genes = {"A": {"partial": True, "low_quality": False, "internal_stops": 0, "protein_length": 50, "description": "x"},
             "B": {"partial": False, "low_quality": True, "internal_stops": 2, "protein_length": 400, "description": "y"}}
    rows = [{"gene_id": "A", "lwr": "0.2", "assign_flags": "possible_fusion,low_family_coverage", "b_label": "cichlid_specific_duplication", "b_cichlid_with2p": 1},
            {"gene_id": "B", "lwr": "0.9", "assign_flags": "", "b_label": "single", "b_cichlid_with2p": 0},
            {"gene_id": "Z", "lwr": "", "assign_flags": "", "b_label": "single", "b_cichlid_with2p": 0}]
    attach_quality(rows, genes)
    assert rows[0]["risk_flags"] == "low_placement_support;fusion_candidate;partial_gene_model;one_species_only"
    assert rows[1]["risk_flags"] == "low_quality_model;internal_stops" and rows[1]["n_risk_flags"] == 2
    assert rows[2]["n_risk_flags"] == 0


def test_own_genes_reads_the_table(tmp_path):
    from types import SimpleNamespace
    p = tmp_path / "g.tsv"
    p.write_text("gene_id\tsymbol\tseqid\tstart\tend\tstrand\tpartial\tlow_quality\tinternal_stops\tprotein_length\t"
                 "description\tchosen_protein\tmitochondrial\n"
                 "X_1\ta\tc1\t10\t20\t+\tTrue\tFalse\t0\t100\td\tp1\tFalse\n"
                 "X_2\tb\tMT\t10\t20\t+\tFalse\tFalse\t0\t100\td\tp2\tTrue\n")
    g = own_genes(SimpleNamespace(genes_table=lambda prefix: SimpleNamespace(local=p)), "X")
    assert g["X_1"]["partial"] is True and g["X_1"]["in_proteome"] is True and g["X_2"]["in_proteome"] is False
    assert truthy("True") and not truthy("False")
