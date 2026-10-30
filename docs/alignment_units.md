# Alignment units (`cichlid-homology alignment-units`, 0.24.0)

A unit is the set of genes to look at together in the protein alignment of a family, for a gene or for a pair of genes of ours.

## The rule

The smallest clade of the rooted family tree that holds the gene (both genes of a pair) and

1. a **human** gene; if the family has none,
2. a **tetrapod** gene (a sarcopterygian that is not the coelacanth); if the family has none,
3. the **whole family**.

The clade is a node of the family tree with our genes in it (`families/<family>/tree.nwk` of the method), rooted like the reference tree where there is one, otherwise by reconciliation with the master species tree
(the same rooting as the clade split and the figure pilot). It is not cut at ancient duplications: for a gene whose clade of the clade table has no human gene because the human gene is across an ancient duplication,
the unit climbs past the duplication to the smallest clade that has one.

## Files

- `alignment_units.tsv`: one row per unit. `unit_id` is `<family>:<H|T|F><k>`: the rule (H, T, F) and the rank by size within the family and rule (1 is the smallest). Also the family, `whole_family`, tree mode, rooting,
  `n_genes`, `n_species`, `n_family_genes`, `share_of_family`, the numbers of human, tetrapod, focal-species, project-species and reference genes, the reference leaves of shadowed species, how many genes and pairs
  asked for the unit, the focal and human gene lists, and the paths of `combined.aln` and `tree.nwk`.
- `alignment_unit_members.tsv.gz`: one row per gene of each unit: `unit_id`, `gene_id` (the leaf name, which is the row name in `combined.aln`), `species`, `kind` (focal, project or reference), `ladder_level`,
  `is_human`, `is_tetrapod`, `shadowed_reference` (the Ensembl copy of a genome that we annotate too: drop it, as the tree calls do).
- `gene_to_unit.tsv`: every gene of the gene table, with its unit, rule, the size of the unit and a status (`ok`, `no family`, `gene not in the family tree`, `family not read`).
- `pair_to_unit.tsv`: every pair, with the unit that holds both genes and a human gene, and the unit of each gene (`unit_a`, `unit_b`). A pair whose genes are in different families has no unit.
- `alignment_units.json` and `.md`: counts by rule, the size of the units, the genes whose unit is the whole family.

## Using a unit

```bash
U=C116_000070:H1; F=<the alignment_file of the unit in alignment_units.tsv>
zcat alignment_unit_members.tsv.gz | awk -F'\t' -v u="$U" '$1==u && $8=="no" {print $2}' > ids.txt     # without the shadowed reference leaves
seqkit grep -f ids.txt "$F" > unit.aln.fasta
```

This only says which rows. The placement methods trim the columns that only the query sequences fill (`--keeplength`), so `combined.aln` alone does not map every residue of a placed gene to a column: that map is
project 4's (the residue step).
