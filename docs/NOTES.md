# Working notes

## Decided by the coordinator (v0.13.0)
- Leaf to species: `LeafSpecies.load(layout, part)`; `species(leaf)`, `of_leaves`, `is_project_gene`,
  `shadowed()`.
- Shadowed species (reference leaves that are also project species): TreeFam: Gacu, Onil. Compara:
  Bspl, Gacu. OrthoFinder: none. Rule: drop their reference leaves before species overlap; placed
  genes stand for the species. Record in run.json. The owner may overrule.
- Pundamilia nyererei and Amphilophus citrinellus are reference-only species, not shadowed.
- Mouse is out of scope for release 1 (an outgroup leaf only). Human and zebrafish only.
- Rooting: reconcile with the master tree, no midpoint rooting, de novo trees included.
- De novo families, release 1: TreeFam 402 of 9,090; Compara 294 of 15,045.

## Open
- Identifier mappings: TreeFam's ENS2allIDs understood; Compara's current Ensembl IDs still need a
  source. Until then compare by gene symbol and say so.
- Project 1 to answer: meaning of `up-to-date` in gene_status.tsv; definition of
  `low_family_coverage`; whether fusion_partner is filled whenever possible_fusion is set.
- Thresholds to settle in the pilot: species overlap, "mostly one-to-one".
- `docs/homology_repo_spec.md` is the pre-v0.13.0 copy: replace it with the corrected sheet.
