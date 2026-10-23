# Working notes

## Decided by the coordinator (v0.13.0)
- Leaf to species: `LeafSpecies.load(layout, part)`; `species(leaf)`, `of_leaves`, `is_project_gene`, `shadowed()`.
- Shadowed species: TreeFam Gacu, Onil; Compara Bspl, Gacu; OrthoFinder none. Rule: drop their reference leaves
  before species overlap; placed genes stand for the species. Recorded in `species.SHADOW_RULE` and run.json.
- Pundamilia nyererei and Amphilophus citrinellus are reference-only, not shadowed.
- Mouse out of scope for release 1. Human and zebrafish only. `none` = no such gene in this method's family.
- Rooting: reconcile with the master tree, no midpoint rooting, de novo trees included.
- De novo families, release 1: TreeFam 402 of 9,090; Compara 294 of 15,045.
- `low_family_coverage` means a different coverage in each method; compare only within a method.
- gene_status.tsv: read `up-to-date` and `aligned` as one outcome until project 1 normalises it.

## Corrections to my own earlier statements
- `docs/homology_repo_spec.md` in 0.1.0 to 0.2.0 was the v0.12.0 sheet ("Written against cichlid-coordinator
  v0.12.0", pin `>=0.12,<0.13`, no mouse scope line, no shadowed-species table). I wrongly said in the 0.2.0
  notes that it was current. 0.2.1 replaces it with the sheet from Dispatch 6 (2026-10-22, coordinator v0.13.0).
- What a gene with `no_informative_sites` looks like is documented in the sheet (section 4) and in project 1's
  OUTPUTS.md: `mode` `place`, blank `lwr` and `n_placements`, no leaf in `tree.nwk`, no row in `placements.tsv`,
  an all-gaps row in `combined.aln` and an all-`-` row in `residue_map.tsv`; no relationship from that method.

## Dispatch 7 (2026-10-22, coordinator v0.14.1) in force
- Aims A1, A2, B1, B2 as in the sheet's section 1; equal standing for the three methods on A1 and A2; B1 from
  Compara and OrthoFinder with TreeFam only as a flagged fallback; B2 with OrthoFinder as the authority.
- `ReferenceNames` gives each TreeFam and Compara reference protein a symbol and official ID; a call resting on a
  `replaced`, `legacy_symbol` or `unmapped` protein is lower confidence. OrthoFinder's HGNC and ZFIN IDs come later.
- Paralog dating is not a goal; the benchmark is optional; mouse out of scope; no comparison with 2018.
- Pin `>=0.14,<0.15`.

## Dispatch 6 (2026-10-22) in force
- One pack replaces every older copy. OUTPUTS.md is the copy from project 1's repository (updated 2026-10-21).
- Fusion columns: the four `fusion_partner*` columns are filled if and only if `possible_fusion` is set.
- `gene_status.tsv` says `aligned`, `skipped` or `NoHit`; `up-to-date` is gone (release 1 was regenerated).
- OrthoFinder: family trees expected to be unrooted IQ-TREE trees, branch lengths, no support values, one per
  orthogroup; which tree is used is the owner's decision, not yet made; results not ready.
- Trees sheet section 4: flags, statuses, multifurcations, run-record rule as project 1 answered.
- Pending with the coordinator: human and zebrafish identifier mappings (compare by symbol, say so); schemas of the
  per-gene tables and of my call tables.

## Open
- Identifier mappings: solved for TreeFam and Compara by `ReferenceNames` (Dispatch 7); OrthoFinder's IDs are later.
- Thresholds: see `pilot_plan.md`; they are settled only after the pilot is run on real data.
- Return type of `load_master_tree()` and import path of `LeafSpecies` are assumed; `cichlid-homology probe` checks them.
