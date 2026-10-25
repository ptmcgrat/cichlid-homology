# Changelog

## 0.6.0
- Goal B labels revised after the 0.5.0 run (`profile.py`): non-cichlid teleosts are one group (TreeFam has only 2 to 4
  species per level); a species carries an extra copy at 2 or more copies (cichlid expansions with 3+ now count);
  states `kept` (share >= 0.5, the valley in both methods' histograms), `lost` (< 0.3), `partial`. Labels: `single`,
  `retained_pair`, `cichlid_specific_duplication`, `cichlid_lost_one_copy`, `partial_retention`, `multi_copy`.
  New columns: `b_pattern` (one character per level), `b_other_*`, `b_cichlid_share2p`, and the B2 counts
  `b_cichlid_with2p`, `b_cichlid_extra_species`. Options `--kept` and `--lost`.
- `expanded_family` flag (10 or more M. zebra genes in the family); it replaces only the ambiguous labels
  (`partial_retention`, `multi_copy`); the original is in `b_label_raw`.
- `symbols.py`: the closest human and zebrafish symbols are compared with the gene's own NCBI symbol (genes table).
  Grades exact, stem (a/b suffix), none, unnamed, no_reference_symbol. The summary reports rates by relationship and
  adjudicates the genes where the two methods name different reference genes.
- Outputs default to `~/share/<name>_<version>_<timestamp>` (`pilot`, `survey`, `survey2`); `--out` overrides.
- 64 offline tests.

## 0.5.0
- Pin `cichlid-coordinator>=0.14,<0.15` (tested against 0.14.1: `cichlid-data check --release v1 --part Homology`
  gives 9 of 9; 32 offline tests passed on the 0.3.0 tree under 0.14.1).
- `docs/homology_repo_spec.md` is now the Dispatch 7 sheet (it was the Dispatch 6 copy).
- `names.py`: uses the coordinator's `ReferenceNames`. Each closest human and zebrafish gene gets its symbol,
  official ID (HGNC / ZFIN), name status and `name_confidence` (`low` for replaced, legacy_symbol, unmapped or
  missing names; the call is kept). All reference orthologs get capped ID and symbol lists.
- The pilot compares the reference gene each method names by official ID (reliable names only) and writes
  `pilot_identity_differences.tsv`; the summary has `identity_of_the_reference_gene` and `reference_names`.
- `probe` also prints the shape of a `ReferenceNames` record. The field names of that record are not in the
  sheet, so `names.info_of` accepts attribute- and dict-style records; check `probe` on first use.
- Kept as they were: the 0.4.0 clade copy-number profile (goal B1) and the older node-based dating columns.
- 52 offline tests.

## 0.4.0
- Goal B by copy-number profile (`profile.py`): the smallest clade containing the gene and an orthologous
  non-teleost anchor (human, else another non-teleost jawed vertebrate); copies per teleost species in it,
  grouped by ladder level; provisional label `single`, `retained_pair`, `cichlid_specific_duplication`,
  `cichlid_lost_one_copy`, `mixed`, `insufficient_species`, `no_orthologous_anchor`. Columns `b_*`.
- Pilot summary has `goal_B` (label counts, by family size, by human relationship, histograms of the share of
  species with two copies) and the agreement of the goal B label between TreeFam and Compara.
- The node-based dating columns stay in the output but are no longer the goal B answer.
- 45 offline tests.

## 0.3.1
- `explain` prints the ladder level of each ancestor and the species that set it, so a date can be judged by hand.
- Alternative levels to compare with the sheet's rule: `paralog_level_overlap` and `dup_level_overlap` (level of the
  species held by both sides only), alongside the all-species level.
- `dup_ref_levels_sides`: teleost levels on each side counting reference species only (our own placed genes fill the
  levels by construction, so the 0.3.0 count was saturated).
- `family_mzeb_genes` column; summary breaks relationships down by family size and reports agreement for small
  families separately.
- `docs/decisions.md`: the open decisions with the evidence so far.
- 36 offline tests.

## 0.3.0 (changes made after reading the first pilot run, 0.2.2)
- Pilot samples genes, not families, so the same genes are called in TreeFam and Compara (the 0.2.2 samples
  shared only 27 genes). Extra genes from de novo families are added per method (`--n-denovo`).
- Rooting: fewest duplications, then fewest losses (0.2.2 left most families with several equal roots).
  `n_tied_roots` is now after the loss tie-break; `n_tied_dup_only` keeps the old count.
- Outgroup (gar, reedfish) copies are counted in the nearest clade outside the duplication, not in the whole tree.
- Teleost pattern: adds the fraction of all teleost species with exactly one copy, and a coverage figure.
- Summary reports relationships with and without expanded families (n_ref or n_focal at least 10).
- `cichlid-homology explain GENE --part TreeFam|Compara` prints the evidence behind a call.
- Review list gains the closest gene and the tie counts of each method.
- 32 offline tests.

## 0.2.2
- `master.py`: read the master tree from the coordinator's `Node` object (children/name), which is what
  `load_master_tree()` returns; the 0.2.0 adapter only accepted Newick text and failed on it. Two tests added.

## 0.2.1
- Documentation only: `docs/homology_repo_spec.md` replaced with the Dispatch 6 sheet (the copy in 0.1.0 to 0.2.0
  was the v0.12.0 sheet); `docs/NOTES.md` corrected and updated with the Dispatch 6 changes. No code change.

## 0.2.0 (pilot code; not yet run on real data)
- `tree.py`: iterative Newick parser, pruning (no recursion limit on 1,500-leaf trees).
- `speciestree.py`, `rooting.py`: root a tree by reconciliation with the master species tree (fewest duplications,
  all rootings scored in one pass; ties reported; up to k tied rootings returned). Never midpoint.
- `orthology.py`: species-overlap duplication/speciation (pairwise between children, so multifurcations work),
  closest homolog and relationship (`1:1`, `1:many`, `many:1`, `many:many`, `paralog`, `none`), teleost pattern.
- `dating.py`: ladder level of a duplication, retained teleost ohnolog versus recent duplication, place-mode
  duplications flagged unresolved, the disputed cyclostome level never used as a date.
- `consensus.py`: agree / majority / split / single / none.
- `pilot.py`, `cichlid-homology pilot`, `cichlid-homology probe`; `master.py` adapter for the master tree.
- `docs/call_table_schema.md` (draft), `docs/pilot_plan.md`; corrected `docs/NOTES.md`.
- 25 offline tests (constructed trees for every case in the sheet).

## 0.1.1
- Pin cichlid-coordinator >=0.13,<0.14; shadow rule; `survey2`.

## 0.1.0
- Scaffold, `read.py`, `survey`.
