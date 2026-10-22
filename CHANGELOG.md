# Changelog

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
