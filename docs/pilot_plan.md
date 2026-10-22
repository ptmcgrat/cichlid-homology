# Pilot (0.2.0): what to run and what it decides

    cichlid-homology probe                       # shapes of the coordinator objects; send the output if anything fails
    cichlid-homology pilot --species Mzeb --n 300 --out ~/pilot

Writes `pilot_genes.TreeFam.tsv`, `pilot_genes.Compara.tsv`, `pilot_summary.json`, `pilot_review.tsv`
(scratch; nothing is pushed). Thresholds are starting values, not decisions.

| To settle | Starting value | Decided by |
|---|---|---|
| species-overlap threshold `t` | 0 (any shared species is a duplication) | `threshold_changes` in the summary at t = 0.1 and 0.25: how many relationships change |
| "mostly one-to-one" | 0.8 of teleost species with exactly one copy | the deciles of `tel_frac_one_of_present` for human 1:1 genes, and the shares at 0.7, 0.8, 0.9 |
| "several teleost levels" for a retained ohnolog | 2 levels on each side | the class counts for genes with a human 1:2, and a hand check of a few |
| tied rootings | check 3, report `n_tied_roots` | `stable_across_roots` |

Rooting: fewest duplications under reconciliation with the master tree; never midpoint.
Shadowed reference leaves are dropped before rooting and before every overlap.
