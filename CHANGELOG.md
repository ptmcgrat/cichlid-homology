# Changelog

## 0.11.2 (published-pair files: both name styles, and a loud failure)
- The OHNOLOGS loader missed the files on mzebra (`drerio.Pairs.Strict.3R.txt`) because it expected underscores
  (`drerio_Pairs_Strict_3R.txt`): nothing was loaded and nothing was said, so the 0.11.0 run had no published-pair evidence.
  `_` or `.` may now separate the parts.
- `--ohnologs` with a folder that holds no recognisable pair file, or that does not exist, now prints a warning, records
  `ohnologs_error` in the summary and puts a warning at the top of the report. On success the run prints which species were
  loaded and the report says so. 140 offline tests.

## 0.11.1 (tree dating: shared species make a pair; the gene's own teleost clade)
- First run on Compara (0.11.0) found two faults. (1) A duplication by reconciliation is not evidence of a pair: the gene
  tree often disagrees with the master species tree (a polytomy among the early teleosts), and 98% of genes that counts and
  parsimony agree are lost duplicates came out as pairs, 93% of those with lineages that share no species. (2) A quarter of
  genes could not be dated because the lowest ancestor of all teleost genes of a large anchor clade also held tetrapod genes
  of other paralogs.
- A pair now needs the two largest lineages to share at least half of their species (`--pair-overlap`, 0.5; the overlap
  histogram has a valley between 0.3 and 0.5). `d_root` is `pair`, `thin` (some shared species, less than half), `disjoint`
  (none: gene-tree incongruence or reciprocal losses; counts cannot tell these apart) or `single`; `d_call` adds `thin_pair`.
- The teleost root is the highest ancestor of the gene inside the anchor clade that holds no non-teleost jawed vertebrate
  (the gene's own teleost clade); the `teleosts_not_monophyletic` status is gone.
- Tiers: a gene tree that shows a pair (`d_call` = pair) is weak independent evidence and gives the `possible` tier. On the
  Compara 0.11.0 table, among the genes where the counts disagree, a tree pair with overlap of 0.5 or more was in a published
  pair 39.6% of the time against 15.6% for genes with no evidence (likelihood ratio 2.5; copies on different sequences 1.4).
- Report: root classes and the overlap histogram. 138 offline tests.

## 0.11.0 (tree-based dating of duplications)
- `treedate.py`: a third line of evidence for B1 and B2, read from the shape of the rooted gene tree. The tree is
  reconciled with the master species tree: a node is a duplication when the lowest common ancestor (in the species
  tree) of its species equals that of one of its children. No share thresholds. For each gene, in its anchor
  clade: `d_root` (`pair` if the teleost genes form two or more lineages, `single` if one), `d_root_level` (a pair at
  the `teleost` level has early teleosts in it; at `percomorph` or lower the tree cannot tell a teleost pair lost in
  the early teleosts from a later duplicate), `d_support` (early teleosts on both sides, one, none), `d_overlap` (the
  species the two largest lineages share as a share of all their species: near 1 for a fully retained pair, 0 for
  reciprocal losses, small for a thin pair or a stray gene), `d_gene_lineage`, `d_focal_lineages`, `d_later_levels`,
  `d_first_later_level`, `d_species_dups`, `d_carriers`/`d_carrier_n` (denovo trees: the cichlid species that hold both
  copies below the deepest later duplication: B2), and `d_call` (pair, lost_duplicate, later_duplicate_in_cichlids,
  later_duplicate_before_cichlids, pair_unresolved_depth).
- Placement trees (TreeFam, Compara): the placed project genes are not used to find duplications (genes on one branch
  cannot be ordered); the structure comes from the reference leaves and each project gene is assigned to the lineage
  it was placed in. Denovo trees (OrthoFinder): every leaf counts.
- Known weakness, tested: one stray gene in the wrong clade can create a duplication that is not there; `d_overlap`
  is small in that case.
- Report: a tree-based dating section with the tree call by B1 category. 135 offline tests; 0.2 ms per gene.

## 0.10.1 (simplified B1 categories, later duplicates before the cichlids)
- `t_category` groups the tiers into six categories: Lost duplicate; Both copies kept (best guess) (R2, R3, probable
  and possible pairs; a zebrafish-only loss is no longer its own category); Later duplicate, in cichlids; Later
  duplicate, before cichlids; Disputed; Unclassified (one group; the reason stays in `t_call`). `t_tier` is the
  confidence (core, probable, possible).
- New: Later duplicate, before cichlids, with `t_when` (teleost, percomorph or ovalentaria level). The B1 definition
  has three outcomes, the third being a duplicate unrelated to the teleost genome duplication, and that includes
  duplicates between the teleost root and the cichlids. A gene qualifies when parsimony places a clade-level gain at
  those levels, the gene has two or more copies, the copies are tandem, and no published pair holds it. The rule is
  strict because counts cannot tell such a duplicate from a teleost-duplication pair kept in the percomorphs (on the
  Compara table of 0.9.1 it finds 127 genes outside the other groups; 45 more already fall in cichlid duplicates
  and 52 in many-copy genes).
- Report: a simplified-category table per method. 123 offline tests.

## 0.10.0 (calibrated parsimony and confidence tiers for B1)
- Parsimony costs are now a gain of 2 and a loss of 1 (`--gain-cost`, `--loss-cost`), chosen in 0.9.1's grid against
  published teleost-duplication pairs (OHNOLOGS v2): recall of known pairs 82% (63.6% with the old robust unit costs),
  false positives 5.8%. `p_robust` now means the label holds at half and twice the gain-to-loss ratio and no longer
  turns the label into `ambiguous`; `p_label_unit` is the label under equal costs, for comparison.
- `ohnologs.py` and `--ohnologs FOLDER`: reads the OHNOLOGS v2 3R pair files (protein-coding, fish genome
  duplication) as independent evidence and adds `t_pairs_relaxed_n`, `t_pairs_strict_n`, `t_pairs_species`. The
  files are given by the owner; nothing is downloaded. Optional: without it the tiers use positions only.
- `tiers.py`: one row per gene in the confidence table for B1: `t_row` (R1 lost early; R2 both kept; R3 both kept,
  zebrafish lost one; R2/3/4 a pair whose depth the counts cannot fix; R5 later cichlid duplicate), `t_call`,
  `t_tier` (core: the share-based and parsimony calls agree; probable: one calls a pair and a published pair
  agrees; possible: one calls a pair and the copies are on different sequences, weak evidence; disputed; unclassified:
  not scored, expanded family or many copies), `t_basis` and `t_positions`. No new cutoffs.
- Report: a tiers section per method. Previewed on the Compara table of 0.9.1 (26,608 genes): 18,575 placed in a
  row: core 15,880 (85%), probable 356, possible 597, disputed 1,742 (9%); core lost duplicates are in a published
  pair for 1.8%, core pairs for 45%, disputed genes for 4.7%.
- 120 offline tests.

## 0.9.1 (calibration grid for the parsimony costs)
- `cichlid-homology pilot ... --parsimony-grid` also labels every gene by parsimony under twelve cost settings:
  a gain costs 1, 2, 3, 4, 6 or 8 times a loss, each with and without a free teleost duplication (one extra copy
  on the branch above the teleost ancestor costs nothing: the genome duplication doubles every gene). New columns
  `pg_g<gain>l1[f]` (label) and `pgm_g<gain>l1[f]` (margin); `f` marks the free-duplication settings.
  Off by default; the existing columns are unchanged.
- Purpose: choose the costs that best reproduce independent pairs (OHNOLOGS v2) instead of fixing them by hand.
- Tests: the free-duplication matrix, the grid columns, the zebrafish tie resolved by a costlier gain or a free
  duplication, and a percomorph-only pair that stays `lost_duplicate` until a free duplication plus a costly gain.
  Cross-checked against a separate implementation on 4,800 random cases (no mismatches).

## 0.9.0 (pilot: lost duplicate or retained ohnolog by parsimony, no share thresholds)
- `parsimony.py`: for the gene's anchor clade (the same clade as the B1 profile), the copy numbers of the
  jawed-vertebrate species are explained by gains and losses on the pruned master species tree (no gene is born
  from nothing; the deepest ancestor has one gene). The teleost ancestor having 1, 2 or 3 copies is compared by
  the fewest events: `lost_duplicate`, `retained_ohnolog`, `multi_copy`. The margin between the best and
  second-best history is the confidence; a tie, or a call that changes when a loss costs 2 gains or the reverse,
  is `ambiguous`. Species with no gene in the clade are scored both as losses and as missing data
  (`p_label`, `p_label_missing`, `p_zero_modes_agree`).
- The events of the best history are reported: gains by ladder level (a later duplicate is a gain below the teleost
  ancestor), losses by species (`p_loss_species`), so a zebrafish loss of a retained ohnolog is explicit.
- Scope: species that are jawed vertebrates in the clade table; the lamprey and invertebrates are left out.
- `summary_p.py` and a report section: label counts, margins, robustness, comparison with the threshold-based B1
  labels, two checks that use nothing from the labels (zebrafish a/b names, positions of the copies), later gains,
  most-lost species, and agreement between TreeFam and Compara.
- Not yet: the scan of every internal node of every tree. This pilot scores the same clades as B1 and tests the
  scoring first.
- 105 offline tests; about 10 ms per unique 60-species clade.

## 0.8.0
- Whole-genome mode: `cichlid-homology pilot --all` calls every gene of the species that has a finished family in a
  method (the 0.7.0 pilot called a biased sample); `--species Mzeb,Acal` runs several species (a folder each).
  Families are streamed one at a time; progress goes to stderr; `--all` skips the old rooting-stability and
  overlap-sensitivity checks.
- Category groups as the owner asked: human OneToOne, OneToTwo, OneToMany, Other (= ManyToOne + ManyToMany), None;
  zebrafish OneToOne, OneToMany (= 1:2 and 1:many), Other, None. The fine category stays in `a_<ref>_cat`; the group
  is `a_<ref>_group`; the report and summaries use the groups.
- Where the extra copies sit (`origin.py`, `owngenes.py`): for genes with 2+ focal copies in the reference clade
  (human clade, else zebrafish), the nearest other copy is tandem_adjacent, tandem_near (at most `--tandem-gap`
  genes between, default 10), same_seqid_distant, or different_seqid, from the genes table positions. Columns
  `orig_*`.
- Artifact-risk flags (`quality.py`): low placement support (lwr < 0.5), fusion candidate, partial or low-quality
  gene model, internal stops, extra copy in one cichlid species only. Columns `risk_flags`, `n_risk_flags`, `gt_*`.
- Candidate tables `pilot_candidates.<method>.tsv` for `cichlid_specific_duplication` genes, and report sections
  for the origin of extra copies and the candidates (carriers, origin, flags, genes in both methods).
- 92 offline tests, including the whole pilot on a fake layout in sample and all-genes mode.

## 0.7.0 (pilot of a simpler A1 and A2; not agreed with the coordinator)
- `clade_counts.py`: for each focal-species gene and each reference species (human, zebrafish), the smallest clade
  holding the gene and a reference gene; counts of reference genes and of focal-species genes in it (other cichlids
  are not counted); category OneToOne, OneToTwo, OneToMany, ManyToOne, ManyToMany, or None; the best reference gene
  is the shortest path (ties flagged). No duplication/speciation rule, so the `paralog` category and the overlap
  threshold are not used for A1 and A2.
- Goal B (`profile.clade_profile_v2`) now uses that human clade (or, with no human gene, the smallest clade with
  another non-teleost jawed vertebrate), with no orthology test. The teleost distribution is summarised for the
  human OneToOne and OneToTwo genes.
- Names, official IDs and the NCBI-symbol check apply to the new best genes (`names.attach_names_a`,
  `symbols.attach_symbol_check_a`).
- `summary_a.py` and `report.py`: overall numbers per method, agreement between TreeFam and Compara (category, the
  reference-gene count, the focal-gene count, the best gene by official ID, adjudication by the gene's own symbol),
  and goal B. The pilot writes `pilot_report.md` and prints it.
- New files: `pilot_category_splits.tsv`, `pilot_best_gene_differences.<ref>.tsv`.
- The old columns (`human_relationship`, paralog levels, dup classes, overlap) are still computed and written, under
  their old names, and are not used by the new summary. They go once the change is agreed.
- 82 offline tests.

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
