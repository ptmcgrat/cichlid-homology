# Decisions to settle (pilot 0.3.0 data; Mzeb, 600 genes sampled in both methods + 100 de novo per method)

Status: OPEN until the owner fills "Decision". Settled since by Dispatch 7: equal standing on A1 and A2 (row 8),
paralog dating is not a goal (rows 5 and 10 partly), official IDs are the common key (row 7). Evidence is descriptive; there is no ground truth yet.

| # | Decision | Evidence so far | How to settle | Proposal | Decision |
|---|---|---|---|---|---|
| 1 | Species-overlap threshold t | t 0 to 0.1 changes human relationship for 3.5% (TreeFam) / 1.3% (Compara) of genes, zebrafish 0.1% / 0.1%; t 0.25: 5.7% / 3.6% human, 1.0% / 5.1% zebrafish | Benchmark (symbols) once mappings exist; until then look at which calls flip | 0 | |
| 2 | "Mostly one-to-one" (human 1:1 genes) | fraction of ALL teleost species with exactly one copy: >=0.7 96%/96%, >=0.8 91%/88%, >=0.9 73%/77% (TreeFam/Compara); 23% / 58% have some species with 2+ copies | A convention; look for a natural break in the distribution | >=0.8 of all teleost species, coverage >=0.5, and report copies-2+ and absent counts | |
| 3 | Retained-ohnolog "several levels" | weaker side has 4 of 4 levels for 80% (TreeFam) / 97% (Compara), >=3 for 96% / 99%: the test is saturated because our own placed species fill the levels | Count reference species only (columns `dup_ref_levels_sides`) | count reference-only levels; level `teleost` stays the main test | |
| 4 | Expanded families (>=10 Mzeb genes) | 125 TreeFam / 66 Compara families hold 73% / 61% of the genes; paralog share rises from 2% (1 gene) to 40% (10+) in TreeFam | Owner policy | call as normal, flag `expanded_family`, report statistics by family size | |
| 5 | Dating a duplication or paralog | TreeFam human paralogs dated outside_vertebrates: 28% (Compara 13%) | Hand check of ~10 genes with `explain`; compare `*_level` with `*_level_overlap` | decide after the hand check | |
| 6 | Place-mode duplications among our genes | of genes with a duplication class, unresolved in place mode: human 42% (TreeFam) / 42% (Compara), zebrafish 86% / 67% | Owner policy; OrthoFinder is the authority | report "unresolved", date from OrthoFinder when it exists | |
| 7 | Ties for the closest gene | closest_tie on human calls: 27% TreeFam, 3.6% Compara; zebrafish 41% / 36% | Symbol agreement with the gene's own NCBI symbol, once mappings exist | report all tied genes; pick by distance; flag | |
| 8 | Which method wins a split | agree when both answer: human 70% (79% small families), zebrafish 63% (77%) | Benchmark; or majority once OrthoFinder is in | no preference yet; keep every call, review list | |
| 9 | Low placement support | lwr<0.5: paralog 47% vs 34% overall (TreeFam), 25% vs 17% (Compara) | Owner policy | mark `low_support`; never `high` confidence | |
| 10 | Disputed cyclostome level | 190 + 43 TreeFam, 45 + 79 Compara calls hit it | Settled by the sheet | report `cyclostome_disputed`, never a date | settled |


## Goal B (added in 0.4.0; replaces rows 2, 3 and 5 for goal B)
Counts per teleost species in the smallest clade holding the gene and an orthologous non-teleost anchor.
Expected: clusters at one copy (duplicate lost) and two (both kept); the middle ground is read from which ladder
levels carry two copies (cichlids only = recent). Thresholds (`frac` 0.7, at least 3 non-cichlid species) are
provisional: set them from the histograms in `pilot_summary.json` (`goal_B`).


## Settled by the 0.5.0 pilot run (Mzeb, 600 genes + 100 de novo per method)
- Goal B threshold: the share of non-cichlid species with 2 copies has two clusters and a valley near 0.5 in both
  methods (TreeFam 0.4 to 0.5: 28 genes; Compara 0.5 to 0.6: 18 genes). 0.5 is used for "kept", 0.3 for "lost".
- Reference gene identity, TreeFam against Compara: 98.5% (human) and 98.3% (zebrafish) the same gene when each
  method has exactly one reference ortholog; 54% and 56% when either has several (even counting any shared ID).
- Both methods give the same goal B label for 89% of the 1,270 genes labelled in both; only 23 conflicts do not
  involve `mixed`.
- Open: how many genes the NCBI-symbol check can test (many cichlid genes have only LOC names).


## 0.7.0 pilot: A1 and A2 as counts (owner's proposal; not yet agreed with the coordinator)
Replaces the overlap-based relationship and the `paralog` category. The pilot's report shows whether the two methods
agree on the counts and on the best reference gene. To decide after reading it: keep the counts, how to treat very
large clades (expanded families), and whether to tell the coordinator that the sheet's `paralog` category is replaced.


## 0.8.0: whole-genome run, groups, origin of extra copies
- Category groups and `--all` as above. Tandem threshold: at most 10 genes between two copies (`--tandem-gap`).
- Open: what counts as a real cichlid-specific duplication (the 0.7.0 pilot found 4 to 12% of extra-copy genes had
  extra copies only in cichlids); how to use the risk flags; which copy is the older one needs OrthoFinder.


## 0.9.0: parsimony instead of share thresholds for lost duplicate / retained ohnolog
- Relative cost of a gain and a loss: 1:1 by default; a call must hold at 1:2 and 2:1 or it is `ambiguous`.
- Zero counts: scored both as losses and as missing data; the report says how often they differ. OPEN: which to use.
- Next, if the scoring works: scan every internal node of every tree for clades with these patterns.


## 0.9.1: calibrating the parsimony costs
- The costs of a gain and a loss are priors. Counts cannot separate "a pair at the teleost root with early losses"
  from "a new duplicate later", so no cost setting resolves every gene; the grid shows how far costs go.
- Next: score the grid against the published pairs (OHNOLOGS v2, zebrafish; medaka, stickleback and Tetraodon once
  their genes are mapped) and pick the setting that agrees best; use gene-level evidence (positions, published
  pairs) for the cases that remain.


## 0.10.0: calibration result and the tiers
- Calibration (0.9.1 grid, Compara, 26,608 genes; positives: 2,637 genes in observed cichlid pairs or in published
  pairs of two or more percomorphs; negatives: 8,192 genes with one human and one zebrafish gene and no published
  pair). Recall / false positives: share-based thresholds (pair at the teleost root) 84.0% / 0.8%; parsimony gain 2
  loss 1: 82.4% / 5.8%; gain 1: 74.2% / 3.0%; robust unit costs: 63.6% / 1.3%; both calls: 79.0% / 0.7%; either:
  87.4% / 5.9%. Higher gain costs overshoot (gain 6-8 with a free teleost duplication: 40-44% false positives).
- The share cutoff (0.5) sits on a flat optimum: cutoffs 0.3 to 0.6 score within 1.2 points of each other.
- A free teleost duplication did not improve the overall score and is not used; the option stays in the grid.
- The percomorph-only pattern (all early teleosts single, percomorphs double) is the same under "a pair at the
  teleost root with early losses" and "a new duplicate in the percomorphs"; no cost setting separates them (a
  gain would have to cost about eight times a loss). Independent evidence (published pairs) is used for these.
- Positions alone are weak evidence for a pair (45 to 58% of genes with two copies that are not pairs are also on
  different sequences), so they give `possible`, not `probable`.
- The published lists are incomplete (about 40% of zebrafish genes named with an a or b are in them): a gene in a
  pair is evidence for one; a gene in no pair is not evidence against.
- For the pipeline the OHNOLOGS files belong in the coordinator's inputs; until then they are given on the command line.


## 0.10.1: simplified categories
- Owner decisions: a zebrafish-only loss is not its own category; pair-like calls are one group, "Both copies kept
  (best guess)", with the confidence kept in `t_tier`; everything not classified is one group, "Unclassified"; a category
  is added for duplicates after the teleost root and before the cichlids that are not ohnologs.
- The later-duplicate-before-cichlids rule needs tandem copies and no published pair because dispersed copies
  (67% of the 909 candidate genes with two or more copies) look like a teleost-duplication pair kept in the percomorphs.


## 0.11.0: tree-based dating
- Counts cannot separate a teleost pair lost in the early teleosts from a new duplicate in the percomorphs, or
  reciprocal losses from a lost duplicate; the tree shape can (a duplication node whose species lie at the teleost
  level against one inside the percomorphs; two lineages with no species in common). Reconciliation with the master
  species tree is used, not species overlap, so reciprocal losses show as a pair.
- The tree-based call is added beside the counts and parsimony, not instead of them; the tiers are not changed yet.
  First look on mzebra: how the tree call falls across the B1 categories, how often a pair has early teleosts on both
  sides, and the overlap of the disputed genes.
- When OrthoFinder arrives its trees are full trees of our own genes, so the same code also dates duplications
  among project genes (B2 carriers).


## 0.11.1: what the first tree-dating run showed
- On Compara (place mode) the tree call is weaker than counts on its own: among dated genes, thresholds recall 89.6% of
  published-pair positives at 0.5% false positives; a tree pair with overlap of 0.5 or more recalls 63.7% at 0.4%. Reconciliation
  alone ("any pair") called 99.8% of positives and 99.7% of negatives a pair, so shared species are required.
- Its value is for the genes where counts disagree: there a tree pair is the best independent signal (likelihood ratio 2.5).
  "Tree OR thresholds" scored J = 85.0 against 83.2 for thresholds alone.
- Most of the information is in the reference species. The real gain is expected on OrthoFinder, where the tree holds our
  own genes (duplications among project genes, carriers for B2) and all ten cichlids.
- Check on the 0.11.0 run: `t_pairs_*` were empty, so the published pairs were not loaded (the run was probably made without
  `--ohnologs`, or the folder did not hold the files); the categories of that run lack the probable tier.


## 0.11.3: the tree call after the clade fix
- Compara, 0.11.2 run (benchmark: 2,637 published-pair positives, 8,192 negatives): tree pair 56.5% / 0.3% (J 56.2);
  thresholds 84.0% / 0.8% (J 83.2); thresholds or tree pair 89.3% / 0.8% (J 88.5); either count call or tree pair 91.4% / 5.9%.
  The 140 positives the thresholds miss and the tree finds are in families the thresholds call many-copy (138 of 140).
- On genes where the counts disagree, a tree pair was in a published pair 47.3% of the time against 18.7% for genes with no
  evidence (likelihood ratio 2.5; copies on different sequences 1.4).
- Overlap histogram: a large mass at 0 (the gene tree disagrees with the species tree, or reciprocal losses), a second
  mode at 0.8 to 1.0, a valley at 0.3 to 0.5; 0.5 is kept.
- Many-copy families: a tree pair is evidence of a pair (39.0% in a published pair, against 8.8% without shared species and
  43.6% for ordinary pairs); in expanded families it is weak (14.9% against 6.0%) and they stay unclassified.


## 0.12.0: where the pair is in the tree
- Compara, 0.11.3 run, 2,637 published-pair positives and 8,192 negatives: the pair at the top of the gene's teleost clade only:
  56.5% / 0.3%; at the deepest supported duplication: 82.4% / 0.5%; thresholds 84.0% / 0.8%; thresholds or the tree 91.0% / 0.9%.
- The tree-based `later_duplicate_before_cichlids` calls of 0.11.3 were mostly the teleost pair below an incongruent top; the
  genuine ones (percomorph or ovalentarian level) were 53 genes in Compara, because in a placement tree duplications among
  project genes are not visible. The tiers' own category (parsimony events, tandem copies, no published pair) found 132.
- To check on the next run: how many genes are now `pair` and whether they look like pairs (different sequences, published
  pairs), how many `later_duplicate_before_cichlids` remain and whether they look like local duplicates (tandem copies, no
  published pair), and the `shared_early` / `early_within` split.


## 0.13.0: confidence by agreement of three calls
- Owner decision: tiers from the votes of thresholds, parsimony and the tree; core all three agree, probable two of three, no
  majority is Disputed; the published pairs are an outside check, not part of the calls.
- The tree is the call most often outvoted when exactly two agree (1,191 of 1,963 on Compara, against parsimony 599 and
  thresholds 173): it misses pairs whose lineages share no species (a disjoint root: 12.5% of the positives), where only
  counts see them.
- Expected on the next run: the same categories as the preview above; a fourth voice (OrthoFinder) would be a second tree.
