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
