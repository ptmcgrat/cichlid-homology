# Changelog

## 0.22.0 (a teleost duplication is not an ancient one just because the tetrapods sit in one lineage)
- `clade-split`: a duplication is cut as ancient only when each of its two biggest lineages carries at least `--min-outgroup-species` (default 1) species that are not teleosts (tetrapods,
  coelacanth, shark, gar, reedfish). Why: the figure pilot showed that when the tree nests all the non-teleost genes in one lineage (gnb1, stag2: the other lineage is teleost-only), the old rule
  cut a teleost duplication in two, so the ohnologs became separate 1:1 clades and the pair was missed, or an out-paralog was called where there is a teleost pair (wnt7ba, pi4kab, loxl2a). A
  duplication whose tetrapod-free lineage still holds gar or reedfish is older than the gar split (pdap1) and is still cut; an ancient paralogy with outgroups on both sides (stk17, dagl) too.
  `--min-outgroup-species 0` gives the rule of 0.16 to 0.21; 2 is the sensitivity run. `figure-pilot` takes the same option so that the clade ids of the gene table apply.
- Deep order: each of the four known relationships is tested strictly and tolerantly (some branch without the outgroup holds 90% of the genes of the group). The strict test failed in half of the clades
  because one stray teleost gene is enough; the tolerant one never passes fewer clades. The score, the flag `deep_order_inconsistent` and the figure use the tolerant one; the record keeps both.
- The pilot tally counts clades, not genes (the "72% attach to the slower lineage" figure was 126 genes from 9 clades; by clade it was 5 of 9). New `strata` summary: per stratum, the clades with a
  teleost split, a later split and no split.
- A gene placed on the branch that leads to its clade (not inside it) is drawn on the stem with a marker.
- `compare-tables` has a section for the genes named a and b: how many are partners, in the same clade, in the same family or in different families, before and after.
  229 offline tests.

## 0.21.1 (is the deep order of the tree right?)
- Each clade record carries `deep_order`: whether the tree recovers the four relationships we know (shark outside the bony fishes; coelacanth with the tetrapods; reedfish outside gar and the
  teleosts; gar outside the teleosts), tested only where the species are present, with a score. Why: the position of a single gar or reedfish gene is unreliable (short internal branches, long
  tips), so it should not date anything; the share of known relationships a tree recovers is a quality score for its deep nodes, with the coelacanth as the control. On the real irx1 tree it is 2 of 4.
  The figure says which failed (flag `deep_order_inconsistent`), and `figure_index.tsv` has the score.
- The pilot summary has a `tally`: how often each relationship is recovered, the scores, the copies of gar and reedfish in clades that have a split, where the non-teleost genes attach by species,
  and whether they attach to the slower lineage (the one with the shorter stem), which tells whether the position of an outgroup says anything about rate asymmetry between the ohnologs.
  222 offline tests.

## 0.21.0 (the figure pilot on the Compara data)
- New command `figure-pilot`: for a stratified sample of genes (or all of them), one summary record and one figure per clade; records are the product (`clade_records.jsonl`, schema in
  `docs/clade_record_schema.md`), figures are for looking (`figures/`, `figure_index.html` as a contact sheet, `figure_index.tsv` with the triage numbers). Compara only; TreeFam is not used.
- `nhx.py`: reads the raw Compara dump (EMF, plain or gz) and its NHX tags: the bootstrap (`B`) and the duplication flag (`D`, `DD`) that the prepared trees lost. `suppmap.py` puts them on our trees
  by the reference proteins below each node, and roots a placed tree where the reference tree is rooted (`reroot_like`), so the clade cuts and the supports line up.
- `lineages.py`: the lineages of a clade from the teleost genes alone. The highest node whose two biggest children share at least half of their teleost species (and at least two) is the split;
  every non-teleost gene is the ancestral copy wherever the tree attaches it, and where it attaches and how well supported that is are reported. This replaces the rule that wanted both sides to be
  teleost-only, which threw away the real duplication of irx1a and irx1b because a gar and a reedfish gene hang inside one lineage (bootstrap 0 and 25). Lineage letters are ordered by the zebrafish gene
  of each lineage, so they mean the same in every method.
- `figrecord.py` builds the record (collapsed tree, lineages, supports, copies per species, our genes with their lineage and placement weight, triage numbers, flags); `figrender.py` draws it
  (parallel with `--workers`). `cladesplit.split_family` gained `group_leaves` so that the figure finds a clade of the clade table in the placed tree.
- `figure-pilot` needs matplotlib (`pip install cichlid-homology[figures]` or `pip install matplotlib`); it says so before it starts if it is missing.
- Sample strata: pairs that match a published pair, genes named a/b but called single copy, other pairs, 1:1 and 1:2 with human, out-paralogs, later duplicates, expansions, ambiguous partners,
  clades with no human gene, placement weight under 0.5, genes in de novo trees. 220 offline tests.

## 0.20.0 (species profiles of the ohnolog clades)
- `clade-split` now also writes, for every clade that holds a teleost-level duplication node (the highest such node is the pair node), the species profile
  of the pair: for each species the number of genes outside the node (the ancestral copy), in lineage A and in lineage B (`clade_species.<part>.clades.json`,
  counts as an `outside.A.B` string per species in the order of `clade_species.<part>.species.json`), the later duplications inside each lineage that cover
  two or more species (lineage, level, species), and the zebrafish genes of each lineage (to assign a focal-species gene to a lineage). The species table is in
  species-tree order with a group (tetrapod, coelacanth, cartilaginous fish, gar or reedfish, early teleost, percomorph, ovalentarian, cichlid) and a mark for
  the salmonid, carp and cyprinid-specific genome duplications; `clade_species.<part>.tree.json` is the species tree pruned to those species. `--no-species-profiles`
  skips it. Why: to draw which tetrapods and teleosts carry each ohnolog, which species miss one member, and which duplicated one of them. 196 offline tests.

## 0.19.0 (synteny in the gene table, ambiguous clades resolved, chromosome tables)
- First real table (0.18.0, with positions): the 2,504 paired genes on different chromosomes fell on only 111 chromosome pairs, the largest carrying 129 to
  220 pairs; 96% of the pairs sat on a chromosome pair holding at least 3 other pairs (published pairs 98%, pairs with no published support 94%), and 191 of
  the 227 low-confidence pairs on a chromosome pair with 5 or more others: the signature of the teleost genome duplication, which the confidence had not used.
- New column `ev_synteny`: the number of OTHER firm pairs (status paired, high or medium before synteny, different chromosomes) on the pair's chromosome pair.
  The confidence now counts four signals (lineages, chromosomes, zebrafish split, synteny of at least `synteny_min_pairs` = 5); against: same lineage, tandem
  copies, a chromosome pair no other pair shares. High: three or more and none against; medium: two and none against, or three with one against; low: otherwise.
- A possible pair (1:1 clade, two focal genes) on a supported chromosome pair is now `pair retained in few teleosts (synteny-supported)`.
- Ambiguous clades (1:2 clade, three or more focal genes): a gene is paired with the gene in the other lineage on a partner chromosome (a chromosome pair with
  at least 5 firm pairs) when each is the other's only candidate: `paired (resolved by synteny)` (its confidence does not count synteny again); otherwise the
  candidates are listed in the new column `ohnolog_candidates`, and a gene with none is `partner not found on a partner chromosome`. `ohnolog_lineage` is now
  filled for every gene of a 1:2 clade.
- `gene-table --positions` also writes `chromosome_pairs.<species>.tsv` (per chromosome pair: pairs by kind and confidence, shares of each chromosome's pairs,
  mutual top partners) and `chromosomes.<species>.tsv` (genes, genes in pairs, top partner, partner chromosomes with 5 or more pairs). 68 columns. 189 offline tests.

## 0.18.0 (the gene table)
- New read-only command `cichlid-homology gene-table`: one row per focal-species gene, 66 columns in 11 groups (identity; family and clade; closest human gene;
  closest zebrafish gene; teleost trend and relationships; human-side cleanliness; B1 call; ohnolog; later duplicates and out-paralogs; flags; provenance),
  from the pilot gene tables of Compara (and TreeFam, for agreement), the clade table of `clade-split`, the published pairs (optional) and the gene positions
  (`--positions`, through the layout). Writes `gene_table.<species>.tsv`, a JSON summary and `gene_table_columns.md`.
- A gene's clade is the clade of its closest zebrafish protein (else that of its closest human protein if it has teleosts). The relationships (`rel_human`,
  `rel_zebrafish`) are the clade's human or zebrafish gene count against the teleost median, not the focal copy number, which is compared in `copy_vs_trend`.
- Ohnolog: a 1:2 clade with two focal genes is a pair (partner gene, symbol, lineage, matching zebrafish gene, confidence high/medium/low from three signals:
  different lineages, different chromosomes, two different zebrafish genes; the published pairs are an outside check); one focal gene means the partner was lost
  (the other zebrafish genes are listed); three or more are ambiguous until OrthoFinder. A 1:1 clade with two focal genes on different chromosomes is a possible pair,
  otherwise a later duplicate. `ohnologs.load_pair_partners` reads the published partners. 182 offline tests.

## 0.17.0 (teleost trend of each clade)
- The clade table now carries the teleost trend: the median number of genes per teleost species present (`teleost_median_genes`), the shares of
  teleost species with one, two and three or more genes, the median over the early teleosts and over the percomorphs and below, the number of zebrafish genes
  in the clade, and `teleost_trend`: 1:1 (median 1), 1:2 (median 2), 1:3+ (median 3 or more), or `mixed 1-2` / `mixed 2-3` (a half-integer median: the species
  are split). The median needs no cutoff.
- Why: for the 1:1 and 1:2 relationships the owner wants the trend across the teleosts (the median gene number close to 1 or to 2), not the copy number of one
  species. On the 0.16.1 table, with the mean as a stand-in, the 13,184 clades with one human gene and teleosts were bimodal: 73% about 1 (mean up to 1.25),
  8% in the valley (1.25 to 1.75), 18% about 2 (1.75 to 2.5), 1% beyond; in the about-2 clades zebrafish kept both copies in 1,545 of 2,316.
- The report tabulates the trend against the number of zebrafish genes in the clade for clades with one and with two human genes. 168 offline tests.

## 0.16.1 (protein IDs in the clade table)
- First full run of 0.16.0 on Compara: 15,227 families with a human gene split into 20,121 clades; every multi-member family used its own recorded root
  (14,942), none needed reconciliation. Clades with exactly one human gene: 16,421 (13,093 from one-human families, 3,328 from split families, as clean
  as the rest: 69% against 72% with every tetrapod species single-copy). Two-human clades are not clean (64% under half the tetrapod species single-copy).
- The clade table now lists the human and zebrafish protein IDs of each clade (`human_proteins`, `zebrafish_proteins`) so that an M. zebra gene joins to
  its clade through its closest human or zebrafish gene. 166 offline tests.

## 0.16.0 (clades of the reference families, split at their ancient duplications)
- New read-only command `cichlid-homology clade-split --part Compara`: for every reference family with a human gene (`--all-families` for the rest),
  root the tree (the tree's own root when it has two children, else by reconciliation with the master species tree; `--rooting reconcile` always
  reconciles), find the ancient duplications and cut there, and score each clade left: a duplication is ancient when it is a reconciliation
  duplication whose two largest children share at least `--pair-overlap` (0.5) of their species, tetrapods lie below it, and either teleosts do
  too or its species-tree node is the amniote ancestor or older (a paralog lost in every fish is still cut). Duplications younger inside the
  tetrapods (human-specific duplicates) or inside the teleosts only (the teleost genome duplication) are not cut but counted.
- Per clade: human genes; tetrapod species, genes, mean genes per species and the share of tetrapod species with exactly one gene; teleost species,
  mean genes per species and the share of the family's teleost species in the clade; gar and reedfish; supported duplications inside it (tetrapod
  clade, tetrapod species, teleost clade, mixed); root-to-tip distances of the tetrapods from the clade's root (median, spread) and the human
  distance relative to the other tetrapods (with at least three of them). Tetrapods are the sarcopterygians without the coelacanth. No cutoffs.
- Outputs: `clade_split.<part>.families.tsv`, `.groups.tsv`, a JSON summary and a report (clades per family by number of human genes; profiles of the
  clades with exactly one and two human genes, split by whether the family itself has one human gene; examples of families with four). 165 offline tests.

## 0.15.0 (mean genes per species, by number of human genes)
- The family survey now counts tetrapod leaves and species (`n_tetrapod_leaves`, `n_tetrapod_species`: the sarcopterygians without the
  coelacanth) and reports, for families with exactly 1 and exactly 2 human genes, the distribution of the mean number of genes per species
  among the tetrapod species present and among the reference teleost species present (genes / species per family): bins (1, 1 to 1.25,
  1.25 to 1.5, 1.5 to 2, 2 to 3, 3 to 5, over 5), quartiles, the mean of the means, the median number of species present and the families
  with no species of the group. Both sources (`--source built` and `--source reference`). 155 offline tests.

## 0.14.3 (one-member families in the reference survey)
- Full run of 0.14.2 on the coordinator's Compara116 input: 41,193 families (40,860 with a tree, 333 alignment only), 1,207,777 leaves, none
  without a species. 7,715 families came back with no leaves: a one-member family has a tree with no parentheses (`Hsap_1;`) that the leaf
  parser could not see, so single human genes were counted as zero. Now a bare name is read as a one-member tree, and a family whose tree
  yields nothing is read from its alignment (`leaf_source` `alignment_because_tree_empty`).
- The report describes the families with human and teleost genes that no project species was assigned to: audit statuses, median leaves and
  teleost species, those with one or two teleost species, those with exactly one human gene, those with gar or reedfish. 153 offline tests.

## 0.14.2 (reference survey: alignments as well as trees)
- First run of 0.14.1 on the real input: `prepared/trees/` holds alignments (`.aln`, 41,193 families), not only trees; reading them as
  Newick gave zero leaves and every count was zero. The audit file was read correctly: 54,308 families (`ok` 23,864, `tiny_tree` 16,996,
  `no_panel_members` 13,115, `cichlid_only` 333); the 13,115 `no_panel_members` have no file.
- The reference reader now keeps every file of a family and uses the tree when there is one, otherwise the sequence names of the alignment
  (FASTA headers or Stockholm lines, a trailing `/start-end` removed). New columns `leaf_source` (tree, alignment, none) and, in the
  report, files by extension, families with a tree file, families with an alignment only, families with no leaves. 150 offline tests.

## 0.14.1 (family survey of the coordinator's prepared reference trees)
- `cichlid-homology family-survey --part Compara --source reference [--prepared PATH]` reads the coordinator's prepared Compara116 trees
  (`Inputs/Families/Compara116/prepared/trees/`), counts human, mouse, zebrafish, gar and reedfish genes, reference teleost genes and species,
  and joins each family with its `family_audit.tsv` status, with whether the project built it (a folder under Results) and with whether any
  project species (or the focal species) was assigned to it. The report gives the histogram of human genes per reference family, the families
  with human genes and no teleost, the families with human and teleost genes that no project species hit, and the audit statuses.
- The file formats of `prepared/` are not documented in the pack, so the reader is adaptive and reports what it finds: file extensions, the
  family id taken from the file name up to the first dot, the audit columns (`family` and `status`, else the first two), audited families with
  no tree file, leaves with no species, unreadable trees. Bracketed comments (NHX) in a tree are ignored.
- Output file names now carry the source: `family_survey.<part>.<built|reference>.tsv|json|md`. 150 offline tests.

## 0.14.0 (family survey: human genes and teleosts in every family)
- New read-only command `cichlid-homology family-survey --part Compara`: for every finished family of a method it counts human, mouse,
  zebrafish, gar and reedfish genes, reference teleost genes and species, other jawed vertebrates, tetrapod leaves, unknown leaves, and our
  own placed genes, and writes `family_survey.<part>.tsv` (one row per family), a JSON summary and a short report: the histogram of human
  genes per family, and how many families hold human genes and no teleosts (reference teleost leaves only: a placed M. zebra gene in a
  tetrapod tree does not count as a teleost). `--limit N` for a quick test, `--species` marks the families that hold a focal-species gene.
- Why: the ohnolog analysis is simplest where the human protein is a single gene with clear mammalian orthologs. Compara splits a family
  into a fish tree and a tetrapod tree, so many families hold human genes and no teleosts and are never reached from the M. zebra gene table.
  This is the first of the measures for the human-side criterion (docs/decisions.md). 147 offline tests.

## 0.13.0 (categories and confidence from the votes of three calls)
- Run of 0.12.0 on Compara: the tree alone finds 82.4% of the 2,637 published-pair positives at 0.5% of the 8,192 negatives
  (thresholds 84.0% / 0.8%; thresholds or tree 91.0% / 0.9%). The genes called pairs look like pairs (85% on different sequences,
  45% in a published pair, 71% human OneToTwo; the 2,092 non-expanded genes with the pair node below the top of the teleost clade
  83%, 41% and 72%); the later duplicates look local (57 to 88% tandem, 5 to 10% in a published pair).
- The three calls (thresholds, parsimony, tree) form a confidence ladder. Votes for a pair, positives : negatives: 3 of 3
  1,877 : 34; 2 of 3 365 : 27; 1 of 3 203 : 427; none 192 : 7,704.
- `tiers.py` rewritten: each call votes `pair`, `lost` or abstains (`t_vote_thr`, `t_vote_par`, `t_vote_tree`); a majority sets the
  category (Both copies kept (best guess), Lost duplicate), `core` when all three agree and `probable` when two do; no
  majority is Disputed. The published pairs are no longer part of the calls (they were `probable`/`possible`): they stay in
  the table (`t_pairs_*`) as an outside check; they still exclude a later duplicate before the cichlids (an ohnolog). A later
  duplicate in the cichlids needs thresholds and parsimony (`core`) or the tree and one of them (`probable`); before the
  cichlids needs parsimony or the tree, two or more tandem copies and no published pair (`possible`, or `probable` with both).
  A pair in a many-copy family on the tree alone is `possible`. `t_row` is R1 (lost), R2/3/4 (pair), R4 (later before cichlids),
  R5 (later in cichlids), other (no majority); R2 and R3 are no longer separate.
- On the 0.12.0 Compara table: Lost duplicate 10,756 (core 9,278, probable 1,478); Both copies kept 6,800 (core 4,222, probable
  1,747, possible 831); Later duplicate in cichlids 671; before cichlids 112; Disputed 1,067 (was 1,828); Unclassified 7,202.
  90.2% of the published-pair positives land in Both copies kept and 94.8% of the negatives in Lost duplicate.
- Report: the tiers section is rewritten for the votes. 143 offline tests.

## 0.12.0 (tree dating: the pair is the deepest supported duplication)
- Run of 0.11.3 on Compara: the later-duplicate calls fell from 15,304 to 2,418 and the categories matched the preview. But
  1,328 of the 1,381 genes the tree called `later_duplicate_before_cichlids` were at the teleost level with a `disjoint` root and
  looked like pairs (80% on different sequences, 35% in a published pair, 71% human OneToTwo): the top of the gene's teleost
  clade is where the tree disagrees with the species tree, and the real teleost duplication is a node below it.
- The pair node is now the DEEPEST duplication with at least `--pair-overlap` (0.5) shared species on the path from the gene up
  through its teleost clade (new columns `d_dup_level`, `d_dup_at_root`; `d_support`, `d_lineages`, `d_overlap` and
  `d_gene_lineage` describe that node). On the 0.11.3 table, calling a pair where the deepest supported duplication is at the
  teleost level gives 82.4% of the 2,637 published-pair positives at 0.5% of the 8,192 negatives (J 81.9); thresholds 84.0% /
  0.8% (J 83.2); thresholds or tree 91.0% / 0.9% (J 90.2).
- The level of a duplication is `teleost` only if early teleosts are in two lineages (`both_early`), or an early teleost has
  genes on both sides (`shared_early`), or every early-teleost gene of the clade is inside the node (`early_within`: a pair whose
  second copy was lost in all the early teleosts; a later duplicate would leave the early genes outside). Otherwise it is the
  level of the species the two sides share. A stray early-teleost gene inside one percomorph lineage therefore does not make a
  percomorph duplicate a teleost pair.
- `d_call`: `pair` (deepest supported duplication at the teleost level), `pair_unresolved_depth` (the gene's whole teleost clade is
  a duplication with no early teleosts), `later_duplicate_before_cichlids`, `later_duplicate_in_cichlids`, `thin_pair`,
  `lost_duplicate`. Tiers: the tree level also feeds the later-duplicate-before-cichlids rule. 144 offline tests.

## 0.11.3 (later duplications need shared species; many-copy pairs)
- Run of 0.11.2 on Compara: the teleost clade fix dated 90% of genes (65% before) and the pair call works: 56.5% of
  published-pair positives are tree pairs against 0.3% of negatives, none of the 9,541 core lost duplicates is a pair, and
  "thresholds or tree pair" recalls 89.3% at 0.8% false positives (thresholds alone 84.0%).
- Fault: 15,304 genes, 8,768 of the core lost duplicates among them, were called `later_duplicate_before_cichlids`. On the path to
  the teleost root of a gene with one copy in every species there were 5 to 9 nested nodes where reconciliation said duplication
  with 18 to 34 species and none in common: the gene tree disagrees with the species tree. A later duplication now needs the same
  shared species as a pair (the two largest children share at least `--pair-overlap`, 0.5).
- Tiers: a gene in a family with many copies (thresholds or parsimony say three or more copies in many teleosts) whose gene tree
  shows a pair at the teleost root is a pair: `possible`, or `probable` with a published pair. On Compara 510 genes were in a
  published pair 39.0% of the time against 8.8% for many-copy genes with no shared species. Expanded families stay unclassified
  (14.9% against 6.0%). Preview on the 0.11.2 table: Both copies kept 6,412 -> 6,922 (core 5,413, probable 555, possible 954);
  Unclassified 8,033 -> 7,523. 142 offline tests.

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
