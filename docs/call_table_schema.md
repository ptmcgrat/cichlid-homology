# Draft schema of the call tables (proposal for the coordinator; to be finalised after the pilot is run)

Written to `Results/v1/Homology/native/calls/` until the coordinator adopts it. Human and zebrafish only.
Release 1; one set of files per focal species `<prefix>`.

## A. `<prefix>.method_calls.tsv` : one row per gene per method (the evidence)
The columns the pilot writes today (`pilot_genes.<method>.tsv`):

| Column | Meaning |
|---|---|
| `gene_id`, `method`, `family`, `tree_mode` | the gene, TreeFam / Compara / OrthoFinder, its family, `place` or `denovo` |
| `placement_status`, `lwr` | from `gene_placement.tsv` (empty for OrthoFinder) |
| `assign_flags`, `fusion_partner` | from the assignment table, kept as evidence |
| `rooting`, `root_cost`, `n_tied_roots`, `roots_checked` | `recorded` or `reconciliation`; duplications at the chosen root; how many rootings tie; how many were checked |
| `leaves_after_prune`, `shadowed_dropped` | leaves left after the shadow rule; reference leaves dropped |
| `parent_children`, `species_in_tree` | children of the gene's parent node (multifurcation size); species left |
| `<ref>_relationship` | `1:1`, `1:many`, `many:1`, `many:many`, `paralog`, `none` (`<ref>` = `human`, `zebrafish`); first number counts reference genes, second counts focal co-orthologs |
| `<ref>_n_ref`, `<ref>_n_focal` | the two counts |
| `<ref>_closest`, `<ref>_closest_is_ortholog`, `<ref>_closest_tie` | closest reference leaf (smallest clade; ties by path length), whether it is an ortholog, whether several were tied |
| `<ref>_n_ref_in_family` | reference genes in this method's family; 0 means `none` |
| `<ref>_stable_across_roots` | the relationship is the same under all tied rootings checked |
| `<ref>_paralog_level` | ladder level of the duplication separating gene and closest reference gene (paralogs); `cyclostome_disputed` is never a date |
| `<ref>_dup_level`, `<ref>_dup_class`, `<ref>_dup_resolved`, `<ref>_dup_levels_sides`, `<ref>_dup_outgroup_single` | for genes with several focal co-orthologs: ladder level of the nearest duplication, class (`teleost_ohnolog`, `teleost_ohnolog_no_outgroup`, `teleost_level_not_spread`, `teleost_level_outgroup_multi`, `lineage_specific:<level>`, `older_than_teleost:<level>`, `cyclostome_disputed`, `unresolved_place_mode`), whether the tree can resolve it, teleost levels on each side, whether gar/reedfish hold one copy |
| `teleost_tel_species`, `teleost_tel_one`, `teleost_tel_two_plus`, `teleost_tel_zero`, `teleost_tel_frac_one_of_present` | copies of the gene's orthologs per teleost species |

To add before the full run: `<ref>_closest_symbol` (when the identifier mappings exist), `reason` (a code), `confidence`.

## B. `<prefix>.calls.tsv` : the decision, one row per gene per question
`gene_id, question, decision, agreement, n_answering, methods, treefam, compara, orthofinder, families, confidence, reason, review`

- `question`: `fusion`, `human_relationship`, `human_closest`, `zebrafish_relationship`, `zebrafish_closest`, `teleost_pattern`, `duplication_class`, `duplication_level`, `paralog_level`.
- `agreement`: `agree`, `majority`, `split`, `single`, `none`. A method whose family holds no gene of the reference species answers nothing (`none` is not a value).
- `decision` is empty for `split` unless a rule prefers one method; every competing value stays in the method columns.
- `reason`: a code: `no_ref_in_family`, `closest_not_ortholog`, `closest_tie`, `unstable_root`, `place_mode_unresolved`, `low_placement_support`, `denovo_family`, `gene_not_a_leaf`, `no_family`, `family_skipped`, `cyclostome_disputed`.
- `review`: True for every `split`, and for `single` where the one answer is low confidence.

## C. `<prefix>.review.tsv` : the rows of B with `review = True`.
## D. `<prefix>.missing_annotation.tsv` : families conserved in other cichlids/teleosts with no gene of the focal species (candidates only): `family, method, n_species_with_gene, n_cichlids_with_gene, species_with_gene`.


## Draft additions in 0.10.0 (B1 tiers)
`t_row` (R1, R2, R3, R2/3/4, R5, other, outside), `t_call` (lost_duplicate, retained_pair, cichlid_specific_duplicate,
unresolved, no_anchor, expanded_family, many_copies), `t_tier` (core, probable, possible, disputed, outside), `t_basis`
(a short reason), `t_positions` (dispersed, tandem or empty), and, when published pairs are given, `t_pairs_relaxed_n`,
`t_pairs_strict_n`, `t_pairs_species`. Parsimony: `p_label` (gain 2, loss 1), `p_label_unit` (equal costs), `p_robust`,
`p_margin`, `p_cost1` to `p_cost3`, and the event columns.
