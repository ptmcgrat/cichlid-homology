# The clade record (cichlid-homology/clade-record/1)

One JSON object per method and clade, written one per line to `clade_records.jsonl` by `figure-pilot`. The figure is drawn from the record alone, and the site can draw the same record.
The schema is for the pilot: it is to be designed with the coordinator after the pilot, like the call tables.

| key | meaning |
|---|---|
| `schema`, `method`, `family`, `clade`, `tree_mode`, `root_source`, `focal_species` | provenance. `root_source` is "like the reference root" when the placed tree was rooted where the reference tree is rooted |
| `n_leaves`, `n_own` | leaves of the clade, and our own genes of the focal species among them |
| `mode` | `split` (two or more lineages), `single` (teleost genes, one lineage) or `none` (no teleost gene) |
| `split` | null, or `{kind: teleost or later, level, overlap, n_lineages}`. A teleost split has early teleosts on both sides; a later split is a duplication among the teleosts only |
| `letters`, `lineages` | one entry per lineage: `id`, `species`, `genes`, `stem` and `depth` (branch lengths from the split to the crown, and inside the crown), `crown_support` (bootstrap of the crown). The letters are ordered by the zebrafish gene of each lineage (irx1a before irx1b), then by size |
| `outgroup` | `inside`: non-teleost genes whose place in the tree is inside a lineage, with `attach_support` (the bootstrap of the node that joins them to the teleosts); `outside`: the others. They are the ancestral copy either way |
| `strip` | `species`, `labels`, `groups`, `counts`: for every reference species `[outside, lineage A, lineage B, ..., stray]` copies. A stray teleost gene is a child of the split that shares no species with the lineages |
| `focal_strip` | the same for our own species (species, labels, own, counts) |
| `genes` | our own genes of the focal species: `gene`, `lineage` (A, B, stray, outside, or "outside the clade"), `lwr` (the weight of the placement), `status`, `n_placements`, `calls` (columns of the gene table) |
| `anchors` | the zebrafish and human genes of each lineage (and `outside`), with their symbols: the names that tie a lineage to the same lineage in another method |
| `tree` | the collapsed tree: nodes `{k: node, support, dup, kids, split}`, blocks `{k: block, label, n_genes, n_species, support, dup, own, other_focal}`, and single leaves `{k: leaf, name, species, group, symbol}` (the non-teleost fishes) |
| `support_items`, `support_available` | the supports of the nodes that matter, for the second panel; false when the tree carries no support values |
| `deep_order` | the four relationships we know, each true, false or null when its species are not in the clade: `shark_outside`, `coelacanth_with_tetrapods`, `reedfish_outside_gar_teleosts`, `gar_outside_teleosts`; with `tested`, `passed`, `score`. A quality score for the deep nodes of the tree: the coelacanth is the control, gar and reedfish are the ones that wander |
| `triage` | `overlap`, `shared_species`, `split_support`, `weakest_crown`, `outgroup_inside`, `stem_ratio`, `deep_order_score`: numbers to rank clades by |
| `flags` | `no_split`, `no_teleost_genes`, `deep_order_inconsistent` (at least one tested relationship not recovered), `outgroup_inside_lineage`, `weak_crown_support` (under 50), `stray_teleost_genes`, `no_support_values`, `clade_not_exact_in_placed_tree`, `whole_tree`, `no_reference_clades`, `gene_outside_clade_node` |

Support values are Ensembl's bootstrap (0 to 100) and its duplication flag. They exist for the Compara reference trees only (read from the raw dump); a placed gene carries the weight of its placement instead.
OrthoFinder trees carry no support yet; a record built from one has `support_available` false until they do.
