# Changelog

## 0.1.1
- Pinned to cichlid-coordinator >=0.13,<0.14.
- `species.py`: the decided rule for reference leaves that are the same genome as a placed species
  (`drop_shadowed`, `SHADOW_RULE`), and a loader for the coordinator's `LeafSpecies`.
- `survey2.py` and `cichlid-homology survey2`: species of leaves, shadowed leaves, coverage of human,
  zebrafish and non-teleost fish, placement support, fusion columns, and a targeted de novo sample.
- Removed `scripts/survey_families.py` (replaced by `cichlid-homology survey`).
- Offline tests for the shadow rule, placement support and fusion counts.

## 0.1.0
- Repository scaffold, `read.py` (`Family`), `cichlid-homology survey`, offline tests.
