# cichlid-homology

Project 3 of the cichlid protein phylogeny resource: every orthology call, for all three methods
(TreeFam, Compara, OrthoFinder), with one rule. See `docs/homology_repo_spec.md`.

Depends on `cichlid-coordinator>=0.13,<0.14` and uses only its public API (`cichlid_data`) and
the `cichlid-data` command line. No hand-built data paths, no rclone, no downloads.

## Install (development)

    python -m pip install -e ../cichlid-coordinator     # once, if not already in the environment
    python -m pip install -e ".[dev]"

## Check and test

    cichlid-data check --release v1 --part Homology
    pytest                                               # offline

## Survey (tasks 1-2, read-only)

    cichlid-homology survey --species Mzeb --n 300 --out ~/survey.json

## Pilot

    cichlid-homology probe
    cichlid-homology pilot --species Mzeb --n 600 --n-denovo 100 --out ~/pilot
    cichlid-homology explain Mzeb_101470768 --part TreeFam

See `docs/pilot_plan.md`.

## Status

0.2.x: pilot code (rooting, relationships, dating, method comparison), tested offline, not yet run on
real data. The call-table schema is a draft (`docs/call_table_schema.md`); nothing is written to
`Results/v1/Homology/` except scratch files under `native/calls/pilot`, which are not pushed.
