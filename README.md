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

## Status

0.1.x: reading layer and survey only. The pilot (rooting, relationships, dating, method
comparison) has not started; the call-table schema is not settled and nothing is written to
`Results/v1/Homology/` yet.
