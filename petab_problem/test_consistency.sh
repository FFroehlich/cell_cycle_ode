#!/usr/bin/env bash
petablint -s petab/cell_cycle_no_M.sbml -m petab/BT20_measurements.tsv -c petab/BT20_conditions.tsv -p petab/BT20_parameters.tsv
petablint -y petab/BT20.yaml