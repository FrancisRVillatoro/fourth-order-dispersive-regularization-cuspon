#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
mkdir -p logs results figures

python code/rerun_regularized_bvp.py | tee logs/reproduce_bvp.log
python code/check_results.py | tee logs/check_results.log
python code/make_cuspon_figures.py | tee logs/reproduce_figures.log
python code/y0_shoot_revised.py | tee logs/y0_shoot_revised.log
python code/y12.py | tee logs/y12.log

echo "Reproduction finished successfully."
