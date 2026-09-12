# Fourth-Order Dispersive Regularization of a Cuspon — reproducibility package

Reproducibility code, numerical data, and figure sources for

**Francisco R. Villatoro, “Fourth-Order Dispersive Regularization of a Cuspon.”**

This repository contains **no manuscript source or article PDF**. It is restricted to scientific-computation code, frozen numerical inputs/outputs, reproducibility logs, citation metadata, licenses, and figures used to document the computations.

## What is reproduced

The main deterministic calculation solves, for `omega = 1`,

```text
epsilon U'''' - (1-U^2) U'' + U = 0
```

on `0 <= xi <= Xmax`, with the reversible center conditions

```text
U'(0) = U'''(0) = 0
```

and two linear stable-tail boundary conditions at `Xmax`. The production rerun uses `Xmax = 20`, `scipy.integrate.solve_bvp`, residual tolerance `1e-6`, boundary-condition tolerance `1e-10`, and at most 60000 collocation nodes.

The repository reproduces:

- the full BVP validation table for `epsilon = 1e-4,...,1e-8`;
- comparison against the archived exploratory BVP table;
- the domain-refinement check `Xmax = 12,16,20` at `epsilon = 1e-8`;
- the scaled-observable validation figure;
- regularized solitary-wave core profiles;
- the `(U,U_xi)` phase-plane projection of the fourth-order homoclinics together with the singular reduced cuspon.

The numerical rerun is **diagnostic**: no existence theorem in the article depends on the BVP solver.

## Repository layout

```text
code/
  rerun_regularized_bvp.py   deterministic full-BVP rerun
  make_cuspon_figures.py     core-profile and phase-plane figures
  y0_shoot_revised.py        diagnostic universal-separatrix shooting
  y12.py                     diagnostic Y1/Y2 correction study
  check_results.py           numerical consistency checks

data/
  regularized_bvp_results_archived.csv
  cuspon_inner_corrections_O1_O2.csv

results/
  bvp_rerun.csv
  numerical_validation_table_rerun.csv
  bvp_vs_archived.csv
  bvp_domain_refinement.csv
  bvp_environment.json
  regularized_bvp_profiles_data.csv
  cuspon_phase_regularized_data.csv
  cuspon_phase_singular_data.csv

figures/
  fig_observables_validation.pdf
  regularized_bvp_profiles.{pdf,png}
  cuspon_phase_portrait.{pdf,png}

logs/
  reproducibility logs for release v1.0.0
```

## Environment

The archived release rerun used:

- Python 3.13.5
- NumPy 2.3.5
- SciPy 1.17.0
- Matplotlib 3.10.8

Exact pinned Python dependencies are in `requirements.txt` and the recorded runtime metadata are in `results/bvp_environment.json`.

## Reproduce in WSL/Linux

A clean virtual environment is recommended:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
bash reproduce_all.sh
```

A complete run performs the BVP calculation twice independently: once for the quantitative validation table and once for the two geometric/profile figures. This is intentional and provides a useful cross-check. On the archived environment the BVP validation rerun itself takes roughly 30 seconds.

To run only the quantitative BVP validation:

```bash
python code/rerun_regularized_bvp.py
python code/check_results.py
```

To regenerate the two additional figures:

```bash
python code/make_cuspon_figures.py
```

## Release-v1.0.0 numerical checks

The rerun agrees with the archived exploratory table to maximum relative differences of approximately

- `5.9e-10` in `U(0)-1`,
- `1.3e-9` in the interface location,
- `1.3e-9` in `U''(0)`.

For `epsilon = 1e-8`, changing the truncation from `Xmax=20` to `12` or `16` changes the three reported observables by at most about `1.2e-10` relative to the `Xmax=20` values.

At `epsilon = 1e-8`, the second-order relative errors of the asymptotic formulas are about `1.02%`, `0.14%`, and `0.20%` for overshoot, interface location, and central curvature. The overshoot data do **not** independently resolve the coefficient `A2`; this is why the article describes the computation as consistent with the expansion rather than as a numerical proof of that coefficient.

## Universal-inner diagnostics

`code/y0_shoot_revised.py` and `code/y12.py` are retained as transparent diagnostics for the universal inner profile and its first two correction problems. They are **not interval-arithmetic certificates**. In particular, the finite-cutoff `Y2` computation is used as a convergence diagnostic, not as part of the rigorous proof.

## Integrity

`SHA256SUMS.txt` contains SHA-256 hashes for the release files. `release_manifest.tsv` records path, size, and hash. Run

```bash
sha256sum -c SHA256SUMS.txt
```

to verify the downloaded release tree (excluding the checksum/manifest files themselves).

## Citation

Citation metadata are provided in `CITATION.cff`. After the GitHub `v1.0.0` release is archived by Zenodo, the version DOI can be used to cite the exact reproducibility snapshot.

## Licenses

- Original source code and shell scripts: MIT License (`LICENSE`).
- Numerical data, generated figures, and repository documentation: CC BY 4.0 (`DATA_LICENSE.md`).
