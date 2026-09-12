# v1.0.0 — first reproducibility release

This release freezes the numerical artifacts used with **“Fourth-Order Dispersive Regularization of a Cuspon.”**

## Included

- deterministic SciPy `solve_bvp` rerun of the full fourth-order boundary-value problem;
- the validation table for the second-order asymptotic observables;
- comparison with the earlier exploratory BVP table;
- `Xmax=12,16,20` domain-refinement check at `epsilon=1e-8`;
- the observable-validation figure;
- regularized-core profile figure;
- projected phase portrait of the regularized homoclinics and singular reduced cuspon;
- diagnostic universal-inner scripts (`Y0`, `Y1`, `Y2`);
- pinned Python dependencies, logs, hashes, and citation metadata.

The numerical computations are diagnostic and are not used to establish the analytical existence theorem.
