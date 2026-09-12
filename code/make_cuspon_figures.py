#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Generate two figures for the fourth-order regularization of the IDD cuspon.

1) regularized_bvp_profiles.{png,pdf}
   Regularized solitary-wave cores for epsilon = 1e-2, 1e-4, 1e-6, 1e-8.

2) cuspon_phase_portrait.{png,pdf}
   Reduced singular cuspon in the (U,U_xi) phase plane together with the
   projections of the fourth-order homoclinic orbits onto the same plane.

The regularized profiles solve, for omega=1,

    epsilon * U^(4) - (1-U^2) * U'' + U = 0,

on the half-line [0,XMAX], with even-center conditions and exact linear
stable-subspace boundary conditions at XMAX. The negative half-line is
obtained by reversibility.

Dependencies: numpy, scipy, matplotlib.
'''

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_bvp
from scipy.special import lambertw


EPSILONS = (1e-2, 1e-4, 1e-6, 1e-8)
XMAX = 20.0
TOL = 1.0e-6
BC_TOL = 1.0e-10
MAX_NODES = 60000

# Second-order asymptotic coefficients used only to construct a robust
# initial guess for solve_bvp.
A0, A1, A2 = 0.3933180528, 1.1065076820, 0.0131925717
K0, K1, K2 = 0.6898804165, 0.3022408770, 0.6573620821


def L_eps(eps: float) -> float:
    return float(lambertw(3.0 / (8.0 * eps)).real / 3.0)


def asymptotic_seed_scales(eps: float) -> tuple[float, float]:
    L = L_eps(eps)
    overshoot_scaled = A0 + A1 / L + A2 / L**2
    curvature_scaled = K0 + K1 / L + K2 / L**2
    overshoot = (eps * L)**(1.0 / 3.0) * overshoot_scaled
    curvature = eps**(-1.0 / 3.0) * L**(2.0 / 3.0) * curvature_scaled
    return overshoot, curvature


def stable_rates(eps: float) -> tuple[float, float]:
    disc = math.sqrt(1.0 - 4.0 * eps)
    lam_s = math.sqrt(2.0 / (1.0 + disc))
    lam_f = math.sqrt((1.0 + disc) / (2.0 * eps))
    return lam_s, lam_f


def ode_and_jacobian(eps: float):
    def rhs(x, y):
        U, Up, Upp, Uppp = y
        return np.vstack(
            (
                Up,
                Upp,
                Uppp,
                ((1.0 - U * U) * Upp - U) / eps,
            )
        )

    def jac(x, y):
        U, Up, Upp, Uppp = y
        J = np.zeros((4, 4, x.size))
        J[0, 1] = 1.0
        J[1, 2] = 1.0
        J[2, 3] = 1.0
        J[3, 0] = (-2.0 * U * Upp - 1.0) / eps
        J[3, 2] = (1.0 - U * U) / eps
        return J

    return rhs, jac


def boundary_conditions(eps: float):
    lam_s, lam_f = stable_rates(eps)
    s = lam_s + lam_f
    p = lam_s * lam_f

    def bc(ya, yb):
        return np.array(
            [
                ya[1],
                ya[3],
                yb[2] + s * yb[1] + p * yb[0],
                yb[3] + s * yb[2] + p * yb[1],
            ]
        )

    def bc_jac(ya, yb):
        A = np.zeros((4, 4))
        B = np.zeros((4, 4))
        A[0, 1] = 1.0
        A[1, 3] = 1.0
        B[2] = [p, s, 1.0, 0.0]
        B[3] = [0.0, p, s, 1.0]
        return A, B

    return bc, bc_jac


def initial_mesh(eps: float, xmax: float) -> np.ndarray:
    L = L_eps(eps)
    delta = eps**(1.0 / 3.0) * L**(-1.0 / 6.0)
    fast = math.sqrt(eps)

    core_end = min(0.8, max(12.0 * delta, 40.0 * fast, 0.02))
    h_core = min(fast / 5.0, delta / 30.0, 2.0e-3)
    n_core = max(
        220,
        min(3200, int(math.ceil(core_end / h_core)) + 1),
    )

    pieces = [np.linspace(0.0, core_end, n_core)]
    if core_end < 2.0:
        pieces.append(np.linspace(core_end, 2.0, 650))
    if xmax > 2.0:
        pieces.append(np.linspace(2.0, xmax, 550))
    return np.unique(np.concatenate(pieces))


def initial_guess(eps: float, x: np.ndarray) -> np.ndarray:
    overshoot, curvature = asymptotic_seed_scales(eps)
    amp = 1.0 + overshoot

    radius = amp / curvature
    r = np.sqrt(x * x + radius * radius)
    U = amp * np.exp(radius - r)

    q = -x / r
    qp = -radius * radius / r**3
    qpp = 3.0 * radius * radius * x / r**5

    return np.vstack(
        (
            U,
            U * q,
            U * (q * q + qp),
            U * (q**3 + 3.0 * q * qp + qpp),
        )
    )


def solve_regularized(eps: float, xmax: float = XMAX):
    x = initial_mesh(eps, xmax)
    y = initial_guess(eps, x)
    rhs, jac = ode_and_jacobian(eps)
    bc, bc_jac = boundary_conditions(eps)

    sol = solve_bvp(
        rhs,
        bc,
        x,
        y,
        fun_jac=jac,
        bc_jac=bc_jac,
        tol=TOL,
        bc_tol=BC_TOL,
        max_nodes=MAX_NODES,
        verbose=0,
    )

    if sol.status != 0:
        raise RuntimeError(
            f"solve_bvp failed for eps={eps:g}: {sol.message}; "
            f"nodes={sol.x.size}; "
            f"max_rms={np.max(sol.rms_residuals):.3e}"
        )
    return sol


def save_figure(fig, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_profile_figure(solutions: dict[float, object], figure_dir: Path, data_dir: Path) -> None:
    xi = np.linspace(-0.6, 0.6, 2401)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    profile_rows = {"xi": xi}

    linestyles = ("-", "--", "-.", ":")
    for eps, ls in zip(EPSILONS, linestyles):
        U = solutions[eps].sol(np.abs(xi))[0]
        profile_rows[f"U_eps_{eps:.0e}"] = U
        exponent = int(round(math.log10(eps)))
        ax.plot(
            xi,
            U,
            linestyle=ls,
            linewidth=1.8,
            label=rf"$\varepsilon=10^{{{exponent}}}$",
        )

    ax.axhline(1.0, linestyle=(0, (2, 2)), linewidth=0.9)
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel(r"$U_\varepsilon(\xi)$")
    ax.set_title("Regularized solitary-wave core")
    ax.legend(frameon=True)
    fig.tight_layout()

    save_figure(fig, figure_dir / "regularized_bvp_profiles")

    fields = list(profile_rows)
    with open(data_dir / "regularized_bvp_profiles_data.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(fields)
        for i in range(xi.size):
            wr.writerow([profile_rows[k][i] for k in fields])


def make_phase_portrait(solutions: dict[float, object], figure_dir: Path, data_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.6))

    # Singular reduced cuspon:
    #     (U_xi)^2 = -log(1-U^2),  0<U<1.
    Uc = np.linspace(1.0e-5, 1.0 - 1.0e-10, 8000)
    Pc = np.sqrt(-np.log1p(-Uc * Uc))

    # One Line2D object (with a NaN separator) keeps the two singular
    # branches in the same color without manually fixing any color.
    Uc_both = np.concatenate((Uc, [np.nan], Uc))
    Pc_both = np.concatenate((Pc, [np.nan], -Pc))
    ax.plot(
        Uc_both,
        Pc_both,
        linestyle=(0, (5, 2)),
        linewidth=2.0,
        label="singular cuspon",
    )

    phase_rows = []
    linestyles = ("-", "--", "-.", ":")
    xi_pos = np.linspace(0.0, 10.0, 5001)

    for eps, ls in zip(EPSILONS, linestyles):
        Y = solutions[eps].sol(xi_pos)
        U = Y[0]
        P = Y[1]

        U_full = np.concatenate((U[::-1], U[1:]))
        P_full = np.concatenate((-P[::-1], P[1:]))

        exponent = int(round(math.log10(eps)))
        ax.plot(
            U_full,
            P_full,
            linestyle=ls,
            linewidth=1.6,
            label=rf"regularized $\varepsilon=10^{{{exponent}}}$",
        )

        for uu, pp in zip(U_full, P_full):
            phase_rows.append((eps, uu, pp))

    ax.set_xlabel(r"$U$")
    ax.set_ylabel(r"$U_\xi$")
    ax.set_title(r"Homoclinic orbit: $(U,U_\xi)$ phase-plane projection")
    ax.set_xlim(0.0, 1.30)
    ax.set_ylim(-3.8, 3.8)
    ax.legend(frameon=True, fontsize=8)
    fig.tight_layout()

    save_figure(fig, figure_dir / "cuspon_phase_portrait")

    with open(data_dir / "cuspon_phase_regularized_data.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["epsilon", "U", "U_xi"])
        wr.writerows(phase_rows)

    with open(data_dir / "cuspon_phase_singular_data.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["U", "U_xi_plus", "U_xi_minus"])
        wr.writerows(zip(Uc, Pc, -Pc))


def main() -> None:
    ap = argparse.ArgumentParser()
    root = Path(__file__).resolve().parent.parent
    ap.add_argument(
        "--figure-dir",
        default=str(root / "figures"),
        help="Directory for PNG/PDF figure outputs.",
    )
    ap.add_argument(
        "--data-dir",
        default=str(root / "results"),
        help="Directory for CSV figure-data outputs.",
    )
    args = ap.parse_args()

    figure_dir = Path(args.figure_dir)
    data_dir = Path(args.data_dir)
    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    solutions = {}
    for eps in EPSILONS:
        print(f"Solving epsilon={eps:.0e} ...", flush=True)
        sol = solve_regularized(eps)
        solutions[eps] = sol
        print(
            f"  nodes={sol.x.size}, "
            f"max_rms={np.max(sol.rms_residuals):.3e}, "
            f"U(0)={sol.y[0,0]:.12f}",
            flush=True,
        )

    make_profile_figure(solutions, figure_dir, data_dir)
    make_phase_portrait(solutions, figure_dir, data_dir)

    print("Generated:")
    for name in (
        "regularized_bvp_profiles.png",
        "regularized_bvp_profiles.pdf",
        "cuspon_phase_portrait.png",
        "cuspon_phase_portrait.pdf",
        "regularized_bvp_profiles_data.csv",
        "cuspon_phase_regularized_data.csv",
        "cuspon_phase_singular_data.csv",
    ):
        base = figure_dir if name.endswith((".png", ".pdf")) else data_dir
        print(" ", base / name)


if __name__ == "__main__":
    main()
