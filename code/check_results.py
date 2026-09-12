#!/usr/bin/env python3
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"


def rows(name):
    with (RESULTS / name).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def f(x):
    return float(x)


def main():
    comp = rows("bvp_vs_archived.csv")
    limits = {
        "u0_relative_difference": 1e-8,
        "x_cross_relative_difference": 1e-8,
        "u_xx0_relative_difference": 1e-8,
    }
    for key, lim in limits.items():
        value = max(f(r[key]) for r in comp)
        print(f"{key}: {value:.6e}  limit={lim:.1e}")
        if value > lim:
            raise SystemExit(f"FAIL: {key} exceeds {lim:g}")

    dom = rows("bvp_domain_refinement.csv")
    for key in ("u0_rel_to_xmax20", "x_cross_rel_to_xmax20", "u_xx0_rel_to_xmax20"):
        value = max(f(r[key]) for r in dom)
        print(f"domain {key}: {value:.6e}  limit=1e-8")
        if value > 1e-8:
            raise SystemExit(f"FAIL: {key} exceeds 1e-8")

    val = rows("numerical_validation_table_rerun.csv")
    if len(val) != 5:
        raise SystemExit(f"FAIL: expected 5 validation rows, found {len(val)}")
    eps = [f(r["epsilon"]) for r in val]
    if eps != [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
        raise SystemExit(f"FAIL: unexpected epsilon grid: {eps}")

    print("PASS: numerical release checks")


if __name__ == "__main__":
    main()
