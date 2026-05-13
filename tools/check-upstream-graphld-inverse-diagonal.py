#!/usr/bin/env python3
"""Generate upstream GraphLD inverse-diagonal conformance fixtures."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import polars as pl

SELECT_LIMIT = 40
SELECT_STRIDE = 2
FULL_PROBE_COUNT = 16
SELECTED_PROBE_COUNT = 12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graphld-root", default=".sync/graphld")
    parser.add_argument("--data-dir", default=".sync/graphld/data/test")
    parser.add_argument("--population", default="EUR")
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def first_population_block(metadata: pl.DataFrame) -> dict:
    if metadata.height == 0:
        raise RuntimeError("No metadata rows found for requested population")
    metadata = metadata.sort(["chrom", "chromStart"])
    return metadata.row(0, named=True)


def rademacher_probes(n_rows: int, n_cols: int, seed: int) -> np.ndarray:
    rng = np.random.RandomState(seed)
    return rng.choice([-1.0, 1.0], size=(n_rows, min(n_rows, n_cols)))


def write_matrix(path: Path, matrix: np.ndarray) -> None:
    np.savetxt(path, np.asarray(matrix, dtype=float), delimiter=",")


def main() -> None:
    args = parse_args()

    graphld_root = Path(args.graphld_root).resolve()
    graphld_src = graphld_root / "src"
    if not graphld_src.exists():
        raise FileNotFoundError(f"GraphLD source directory not found: {graphld_src}")
    sys.path.insert(0, os.fspath(graphld_src))

    from graphld.io import load_ldgm, read_ldgm_metadata  # type: ignore

    data_dir = Path(args.data_dir).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata = read_ldgm_metadata(os.fspath(data_dir / "metadata.csv"), populations=args.population)
    row = first_population_block(metadata)
    edgelist_name = row["name"]
    snplist_name = row["snplistName"]
    block_name = Path(edgelist_name).stem

    ldgm = load_ldgm(
        os.fspath(data_dir / edgelist_name),
        snplist_path=os.fspath(data_dir / snplist_name),
    )

    full_n = int(ldgm.shape[0])
    selected_indices_zero = np.arange(0, min(full_n, SELECT_LIMIT), SELECT_STRIDE, dtype=int)
    selected = ldgm[selected_indices_zero]

    full_probes = rademacher_probes(full_n, FULL_PROBE_COUNT, args.seed)
    selected_probes = rademacher_probes(int(selected.shape[0]), SELECTED_PROBE_COUNT, args.seed + 1)

    full_exact = np.asarray(ldgm.inverse_diagonal(method="exact"), dtype=float).reshape(-1)
    full_hutchinson, full_hutchinson_solved = ldgm.inverse_diagonal(
        method="hutchinson",
        initialization=(full_probes, full_probes.copy()),
    )
    full_xdiag, full_xdiag_solved = ldgm.inverse_diagonal(
        method="xdiag",
        initialization=(full_probes, full_probes.copy()),
    )

    selected_exact = np.asarray(selected.inverse_diagonal(method="exact"), dtype=float).reshape(-1)
    selected_hutchinson, selected_hutchinson_solved = selected.inverse_diagonal(
        method="hutchinson",
        initialization=(selected_probes, selected_probes.copy()),
    )
    selected_xdiag, selected_xdiag_solved = selected.inverse_diagonal(
        method="xdiag",
        initialization=(selected_probes, selected_probes.copy()),
    )

    pl.DataFrame(
        {
            "block_name": [block_name],
            "edgelist_name": [edgelist_name],
            "snplist_name": [snplist_name],
            "population": [args.population],
            "full_n": [full_n],
            "selected_n": [int(selected.shape[0])],
            "selected_limit": [int(min(full_n, SELECT_LIMIT))],
            "selected_stride": [SELECT_STRIDE],
            "full_probe_cols": [int(full_probes.shape[1])],
            "selected_probe_cols": [int(selected_probes.shape[1])],
            "seed": [int(args.seed)],
        }
    ).write_csv(out_dir / "summary.csv")

    pl.DataFrame(
        {
            "exact": full_exact,
            "hutchinson": np.asarray(full_hutchinson, dtype=float).reshape(-1),
            "xdiag": np.asarray(full_xdiag, dtype=float).reshape(-1),
        }
    ).write_csv(out_dir / "full_diagonal.csv")
    pl.DataFrame(
        {
            "exact": selected_exact,
            "hutchinson": np.asarray(selected_hutchinson, dtype=float).reshape(-1),
            "xdiag": np.asarray(selected_xdiag, dtype=float).reshape(-1),
        }
    ).write_csv(out_dir / "selected_diagonal.csv")

    write_matrix(out_dir / "full_probes.csv", full_probes)
    write_matrix(out_dir / "selected_probes.csv", selected_probes)
    write_matrix(out_dir / "full_hutchinson_solved.csv", full_hutchinson_solved)
    write_matrix(out_dir / "selected_hutchinson_solved.csv", selected_hutchinson_solved)
    write_matrix(out_dir / "full_xdiag_solved.csv", full_xdiag_solved)
    write_matrix(out_dir / "selected_xdiag_solved.csv", selected_xdiag_solved)


if __name__ == "__main__":
    main()
