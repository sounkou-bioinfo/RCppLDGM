#!/usr/bin/env python3
"""Generate pinned GraphLD BLUP/clump outputs for RcppLDGM conformance checks."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import polars as pl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graphld-root", default=".sync/graphld")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sumstats", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--population", default="EUR")
    parser.add_argument("--sigmasq", type=float, default=0.01)
    parser.add_argument("--sample-size", type=float, required=True)
    parser.add_argument("--rsq-threshold", type=float, default=0.1)
    parser.add_argument("--chisq-threshold", type=float, default=30.0)
    return parser.parse_args()


def require_columns(frame: pl.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def main() -> None:
    args = parse_args()

    graphld_root = Path(args.graphld_root).resolve()
    graphld_src = graphld_root / "src"
    if not graphld_src.exists():
        raise FileNotFoundError(f"GraphLD source directory not found: {graphld_src}")
    sys.path.insert(0, os.fspath(graphld_src))

    import graphld  # type: ignore

    metadata_path = Path(args.metadata).resolve()
    sumstats_path = Path(args.sumstats).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    sumstats = pl.read_csv(sumstats_path, separator="\t")
    required_columns = ["SNP", "CHR", "POS", "A1", "A2", "REF", "ALT", "Z", "N"]
    require_columns(sumstats, required_columns, "prepared summary statistics")

    blup = graphld.run_blup(
        os.fspath(metadata_path),
        sumstats,
        sigmasq=float(args.sigmasq),
        sample_size=float(args.sample_size),
        populations=args.population,
        run_in_serial=True,
        match_by_position=False,
    )
    require_columns(blup, required_columns + ["weight"], "GraphLD BLUP output")
    blup.select(required_columns + ["weight"]).write_csv(out_dir / "blup.csv")

    clump = graphld.run_clump(
        sumstats,
        ldgm_metadata_path=os.fspath(metadata_path),
        rsq_threshold=float(args.rsq_threshold),
        chisq_threshold=float(args.chisq_threshold),
        populations=args.population,
        run_in_serial=True,
        match_by_position=True,
    )
    require_columns(clump, required_columns + ["is_index"], "GraphLD clump output")
    clump.select(required_columns + ["is_index"]).write_csv(out_dir / "clump.csv")


if __name__ == "__main__":
    main()
