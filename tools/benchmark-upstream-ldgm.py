#!/usr/bin/env python3
"""Benchmark pinned upstream Python ldgm operations.

The R wrapper `tools/benchmark-upstream-ldgm.R` calls this script so RcppLDGM
can report apples-to-apples timings against the checked-out upstream Python code
for the examples listed in the generated golden manifest.
"""

from __future__ import annotations

import argparse
import csv
import inspect
import subprocess
import sys
import time
from pathlib import Path
from typing import Any


DEFAULT_EXAMPLES = [
    "single_tree_ts_n2_2_mutations",
    "single_tree_ts_mutation_n3",
    "multiple_snps_branch",
    "figure_one_example",
    "gils_example_tree",
    "single_tree_all_samples_one_mutation_n3",
    "supplementary_example",
    "triangle_example",
    "two_tree_mutation_ts",
    "two_tree_two_mrcas",
]


def require_upstream_modules(repo_root: Path):
    sys.path.insert(0, str(repo_root))
    missing = []
    for module in ("networkx", "msprime", "tskit", "numpy", "pandas"):
        try:
            __import__(module)
        except Exception as exc:  # pragma: no cover - dependency diagnostic path
            missing.append(f"{module}: {exc}")
    if missing:
        msg = "Missing upstream Python dependencies:\n  " + "\n  ".join(missing)
        msg += "\nRun `make upstream-python` before benchmarking upstream ldgm."
        raise SystemExit(msg)

    import ldgm  # type: ignore
    from tests import utility_functions  # type: ignore

    return ldgm, utility_functions


def git_commit(repo_root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo_root, text=True
        ).strip()
    except Exception:
        return "unknown"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="CSV path for upstream timing rows")
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Repeat each timed operation this many times and report elapsed seconds per call",
    )
    parser.add_argument("--path-weight-threshold", type=float, default=100.0)
    parser.add_argument("--edge-weight-threshold", type=float, default=None)
    parser.add_argument("--make-sibs", action="store_true")
    parser.add_argument("--examples", nargs="*", default=DEFAULT_EXAMPLES)
    return parser.parse_args()


def timed(operation: str, func, batch_size: int) -> tuple[str, Any, float]:
    start = time.perf_counter()
    result = None
    for _ in range(batch_size):
        result = func()
    return operation, result, (time.perf_counter() - start) / batch_size


def main() -> int:
    args = parse_args()
    if args.iterations < 1:
        raise SystemExit("--iterations must be positive")
    if args.batch_size < 1:
        raise SystemExit("--batch-size must be positive")

    repo_root = Path(__file__).resolve().parents[1]
    ldgm, utility_functions = require_upstream_modules(repo_root)
    commit = git_commit(repo_root)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for example in args.examples:
        constructor = getattr(utility_functions, example, None)
        if constructor is None or not callable(constructor) or inspect.getmodule(constructor) is not utility_functions:
            raise SystemExit(f"upstream example constructor not found: {example}")

        for iteration in range(1, args.iterations + 1):
            op, bricked, elapsed = timed(
                "brick_ts",
                lambda: ldgm.brick_ts(
                    constructor(),
                    recombination_freq_threshold=None,
                    progress=False,
                ),
                args.batch_size,
            )
            rows.append(
                {
                    "implementation": "upstream_python",
                    "example": example,
                    "operation": op,
                    "iteration": iteration,
                    "elapsed_sec": elapsed,
                    "batch_size": args.batch_size,
                    "num_edges": int(bricked.num_edges),
                    "num_mutations": int(bricked.num_mutations),
                    "brick_graph_edges": "",
                    "reduced_edges": "",
                    "upstream_commit": commit,
                }
            )

            op, brick_graph, elapsed = timed(
                "brick_haplo_graph",
                lambda: ldgm.brick_haplo_graph(
                    bricked,
                    edge_weight_threshold=args.edge_weight_threshold,
                    make_sibs=args.make_sibs,
                    progress=False,
                ),
                args.batch_size,
            )
            rows.append(
                {
                    "implementation": "upstream_python",
                    "example": example,
                    "operation": op,
                    "iteration": iteration,
                    "elapsed_sec": elapsed,
                    "batch_size": args.batch_size,
                    "num_edges": int(bricked.num_edges),
                    "num_mutations": int(bricked.num_mutations),
                    "brick_graph_edges": int(brick_graph.number_of_edges()),
                    "reduced_edges": "",
                    "upstream_commit": commit,
                }
            )

            op, reduced, elapsed = timed(
                "reduce_graph",
                lambda: ldgm.reduce_graph(
                    brick_graph,
                    bricked,
                    path_weight_threshold=args.path_weight_threshold,
                    num_processes=1,
                    progress=False,
                ),
                args.batch_size,
            )
            rows.append(
                {
                    "implementation": "upstream_python",
                    "example": example,
                    "operation": op,
                    "iteration": iteration,
                    "elapsed_sec": elapsed,
                    "batch_size": args.batch_size,
                    "num_edges": int(bricked.num_edges),
                    "num_mutations": int(bricked.num_mutations),
                    "brick_graph_edges": int(brick_graph.number_of_edges()),
                    "reduced_edges": int(reduced.number_of_edges()),
                    "upstream_commit": commit,
                }
            )

            op, snplist, elapsed = timed(
                "make_snplist",
                lambda: ldgm.make_snplist(bricked),
                args.batch_size,
            )
            rows.append(
                {
                    "implementation": "upstream_python",
                    "example": example,
                    "operation": op,
                    "iteration": iteration,
                    "elapsed_sec": elapsed,
                    "batch_size": args.batch_size,
                    "num_edges": int(bricked.num_edges),
                    "num_mutations": int(bricked.num_mutations),
                    "brick_graph_edges": "",
                    "reduced_edges": "",
                    "upstream_commit": commit,
                }
            )

            op, final_ldgm, elapsed = timed(
                "make_ldgm",
                lambda: ldgm.make_ldgm(
                    constructor(),
                    path_weight_threshold=args.path_weight_threshold,
                    recombination_freq_threshold=None,
                    num_processes=1,
                    progress=False,
                )[0],
                args.batch_size,
            )
            rows.append(
                {
                    "implementation": "upstream_python",
                    "example": example,
                    "operation": op,
                    "iteration": iteration,
                    "elapsed_sec": elapsed,
                    "batch_size": args.batch_size,
                    "num_edges": int(bricked.num_edges),
                    "num_mutations": int(bricked.num_mutations),
                    "brick_graph_edges": "",
                    "reduced_edges": int(final_ldgm.number_of_edges()),
                    "upstream_commit": commit,
                }
            )

    fieldnames = [
        "implementation",
        "example",
        "operation",
        "iteration",
        "elapsed_sec",
        "batch_size",
        "num_edges",
        "num_mutations",
        "brick_graph_edges",
        "reduced_edges",
        "upstream_commit",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
