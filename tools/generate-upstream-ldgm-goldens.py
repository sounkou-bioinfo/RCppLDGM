#!/usr/bin/env python3
"""Generate upstream ldgm conformance artifacts.

This script intentionally uses the pinned upstream Python implementation in this
repository (the `ldgm/` package and `tests/utility_functions.py`) to produce
CSV/JSON artifacts for RcppLDGM conformance checks. These are not hand-authored
fixtures: every output records the upstream commit, constructor name, options,
and command.

Example:
    python tools/generate-upstream-ldgm-goldens.py --out inst/extdata/ldgm-goldens
"""

from __future__ import annotations

import argparse
import csv
import inspect
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


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
        msg += "\nInstall upstream ldgm test dependencies before generating goldens."
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


def edge_rows(graph) -> list[dict[str, Any]]:
    rows = []
    for source, target in graph.edges():
        rows.append(
            {
                "from": int(source),
                "to": int(target),
                "weight": float(graph.get_edge_data(source, target)["weight"]),
            }
        )
    rows.sort(key=lambda row: (row["from"], row["to"], row["weight"]))
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def snplist_input_rows(bricked) -> list[dict[str, Any]]:
    sites = {site.id: site for site in bricked.sites()}
    rows = []
    for mutation in bricked.mutations():
        site = sites[mutation.site]
        rows.append(
            {
                "site": int(site.id),
                "position": float(site.position),
                "ancestral_state": site.ancestral_state,
                "mutation": int(mutation.id),
                "derived_state": mutation.derived_state,
            }
        )
    rows.sort(key=lambda row: row["mutation"])
    return rows


def default_examples(utility_functions) -> list[str]:
    # Default goldens are upstream constructors that generate a valid mutated
    # tree sequence and pass the bricking -> brick graph -> reduction pipeline.
    # Historical draft names `two_tree_ts_n2` and `two_tree_ts_extra_length_n3`
    # are not constructors in the pinned upstream tests; the similarly named
    # `two_tree_ts` and `two_tree_ts_extra_length` contain no mutations and
    # upstream raises `ValueError("Tree sequence must contain mutations")`.
    preferred = [
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
    available = {
        name
        for name, value in vars(utility_functions).items()
        if callable(value) and not name.startswith("_")
    }
    return [name for name in preferred if name in available]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="inst/extdata/ldgm-goldens")
    parser.add_argument("--path-weight-threshold", type=float, default=100.0)
    parser.add_argument("--edge-weight-threshold", type=float, default=None)
    parser.add_argument("--make-sibs", action="store_true")
    parser.add_argument(
        "--examples",
        nargs="*",
        help="Specific tests.utility_functions constructors. Use 'all' for every callable constructor.",
    )
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[1]
    ldgm, utility_functions = require_upstream_modules(repo_root)

    if not args.examples:
        examples = default_examples(utility_functions)
    elif args.examples == ["all"]:
        examples = sorted(
            name
            for name, value in vars(utility_functions).items()
            if callable(value)
            and not name.startswith("_")
            and inspect.getmodule(value) is utility_functions
        )
    else:
        examples = args.examples

    out_dir = (repo_root / args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []

    commit = git_commit(repo_root)
    command = " ".join(sys.argv)

    for name in examples:
        constructor = getattr(utility_functions, name, None)
        if constructor is None or not callable(constructor):
            failures.append({"example": name, "error": "constructor not found"})
            if not args.continue_on_error:
                break
            continue

        try:
            ts = constructor()
            bricked = ldgm.brick_ts(ts, recombination_freq_threshold=None)
            brick_graph = ldgm.brick_haplo_graph(
                bricked,
                edge_weight_threshold=args.edge_weight_threshold,
                make_sibs=args.make_sibs,
                progress=False,
            )
            bricks_to_muts = ldgm.utility.get_mut_edges(bricked)
            reduced = ldgm.reduce_graph(
                brick_graph,
                bricked,
                path_weight_threshold=args.path_weight_threshold,
                num_processes=1,
                progress=False,
            )
            snplist = ldgm.make_snplist(bricked)

            prefix = out_dir / name
            brick_graph_file = prefix.with_suffix(".brick_graph.csv").name
            brick_map_file = prefix.with_suffix(".bricks_to_muts.csv").name
            reduced_file = prefix.with_suffix(".reduced_edgelist.csv").name
            snplist_file = prefix.with_suffix(".snplist.csv").name
            snplist_input_file = prefix.with_suffix(".snplist_input.csv").name
            metadata_file = prefix.with_suffix(".metadata.json").name

            write_csv(out_dir / brick_graph_file, edge_rows(brick_graph), ["from", "to", "weight"])
            write_csv(out_dir / reduced_file, edge_rows(reduced), ["from", "to", "weight"])
            write_csv(
                out_dir / brick_map_file,
                [
                    {"brick": int(brick), "mutation": int(mutations[0]), "mutations": ";".join(map(str, mutations))}
                    for brick, mutations in bricks_to_muts.items()
                ],
                ["brick", "mutation", "mutations"],
            )
            snplist.to_csv(out_dir / snplist_file, index=False)
            write_csv(
                out_dir / snplist_input_file,
                snplist_input_rows(bricked),
                ["site", "position", "ancestral_state", "mutation", "derived_state"],
            )

            metadata = {
                "example": name,
                "upstream_commit": commit,
                "command": command,
                "path_weight_threshold": args.path_weight_threshold,
                "edge_weight_threshold": args.edge_weight_threshold,
                "make_sibs": args.make_sibs,
                "num_sites": int(bricked.num_sites),
                "num_mutations": int(bricked.num_mutations),
                "num_edges": int(bricked.num_edges),
                "brick_graph_edges": int(brick_graph.number_of_edges()),
                "reduced_edges": int(reduced.number_of_edges()),
            }
            (out_dir / metadata_file).write_text(json.dumps(metadata, indent=2) + "\n")

            manifest_rows.append(
                {
                    "example": name,
                    "brick_graph": brick_graph_file,
                    "bricks_to_muts": brick_map_file,
                    "reduced_edgelist": reduced_file,
                    "snplist": snplist_file,
                    "snplist_input": snplist_input_file,
                    "metadata": metadata_file,
                    "path_weight_threshold": args.path_weight_threshold,
                    "edge_weight_threshold": "" if args.edge_weight_threshold is None else args.edge_weight_threshold,
                    "make_sibs": args.make_sibs,
                    "upstream_commit": commit,
                }
            )
            print(f"generated {name}")
        except Exception as exc:  # pragma: no cover - upstream dependent path
            failures.append({"example": name, "error": repr(exc)})
            print(f"failed {name}: {exc}", file=sys.stderr)
            if not args.continue_on_error:
                break

    write_csv(
        out_dir / "manifest.csv",
        manifest_rows,
        [
            "example",
            "brick_graph",
            "bricks_to_muts",
            "reduced_edgelist",
            "snplist",
            "snplist_input",
            "metadata",
            "path_weight_threshold",
            "edge_weight_threshold",
            "make_sibs",
            "upstream_commit",
        ],
    )
    if failures:
        write_csv(out_dir / "failures.csv", failures, ["example", "error"])

    if not manifest_rows:
        raise SystemExit("No goldens generated")
    if failures and not args.continue_on_error:
        raise SystemExit(f"Stopped after failure: {failures[-1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
