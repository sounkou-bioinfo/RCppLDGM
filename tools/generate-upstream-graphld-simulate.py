#!/usr/bin/env python3
"""Generate upstream GraphLD simulation fixtures for RcppLDGM conformance checks.

This script intentionally uses the pinned upstream GraphLD checkout and writes a
small CSV artifact plus a manifest with exact command/context metadata.

Examples:
    python tools/generate-upstream-graphld-simulate.py \
      --graphld-root .sync/graphld \
      --metadata .sync/graphld/data/test/metadata.csv \
      --out .sync/graphld-simulate-goldens
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any


def git_commit(repo_root: Path) -> str:
  try:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo_root, text=True).strip()
  except Exception:
    return "unknown"


def load_graphld_module(graphld_root: Path):
  src = graphld_root / "src"
  if not src.exists():
    raise RuntimeError(f"GraphLD source directory not found at {src}")
  if str(src) not in sys.path:
    sys.path.insert(0, str(src))

  import graphld  # type: ignore
  from graphld.io import read_ldgm_metadata  # type: ignore

  return graphld, read_ldgm_metadata


def parse_args() -> argparse.Namespace:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--graphld-root", default=".sync/graphld")
  parser.add_argument("--metadata", default=None, help="Path to GraphLD metadata CSV")
  parser.add_argument("--out", required=True)
  parser.add_argument("--population", default="EUR")
  parser.add_argument("--chromosomes", nargs="*", default=None, type=int)
  parser.add_argument("--max-blocks", type=int, default=1)
  parser.add_argument("--sample-size", type=float, default=1000.0)
  parser.add_argument("--heritability", type=float, default=0.5)
  parser.add_argument("--component-variance", nargs="*", type=float, default=[1.0])
  parser.add_argument("--component-weight", nargs="*", type=float, default=[1.0])
  parser.add_argument("--alpha-param", type=float, default=-1)
  parser.add_argument("--random-seed", type=int, default=None)
  parser.add_argument("--annotation-columns", nargs="*", default=None)
  parser.add_argument("--name", default="default")
  return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
  import csv

  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


def to_int_list(values: list[int] | None) -> list[int] | None:
  if values is None:
    return None
  return [int(v) for v in values]


def main() -> int:
  args = parse_args()

  graphld_root = Path(args.graphld_root).resolve()
  metadata_path = args.metadata or (graphld_root / "data/test/metadata.csv")
  metadata_path = Path(metadata_path).resolve()

  if not metadata_path.exists():
    raise RuntimeError(f"Metadata path not found: {metadata_path}")

  graphld_module, read_metadata = load_graphld_module(graphld_root)
  command = " ".join(sys.argv)

  metadata = read_metadata(
    str(metadata_path),
    populations=args.population,
    chromosomes=to_int_list(args.chromosomes),
    max_blocks=args.max_blocks,
  )

  # Limit to the requested subset so generated fixtures stay small.
  out_dir = Path(args.out).resolve()
  out_dir.mkdir(parents=True, exist_ok=True)
  workspace_dir = out_dir / f"{args.name}_workspace"
  workspace_dir.mkdir(parents=True, exist_ok=True)
  metadata_dir = metadata_path.parent

  required_resources = set()
  required_resources.update(metadata.get_column("snplistName").to_list())
  required_resources.update(metadata.get_column("name").to_list())

  for resource in sorted(required_resources):
    source = metadata_dir / resource
    target = workspace_dir / resource
    if not source.exists():
      raise FileNotFoundError(f"Missing upstream resource: {source}")
    if not target.exists():
      target.write_bytes(source.read_bytes())

  filtered_metadata_path = workspace_dir / f"{args.name}.metadata.csv"
  metadata.write_csv(filtered_metadata_path)

  result = graphld_module.run_simulate(
    sample_size=args.sample_size,
    heritability=args.heritability,
    component_variance=args.component_variance,
    component_weight=args.component_weight,
    alpha_param=args.alpha_param,
    random_seed=args.random_seed,
    annotation_columns=args.annotation_columns,
    ldgm_metadata_path=str(filtered_metadata_path),
    populations=args.population,
    chromosomes=to_int_list(args.chromosomes),
    run_in_serial=True,
    verbose=False,
  )
  output_csv = out_dir / f"{args.name}.simulate.csv"
  output_json = out_dir / f"{args.name}.simulate.json"

  result.write_csv(output_csv)

  metadata_payload = {
    "scenario": args.name,
    "upstream_commit": git_commit(graphld_root),
    "command": command,
    "sample_size": args.sample_size,
    "heritability": args.heritability,
    "component_variance": args.component_variance,
    "component_weight": args.component_weight,
    "alpha_param": args.alpha_param,
    "random_seed": args.random_seed,
    "annotation_columns": args.annotation_columns,
    "population": args.population,
    "chromosomes": args.chromosomes,
    "max_blocks": args.max_blocks,
    "run_in_serial": True,
    "num_rows": int(result.height),
    "upstream_metadata_path": str(filtered_metadata_path),
  }
  output_json.write_text(__import__("json").dumps(metadata_payload, indent=2), encoding="utf-8")

  write_csv(
    out_dir / "manifest.csv",
    [
      {
        "scenario": args.name,
        "output": output_csv.name,
        "meta": str(filtered_metadata_path.relative_to(out_dir)),
        "meta_json": output_json.name,
        "population": args.population,
        "chromosomes": "" if args.chromosomes is None else ",".join(map(str, args.chromosomes)),
        "max_blocks": args.max_blocks,
      }
    ],
    ["scenario", "output", "meta", "meta_json", "population", "chromosomes", "max_blocks"],
  )

  return 0


if __name__ == "__main__":
  raise SystemExit(main())
