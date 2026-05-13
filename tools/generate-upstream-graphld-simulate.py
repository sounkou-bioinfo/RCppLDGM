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
  parser.add_argument("--populations", nargs="*", default=None)
  parser.add_argument("--chromosomes", nargs="*", default=None, type=int)
  parser.add_argument("--max-blocks", type=int, default=1)
  parser.add_argument("--sample-size", type=float, default=1000.0)
  parser.add_argument("--heritability", type=float, default=0.5)
  parser.add_argument("--component-variance", nargs="*", type=float, default=[1.0])
  parser.add_argument("--component-weight", nargs="*", type=float, default=[1.0])
  parser.add_argument("--alpha-param", type=float, default=-1)
  parser.add_argument("--random-seed", type=int, default=None)
  parser.add_argument("--annotation-columns", nargs="*", default=None)
  parser.add_argument("--fixture", default="default")
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


def normalize_populations(population: str, populations: list[str] | None) -> list[str]:
  values = populations if populations else [population]
  values = [str(value) for value in values if str(value)]
  if not values:
    raise ValueError("at least one population must be provided")
  return values


def populations_arg(populations: list[str]) -> str | list[str]:
  if len(populations) == 1:
    return populations[0]
  return populations


def write_toy_multi_population_fixture(workspace_dir: Path) -> Path:
  edge_lines = ["0,0,2", "1,1,3", "2,2,4", "0,1,0.1", "1,2,0.2"]
  for stem in ("eur_block", "afr_block"):
    population = "EUR" if stem == "eur_block" else "AFR"
    (workspace_dir / f"{stem}.{population}.edgelist").write_text("\n".join(edge_lines) + "\n", encoding="utf-8")

  write_csv(
    workspace_dir / "eur_block.snplist",
    [
      {"index": 0, "anc_alleles": "A", "deriv_alleles": "G", "af": 0.40, "site_ids": "rs1", "position": 10, "swap": "+"},
      {"index": 1, "anc_alleles": "C", "deriv_alleles": "T", "af": 0.20, "site_ids": "rs2", "position": 20, "swap": "+"},
      {"index": 2, "anc_alleles": "T", "deriv_alleles": "C", "af": 0.30, "site_ids": "rs3", "position": 30, "swap": "+"},
    ],
    ["index", "anc_alleles", "deriv_alleles", "af", "site_ids", "position", "swap"],
  )
  write_csv(
    workspace_dir / "afr_block.snplist",
    [
      {"index": 0, "anc_alleles": "A", "deriv_alleles": "G", "af": 0.10, "site_ids": "rs1", "position": 10, "swap": "+"},
      {"index": 1, "anc_alleles": "C", "deriv_alleles": "T", "af": 0.60, "site_ids": "rs2", "position": 20, "swap": "+"},
      {"index": 2, "anc_alleles": "T", "deriv_alleles": "C", "af": 0.80, "site_ids": "rs3", "position": 30, "swap": "+"},
    ],
    ["index", "anc_alleles", "deriv_alleles", "af", "site_ids", "position", "swap"],
  )
  metadata_path = workspace_dir / "multi.metadata.csv"
  write_csv(
    metadata_path,
    [
      {
        "chrom": 1,
        "chromStart": 1,
        "chromEnd": 100,
        "population": "EUR",
        "name": "eur_block.EUR.edgelist",
        "snplistName": "eur_block.snplist",
        "numVariants": 3,
        "numIndices": 3,
        "numEntries": 5,
        "info": "toy_multi_population",
      },
      {
        "chrom": 1,
        "chromStart": 101,
        "chromEnd": 200,
        "population": "AFR",
        "name": "afr_block.AFR.edgelist",
        "snplistName": "afr_block.snplist",
        "numVariants": 3,
        "numIndices": 3,
        "numEntries": 5,
        "info": "toy_multi_population",
      },
    ],
    ["chrom", "chromStart", "chromEnd", "population", "name", "snplistName", "numVariants", "numIndices", "numEntries", "info"],
  )
  return metadata_path


def main() -> int:
  args = parse_args()

  graphld_root = Path(args.graphld_root).resolve()
  graphld_module, read_metadata = load_graphld_module(graphld_root)
  command = " ".join(sys.argv)
  selected_populations = normalize_populations(args.population, args.populations)

  out_dir = Path(args.out).resolve()
  out_dir.mkdir(parents=True, exist_ok=True)
  workspace_dir = out_dir / f"{args.name}_workspace"
  workspace_dir.mkdir(parents=True, exist_ok=True)

  if args.fixture == "default":
    metadata_path = args.metadata or (graphld_root / "data/test/metadata.csv")
    metadata_path = Path(metadata_path).resolve()
    if not metadata_path.exists():
      raise RuntimeError(f"Metadata path not found: {metadata_path}")
    metadata_dir = metadata_path.parent
    metadata = read_metadata(
      str(metadata_path),
      populations=populations_arg(selected_populations),
      chromosomes=to_int_list(args.chromosomes),
      max_blocks=args.max_blocks,
    )

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
  elif args.fixture == "toy_multi_population":
    metadata_path = write_toy_multi_population_fixture(workspace_dir)
    metadata = read_metadata(
      str(metadata_path),
      populations=populations_arg(selected_populations),
      chromosomes=to_int_list(args.chromosomes),
      max_blocks=args.max_blocks,
    )
  else:
    raise ValueError(f"Unknown fixture: {args.fixture}")

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
    populations=populations_arg(selected_populations),
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
    "populations": selected_populations,
    "chromosomes": args.chromosomes,
    "max_blocks": args.max_blocks,
    "fixture": args.fixture,
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
        "populations": ",".join(selected_populations),
        "chromosomes": "" if args.chromosomes is None else ",".join(map(str, args.chromosomes)),
        "max_blocks": args.max_blocks,
        "fixture": args.fixture,
      }
    ],
    ["scenario", "output", "meta", "meta_json", "population", "populations", "chromosomes", "max_blocks", "fixture"],
  )

  return 0


if __name__ == "__main__":
  raise SystemExit(main())
