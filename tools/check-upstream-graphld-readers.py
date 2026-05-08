#!/usr/bin/env python3
"""Generate pinned GraphLD reader outputs for RcppLDGM conformance checks.

This intentionally imports upstream modules by file path so the checks can run
without importing graphld.__init__ and its optional SuiteSparse dependency.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import numpy as np
import polars as pl


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load module spec for {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_gene_sets(gmt_file: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    with gmt_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                gene_sets[fields[0]] = fields[2:]
    return gene_sets


def load_gene_table_like_graphld(gene_table_path: Path, chromosomes: list[int] | None = None) -> pl.DataFrame:
    schema = {
        "gene_id": pl.Utf8,
        "gene_id_version": pl.Utf8,
        "gene_name": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "CHR": pl.Utf8,
    }
    gene_table = (
        pl.scan_csv(gene_table_path, schema=schema, separator="\t", has_header=True)
        .filter(pl.col("CHR").is_in([str(i) for i in range(1, 23)]))
        .filter(pl.col("gene_id").is_not_null())
        .with_columns(pl.col("gene_name").fill_null("NA"))
        .with_columns(((pl.col("start") + pl.col("end")) / 2).alias("midpoint"))
        .with_columns(pl.col("midpoint").alias("POS"))
        .sort(pl.col("CHR").cast(pl.Int64), "midpoint")
        .collect()
    )
    if chromosomes:
        gene_table = gene_table.filter(pl.col("CHR").cast(pl.Int64).is_in(chromosomes))
    return gene_table.with_columns(pl.col("midpoint").cast(pl.Int64).alias("POS"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graphld-root", default=".sync/graphld")
    parser.add_argument("--data-dir", default=".sync/graphld/data/test")
    parser.add_argument("--out", required=True)
    parser.add_argument("--num-vcf-rows", type=int, default=25)
    args = parser.parse_args()

    graphld_root = Path(args.graphld_root).resolve()
    data_dir = Path(args.data_dir).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    io_mod = load_module("graphld_io_file", graphld_root / "src" / "graphld" / "io.py")
    vcf_mod = load_module("graphld_vcf_file", graphld_root / "src" / "graphld" / "vcf_io.py")
    parquet_mod = load_module("graphld_parquet_file", graphld_root / "src" / "graphld" / "parquet_io.py")
    genesets_mod = load_module("graphld_score_genesets_file", graphld_root / "src" / "score_test" / "genesets.py")

    vcf = vcf_mod.read_gwas_vcf(str(data_dir / "example.gwas.vcf"), num_rows=args.num_vcf_rows)
    vcf.select(["CHR", "POS", "ID", "REF", "ALT", "ES", "SE", "LP", "AF", "Z"]).write_csv(out_dir / "vcf.csv")

    parquet_file = data_dir / "example_multi_trait.parquet"
    traits = parquet_mod.get_parquet_traits(parquet_file)
    (out_dir / "parquet_traits.txt").write_text("\n".join(traits) + "\n", encoding="utf-8")
    for trait in traits:
        parquet_mod.read_parquet_sumstats(parquet_file, trait=trait).write_csv(out_dir / f"parquet_{trait}.csv")

    bed = io_mod.read_bed(str(data_dir / "annot" / "test_regions.bed"))
    bed.write_csv(out_dir / "bed.csv")

    annotations = io_mod.load_annotations(
        str(data_dir / "annot"),
        chromosome=1,
        add_positions=False,
        exclude_bed=True,
    )
    annotations.write_csv(out_dir / "annotations.csv")

    gene_table = load_gene_table_like_graphld(graphld_root / "tests" / "score_test_data" / "genes_test.tsv", chromosomes=[22])
    gene_sets = read_gene_sets(graphld_root / "tests" / "score_test_data" / "test_symbols.gmt")
    variant_table = pl.DataFrame(
        {
            "CHR": [22, 22, 22, 22],
            "POS": [15281327, 30265000, 39700000, 50799284],
            "RSID": ["stub1", "stub2", "stub3", "stub4"],
        }
    )
    weights = np.array([1.0, 0.5], dtype=float)
    matrix = genesets_mod.gene_variant_matrix(variant_table, gene_table, weights).toarray()
    pl.DataFrame(matrix, schema=[f"gene_{i}" for i in range(matrix.shape[1])]).write_csv(out_dir / "gene_variant_matrix.csv")
    genesets_mod.convert_gene_set_to_gene_annotations(gene_sets, gene_table).write_csv(out_dir / "gene_set_annotations.csv")


if __name__ == "__main__":
    main()
