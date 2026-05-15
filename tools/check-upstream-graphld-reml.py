#!/usr/bin/env python3
"""Generate pinned upstream GraphLD graphREML outputs for conformance checks."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import polars as pl

FIXTURES = ("default", "synthetic_multiblock")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graphld-root", default=".sync/graphld")
    parser.add_argument("--data-dir", default=".sync/graphld/data/test")
    parser.add_argument("--population", default="EUR")
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--num-iterations", type=int, default=3)
    parser.add_argument("--fixture", choices=FIXTURES, default="default")
    return parser.parse_args()


def load_metadata(data_dir: Path, population: str) -> pl.DataFrame:
    from graphld.io import read_ldgm_metadata  # type: ignore

    metadata = read_ldgm_metadata(os.fspath(data_dir / "metadata.csv"), populations=population)
    if metadata.is_empty():
        raise RuntimeError(f"No metadata rows found for population {population!r}")
    return metadata.sort(["chrom", "chromStart"])


def load_default_inputs(data_dir: Path, population: str) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    import graphld as gld  # type: ignore

    metadata = load_metadata(data_dir, population)
    summary_stats = gld.read_ldsc_sumstats(os.fspath(data_dir / "example.sumstats"))
    annotation_data = gld.load_annotations(
        os.fspath(data_dir / "annot"),
        chromosome=1,
        add_positions=False,
        exclude_bed=True,
    )
    return metadata, summary_stats, annotation_data


def load_synthetic_multiblock_inputs(data_dir: Path, population: str) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    from graphld.io import load_ldgm  # type: ignore

    metadata = load_metadata(data_dir, population)
    frames: list[pl.DataFrame] = []
    for block_index, row in enumerate(metadata.iter_rows(named=True), start=1):
        ldgm = load_ldgm(
            os.fspath(data_dir / row["name"]),
            snplist_path=os.fspath(data_dir / row["snplistName"]),
        )
        variant_info = ldgm.variant_info.head(40)
        if variant_info.is_empty():
            continue
        sign = 1.0 if block_index % 2 == 1 else -1.0
        z_values = np.linspace(0.2, 1.0, variant_info.height, dtype=float) * sign
        frames.append(
            pl.DataFrame(
                {
                    "SNP": variant_info["site_ids"],
                    "CHR": pl.Series([int(row["chrom"])] * variant_info.height, dtype=pl.Int64),
                    "POS": variant_info["position"].cast(pl.Int64),
                    "A1": variant_info["deriv_alleles"],
                    "A2": variant_info["anc_alleles"],
                    "Z": z_values,
                    "N": pl.Series([100000.0] * variant_info.height, dtype=pl.Float64),
                }
            )
        )
    if not frames:
        raise RuntimeError("Synthetic multiblock graphREML fixture produced no summary-stat rows")
    summary_stats = pl.concat(frames, how="vertical")
    annotation_data = (
        summary_stats.select(["SNP", "CHR", "POS"])
        .unique(subset=["SNP"], keep="first")
        .with_columns(pl.lit(1.0).alias("base"))
    )
    return metadata, summary_stats, annotation_data


def load_fixture_inputs(
    fixture: str,
    data_dir: Path,
    population: str,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    if fixture == "default":
        return load_default_inputs(data_dir, population)
    if fixture == "synthetic_multiblock":
        return load_synthetic_multiblock_inputs(data_dir, population)
    raise ValueError(f"Unknown fixture: {fixture}")


def merge_summary_annotations(summary_stats: pl.DataFrame, annotation_data: pl.DataFrame) -> pl.DataFrame:
    merged_data = summary_stats.join(annotation_data, on=["SNP"], how="right")
    if "CHR_right" in merged_data.columns and "POS_right" in merged_data.columns:
        merged_data = (
            merged_data.drop("CHR")
            .drop("POS")
            .rename({"CHR_right": "CHR", "POS_right": "POS"})
        )
    return merged_data.unique(subset=["SNP"], keep="first")


def write_block_metrics(
    out_dir: Path,
    metadata: pl.DataFrame,
    merged_data: pl.DataFrame,
    data_dir: Path,
    seed: int,
    sample_size: float,
) -> None:
    from graphld.heritability import GraphREML  # type: ignore
    from graphld.io import load_ldgm, partition_variants  # type: ignore

    block_data = partition_variants(metadata, merged_data)
    metric_rows: list[dict[str, float | int | str]] = []
    h2_frames: list[pl.DataFrame] = []
    params = np.zeros((1, 1), dtype=float)
    for row, block in zip(metadata.iter_rows(named=True), block_data, strict=False):
        if len(block) == 0:
            continue
        ldgm = load_ldgm(
            os.fspath(data_dir / row["name"]),
            snplist_path=os.fspath(data_dir / row["snplistName"]),
        )
        ldgm, pz = GraphREML._initialize_block_zscores(ldgm, block, ["base"], False, False)
        if pz is None:
            continue
        pz = pz / np.sqrt(sample_size)
        ldgm.times_scalar(1.0 / sample_size)
        annotations_matrix = ldgm.variant_info.select(["base"]).to_numpy()
        old_variant_h2 = np.zeros(len(ldgm.variant_info), dtype=float)
        likelihood, gradient, hessian, per_variant_h2 = GraphREML._compute_block_likelihood(
            ldgm=ldgm,
            Pz=pz,
            annotations=annotations_matrix,
            params=params,
            link_fn_denominator=6e6,
            old_variant_h2=old_variant_h2,
            num_samples=100,
            likelihood_only=False,
            seed=seed,
        )
        block_name = Path(row["name"]).stem
        metric_rows.append(
            {
                "block_name": block_name,
                "sample_size": sample_size,
                "seed": int(seed),
                "likelihood": float(likelihood),
                "gradient": float(np.asarray(gradient, dtype=float).reshape(-1)[0]),
                "hessian": float(np.asarray(hessian, dtype=float).reshape(-1)[0]),
                "n_hessian_rows": int(np.asarray(hessian).shape[0]),
                "n_variant_rows": int(len(per_variant_h2)),
                "n_active_indices": int(len(pz)),
                "per_variant_h2_sum": float(np.asarray(per_variant_h2, dtype=float).sum()),
            }
        )
        h2_frames.append(
            pl.DataFrame(
                {
                    "block_name": [block_name] * len(per_variant_h2),
                    "variant_row": np.arange(1, len(per_variant_h2) + 1, dtype=int),
                    "per_variant_h2": np.asarray(per_variant_h2, dtype=float).reshape(-1),
                }
            )
        )
    if not metric_rows:
        raise RuntimeError("GraphREML fixture produced no non-empty blocks")
    pl.DataFrame(metric_rows).write_csv(out_dir / "reml_metrics.csv")
    pl.concat(h2_frames, how="vertical").write_csv(out_dir / "per_variant_h2.csv")


def main() -> None:
    args = parse_args()

    graphld_root = Path(args.graphld_root).resolve()
    graphld_src = graphld_root / "src"
    if not graphld_src.exists():
        raise FileNotFoundError(f"GraphLD source directory not found: {graphld_src}")
    sys.path.insert(0, os.fspath(graphld_src))

    import graphld as gld  # type: ignore
    from graphld.cli import write_convergence_results, write_results, write_tall_results  # type: ignore
    from graphld.io import partition_variants  # type: ignore

    data_dir = Path(args.data_dir).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata, summary_stats, annotation_data = load_fixture_inputs(args.fixture, data_dir, args.population)
    merged_data = merge_summary_annotations(summary_stats, annotation_data)
    sample_size = float(summary_stats["N"].mean())
    write_block_metrics(out_dir, metadata, merged_data, data_dir, args.seed, sample_size)

    model = gld.ModelOptions()
    method = gld.MethodOptions(
        num_iterations=args.num_iterations,
        run_serial=True,
        verbose=False,
        score_test_hdf5_file_name=os.fspath(out_dir / "reml_score.h5"),
        score_test_hdf5_trait_name="trait",
    )
    summary = gld.run_graphREML(
        model_options=model,
        method_options=method,
        summary_stats=summary_stats,
        annotation_data=annotation_data,
        ldgm_metadata_path=os.fspath(data_dir / "metadata.csv"),
        populations=args.population,
    )
    pl.DataFrame(
        {
            "parameter": [float(np.asarray(summary["parameters"], dtype=float).reshape(-1)[0])],
            "heritability": [float(np.asarray(summary["heritability"], dtype=float).reshape(-1)[0])],
            "enrichment": [float(np.asarray(summary["enrichment"], dtype=float).reshape(-1)[0])],
            "final_likelihood": [float(summary["log"]["final_likelihood"])],
            "trust_region_lambda": [float(np.asarray(summary["log"]["trust_region_lambdas"], dtype=float).reshape(-1)[0])],
            "converged": [bool(summary["log"]["converged"])],
            "num_iterations": [int(summary["log"]["num_iterations"])],
        }
    ).write_csv(out_dir / "reml_summary.csv")
    pl.DataFrame(
        {
            "iteration": np.arange(1, len(summary["likelihood_history"]) + 1),
            "likelihood": np.asarray(summary["likelihood_history"], dtype=float).reshape(-1),
            "trust_region_lambda": np.asarray(summary["log"]["trust_region_lambdas"], dtype=float).reshape(-1),
        }
    ).write_csv(out_dir / "reml_history.csv")
    pl.DataFrame(
        {
            "jackknife_block": np.arange(1, np.asarray(summary["jackknife_params"]).shape[0] + 1),
            "base": np.asarray(summary["jackknife_params"], dtype=float).reshape(-1),
        }
    ).write_csv(out_dir / "reml_jackknife_params.csv")
    pl.DataFrame(
        {
            "jackknife_block": np.arange(1, np.asarray(summary["jackknife_h2"]).shape[0] + 1),
            "base": np.asarray(summary["jackknife_h2"], dtype=float).reshape(-1),
        }
    ).write_csv(out_dir / "reml_jackknife_h2.csv")
    pl.DataFrame(
        {
            "jackknife_block": np.arange(1, np.asarray(summary["jackknife_enrichment"]).shape[0] + 1),
            "base": np.asarray(summary["jackknife_enrichment"], dtype=float).reshape(-1),
        }
    ).write_csv(out_dir / "reml_jackknife_enrichment.csv")

    parameter_values = np.asarray(summary["parameters"], dtype=float).reshape(-1)
    parameter_se = np.asarray(summary["parameters_se"], dtype=float).reshape(-1)
    parameter_log10pval = np.asarray(summary["parameters_log10pval"], dtype=float).reshape(-1)
    heritability_values = np.asarray(summary["heritability"], dtype=float).reshape(-1)
    heritability_se = np.asarray(summary["heritability_se"], dtype=float).reshape(-1)
    heritability_log10pval = np.asarray(summary["heritability_log10pval"], dtype=float).reshape(-1)
    enrichment_values = np.asarray(summary["enrichment"], dtype=float).reshape(-1)
    enrichment_se = np.asarray(summary["enrichment_se"], dtype=float).reshape(-1)
    enrichment_log10pval = np.asarray(summary["enrichment_log10pval"], dtype=float).reshape(-1)

    write_results(
        os.fspath(out_dir / "reml_parameters.csv"),
        parameter_values,
        parameter_se,
        parameter_log10pval,
        model.annotation_columns,
        "trait",
    )
    write_results(
        os.fspath(out_dir / "reml_heritability.csv"),
        heritability_values,
        heritability_se,
        heritability_log10pval,
        model.annotation_columns,
        "trait",
    )
    write_results(
        os.fspath(out_dir / "reml_enrichment.csv"),
        enrichment_values,
        enrichment_se,
        enrichment_log10pval,
        model.annotation_columns,
        "trait",
    )
    for trait_name in ("trait1", "trait2"):
        write_results(
            os.fspath(out_dir / "reml_parameters_multi.csv"),
            parameter_values,
            parameter_se,
            parameter_log10pval,
            model.annotation_columns,
            trait_name,
        )
        write_results(
            os.fspath(out_dir / "reml_heritability_multi.csv"),
            heritability_values,
            heritability_se,
            heritability_log10pval,
            model.annotation_columns,
            trait_name,
        )
        write_results(
            os.fspath(out_dir / "reml_enrichment_multi.csv"),
            enrichment_values,
            enrichment_se,
            enrichment_log10pval,
            model.annotation_columns,
            trait_name,
        )
    write_tall_results(os.fspath(out_dir / "reml_tall.csv"), model, summary)
    write_convergence_results(os.fspath(out_dir / "reml_convergence.csv"), summary)
    write_convergence_results(os.fspath(out_dir / "reml_convergence_multi.csv"), summary)
    write_convergence_results(os.fspath(out_dir / "reml_convergence_multi.csv"), summary)
    try:
        tall_multi_path = out_dir / "reml_tall_multi.csv"
        write_tall_results(os.fspath(tall_multi_path), model, summary)
        write_tall_results(os.fspath(tall_multi_path), model, summary)
    except Exception as exc:  # noqa: BLE001
        (out_dir / "reml_tall_multi_error.txt").write_text(f"{type(exc).__name__}: {exc}\n", encoding="utf-8")
    else:
        raise RuntimeError("expected GraphLD tall multi-write to fail on existing output file")


if __name__ == "__main__":
    main()
