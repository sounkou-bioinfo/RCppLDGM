#!/usr/bin/env python3
"""Generate pinned upstream GraphLD graphREML block metrics for conformance checks."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import polars as pl


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--graphld-root", default=".sync/graphld")
    parser.add_argument("--data-dir", default=".sync/graphld/data/test")
    parser.add_argument("--population", default="EUR")
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--num-iterations", type=int, default=3)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    graphld_root = Path(args.graphld_root).resolve()
    graphld_src = graphld_root / "src"
    if not graphld_src.exists():
        raise FileNotFoundError(f"GraphLD source directory not found: {graphld_src}")
    sys.path.insert(0, os.fspath(graphld_src))

    import graphld as gld  # type: ignore
    from graphld.cli import write_convergence_results, write_results, write_tall_results  # type: ignore
    from graphld.heritability import GraphREML  # type: ignore
    from graphld.io import load_ldgm, partition_variants, read_ldgm_metadata  # type: ignore

    data_dir = Path(args.data_dir).resolve()
    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    sumstats = gld.read_ldsc_sumstats(os.fspath(data_dir / "example.sumstats"))
    annotations = gld.load_annotations(
        os.fspath(data_dir / "annot"),
        chromosome=1,
        add_positions=False,
        exclude_bed=True,
    )

    merged_data = sumstats.join(annotations, on=["SNP"], how="right")
    if "CHR_right" in merged_data.columns and "POS_right" in merged_data.columns:
        merged_data = (
            merged_data
            .drop("CHR")
            .drop("POS")
            .rename({"CHR_right": "CHR", "POS_right": "POS"})
        )
    merged_data = merged_data.unique(subset=["SNP"], keep="first")

    metadata = read_ldgm_metadata(os.fspath(data_dir / "metadata.csv"), populations=args.population)
    block_data = partition_variants(metadata, merged_data)

    selected_block = None
    selected_name = None
    for row, block in zip(metadata.iter_rows(named=True), block_data, strict=False):
        if len(block) > 0:
            selected_block = block
            selected_name = row["name"]
            snplist_name = row["snplistName"]
            break
    if selected_block is None or selected_name is None:
        raise RuntimeError("No non-empty GraphREML block found in upstream test data")

    ldgm = load_ldgm(
        os.fspath(data_dir / selected_name),
        snplist_path=os.fspath(data_dir / snplist_name),
    )
    ldgm, Pz = GraphREML._initialize_block_zscores(ldgm, selected_block, ["base"], False, False)
    if Pz is None:
        raise RuntimeError("GraphREML block initialization unexpectedly returned no matched variants")

    sample_size = float(sumstats["N"].mean())
    Pz = Pz / np.sqrt(sample_size)
    ldgm.times_scalar(1.0 / sample_size)

    annotations_matrix = ldgm.variant_info.select(["base"]).to_numpy()
    params = np.zeros((1, 1), dtype=float)
    old_variant_h2 = np.zeros(len(ldgm.variant_info), dtype=float)
    likelihood, gradient, hessian, per_variant_h2 = GraphREML._compute_block_likelihood(
        ldgm=ldgm,
        Pz=Pz,
        annotations=annotations_matrix,
        params=params,
        link_fn_denominator=6e6,
        old_variant_h2=old_variant_h2,
        num_samples=100,
        likelihood_only=False,
        seed=args.seed,
    )

    pl.DataFrame(
        {
            "block_name": [Path(selected_name).stem],
            "sample_size": [sample_size],
            "seed": [int(args.seed)],
            "likelihood": [float(likelihood)],
            "gradient": [float(np.asarray(gradient, dtype=float).reshape(-1)[0])],
            "hessian": [float(np.asarray(hessian, dtype=float).reshape(-1)[0])],
            "n_hessian_rows": [int(np.asarray(hessian).shape[0])],
            "n_variant_rows": [int(len(per_variant_h2))],
            "n_active_indices": [int(len(Pz))],
            "per_variant_h2_sum": [float(np.asarray(per_variant_h2, dtype=float).sum())],
        }
    ).write_csv(out_dir / "reml_metrics.csv")
    pl.DataFrame({"per_variant_h2": np.asarray(per_variant_h2, dtype=float).reshape(-1)}).write_csv(out_dir / "per_variant_h2.csv")

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
        summary_stats=sumstats,
        annotation_data=annotations,
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
