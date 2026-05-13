#!/usr/bin/env python3
"""Check RcppLDGM score-test HDF5 output with upstream GraphLD Python I/O.

This optional conformance check imports the pinned upstream GraphLD
``score_test_io.py`` from ``.sync/graphld`` and verifies that a small HDF5 file
written by ``ldgm_write_score_test_hdf5()`` can be read through GraphLD's h5py
path. If h5py is not installed the check skips by default; set
``RCPP_LDGM_REQUIRE_H5PY=1`` to make that a hard failure.
"""

from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

import numpy as np


def _as_bool(value: str | None) -> bool:
    return str(value or "").lower() in {"1", "true", "yes", "y"}


def _require_or_skip(message: str) -> None:
    if _as_bool(os.environ.get("RCPP_LDGM_REQUIRE_H5PY")):
        raise RuntimeError(message)
    print(f"SKIP: {message}")
    raise SystemExit(0)


def _score_test_module_dir(repo_root: Path) -> Path:
    return repo_root / ".sync" / "graphld" / "src" / "score_test"


def _load_score_test_io(repo_root: Path):
    module_path = _score_test_module_dir(repo_root) / "score_test_io.py"
    if not module_path.exists():
        _require_or_skip(f"upstream GraphLD score_test_io.py not found at {module_path}")

    spec = importlib.util.spec_from_file_location("graphld_upstream_score_test_io", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load import spec for {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_score_test(repo_root: Path):
    module_dir = _score_test_module_dir(repo_root)
    if not (module_dir / "score_test.py").exists():
        _require_or_skip(f"upstream GraphLD score_test.py not found at {module_dir}")
    sys.path.insert(0, str(module_dir))
    try:
        import score_test  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on optional Python env
        _require_or_skip(f"could not import upstream GraphLD score_test.py: {exc}")
    return score_test


def _load_convert_scores(repo_root: Path):
    module_dir = _score_test_module_dir(repo_root)
    module_path = module_dir / "convert_scores.py"
    if not module_path.exists():
        _require_or_skip(f"upstream GraphLD convert_scores.py not found at {module_path}")
    sys.path.insert(0, str(module_dir))
    spec = importlib.util.spec_from_file_location("graphld_upstream_convert_scores", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load import spec for {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _to_list(series_or_array):
    values = series_or_array.to_list() if hasattr(series_or_array, "to_list") else list(series_or_array)
    decoded = []
    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(value)
    return decoded


def _decode_attr_string(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _decode_groups(groups):
    return {
        key: [value.decode("utf-8") if isinstance(value, bytes) else value for value in values]
        for key, values in groups.items()
    }


def _normalize_list(values):
    normalized = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf-8")
        if isinstance(value, str):
            stripped = value.strip()
            if stripped.isdigit():
                normalized.append(int(stripped))
                continue
            normalized.append(stripped)
            continue
        if isinstance(value, (np.integer, int)):
            normalized.append(int(value))
            continue
        if isinstance(value, (np.floating, float)) and float(value).is_integer():
            normalized.append(int(value))
            continue
        normalized.append(value)
    return normalized


def _assert_gene_conversion(
    score_test_io,
    convert_scores,
    variant_hdf5: Path,
    gene_table_path: Path,
    r_gene_hdf5: Path,
    trait_name: str,
    second_trait: str,
):
    import h5py

    expected_groups = {"body": [trait_name], "combined": [trait_name, second_trait]}

    with tempfile.NamedTemporaryFile(prefix="graphld-upstream-gene-", suffix=".h5", delete=False) as handle:
        upstream_gene_hdf5 = Path(handle.name)
    try:
        convert_scores.convert_hdf5(str(variant_hdf5), str(upstream_gene_hdf5), str(gene_table_path), [1.0])

        r_row_data = score_test_io.load_row_data(str(r_gene_hdf5))
        upstream_row_data = score_test_io.load_row_data(str(upstream_gene_hdf5))
        expected_columns = {"CHR", "POS", "gene_id", "gene_name", "jackknife_blocks"}
        if not expected_columns.issubset(set(r_row_data.columns)):
            raise AssertionError(f"gene row_data columns {r_row_data.columns} do not include {sorted(expected_columns)}")
        if not expected_columns.issubset(set(upstream_row_data.columns)):
            raise AssertionError(
                f"upstream gene row_data columns {upstream_row_data.columns} do not include {sorted(expected_columns)}"
            )

        for column in ["CHR", "POS", "gene_id", "gene_name", "jackknife_blocks"]:
            r_values = _normalize_list(_to_list(r_row_data[column]))
            upstream_values = _normalize_list(_to_list(upstream_row_data[column]))
            if r_values != upstream_values:
                raise AssertionError(
                    f"R gene row_data column {column!r} does not match upstream: "
                    f"{r_values!r} vs {upstream_values!r}"
                )

        r_groups = _decode_groups(score_test_io.get_trait_groups(str(r_gene_hdf5)))
        upstream_groups = _decode_groups(score_test_io.get_trait_groups(str(upstream_gene_hdf5)))
        if r_groups != expected_groups:
            raise AssertionError(f"unexpected R gene trait groups: {r_groups!r}")
        if upstream_groups != expected_groups:
            raise AssertionError(f"unexpected upstream gene trait groups: {upstream_groups!r}")

        r_trait_names = score_test_io.get_trait_names(str(r_gene_hdf5))
        upstream_trait_names = score_test_io.get_trait_names(str(upstream_gene_hdf5))
        if r_trait_names != upstream_trait_names:
            raise AssertionError(f"gene trait names do not match upstream: {r_trait_names!r} vs {upstream_trait_names!r}")

        for name, expected_gradient in {
            trait_name: np.array([-0.1, 0.3]),
            second_trait: np.array([0.6, -0.1]),
        }.items():
            r_trait = score_test_io.load_trait_hdf5(str(r_gene_hdf5), name)
            upstream_trait = score_test_io.load_trait_hdf5(str(upstream_gene_hdf5), name)
            np.testing.assert_allclose(r_trait["gradient"], expected_gradient, rtol=0, atol=1e-12)
            np.testing.assert_allclose(r_trait["gradient"], upstream_trait["gradient"], rtol=0, atol=1e-12)

        r_gene_trait = score_test_io.load_trait_hdf5(str(r_gene_hdf5), trait_name)
        if "hessian" not in r_gene_trait:
            raise AssertionError("R gene trait data does not contain hessian")
        if "parameters" not in r_gene_trait:
            raise AssertionError("R gene trait data does not contain parameters group")
        if "posterior_scale" not in r_gene_trait:
            raise AssertionError("R gene trait data does not contain projected posterior_scale")
        np.testing.assert_allclose(r_gene_trait["hessian"], np.array([-0.03, -0.03]), rtol=0, atol=1e-12)
        np.testing.assert_allclose(r_gene_trait["posterior_scale"], np.array([2.5, 0.5]), rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            r_gene_trait["parameters"]["parameters"], np.array([0.4, -0.5]), rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            r_gene_trait["parameters"]["jackknife_parameters"],
            np.array([[0.41, -0.49], [0.39, -0.51]]),
            rtol=0,
            atol=1e-12,
        )

        with h5py.File(r_gene_hdf5, "r") as handle:
            if _decode_attr_string(handle.attrs["data_type"]) != "gene":
                raise AssertionError(f"unexpected gene data_type attribute: {handle.attrs['data_type']!r}")
            keys = [_decode_attr_string(value) for value in handle.attrs["keys"]]
            if keys != ["gene_id", "gene_name"]:
                raise AssertionError(f"unexpected gene keys attribute: {keys!r}")
    finally:
        upstream_gene_hdf5.unlink(missing_ok=True)


def main(argv: list[str]) -> int:
    if len(argv) not in {3, 5, 7, 9, 11}:
        print(
            "usage: check-graphld-hdf5-python-interop.py <file.h5> <trait> "
            "[expected-score.tsv expected-jackknife.tsv "
            "[surrogate-map.h5 block-name [second-trait expected-meta.tsv [gene-table.tsv gene-file.h5]]]]",
            file=sys.stderr,
        )
        return 2

    try:
        import h5py
    except Exception as exc:  # pragma: no cover - depends on local Python environment
        _require_or_skip(f"Python h5py is unavailable: {exc}")

    try:
        import polars as pl
    except Exception as exc:  # pragma: no cover - depends on local Python environment
        _require_or_skip(f"Python polars is unavailable: {exc}")

    repo_root = Path(__file__).resolve().parents[1]
    score_test_io = _load_score_test_io(repo_root)

    hdf5_path = Path(argv[1])
    trait_name = argv[2]

    with h5py.File(hdf5_path, "r") as handle:
        if _decode_attr_string(handle.attrs["metadata"]) != "":
            raise AssertionError(f"unexpected metadata attribute: {handle.attrs['metadata']!r}")
        if _decode_attr_string(handle.attrs["data_type"]) != "variant":
            raise AssertionError(f"unexpected data_type attribute: {handle.attrs['data_type']!r}")
        keys = [_decode_attr_string(value) for value in handle.attrs["keys"]]
        if keys != ["RSID", "POS"]:
            raise AssertionError(f"unexpected keys attribute: {keys!r}")
        if _decode_attr_string(handle.attrs["source"]) != "graphld-hdf5-interop":
            raise AssertionError(f"unexpected source attribute: {handle.attrs['source']!r}")

    row_data = score_test_io.load_row_data(str(hdf5_path))
    expected_columns = {"CHR", "POS", "RSID", "AF", "jackknife_blocks"}
    if not expected_columns.issubset(set(row_data.columns)):
        raise AssertionError(f"row_data columns {row_data.columns} do not include {sorted(expected_columns)}")

    if _to_list(row_data["CHR"]) != [1, 1, 2]:
        raise AssertionError(f"unexpected CHR values: {_to_list(row_data['CHR'])}")
    if _to_list(row_data["POS"]) != [10, 20, 30]:
        raise AssertionError(f"unexpected POS values: {_to_list(row_data['POS'])}")
    if _to_list(row_data["RSID"]) != ["rs1", "rs2", "rs3"]:
        raise AssertionError(f"unexpected RSID values: {_to_list(row_data['RSID'])}")
    np.testing.assert_allclose(row_data["AF"].to_numpy(), np.array([0.1, 0.2, 0.3]), rtol=0, atol=1e-12)
    if _to_list(row_data["jackknife_blocks"]) != [0, 1, 1]:
        raise AssertionError(f"unexpected jackknife blocks: {_to_list(row_data['jackknife_blocks'])}")

    trait_names = score_test_io.get_trait_names(str(hdf5_path))
    if trait_name not in trait_names:
        raise AssertionError(f"trait {trait_name!r} not found in {trait_names}")

    trait_groups = score_test_io.get_trait_groups(str(hdf5_path))
    decoded_groups = _decode_groups(trait_groups)
    if decoded_groups != {"body": [trait_name], "combined": [trait_name, trait_name + "_b"]}:
        raise AssertionError(f"unexpected trait groups: {trait_groups}")

    trait_data = score_test_io.load_trait_hdf5(str(hdf5_path), trait_name)
    if "gradient" not in trait_data:
        raise AssertionError("trait data does not contain gradient")
    if "hessian" not in trait_data:
        raise AssertionError("trait data does not contain hessian")
    if "parameters" not in trait_data:
        raise AssertionError("trait data does not contain parameters group")
    if "posterior_scale" not in trait_data:
        raise AssertionError("trait data does not contain posterior_scale")
    np.testing.assert_allclose(trait_data["gradient"], np.array([0.1, -0.2, 0.3]), rtol=0, atol=1e-12)
    np.testing.assert_allclose(trait_data["hessian"], np.array([-0.01, -0.02, -0.03]), rtol=0, atol=1e-12)
    np.testing.assert_allclose(trait_data["parameters"]["parameters"], np.array([0.4, -0.5]), rtol=0, atol=1e-12)
    np.testing.assert_allclose(trait_data["posterior_scale"], np.array([1.5, 1.0, 0.5]), rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        trait_data["parameters"]["jackknife_parameters"],
        np.array([[0.41, -0.49], [0.39, -0.51]]),
        rtol=0,
        atol=1e-12,
    )

    point_estimates = None
    jackknife_estimates = None
    score_test = None
    annotations = None
    if len(argv) in {5, 9, 11}:
        score_test = _load_score_test(repo_root)
        expected_score = pl.read_csv(argv[3], separator="\t")
        expected_jackknife = pl.read_csv(argv[4], separator="\t")
        annotations = pl.DataFrame(
            {
                "RSID": ["rs1", "rs2", "rs3"],
                "annot_a": [1.0, 0.0, 1.0],
                "annot_b": [0.0, 1.0, 1.0],
            }
        )
        trait_object = score_test_io.load_trait_data(str(hdf5_path), trait_name, row_data)
        annot_object = score_test.VariantAnnot(annotations, ["annot_a", "annot_b"])
        point_estimates, jackknife_estimates = score_test.run_score_test(trait_object, annot_object)
        np.testing.assert_allclose(
            point_estimates.ravel(), expected_score["score"].to_numpy(), rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(
            jackknife_estimates,
            expected_jackknife.select(["annot_a", "annot_b"]).to_numpy(),
            rtol=0,
            atol=1e-12,
        )
        z_scores = point_estimates.ravel() / np.std(jackknife_estimates, axis=0) / np.sqrt(jackknife_estimates.shape[0] - 1)
        np.testing.assert_allclose(z_scores, expected_score["z"].to_numpy(), rtol=0, atol=1e-12)

    if len(argv) in {7, 9, 11}:
        surrogate_path = Path(argv[5])
        block_name = argv[6]
        with h5py.File(surrogate_path, "r") as handle:
            if block_name not in handle:
                raise AssertionError(f"surrogate block {block_name!r} not found in {surrogate_path}")
            np.testing.assert_array_equal(handle[block_name][:], np.array([0, 2, -1]))

    if len(argv) in {9, 11}:
        if score_test is None or annotations is None or point_estimates is None or jackknife_estimates is None:
            raise AssertionError("meta-analysis check requires score-test estimates from the first trait")
        from meta_analysis import MetaAnalysis  # type: ignore

        second_trait = argv[7]
        expected_meta = pl.read_csv(argv[8], separator="\t")
        second_trait_object = score_test_io.load_trait_data(str(hdf5_path), second_trait, row_data)
        annot_object = score_test.VariantAnnot(annotations, ["annot_a", "annot_b"])
        second_point_estimates, second_jackknife_estimates = score_test.run_score_test(
            second_trait_object, annot_object
        )
        meta = MetaAnalysis()
        meta.update(point_estimates, jackknife_estimates)
        meta.update(second_point_estimates, second_jackknife_estimates)
        np.testing.assert_allclose(
            meta.point_estimates.ravel(), expected_meta["score"].to_numpy(), rtol=0, atol=1e-12
        )
        np.testing.assert_allclose(meta.z_scores.ravel(), expected_meta["z"].to_numpy(), rtol=0, atol=1e-12)

    if len(argv) == 11:
        gene_table_path = Path(argv[9])
        gene_hdf5_path = Path(argv[10])
        convert_scores = _load_convert_scores(repo_root)
        _assert_gene_conversion(
            score_test_io,
            convert_scores,
            hdf5_path,
            gene_table_path,
            gene_hdf5_path,
            trait_name,
            second_trait,
        )

    print("GraphLD Python score_test_io HDF5 interop check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
