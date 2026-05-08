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
from pathlib import Path

import numpy as np


def _as_bool(value: str | None) -> bool:
    return str(value or "").lower() in {"1", "true", "yes", "y"}


def _require_or_skip(message: str) -> None:
    if _as_bool(os.environ.get("RCPP_LDGM_REQUIRE_H5PY")):
        raise RuntimeError(message)
    print(f"SKIP: {message}")
    raise SystemExit(0)


def _load_score_test_io(repo_root: Path):
    module_path = repo_root / ".sync" / "graphld" / "src" / "score_test" / "score_test_io.py"
    if not module_path.exists():
        _require_or_skip(f"upstream GraphLD score_test_io.py not found at {module_path}")

    spec = importlib.util.spec_from_file_location("graphld_upstream_score_test_io", module_path)
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


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("usage: check-graphld-hdf5-python-interop.py <file.h5> <trait>", file=sys.stderr)
        return 2

    try:
        import h5py  # noqa: F401  # imported to produce a clear optional-dependency check
    except Exception as exc:  # pragma: no cover - depends on local Python environment
        _require_or_skip(f"Python h5py is unavailable: {exc}")

    try:
        import polars  # noqa: F401
    except Exception as exc:  # pragma: no cover - depends on local Python environment
        _require_or_skip(f"Python polars is unavailable: {exc}")

    repo_root = Path(__file__).resolve().parents[1]
    score_test_io = _load_score_test_io(repo_root)

    hdf5_path = Path(argv[1])
    trait_name = argv[2]

    row_data = score_test_io.load_row_data(str(hdf5_path))
    expected_columns = {"CHR", "POS", "RSID", "jackknife_blocks"}
    if not expected_columns.issubset(set(row_data.columns)):
        raise AssertionError(f"row_data columns {row_data.columns} do not include {sorted(expected_columns)}")

    if _to_list(row_data["CHR"]) != [1, 1, 2]:
        raise AssertionError(f"unexpected CHR values: {_to_list(row_data['CHR'])}")
    if _to_list(row_data["POS"]) != [10, 20, 30]:
        raise AssertionError(f"unexpected POS values: {_to_list(row_data['POS'])}")
    if _to_list(row_data["RSID"]) != ["rs1", "rs2", "rs3"]:
        raise AssertionError(f"unexpected RSID values: {_to_list(row_data['RSID'])}")
    if _to_list(row_data["jackknife_blocks"]) != [0, 1, 1]:
        raise AssertionError(f"unexpected jackknife blocks: {_to_list(row_data['jackknife_blocks'])}")

    trait_names = score_test_io.get_trait_names(str(hdf5_path))
    if trait_name not in trait_names:
        raise AssertionError(f"trait {trait_name!r} not found in {trait_names}")

    trait_data = score_test_io.load_trait_hdf5(str(hdf5_path), trait_name)
    if "gradient" not in trait_data:
        raise AssertionError("trait data does not contain gradient")
    if "hessian" not in trait_data:
        raise AssertionError("trait data does not contain hessian")
    if "parameters" not in trait_data:
        raise AssertionError("trait data does not contain parameters group")
    np.testing.assert_allclose(trait_data["gradient"], np.array([0.1, -0.2, 0.3]), rtol=0, atol=1e-12)
    np.testing.assert_allclose(trait_data["hessian"], np.array([-0.01, -0.02, -0.03]), rtol=0, atol=1e-12)
    np.testing.assert_allclose(trait_data["parameters"]["parameters"], np.array([0.4, -0.5]), rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        trait_data["parameters"]["jackknife_parameters"],
        np.array([[0.41, -0.49], [0.39, -0.51]]),
        rtol=0,
        atol=1e-12,
    )

    print("GraphLD Python score_test_io HDF5 interop check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
