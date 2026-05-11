#!/usr/bin/env bash
set -euo pipefail

# Create an isolated Python environment for running pinned upstream ldgm
# conformance generators. Keep it under .sync/ so it is not part of the R
# package source and does not affect the system Python installation.
VENV_DIR="${RCPP_LDGM_UPSTREAM_PYTHON_VENV:-.sync/ldgm-python}"
PYTHON_BIN="${PYTHON:-python}"
SCIPYSPARSE_PIN="${RCPP_LDGM_SCIPYSPARSE_PIN:-scikit-sparse>=0.4.12,<0.5.0}"

if [ -d "${VENV_DIR}" ]; then
  rm -rf "${VENV_DIR}"
fi

if command -v uv >/dev/null 2>&1; then
  uv venv "${VENV_DIR}" --python "${PYTHON_BIN}"
  uv pip install --python "${VENV_DIR}/bin/python" \
    networkx \
    msprime \
    tskit \
    scipy \
    polars \
    pandas \
    tqdm \
    click \
    h5py \
    filelock \
    "${SCIPYSPARSE_PIN}"
else
  "${PYTHON_BIN}" -m venv "${VENV_DIR}"
  "${VENV_DIR}/bin/python" -m pip install --upgrade pip
  "${VENV_DIR}/bin/python" -m pip install \
    networkx \
    msprime \
    tskit \
    scipy \
    polars \
    pandas \
    tqdm \
    click \
    h5py \
    filelock \
    "${SCIPYSPARSE_PIN}"
fi

"${VENV_DIR}/bin/python" - <<'PY'
modules = ['networkx', 'msprime', 'tskit', 'scipy', 'polars', 'pandas', 'tqdm', 'click', 'h5py', 'filelock', 'numpy', 'sksparse']
for name in modules:
    module = __import__(name)
    print(f"{name} {getattr(module, '__version__', 'OK')}")
PY

printf '\nUpstream Python environment ready: %s\n' "${VENV_DIR}/bin/python"
