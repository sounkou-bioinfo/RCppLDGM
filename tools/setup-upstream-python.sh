#!/usr/bin/env bash
set -euo pipefail

# Create an isolated Python environment for running pinned upstream ldgm
# conformance generators. Keep it under .sync/ so it is not part of the R
# package source and does not affect the system Python installation.
VENV_DIR="${RCPP_LDGM_UPSTREAM_PYTHON_VENV:-.sync/ldgm-python}"
PYTHON_BIN="${PYTHON:-python}"

"${PYTHON_BIN}" -m venv "${VENV_DIR}"
"${VENV_DIR}/bin/python" -m pip install --upgrade pip
"${VENV_DIR}/bin/python" -m pip install \
  networkx \
  msprime \
  tskit \
  scipy \
  polars \
  pandas \
  tqdm

"${VENV_DIR}/bin/python" - <<'PY'
modules = ['networkx', 'msprime', 'tskit', 'numpy', 'pandas', 'scipy', 'polars', 'tqdm']
for name in modules:
    module = __import__(name)
    print(f"{name} {getattr(module, '__version__', 'OK')}")
PY

printf '\nUpstream Python environment ready: %s\n' "${VENV_DIR}/bin/python"
