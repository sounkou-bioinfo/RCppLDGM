# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
PKGNAME := $(shell sed -n 's/Package: *\([^ ]*\)/\1/p' DESCRIPTION)
PKGVERS := $(shell sed -n 's/Version: *\([^ ]*\)/\1/p' DESCRIPTION)


all: check



rd:
	R -e 'roxygen2::roxygenize(load_code = "source")'
vig:
	R -e "tools::buildVignettes(dir = '.')"
vig-md:
	R -e "for (f in Sys.glob('vignettes/*.Rmd')) { out <- sub('\\\\.Rmd$$', '.md', f); rmarkdown::render(f, output_format = rmarkdown::md_document(variant = 'gfm'), output_file = basename(out), output_dir = dirname(out), quiet = FALSE, envir = new.env(parent = globalenv())) }"
build:  install_deps
	R CMD build .

check: build
	R CMD check --as-cran --no-manual $(PKGNAME)_$(PKGVERS).tar.gz

install_deps:
	R \
	-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
	-e 'remotes::install_deps(dependencies = TRUE)'

install: build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz
install2:
	R CMD INSTALL --no-configure .

install3:
	R CMD INSTALL .
clean:
	@rm -rf $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck

UPSTREAM_PYTHON := $(CURDIR)/.sync/ldgm-python/bin/python

# Development targets
dev-install:
	R CMD INSTALL --preclean .
dev-install-debug-win:
	RTINYCC_DEBUG_BUILD=1 R CMD INSTALL --preclean .

test1:
	R -e "tinytest::test_package('$(PKGNAME)', testdir = 'inst/tinytest', ncpu=1L)"
test2:
	R -e "tinytest::test_package('$(PKGNAME)', testdir = 'inst/tinytest', ncpu=2L)"
test0:
	R -e "tinytest::test_package('$(PKGNAME)', testdir = 'inst/tinytest')"
test: install
	R -e "tinytest::test_package('$(PKGNAME)', testdir = 'inst/tinytest')"

warn-test: dev-install
	R -e "options(warn = 2); tinytest::test_package('$(PKGNAME)', testdir = 'inst/tinytest')"

upstream-python:
	tools/setup-upstream-python.sh

upstream-ldgm-goldens: upstream-python
	.sync/ldgm-python/bin/python tools/generate-upstream-ldgm-goldens.py --out .sync/ldgm-goldens

upstream-ldgm-conformance: dev-install upstream-ldgm-goldens
	RCPP_LDGM_UPSTREAM_GOLDENS=.sync/ldgm-goldens Rscript tools/check-upstream-ldgm-goldens.R

upstream-graphld-smoke:
	RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test Rscript tools/check-upstream-graphld-data.R

upstream-graphld-simulate-conformance: dev-install upstream-python
	RCPP_LDGM_PYTHON=$(UPSTREAM_PYTHON) \
	RCPP_LDGM_GRAPHLD_ROOT=.sync/graphld \
	RCPP_LDGM_GRAPHLD_SIM_METADATA=.sync/graphld/data/test/metadata.csv \
	RCPP_LDGM_GRAPHLD_POP=EUR \
	RCPP_LDGM_GRAPHLD_MAX_BLOCKS=2 \
	RCPP_LDGM_SIM_RANDOM_SEED=42 \
	Rscript tools/check-upstream-graphld-simulate.R

upstream-graphld-reader-conformance: dev-install upstream-python
	RCPP_LDGM_PYTHON=$(UPSTREAM_PYTHON) \
	RCPP_LDGM_GRAPHLD_DATA=.sync/graphld/data/test Rscript tools/check-upstream-graphld-readers.R

upstream-graphld-hdf5-interop: dev-install upstream-python
	RCPP_LDGM_PYTHON=$(UPSTREAM_PYTHON) \
	Rscript tools/check-graphld-hdf5-python-interop.R

upstream-conformance:
	make upstream-ldgm-conformance upstream-graphld-smoke upstream-graphld-reader-conformance upstream-graphld-hdf5-interop upstream-graphld-simulate-conformance

benchmark-precision: dev-install
	Rscript tools/benchmark-precision.R

benchmark-upstream-ldgm: dev-install upstream-ldgm-goldens
	RCPP_LDGM_UPSTREAM_GOLDENS=.sync/ldgm-goldens Rscript tools/benchmark-upstream-ldgm.R

rdm: install
	R -e "rmarkdown::render('README.Rmd')"
	perl -pi -e 's/[ \t]+$$//' README.md
.PHONY: all rd vig vig-md build check install_deps install clean dev-install dev-install-debug-win test0 test1 test2 test warn-test upstream-python upstream-ldgm-goldens upstream-ldgm-conformance upstream-graphld-simulate-conformance upstream-graphld-smoke upstream-graphld-reader-conformance upstream-graphld-hdf5-interop upstream-conformance benchmark-precision benchmark-upstream-ldgm rdm
