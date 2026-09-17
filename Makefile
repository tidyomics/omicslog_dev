SHELL := bash

.PHONY: help r-restore r-snapshot r-init r-check py-install py-lock py-run py-check

help:
	@echo "Common commands:"
	@echo "  make r-init        - Initialize renv (run inside R container)"
	@echo "  make r-restore     - Restore R deps from renv.lock"
	@echo "  make r-snapshot    - Snapshot R deps to renv.lock"
	@echo "  make r-check       - Run R smoke test"
	@echo "  make py-install    - Poetry install (sync venv)"
	@echo "  make py-lock       - Poetry lock"
	@echo "  make py-run ARGS='...' - Run python with Poetry"
	@echo "  make py-check      - Run Python smoke test"

# -------- R --------

r-init:
	R -q -e "renv::init(bare = TRUE)"

r-restore:
	R -q -e "renv::restore(prompt = FALSE)"

r-snapshot:
	R -q -e "renv::snapshot(prompt = FALSE)"

r-check:
	R -q -f analysis/r/smoke.R

# -------- Python --------

py-install:
	poetry install

py-lock:
	poetry lock

py-run:
	poetry run python $(ARGS)

py-check:
	poetry run python analysis/py/smoke.py

.PHONY: q-render q-preview q-check rmd-render

q-check:
	quarto check

q-render:
	quarto render

q-preview:
	quarto preview --no-watch-inputs

# If you still render legacy .Rmd directly (optional)
rmd-render:
	R -q -e "rmarkdown::render(commandArgs(TRUE)[1])"

