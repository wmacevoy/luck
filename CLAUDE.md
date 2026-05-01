# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

Research notes / book titled "Luck" by Warren D. MacEvoy on a probability-theory concept (a "luck" statistic built from a Mahalanobis-style radius for normal / chi-square / multinomial models). The repo contains:

1. A LaTeX (Tufte-book) manuscript at the repo root (`main.tex` and `\include`d chapters).
2. Three parallel numerical implementations of the same mathematical objects (`zluck`, `mnluck`, `chi2luck`, `mnprobln`, `mnsamp`, ...) in **Scilab** (`scilab/`), **R** (`R/`), and **Python** (`python/`), plus symbolic derivations in **Maxima** (`maxima/`).

Changes to a definition usually need to be reflected in more than one of those backends to keep them consistent — check the others before declaring a change done.

## Building the manuscript

`make all` (top-level `Makefile`) runs `makeindex main` (if `main.idx` exists) then `pdflatex main`. It does **not** loop pdflatex automatically — run `make all` two or three times for cross-references and the index to settle. Asymptote figures live in `graphics/*.asy` and build to `graphics/*.pdf` via the `graphs` target. `make clean` removes `main.pdf`, `main.log`, and editor backup files.

The chapters included from `main.tex` are: `introduction`, `normal`, `chi2`, `multinomial`, `computation`, `randomness`, `conclusion`, `proofs`, `cliffnotes`. `cliffnotes.tex` is referenced by `\include` but may not always be present — be aware when editing the include list.

A devcontainer (`.devcontainer/`) provides a Debian image with the TeX Live + Asymptote toolchain. `setup.sh` at the root is the apt-based install list for a bare Linux host (texlive-*, asymptote, scilab, wxmaxima).

## R / Jupyter environment (under `R/`)

The R directory is self-contained and uses a **local conda env at `R/.venv`** built from `R/environment.yml` (R 4.4 + tidyverse + IRkernel + jupyterlab + tectonic/latexmk). All shell helpers `cd` into their own directory first, so invoke them by path rather than copying commands into other directories.

- `R/setup.sh` — create or update `.venv`. Flags: `--restart` (recreate env), `--reset` (clear `data/`), `--debug` (also clear `debug/`, set `config.json.debug=true`).
- `R/run.sh` — launches `jupyter lab` inside the env (notebooks: `liars.ipynb`, `lucknorm.ipynb`, `uranium_eda.ipynb`).
- `R/R.sh` / `R/python.sh` / `R/jupyter.sh` — run `R` / `python` / `jupyter` inside `.venv` (thin wrappers that source `context.sh`).
- `R/testthat.sh` — runs `testthat::test_dir("tests/testthat", stop_on_failure=TRUE, stop_on_warning=TRUE)`. To run a single test file: `R/R.sh -e 'testthat::test_file("tests/testthat/test-mnluck.R")'`.
- `R/context.sh` — defines `conda_exe`, `conda_venv`, `python_exe`, `R_exe`, `jupyter_exe`, `config_json`. Source it (don't re-implement) when adding new helper scripts; it auto-discovers `conda` from `$HOME/miniforge3` or `$HOME/miniconda3`.

`R/bibtex` and `R/xelatex` are local shim binaries copied into `.venv/bin/` by `setup.sh` so notebooks/quarto can find a working TeX toolchain inside the env.

## Scilab sources (`scilab/`)

Files starting with `main_` (e.g. `main_chi2.sci`, `main_mn1.sci`) are runnable demo/figure-generation scripts that `exec()` the library `.sci` files; many write PDFs into `../img/`. Files starting with `test_` are test functions (return `%T` on success) but there is no test runner — execute them by `exec`'ing the file and calling the function in a Scilab session. `q-*.sci` files are exam/quiz problem scripts.

## Cross-language correspondences

When editing a numerical routine, the same name typically exists in multiple backends:

| concept | scilab | R | python | maxima |
|---|---|---|---|---|
| luck z-score | `zluck.sci` | `R/zluck.R` | (in `normal.py`) | — |
| multinomial luck | `mnluck.sci`, `mulluck.sci` | `R/mnluck.R` | — | `multinomial.maxima` |
| chi² luck | `chi2luck.sci` | (via `pchisq` in `mnluck.R`) | — | `chi2luck.wxm` |
| log-prob | `mnprobln.sci`, `mulprobln.sci` | `R/mnprobln.R` | — | — |
| sampler | `mnsamp.sci`, `mulsamp.sci` | `R/mnsamp.R` | — | — |

Keep formulas in sync across backends and across `proofs.tex` / `multinomial.tex` / `normal.tex` / `chi2.tex`.
