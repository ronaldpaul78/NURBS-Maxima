# NURBS-Maxima — SoftwareX manuscript (LaTeX, Overleaf-ready)

## Overleaf
1. Zip this folder (`docs/softwarex/`) and upload it as a new project (New Project → Upload Project).
2. Compiler: pdfLaTeX. Main document: `main.tex`. `elsarticle.cls` and `elsarticle-num.bst` are already in Overleaf.
3. `main.pdf` in this folder is a local build for preview only.

## Files
| file | content |
|---|---|
| `main.tex` | manuscript (SoftwareX Original Software Publication structure) |
| `references.bib` | bibliography (every entry checked against Crossref, arXiv or a publisher catalogue on 2026-11-25) |
| `bench_table.tex` | Table 4 (benchmarks), from `tests/benchmark.mac` + `tests/benchmark_geomdl.py` |
| `figures/make_figures.mac` | regenerates every figure and figure-data file |
| `figures/torus.pdf`, `figures/vase.pdf` | Fig. 2 |
| `figures/pole_sweep.csv` | Fig. 3a |
| `figures/vase_meridian.csv` | Fig. 3b |
| `figures/weight_study.csv` | Fig. 4 |
| `figures/benchmarks_run.md` | raw output of the benchmark run used in Table 4 |

Regenerate figures: `cd docs/softwarex/figures && maxima -b make_figures.mac` (~15 s).

## Reproducing the numbers of the paper
Run from the repository root (`NURBS-Maxima/`). On Windows use the `maxima` command installed with
wxMaxima, or open each file in wxMaxima and run it with `batch("...")`.

| paper item | command | expected result |
|---|---|---|
| 95 regression tests (Section 3, Validation) | `maxima -b tests/run_tests.mac` | `95/95 tests passed` |
| Table 3 (accuracy) | `maxima -b tests/accuracy_table.mac` | the eight errors of Table 3 |
| Cross-validation with geomdl (2.3e-13, degree-elevation finding) | `pip install geomdl==5.5.1`, then `maxima -b tests/crossval_export.mac` and `python tests/crossval_geomdl.py` | `RESULT: PASSED`, report in `tests/crossval/report.md` |
| Example 1 (sphere pole, Gauss-Bonnet) | `maxima -b examples/04_torus.mac` | values quoted in Example 1 |
| Example 2 (weight w) | `maxima -b examples/06_symbolic_weight_study.mac` | `kappa(1/2) = w`, `K = w^2/4` |
| Example 3 (vase) | `maxima -b examples/07_scientific_workflow.mac` | `int K dA = 1.1192965` vs boundary term `1.1192980` |
| Figures 2-4 | `cd docs/softwarex/figures` and `maxima -b make_figures.mac` | PDF and CSV files of the figures |
| Table 4 (benchmarks) | `maxima -b tests/benchmark.mac` and `python tests/benchmark_geomdl.py` | `docs/benchmarks.md` (times depend on the machine) |

## Before submission
* Check the current SoftwareX Guide for Authors / OSP template (the official page could not be read automatically).
* Delete the `\todo` macro definition at the top of `main.tex`: any placeholder left will then stop compilation.
