# NURBS-Maxima — SoftwareX manuscript (LaTeX, Overleaf-ready)

## Overleaf
1. Zip this folder (`docs/softwarex/`) and upload it as a new project (New Project → Upload Project).
2. Compiler: pdfLaTeX. Main document: `main.tex`. `elsarticle.cls` and `elsarticle-num.bst` are already in Overleaf.
3. `main.pdf` in this folder is a local build for preview only.

## Files
| file | content |
|---|---|
| `main.tex` | manuscript (SoftwareX Original Software Publication structure) |
| `references.bib` | bibliography; entries marked `TO VERIFY` print a red note |
| `bench_table.tex` | Table 4 (benchmarks), from `tests/benchmark.mac` + `tests/benchmark_geomdl.py` |
| `figures/make_figures.mac` | regenerates every figure and figure-data file |
| `figures/torus.pdf`, `figures/vase.pdf` | Fig. 2 |
| `figures/pole_sweep.csv` | Fig. 3a |
| `figures/vase_meridian.csv` | Fig. 3b |
| `figures/weight_study.csv` | Fig. 4 |
| `figures/benchmarks_run.md` | raw output of the benchmark run used in Table 4 |

Regenerate figures: `cd docs/softwarex/figures && maxima -b make_figures.mac` (~15 s).

## Before submission
* Replace every red `[...]` placeholder (authors, affiliations, e-mail, repository URL/DOI, CRediT, AI declaration, acknowledgements, impact/users).
* Verify the cells marked † in Table 2 and the `TO VERIFY` references.
* Check the current SoftwareX Guide for Authors / OSP template (the official page could not be read automatically).
* Repeat the benchmarks on the final machine (Table 4).
* Delete the `\todo` macro definition at the top of `main.tex`: any placeholder left will then stop compilation.
