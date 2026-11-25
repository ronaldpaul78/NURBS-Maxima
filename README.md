# NURBS-Maxima

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24523237.svg)](https://doi.org/10.5281/zenodo.24523237)

**NURBS curves and surfaces for the computer algebra system [Maxima](https://maxima.sourceforge.io/)
and its interface [wxMaxima](https://wxmaxima-developers.github.io/wxmaxima/).**

NURBS-Maxima implements the classical algorithms of Piegl and Tiller (*The NURBS Book*) in pure
Maxima. The same code runs in **exact arithmetic** (rationals, radicals, symbolic weights), in
**double precision** and in **arbitrary precision** (bigfloats). This makes it possible to prove
identities symbolically (the NURBS circle satisfies `x² + y² = 1` exactly), to obtain curvatures
as closed-form functions of a design weight `w`, and to evaluate reliably near singular points
where double precision breaks down.

Version 2.0.0 · MIT licence · tested with Maxima 5.46.0 (GCL) and geomdl 5.5.1.

## Features by stage

| Stage | Content |
|---|---|
| 1. Mathematical core | Knot vectors (clamped, uniform, validation, multiplicity, continuity), data parametrisation, B-spline bases (A2.1, A2.2, A2.3, Cox–de Boor) and exact piecewise polynomials. |
| 2. NURBS curves | Evaluation (A4.1), rational derivatives (A4.2), exact piecewise expressions, weights (also symbolic), affine and projective transformations, knot insertion (A5.1), refinement (A5.4), degree elevation (A5.9), Bézier decomposition, splitting, reversal, global interpolation (A9.1, LU solver) and point projection (Newton). |
| 3. NURBS surfaces | Evaluation (A4.3), partial derivatives (A4.4), normals, isocurves, grids and triangle meshes; knot insertion (A5.3), refinement (A5.5), degree elevation and Bézier decomposition along u or v; point projection. Primitives: Bézier, bilinear, extrusion, ruled, revolution, torus, sphere, cylinder. |
| 4. Differential geometry | Curvature, torsion, Frenet frame and arc length of curves. One shared helper for fundamental forms, Gaussian, mean and principal curvatures; principal directions, point classification, singular points, area, volume, `∬K dA` (Gauss–Bonnet) and exact symbolic curvatures. |
| 5. Visualisation and export | `nurbs_plot` (inline in wxMaxima), composable `draw` objects, PNG/SVG/PDF/EPS images, OBJ, STL, VTK (with curvature fields for ParaView), CSV and geomdl-compatible JSON. |
| 6. Validation | rtest suite (95 problems), geomdl reference values, JSON cross-validation with geomdl, benchmarks, seven examples, API reference and article draft. |

## Layout

```
NURBS-Maxima/
├── nurbs_maxima.mac                the package (all stages)
├── examples/
│   ├── 01_bezier_curve.mac         basic: polynomial and rational Bézier, conics, degree elevation
│   ├── 02_circle.mac               basic: exact circle, weights, transformations
│   ├── 03_bezier_surface.mac       intermediate: bicubic patch, curvatures, refinement, projection
│   ├── 04_torus.mac                case study 1: Gauss–Bonnet, singular poles, float vs bfloat
│   ├── 05_modelling_operations.mac intermediate: insertion, refinement, elevation, decomposition
│   ├── 06_symbolic_weight_study.mac case study 2: curvature and area as functions of a weight w
│   └── 07_scientific_workflow.mac  case study 3: interpolation → revolution → K/H fields → VTK/STL
├── tests/
│   ├── rtest_nurbs.mac             95 problems, standard rtest format
│   ├── run_tests.mac               loads the package and runs the rtest suite
│   ├── accuracy_table.mac          recomputes the accuracy table of the paper
│   ├── test_nurbs.mac              rtest suite + comparison with geomdl references
│   ├── compare_geomdl.py           generates reference_geomdl.mac with geomdl
│   ├── reference_geomdl.mac        generated reference values
│   ├── crossval_export.mac         exports test objects to JSON
│   ├── crossval_geomdl.py          re-evaluates them in geomdl → crossval/report.md
│   ├── benchmark.mac               Maxima timings → benchmark_maxima.csv
│   └── benchmark_geomdl.py         geomdl timings → docs/benchmarks.md
├── docs/
│   ├── API.md                      function reference
│   ├── benchmarks.md               measured timings and application range
│   ├── paper.md                    article draft (English)
│   ├── articulo.md                 earlier Spanish draft (v1, superseded by paper.md)
│   └── programacion_nurbs.wxm      Spanish teaching notebook: how the package is programmed
├── README.md
└── LICENSE
```

## Installation and usage

Requirements: Maxima 5.40 or later and gnuplot for graphics (wxMaxima ships both). Python 3 with
`geomdl` is only needed to regenerate the validation data.

In wxMaxima:

```maxima
load("C:/path/NURBS-Maxima/nurbs_maxima.mac")$

C : nurbs_circle(1)$
ratsimp(nurbs_curve_point(C, 1/8)), algebraic;           /* exact: [1/sqrt(2), 1/sqrt(2)] */
X : nurbs_curve_span_expr(C, 2, t)$
ratsimp(X . X), algebraic;                               /* 1 */
nurbs_plot(C)$                                           /* inline plot */

T : nurbs_torus(2, 1/2)$
nurbs_gaussian_curvature(T, 0.1, 0.0);                   /* 0.8 (= 0.7999999999999985) */
nurbs_surface_area(T);                                   /* 4 pi^2 R r */
nurbs_export_vtk(T, "C:/temp/torus.vtk", 40, 20)$        /* open in ParaView */

/* symbolic weight */
assume(w > 0)$
Q : nurbs_bezier_curve([[-1,0],[0,1],[1,0]], [1, w, 1])$
ratsimp(nurbs_curve_curvature(Q, 1/2));                  /* w */

/* arbitrary precision */
fpprec : 50$
nurbs_curve_point(nurbs_bfloat(C), bfloat(1/7));
```

Run an example with `batch("C:/path/NURBS-Maxima/examples/04_torus.mac")$`. Examples save their
images and meshes next to the `.mac` file. Set `nurbs_show_plots : false$` before `batch` to skip
the plot windows.

For dense sampling convert first with `nurbs_float(obj)`; `nurbs_compile()` compiles the
basis-function kernels for a further 1.3–2.9x speed-up.

## Tests and validation

From the `NURBS-Maxima` folder:

```bash
maxima -b tests/run_tests.mac        # rtest suite: 95/95 pass (about 28 s)
maxima -b tests/test_nurbs.mac       # rtest suite + 14 checks against geomdl reference values
maxima -b tests/crossval_export.mac && python3 tests/crossval_geomdl.py   # JSON cross-validation
maxima -b tests/benchmark.mac && python3 tests/benchmark_geomdl.py       # timings
```

Inside wxMaxima the suite can also be run directly after loading the package:
`batch("C:/path/NURBS-Maxima/tests/rtest_nurbs.mac", test)$`.

Cross-validation with geomdl 5.5.1 over 16 curves and surfaces (points, first and second
derivatives, normals) gives a maximum absolute error of 2.3·10⁻¹³; knot insertion, refinement,
degree elevation and Bézier decomposition agree to 1.8·10⁻¹⁵ (see `tests/crossval/report.md`).

## Application range

NURBS-Maxima is **not** a meshing or CAD kernel. In float arithmetic it is 30–50x slower than
geomdl on CPython (`docs/benchmarks.md`). It is meant for:

* exact and symbolic analysis (closed-form curvature, parametric studies in a weight `w`);
* arbitrary-precision evaluation near singular or ill-conditioned points;
* verification of algorithms and reference values for other libraries;
* teaching and mathematical prototyping in wxMaxima;
* moderate meshes (10³–10⁴ points) exported to ParaView, slicers or Python.

## Conventions

Mathematical indices start at 0, as in the literature: `P_i` is `C@ctrlpts[i+1]`. The normal is
the normalised `S_u × S_v`; with the outward normal a convex surface has `H < 0`. One tolerance,
`nurbs_tol` (default 1e-10), governs every floating-point zero test. Details are in
[`docs/API.md`](docs/API.md).

## Known limitations

* Interpreted code: see the application range above.
* `nurbs_curvature_sym` with irrational weights (`√2/2`) produces very large expressions.
* No trimmed surfaces, surface–surface intersection, knot removal or surface fitting yet.

## Licence

MIT. See [LICENSE](LICENSE).

## References

* L. Piegl, W. Tiller. *The NURBS Book*, 2nd ed. Springer, 1997.
* O. R. Bingol, A. Krishnamurthy. NURBS-Python: An open-source object-oriented NURBS modeling
  framework in Python. *SoftwareX* 9 (2019) 85–94.

## How to cite

Santamaria Silupu, R. P., Castillo Maza, Y. A., Jimenez Vilcherrez, J. K., Urbina Guzman, R. T.,
Sosa Gonzales, V. B. *NURBS-Maxima: exact, symbolic and arbitrary-precision NURBS curves, surfaces
and differential geometry for Maxima* (v2.0.0). Zenodo, 2026.
https://doi.org/10.5281/zenodo.24523238

The DOI 10.5281/zenodo.24523237 always points to the latest version. See also `CITATION.cff`.
