# NURBS-Maxima: exact, symbolic and arbitrary-precision NURBS geometry and differential geometry in a free computer algebra system

**Draft v0.2 (2026-11-25).** Authors: *[to complete]*. Affiliation: *[to complete]*. Corresponding author: *[to complete]*.

> Notes for the authors (remove before submission)
> * Every number in Sections 5 and 6 is produced by a file in the repository; the file is named next to each table. They were measured with Maxima 5.46.0 on GCL 2.6.14, geomdl 5.5.1 on CPython 3.13.18, Intel Xeon @ 2.30 GHz (cloud container, 4 cores). Repeat the runs on the final machine and, ideally, on Maxima/SBCL.
> * Cells marked † in Table 1 and the references marked *[to verify]* describe third-party products from general knowledge. Check them against current documentation before submission.
> * The structure follows a software paper for a geometric-modelling or scientific-software journal (e.g. *Computer-Aided Design*, *Computer Aided Geometric Design*, *SoftwareX*, *ACM TOMS*). Adapt length and sections to the target journal.

## Abstract

Non-uniform rational B-splines (NURBS) are the standard representation of geometry in computer-aided design and isogeometric analysis. Open libraries implement them in floating-point arithmetic only, and the computer algebra systems that offer splines are proprietary or lack rational surfaces and their differential geometry. We present **NURBS-Maxima**, an open-source (MIT) package for the free computer algebra system Maxima that implements the algorithms of Piegl and Tiller for curves and surfaces—evaluation, rational derivatives, knot insertion and refinement, degree elevation, Bézier decomposition, global interpolation and point inversion—together with an integrated differential-geometry layer (fundamental forms, Gaussian, mean and principal curvatures, singular points, area, volume and total curvature) and export to VTK, STL, OBJ, CSV and geomdl JSON. A single code base runs in exact rational/radical arithmetic, with symbolic weights, in double precision and in bigfloat arithmetic; the parameter type selects the arithmetic. We validate the package with 95 regression problems, analytical identities proved symbolically, 14 comparisons with reference values from NURBS-Python (geomdl) and a JSON round-trip cross-validation over 16 objects with a maximum absolute discrepancy of 2.3·10⁻¹³ in points, first and second derivatives and normals. Three case studies show what the combination enables: Gauss–Bonnet verification and the behaviour of curvature at the singular poles of a NURBS sphere, where double precision fails (error 38 at 10⁻⁹ from the pole) while 50-digit arithmetic keeps the error below 10⁻³³; closed-form curvature of conics and patches as functions of a design weight `w`; and a design-to-visualisation workflow ending in ParaView and 3D-printing files. Float evaluation is 35–50 times slower than geomdl, so we state the intended range explicitly: exact and symbolic analysis, high-precision reference computations, verification and teaching, and meshes up to about 10⁴ points.

**Keywords:** NURBS; B-splines; computer algebra; Maxima; differential geometry; Gaussian curvature; arbitrary precision; reproducible software.

## 1. Introduction

NURBS combine the local control of B-splines with the exact representation of conics and quadrics through rational weights [1, 2]. They are the geometric basis of CAD exchange standards and of isogeometric analysis [6]. Mature open implementations—openNURBS, SISL, OpenCASCADE, verb and NURBS-Python/geomdl [3]—are numerical libraries: every operation is carried out in double precision, and quantities such as curvature are evaluated pointwise.

Three tasks need more than that.

1. **Symbolic verification.** Proving an identity—that the quadratic rational circle has constant curvature, that a refined surface is the same surface, that a patch lies on a quadric—requires exact arithmetic rather than agreement to machine precision.
2. **Closed-form design sensitivities.** A designer who varies a weight `w` wants the curvature or the area as a function of `w`, not a table of samples. This needs the NURBS algorithms to run on symbolic data.
3. **Ill-conditioned evaluation.** Near degenerate points (the poles of a sphere of revolution, collapsed edges, cusps) the formulas for `K` and `H` divide quantities that both tend to zero. Double precision loses all significant digits; arbitrary precision recovers them.

Computer algebra systems are the natural environment for these tasks, but their spline support is either proprietary (Mathematica, Maple, MATLAB toolboxes) or, in the case of free systems, absent: to our knowledge *[to verify with a systematic search]* there is no complete NURBS package for Maxima [7], a free descendant of Macsyma widely used in teaching through its notebook interface wxMaxima.

The contributions of this work are:

* a pure-Maxima implementation of the core NURBS algorithms for curves and surfaces, including the modelling operations (A5.1, A5.3–A5.6, A5.9), global interpolation with an LU solver and Newton point inversion, that works unchanged in exact, symbolic, double and bigfloat arithmetic;
* an integrated differential-geometry layer built on one shared evaluation of the second-order derivatives, with consistent singular-point handling and global integrals;
* a reproducible validation protocol: Maxima's standard `rtest` suite, symbolic identities, geometric invariants and a two-way cross-validation with geomdl through the package's own JSON export;
* honest performance measurements, including Maxima's `compile`, and an explicit statement of the package's range of application;
* three case studies that exercise the combination of symbolic, high-precision and visual output.

## 2. State of the art

Table 1 compares NURBS-Maxima with the tools that a researcher would consider for the same tasks. "Symbolic" means that the NURBS evaluation itself accepts symbolic data (weights, control points or parameters) and returns closed-form expressions. "Differential geometry" means built-in Gaussian/mean curvature of NURBS surfaces.

**Table 1.** Comparison matrix. † = from general knowledge of the product, to verify against current documentation.

| | NURBS-Maxima | geomdl (NURBS-Python) [3] | MATLAB NURBS Toolbox [10] / Curve Fitting Toolbox | Octave `nurbs` package [11] | Mathematica `BSplineFunction` [12] | Maple `CurveFitting` [13] |
|---|---|---|---|---|---|---|
| Licence | MIT, free host (Maxima, GPL) | MIT | toolbox free/GPL†, host proprietary | GPL-3† | proprietary | proprietary |
| Rational curves and surfaces | yes | yes | yes | yes | yes (weights)† | B-spline basis only, no NURBS objects† |
| Exact arithmetic (rationals, radicals) | yes, all algorithms | no | no (numeric)† | no | partial: basis functions symbolic, spline objects numeric† | basis functions symbolic† |
| Symbolic weights / closed-form curvature | yes | no | no | no | no built-in† | no built-in† |
| Arbitrary precision | yes (bigfloat) | no | no (Symbolic Toolbox separately)† | no | yes (general)† | yes (general)† |
| Knot insertion / refinement / degree elevation | yes / yes / yes | yes / yes / yes (elevation intended for Bézier, see §5.3) | yes† | yes / yes / yes† | not exposed† | no† |
| Bézier decomposition, point inversion | yes, yes | yes, no† | partial† | partial† | no† | no† |
| Integrated differential geometry (K, H, principal, singular points, ∬K dA) | yes | tangents/normals/curvature of curves only† | no† | no (GeoPDEs builds on it)† | general tools, not NURBS-specific† | general tools† |
| Export | VTK (+K, H fields), STL, OBJ, CSV, JSON (geomdl), PNG/SVG/PDF | OBJ, STL, OFF, VTK, JSON, YAML, CSV, CFG† | via MATLAB (STL)† | IGES import, mesh via Octave† | STL, OBJ, VTK, PLY… via `Export`† | STL, PLY via `Export`† |
| Teaching transparency (algorithms readable in the host language) | yes | yes | yes | yes | no (built-in kernel) | no (built-in kernel) |

The gap that NURBS-Maxima fills is the intersection of three columns: a free licence, NURBS objects that accept exact and symbolic data, and differential geometry of NURBS surfaces built into the package. geomdl is the reference for float performance and file formats, and we use it as an independent oracle (§5).

## 3. Background

### 3.1 B-spline bases

For a non-decreasing knot vector `U = {u_0, …, u_m}` the degree-`p` basis functions are given by the Cox–de Boor recursion [8, 9]

```
N_{i,0}(u) = 1 if u_i ≤ u < u_{i+1}, 0 otherwise
N_{i,p}(u) = (u − u_i)/(u_{i+p} − u_i) N_{i,p−1}(u) + (u_{i+p+1} − u)/(u_{i+p+1} − u_{i+1}) N_{i+1,p−1}(u).
```

The package implements span search (A2.1), the `p+1` non-zero functions (A2.2) and their derivatives (A2.3) [1]. A2.2 uses only field operations on `u` and the knots, so a symbolic `u` produces the exact polynomial pieces with no extra code.

### 3.2 NURBS curves and surfaces

A NURBS curve of degree `p` is `C(u) = Σ N_{i,p}(u) w_i P_i / Σ N_{i,p}(u) w_i`; a surface is the tensor product with weights `w_{ij}`. The implementation works with homogeneous points `P_i^w = (w_i P_i, w_i)`: it evaluates a polynomial B-spline in R⁴ and projects. Rational derivatives follow from the Leibniz rule (A4.2, A4.4):

```
C^(k) = ( A^(k) − Σ_{i=1..k} binom(k,i) w^(i) C^(k−i) ) / w.
```

Projective invariance lets any 4×4 homogeneous matrix act on the weighted control points, so rotations, scalings and perspective maps are exact.

### 3.3 Differential geometry

For a regular surface with normal `n = S_u × S_v / |S_u × S_v|`, `E = S_u·S_u`, `F = S_u·S_v`, `G = S_v·S_v`, `L = S_uu·n`, `M = S_uv·n`, `N = S_vv·n`, and

```
K = (LN − M²)/(EG − F²),   H = (EN − 2FM + GL)/(2(EG − F²)),   κ_{1,2} = H ± √(H² − K)   [4].
```

A point is singular when `S_u × S_v = 0`, as at the poles of a surface of revolution whose profile meets the axis. The Gauss–Bonnet theorem gives `∬K dA = 2πχ` on a closed surface: 4π for the sphere and every ellipsoid, 0 for the torus.

## 4. Design and implementation

### 4.1 Representation

Curves and surfaces are Maxima structures, `nurbscrv(degree, knots, ctrlpts, weights)` and `nurbssrf(degree_u, degree_v, knots_u, knots_v, ctrlpts, weights)`. Mathematical indices start at 0, as in [1]. Constructors validate knot vectors (length, monotonicity, multiplicity) and weights (positive when numeric; symbolic weights are accepted). The source (about 1 350 lines) is organised in the six stages of Table 2.

**Table 2.** Package organisation.

| Stage | Main functions |
|---|---|
| 1 Core | `nurbs_find_span`, `nurbs_basis_funs`, `nurbs_ders_basis_funs`, `bspline_basis_piecewise`, `knots_*`, `params_*` |
| 2 Curves | `nurbs_curve_point/derivs`, `nurbs_transform`, `nurbs_curve_insert_knot/refine/elevate_degree/decompose/split`, `nurbs_interpolate_curve`, `nurbs_curve_project` |
| 3 Surfaces | `nurbs_surface_point/derivs/normal`, `nurbs_surface_insert_knot/refine/elevate_degree/decompose/project`, `nurbs_revolve`, `nurbs_torus`, `nurbs_sphere`, `nurbs_surface_mesh` |
| 4 Differential geometry | `nurbs_surface_geometry` (shared helper), `nurbs_gaussian/mean_curvature`, `nurbs_principal_*`, `nurbs_classify_point`, `nurbs_singular_points`, `nurbs_curvature_sym`, `nurbs_surface_area/volume`, `nurbs_total_gaussian_curvature` |
| 5 Output | `nurbs_plot`, `nurbs_draw_*`, `nurbs_export_{plot,obj,stl,vtk,json,curve_csv}` |
| 6 Validation | `tests/` (rtest suite, geomdl references, JSON cross-validation, benchmarks) |

### 4.2 One code base, four arithmetics

No routine forces a conversion. The *kind* of the evaluation parameter decides the result: exact in, exact out; float in, float out; bigfloat in, bigfloat out at precision `fpprec`. `nurbs_float` and `nurbs_bfloat` convert whole objects for dense sampling. Symbolic weights flow through every algorithm because none of them branches on the sign or size of a weight; only the constructor checks positivity, and only when the weight is a number.

### 4.3 Tolerance policy

A single global tolerance `nurbs_tol` (default 10⁻¹⁰) governs every floating-point decision. All zero tests go through one predicate that distinguishes three cases and always returns `true` or `false`: a floating-point value is zero when `|x| ≤ nurbs_tol`; an exact constant when `radcan(x) = 0`; a symbolic expression when `ratsimp(x) = 0`. Singular, umbilic and planar points use relative versions of the same tolerance, e.g. `|S_u × S_v| ≤ nurbs_tol · max(1, |S_u||S_v|)`. The umbilic test uses `H² − K` instead of `κ_1 − κ_2` to avoid the square root, which would amplify round-off from 10⁻¹⁶ to 10⁻⁸.

### 4.4 Shared differential-geometry helper

`nurbs_surface_geometry(S, u, v)` performs one call to the A4.4 derivative routine (order 2) and returns `[S, S_u, S_v, n, (E,F,G), (L,M,N), K, H, (κ_1, κ_2), singular]`. Every curvature function, the VTK exporter (which writes K and H fields) and the integrand of `∬K dA` are thin wrappers around it, so a quantity is never recomputed and all of them agree on what is singular. Symbolic curvature on a whole patch (`nurbs_curvature_sym`) avoids the square root in `n` by using the unnormalised normal, which keeps `K(s, t; w)` rational.

### 4.5 Modelling operations

Curve algorithms follow [1]: knot insertion A5.1, refinement A5.4, degree elevation A5.9 (Bézier extraction, elevation and knot removal of the inserted copies), and Bézier decomposition by refinement to full multiplicity. Surface versions (A5.3, A5.5, and elevation along one direction) apply the curve algorithm to each column (direction u) or row (direction v) of the control net, which is how [1] derives them. Global interpolation (A9.1) solves the banded collocation system by LU factorisation with partial pivoting (`linsolve_by_lu`) in the float, bigfloat or exact rational field, after checking that there are at least `p+1` points and no coincident consecutive points. Point inversion uses Newton–Raphson on `C'(u)·(C(u) − P) = 0` (curves) and the 2×2 system of [1, eq. 6.6] (surfaces), started from the closest sample and with the stopping criteria of [1, §6.1].

### 4.6 Visualisation and exchange

Plots use Maxima's `draw` package and gnuplot, inline in wxMaxima. Meshes go to OBJ (with normals), ASCII STL and legacy VTK with `gaussian_curvature`, `mean_curvature` and normal point data for ParaView. The JSON exporter writes geomdl's schema with 17 significant digits and a `rational` flag that reflects the actual weights, so geometry can move between both systems without loss.

## 5. Validation

### 5.1 Protocol

Four levels, all runnable with one command each (README):

1. **Regression suite** (`tests/rtest_nurbs.mac`, standard `rtest` format, 95 problems, ≈28 s): textbook values (Example 2.3 of [1]), symbolic identities (the circle satisfies `x² + y² = 1` on all four spans; Bézier bases equal Bernstein polynomials; `K = −1/(1+x²+y²)²` for `z = xy`), invariants (partition of unity; invariance under knot insertion, refinement, degree elevation, weight scaling, affine and projective maps; Gauss–Bonnet on sphere, torus and ellipsoid; torus area and volume), symbolic-weight, bigfloat and export checks, and a robustness check that user-defined functions named `f`, `g` or `id` do not alter results (Maxima resolves `f(x)` to a global function even when a local variable `f` holds another one, so the package calls functional arguments through `apply`).
2. **Reference values from an independent library** (`tests/test_nurbs.mac`, Part B, 14 checks against `tests/reference_geomdl.mac` generated by `tests/compare_geomdl.py`).
3. **JSON round-trip cross-validation** (`tests/crossval_export.mac` + `tests/crossval_geomdl.py`): Maxima exports objects with `nurbs_export_json`, geomdl imports and re-evaluates them, and the script reports the maximum absolute error of points, first and second derivatives and normals, plus the control nets produced by the modelling operations on identical input.
4. **Examples** with known answers (`examples/`), whose printed values are quoted in Section 6.

### 5.2 Accuracy

**Table 3.** Representative accuracy (`tests/accuracy_table.mac`; the rtest suite checks the same quantities against tolerances, and all 95 problems pass).

| Check | Max. abs. error |
|---|---|
| Circle `x² + y² = 1`, symbolic, 4 spans | 0 (exact) |
| Knot insertion, refinement, degree elevation: curve unchanged | 0 (exact arithmetic) |
| Circle `r = 3`: curvature `1/3` (14 points) / length `6π` | 1.7·10⁻¹⁶ / 3.6·10⁻¹⁵ |
| Projective invariance, rational cubic (11 points) | 8.9·10⁻¹⁶ |
| Torus `K` vs analytic value (8×10 grid) | 1.8·10⁻¹⁵ |
| Torus area `4π²Rr` / volume `2π²Rr²` | 1.4·10⁻¹² / 4.5·10⁻¹³ |
| Sphere `∬K dA = 4π` | 3.7·10⁻¹² |
| Torus `∬K dA = 0` | 6.4·10⁻¹⁶ |

**Table 4.** Cross-validation with geomdl 5.5.1 (`tests/crossval/report.md`; 11 curve parameters, 7×7 surface parameters).

| Case | Point | 1st derivatives | 2nd derivatives | Normal |
|---|---|---|---|---|
| Bézier, B-spline 3D, interpolated curves | ≤ 8.9·10⁻¹⁶ | ≤ 1.8·10⁻¹⁵ | ≤ 2.8·10⁻¹⁴ | – |
| Rational curves (original, inserted, refined, elevated) | ≤ 8.9·10⁻¹⁶ | ≤ 3.6·10⁻¹⁵ | ≤ 2.3·10⁻¹³ | – |
| Circle and elevated circle | 1.1·10⁻¹⁶ | ≤ 1.8·10⁻¹⁵ | 7.1·10⁻¹⁵ | – |
| Bézier patch, cylinder, sphere | ≤ 2.2·10⁻¹⁶ | ≤ 8.9·10⁻¹⁶ | ≤ 8.3·10⁻¹⁵ | ≤ 2.2·10⁻¹⁶ |
| Torus (original, inserted in u, refined in v, elevated in v) | ≤ 4.4·10⁻¹⁶ | 1.8·10⁻¹⁵ | ≤ 1.8·10⁻¹⁴ | ≤ 2.2·10⁻¹⁶ |

**Table 5.** Modelling operations on identical input (homogeneous control points; knot vectors agree exactly).

| Operation | Max. abs. difference |
|---|---|
| Knot insertion A5.1 (u = 0.4, twice) | 2.2·10⁻¹⁶ |
| Refinement A5.4 vs repeated insertion in geomdl | 1.8·10⁻¹⁵ |
| Degree elevation A5.9, 3 → 5, compared in Bézier form | 1.8·10⁻¹⁵ |
| Surface knot insertion A5.3 (u = 0.3) | 4.4·10⁻¹⁶ |
| Bézier decomposition: curve (3 segments) / torus (16 patches) | 8.9·10⁻¹⁶ / 0 |

The only visible discrepancy in Part B is the arc length (1.75·10⁻³): geomdl approximates it by a polyline through evaluated points, while NURBS-Maxima uses Gauss–Legendre quadrature, which matches `2π` for the unit circle to 10⁻¹⁰. The difference is the reference's discretisation error.

### 5.3 A finding about degree elevation

Applying `geomdl.operations.degree_operations` directly to a non-Bézier rational cubic with interior knots changed the curve by up to 0.336 over 51 samples in our test, whereas A5.9 leaves it unchanged (Table 5 compares both after Bézier decomposition, where they agree). The geomdl routine appears to be intended for Bézier input *[confirm in the geomdl documentation]*; the point for users is that degree elevation of general NURBS should be validated, which the exact arithmetic of NURBS-Maxima makes straightforward: example 05 proves `C_elevated(t) − C(t) ≡ 0` symbolically on an interior span.

### 5.4 Performance and range of application

**Table 6.** Mean time per call in ms (`docs/benchmarks.md`; each measurement loops for at least 1 s).

| Operation | exact | float | bfloat (32 digits) | float + `compile` | geomdl | float / geomdl |
|---|---|---|---|---|---|---|
| Curve point (circle) | 0.32 | 0.25 | 2.27 | 0.087 | 0.0050 | 50× |
| Curve derivatives, order 2 | 0.94 | 0.62 | 5.50 | 0.45 | 0.014 | 44× |
| Surface point (torus) | 0.89 | 0.51 | 5.05 | 0.27 | 0.011 | 46× |
| Surface derivatives, order 2 | 4.42 | 2.65 | 18.8 | 1.67 | 0.074 | 36× |
| K, H, κ₁, κ₂ (shared helper) | 80.7 | 2.58 | 26.5 | 1.94 | n/a | – |
| 31×31 torus mesh with normals | – | 1.23 s | – | – | – | – |

Three observations follow. (i) Converting to floats with `nurbs_float` is the decisive optimisation: exact evaluation of the torus keeps `√2` symbolic and is 31 times slower for curvature. (ii) `compile()` of the span-search and basis kernels gives 1.3–2.9×; the gain is largest where those kernels dominate (curve points) and small where the cost is generic list arithmetic, and compiled code still returns exact results for exact input. (iii) bigfloat costs 5–10× float, which is affordable for the few evaluations where precision matters (§6.1).

NURBS-Maxima is therefore **not** a replacement for a meshing or CAD kernel. Its intended range is: exact and symbolic analysis and parametric studies; high-precision reference values for testing other libraries; verification of algorithms; teaching; and moderate meshes (10³–10⁴ points) exported for visualisation or fabrication.

## 6. Case studies

### 6.1 Gauss–Bonnet and the singular poles of the sphere (`examples/04_torus.mac`)

The exact NURBS torus (`R = 2`, `r = 1/2`, 9×9 net, degree 2×2) satisfies its implicit equation to 4.4·10⁻¹⁶. Along a meridian the computed `K` agrees with the analytic value to all printed digits, and the point classification follows the sign of `K` (Table 7).

**Table 7.** Torus curvature along `u = 0` (outward normal, so `H < 0` on the outer, convex part).

| v | K | K exact | H | type |
|---|---|---|---|---|
| 0 | 0.80000000 | 0.80000000 | −1.20000000 | elliptic |
| 0.125 | 0.60088442 | 0.60088442 | −1.15022110 | elliptic |
| 0.25 | 0.00000000 | 0.00000000 | −1.00000000 | parabolic |
| 0.375 | −0.85894894 | −0.85894894 | −0.78526277 | hyperbolic |
| 0.5 | −1.33333333 | −1.33333333 | −0.66666667 | hyperbolic |

Area 39.478417604356 (exact 4π²Rr = 39.478417604357), volume 9.869604401089 (2π²Rr²), and `∬K dA = −6.4·10⁻¹⁶` (Gauss–Bonnet: 0).

For closed surfaces of genus 0 the sphere gives `∬K dA = 12.566370614356` (4π to 3.7·10⁻¹²). The ellipsoid with semi-axes 3, 2, 1, obtained *exactly* by scaling the NURBS sphere, converges more slowly because its curvature varies by two orders of magnitude (0.028 near the poles of the revolution parametrisation versus 2.25 at the equator): the error is 7.4·10⁻⁴ with 2×2 Gauss cells per span and −3.0·10⁻⁶ with 4×4, consistent with the theorem.

The sphere of revolution has singular points at its poles: the grid search finds 18 singular nodes, all at `z = ±1`, and `nurbs_classify_point` returns `singular` there. Near a pole, `K` is a quotient of two quantities that vanish together. Table 8 evaluates `K − 1` on the unit sphere at distance `ε` (in parameter) from the pole, in double precision and in 50-digit bigfloat (with `nurbs_tol` lowered locally to 10⁻⁴⁰).

**Table 8.** `K − 1` near the pole of the unit NURBS sphere.

| ε | double | bigfloat, `fpprec = 50` |
|---|---|---|
| 10⁻⁴ | 1.4·10⁻⁹ | 1.5·10⁻⁴⁴ |
| 10⁻⁶ | −1.4·10⁻⁵ | 1.5·10⁻⁴¹ |
| 10⁻⁸ | 9.8·10⁻³ | 1.7·10⁻³⁶ |
| 10⁻⁹ | −37.9 | −1.7·10⁻³⁴ |
| 10⁻¹² | 9.8·10⁻⁷ | 6.7·10⁻⁴¹ |

Double precision loses accuracy as `ε²` (catastrophic cancellation in `EG − F²` and `LN − M²`) and returns a curvature of −37 instead of 1 at `ε = 10⁻⁹`; the erratic value at 10⁻¹² shows that the double results in this region are noise. The same code in bigfloat keeps at least 33 correct digits. This is the kind of reference computation the package is designed for.

### 6.2 Symbolic study of a design weight (`examples/06_symbolic_weight_study.mac`)

*Curve.* For the quadratic rational Bézier curve with control points (−1, 0), (0, 1), (1, 0) and weights 1, `w`, 1, the package returns the rational parametrisation and, after implicitisation, a conic whose discriminant is `64(w − 1)w²(w + 1)`: an ellipse for `w < 1`, a parabola for `w = 1` and a hyperbola for `w > 1`, the classical result [1, §7.2], here derived rather than quoted. The curvature is obtained in closed form: at the apex `κ(1/2) = w`, and at the end points `κ(0) = 1/(2^{3/2} w²)`. For `w = 1/√2` the package proves symbolically that the curve is the circular arc `x² + (y + 1)² = 2`. Table 9 adds the arc length, computed numerically, which has no elementary closed form.

**Table 9.** Apex curvature (exactly `w`) and arc length versus `w`.

| w | 0.25 | 0.5 | 0.707 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|
| κ(1/2) | 0.25 | 0.5 | 0.7071 | 1 | 2 | 4 |
| length | 2.0703 | 2.1580 | 2.2214 | 2.2956 | 2.4564 | 2.5981 |

*Surface.* A biquadratic Bézier patch over `[0, 2]²` with the centre control point raised to height 1 and given weight `w` has, by `nurbs_curvature_sym`, centre height `w/(w + 3)`, Gaussian curvature at the centre `K = w²/4` exactly, and `dK/dw = w/2 > 0`: increasing the weight sharpens the bump monotonically. The mean curvature at the centre is also closed-form (a radical expression in `w`). Table 10 shows the values and the area, and the symbolic `K(w = 3)` matches the numeric evaluation of the corresponding patch to 0.

**Table 10.** Centre curvature (closed form) and area of the weighted patch.

| w | 0.5 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|
| K (= w²/4) | 0.0625 | 0.25 | 1 | 4 | 16 |
| H | −0.25 | −0.5 | −1 | −2 | −4 |
| area | 4.0619 | 4.1727 | 4.3884 | 4.6851 | 4.9753 |

Numerical libraries can produce Tables 9 and 10 only by sampling `w`; the closed forms `κ = w` and `K = w²/4` give the designer the sensitivity directly.

### 6.3 From measured profile to ParaView and 3D printing (`examples/07_scientific_workflow.mac`)

Eight measured points of a vase profile are interpolated by a cubic (A9.1, centripetal parameters, LU solve; residual 4.4·10⁻¹⁶) and revolved into an exact rational surface (degrees 2×3, 9×8 net). On a 48×40 grid the shared helper produces the curvature fields: `K ∈ [−3.672, 1.018]`, `H ∈ [−1.027, −0.135]`; 53.7 % of the points are elliptic and 46.3 % hyperbolic, and `K` changes sign at `v ≈ 21/40`, the inflection of the profile, which the VTK field shows as a parabolic circle. The lateral area is 21.2322 and `∬K dA = 1.1192965`. Because the surface is an annulus (χ = 0), Gauss–Bonnet with boundary reduces for a surface of revolution to `∬K dA = 2π(r'(s₀) − r'(s₁))`, with `r'` the radial slope of the profile per unit arc length; evaluating the right-hand side from the profile gives 1.1192980, a relative agreement of 1.4·10⁻⁶ that validates the curvature field without any reference library. Refining the surface in `v` leaves it unchanged (1.1·10⁻¹⁶). The workflow ends with `07_vase.vtk` (K, H and normals for ParaView), `07_vase.stl` (3D printing), `07_vase.obj`, `07_vase.json` (geomdl) and a PNG.

## 7. Discussion and limitations

The package trades speed for generality. The cost is visible in Table 6 and is mitigated, not removed, by floats and `compile`. Symbolic curvature with irrational weights (e.g. `√2/2` in quadrics) produces large expressions; rational or symbolic weights are where `nurbs_curvature_sym` shines. Missing features include trimmed surfaces, surface–surface intersection, knot removal as a public operation, and surface fitting. Maxima's dynamic scoping and reserved names required care in the implementation (local helper names are prefixed `nurbs__`). On the positive side, the algorithms remain readable: a Spanish-language teaching notebook (`docs/programacion_nurbs.wxm`) walks through the source.

## 8. Conclusions

NURBS-Maxima brings the standard NURBS algorithms, a coherent differential-geometry layer and the usual exchange formats to a free computer algebra system, with one code base for exact, symbolic, double and arbitrary-precision arithmetic. Cross-validation against geomdl shows agreement to 2.3·10⁻¹³, and the case studies show results that floating-point libraries cannot give directly: curvature at a singular pole to 33 digits, and design sensitivities in closed form.

## Availability and reproducibility

* Source, examples, tests and data: *[repository URL and archived DOI, to complete]*. Licence: MIT.
* Requirements: Maxima ≥ 5.40 (tested 5.46.0), gnuplot; for validation, Python 3 and geomdl 5.x.
* Commands: `maxima -b tests/run_tests.mac` (95/95), `maxima -b tests/test_nurbs.mac` (95 + 14), `maxima -b tests/accuracy_table.mac` (Table 3), `maxima -b tests/crossval_export.mac && python3 tests/crossval_geomdl.py`, `maxima -b tests/benchmark.mac && python3 tests/benchmark_geomdl.py`; each example runs with `batch("examples/0N_*.mac")`.

## References

1. L. Piegl, W. Tiller. *The NURBS Book*, 2nd ed. Springer, Berlin, 1997.
2. G. Farin. *Curves and Surfaces for CAGD: A Practical Guide*, 5th ed. Morgan Kaufmann, 2002.
3. O. R. Bingol, A. Krishnamurthy. NURBS-Python: An open-source object-oriented NURBS modeling framework in Python. *SoftwareX* 9 (2019) 85–94.
4. M. P. do Carmo. *Differential Geometry of Curves and Surfaces*. Prentice-Hall, 1976.
5. D. F. Rogers. *An Introduction to NURBS: With Historical Perspective*. Morgan Kaufmann, 2001.
6. T. J. R. Hughes, J. A. Cottrell, Y. Bazilevs. Isogeometric analysis: CAD, finite elements, NURBS, exact geometry and mesh refinement. *Comput. Methods Appl. Mech. Engrg.* 194 (2005) 4135–4195.
7. Maxima, a Computer Algebra System, version 5.46.0. https://maxima.sourceforge.io/ *[to verify: version/year cited]*
8. C. de Boor. On calculating with B-splines. *J. Approx. Theory* 6 (1972) 50–62.
9. M. G. Cox. The numerical evaluation of B-splines. *IMA J. Appl. Math.* 10 (1972) 134–149.
10. M. Spink. NURBS Toolbox for MATLAB. MATLAB Central File Exchange. *[to verify: URL, year, licence]*
11. M. Spink, D. Claxton, C. de Falco, R. Vázquez. The NURBS toolbox (Octave Forge `nurbs` package). *[to verify: version and URL]*
12. Wolfram Research. `BSplineFunction`, Wolfram Language & System Documentation Center. *[to verify]*
13. Maplesoft. `CurveFitting[BSpline]` and `CurveFitting[BSplineCurve]`, Maple help. *[to verify]*
14. C. de Falco, A. Reali, R. Vázquez. GeoPDEs: a research tool for isogeometric analysis of PDEs. *Adv. Eng. Softw.* 42 (2011) 1020–1034.
