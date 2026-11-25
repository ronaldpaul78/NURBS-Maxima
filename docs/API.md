# NURBS-Maxima API reference

Version 2.0.0. Load the whole package with

```maxima
load("path/to/NURBS-Maxima/nurbs_maxima.mac")$
```

## Conventions

| Topic | Convention |
|---|---|
| Indices | Knot spans, basis functions `N_{i,p}` and control points `P_i` are numbered from **0**, as in Piegl & Tiller, *The NURBS Book* (2nd ed., 1997). Maxima lists start at 1, so `P_i` is `C@ctrlpts[i+1]`. |
| Point | A list `[x, y]` or `[x, y, z]`. |
| Curve | Structure `nurbscrv(degree, knots, ctrlpts, weights)`; fields are read with `@`: `C@degree`, `C@knots`, `C@ctrlpts`, `C@weights`. |
| Surface | Structure `nurbssrf(degree_u, degree_v, knots_u, knots_v, ctrlpts, weights)` with `ctrlpts[i+1][j+1] = P_{i,j}` (rows along u). |
| Arithmetic | Every routine accepts exact numbers (rationals, radicals), **symbolic weights**, double floats or **bigfloats**. The kind of the parameter decides the kind of the result: exact in, exact out; float in, float out; bfloat in, bfloat out (precision `fpprec`). |
| Options | Optional arguments are written `key = value`. Plot options that do not belong to the package are passed to `draw` unchanged. |
| Direction | Surface operations take `dir = 'u` or `'v`. |
| Normal | `n = S_u × S_v / |S_u × S_v|`. For `nurbs_sphere` and `nurbs_torus` it points outwards, so a convex region has `H < 0` (the unit sphere gives `K = 1`, `H = -1`). |

### Global variables

| Variable | Default | Meaning |
|---|---|---|
| `nurbs_version` | `"2.0.0"` | Package version. |
| `nurbs_tol` | `1.0e-10` | The single tolerance of the package (see below). |
| `nurbs_inline` | `auto` | `auto`, `true` or `false`. With `auto`, plots are drawn inline when running inside wxMaxima. |
| `nurbs_show_plots` | `true` | The example files open plot windows only when this is `true`. |

### Tolerance semantics

All zero tests go through one internal predicate, which always returns `true` or `false` (never `unknown`):

* a floating-point value (float or bfloat) is zero when `|x| <= nurbs_tol`;
* an exact constant (rationals, radicals such as `sqrt(2)`) is zero when `radcan(x) = 0`;
* a symbolic expression (for example one that contains a weight `w`) is zero when `ratsimp(x) = 0`.

Relative tests use the same variable: a point is singular when `|S_u × S_v| <= nurbs_tol · max(1, |S_u| |S_v|)`, umbilic when `|H² − K| <= nurbs_tol · max(1, H²)`, and planar when `|K|, |H| <= nurbs_tol`. With bfloats, lower the tolerance locally, e.g. `block([nurbs_tol: 1e-40], ...)`.

---

## Stage 1. Mathematical core

### Knot vectors and parameters

| Function | Description |
|---|---|
| `knots_clamped_uniform(n, p)` | Clamped uniform knot vector for `n` control points and degree `p` (p+1 zeros, equally spaced interior knots, p+1 ones). |
| `knots_uniform(n, p)` | Unclamped uniform vector `0, 1/(n+p), …, 1`. |
| `knots_bezier(p)` | `[0,…,0, 1,…,1]` with p+1 repetitions each. |
| `knots_validate(U, p, n)` | Checks length `n+p+1`, monotonicity and multiplicity `≤ p+1`. Returns `true` or raises an explanatory error. |
| `knots_multiplicity(u, U)` | Multiplicity of `u` in `U`. |
| `knots_unique(U)` | Distinct knot values. |
| `knots_normalize(U)` | Rescales `U` to `[0, 1]`. |
| `knots_is_clamped(U, p)` | `true` if both ends have multiplicity `p+1`. |
| `knots_spans(U)` | Non-empty spans `[[i, U_i, U_{i+1}], …]`. |
| `knots_domain_spans(U, p)` | Non-empty spans inside the valid domain `[U_p, U_{n+1}]`. |
| `knots_continuity(U, p)` | `[[u, k], …]`: the curve is `C^k` at each interior knot (`k = p − multiplicity`). |
| `params_uniform(Q)`, `params_chord(Q)`, `params_centripetal(Q)` | Parameters for the data points `Q`. Chord and centripetal raise an error on coincident consecutive points. |
| `knots_averaging(t, p)` | Knot vector by averaging (P&T eq. 9.8). |

### Basis functions

| Function | Description |
|---|---|
| `nurbs_find_span(n, p, u, U)` | Algorithm A2.1. Index `i` (from 0) with `U_i ≤ u < U_{i+1}`; `n` = number of control points − 1. |
| `nurbs_basis_funs(i, u, p, U)` | Algorithm A2.2. `[N_{i-p,p}(u), …, N_{i,p}(u)]`. With a symbolic `u` it returns polynomials. |
| `nurbs_ders_basis_funs(i, u, p, n, U)` | Algorithm A2.3. Matrix `D` with `D[k+1, j+1] = dᵏ/duᵏ N_{i-p+j,p}(u)`, `k = 0..n`. |
| `bspline_basis(i, p, u, U)` | Cox–de Boor recursion for a single `N_{i,p}(u)`. |
| `bspline_basis_all(p, u, U)` | Full vector `[N_{0,p}(u), …, N_{n,p}(u)]`. |
| `bspline_basis_piecewise(i, p, U, u)` | Exact piecewise polynomials of `N_{i,p}`: `[[[a,b], poly], …]`. |
| `bernstein(i, n, t)` | Bernstein polynomial. |

---

## Stage 2. NURBS curves

### Construction

| Function | Description |
|---|---|
| `nurbs_curve(p, U, P [, W])` | NURBS curve; all weights are 1 if `W` is omitted. Validates the data. Numeric weights must be positive; symbolic weights are accepted. |
| `bspline_curve(p, U, P)` | Non-rational B-spline. |
| `nurbs_bezier_curve(P [, W])` | Bézier curve, rational when weights are given. |
| `nurbs_line(A, B)` | Straight segment. |
| `nurbs_circle(r [, centre])` | Exact 9-point circle, degree 2, weights `1, √2/2`. 2D, or 3D if the centre has 3 coordinates. |
| `nurbs_circle_3d(C, r, e1, e2)` | Circle in the plane spanned by the orthonormal vectors `e1, e2`. |
| `nurbs_interpolate_curve(Q, p [, method])` | Algorithm A9.1, global interpolation; `method` = `'chord` (default), `'centripetal` or `'uniform`. Requires at least `p+1` points and no coincident consecutive points. The system is solved by LU factorisation (`linsolve_by_lu`, float, bigfloat or exact rational field). |

### Evaluation

| Function | Description |
|---|---|
| `nurbs_curve_point(C, u)` | Algorithm A4.1: `C(u)`. |
| `nurbs_curve_derivs(C, u, d)` | Algorithm A4.2: `[C, C', …, C^(d)]` of the rational curve. |
| `nurbs_curve_tangent(C, u)` | Unit tangent. |
| `nurbs_curve_span_expr(C, i, t)` | Exact rational expression of the curve on span `i` as a function of `t`. |
| `nurbs_curve_expr(C, t)` | Piecewise representation `[[[a,b], [x(t), y(t), …]], …]`. |
| `nurbs_curve_domain(C)` | `[U_p, U_{n+1}]`. |
| `nurbs_curve_ctrlptsw(C)` | Homogeneous points `[w x, w y, (w z), w]`. |
| `nurbs_curve_samples(C, n)` | `n+1` points equally spaced in the parameter. |

### Weights

| Function | Description |
|---|---|
| `nurbs_set_weight(C, i, w)` | New curve with `w_i = w` (`w` may be a symbol). |
| `nurbs_scale_weights(C, c)` | Multiplies all weights by `c`; the curve is unchanged. |
| `nurbs_is_rational(obj)` | `true` unless all weights are provably equal. Safe with symbolic weights. |

### Transformations (curves and surfaces)

| Function | Description |
|---|---|
| `nurbs_transform(obj, M)` | Applies a homogeneous 4×4 matrix, affine or **projective**, to the weighted control points. Exact, by projective invariance of NURBS. |
| `nurbs_translate(obj, t)` | Translation by the vector `t`. |
| `nurbs_scale(obj, s)` | Scaling by a scalar `s` or by `[sx, sy, sz]`. |
| `nurbs_rotate(obj, axis, θ)` | Rotation; `axis` is `x`, `y`, `z` or a vector (Rodrigues' formula). |
| `nurbs_translation_matrix(t)`, `nurbs_scaling_matrix(s)`, `nurbs_rotation_matrix(axis, θ)` | 4×4 matrices. |

### Modelling algorithms

| Function | Description |
|---|---|
| `nurbs_curve_insert_knot(C, u, r)` | Algorithm A5.1: inserts `u` r times without changing the curve. `u` must lie in the domain and the final multiplicity must not exceed `p`. |
| `nurbs_curve_refine(C, X)` | Algorithm A5.4: inserts the sorted list of knots `X` at once. |
| `nurbs_curve_elevate_degree(C, t)` | Algorithm A5.9: raises the degree by `t`; the curve is unchanged. |
| `nurbs_curve_decompose(C)` | Bézier decomposition (A5.6, via refinement). List of Bézier segments; each keeps its parameter interval `[a, b]`. |
| `nurbs_curve_split(C, u)` | Splits at `u`, which must be **strictly interior** to the domain. Returns `[C1, C2]`. |
| `nurbs_curve_reverse(C)` | Reverses the parametrisation. |
| `nurbs_curve_project(C, P [, samples=20, maxiter=50])` | Point inversion/projection (P&T §6.1): Newton–Raphson on `C'(u)·(C(u)−P) = 0` from the best of `samples` points per span; handles closed curves. Returns `[u, C(u), distance]` in float (bfloat if `P` is bfloat). |

---

## Stage 3. NURBS surfaces

### Construction

| Function | Description |
|---|---|
| `nurbs_surface(p, q, U, V, P [, W])` | Tensor-product surface; `P` is a list of rows. |
| `bspline_surface(p, q, U, V, P)` | Non-rational surface. |
| `nurbs_bezier_surface(P [, W])` | Bézier patch. |
| `nurbs_bilinear(P00, P10, P01, P11)` | Bilinear patch. |
| `nurbs_extrude(C, d)` | Extrusion of `C` along the vector `d`. |
| `nurbs_ruled(C1, C2)` | Ruled surface between two curves with the same degree and knots. |
| `nurbs_revolve(C)` | 360° surface of revolution about the z axis; the profile lies in the xz plane. |
| `nurbs_torus(R, r)`, `nurbs_sphere(r)`, `nurbs_cylinder(r, h)` | Exact primitives. |

### Evaluation

| Function | Description |
|---|---|
| `nurbs_surface_point(S, u, v)` | Algorithm A4.3. |
| `nurbs_surface_derivs(S, u, v, d)` | Algorithm A4.4. Matrix `SKL` with `SKL[k+1, l+1] = ∂^{k+l}S/∂uᵏ∂vˡ`. |
| `nurbs_surface_tangents(S, u, v)` | `[S_u, S_v]`. |
| `nurbs_surface_normal(S, u, v)` | Unit normal, or `und` at a singular point. |
| `nurbs_surface_span_expr(S, i, j, s, t)` | Exact rational expression of patch `(i, j)`. |
| `nurbs_surface_isocurve(S, dir, val)` | Exact isoparametric curve as a NURBS curve (`dir = 'u` fixes u). |
| `nurbs_surface_grid(S, nu, nv)` | `(nu+1)×(nv+1)` grid of points. |
| `nurbs_surface_mesh(S, nu, nv)` | `[vertices, triangles (1-based indices), normals]`; one derivative evaluation per vertex. |
| `nurbs_surface_domain(S)` | `[[u0,u1],[v0,v1]]`. |
| `nurbs_surface_ctrlptsw(S)` | Homogeneous control net. |

### Modelling algorithms

Each operation applies the curve algorithm to every column (`'u`) or row (`'v`) of the net.

| Function | Description |
|---|---|
| `nurbs_surface_insert_knot(S, dir, val, r)` | Algorithm A5.3: knot insertion along u or v. |
| `nurbs_surface_refine(S, dir, X)` | Algorithm A5.5: knot refinement along u or v. |
| `nurbs_surface_elevate_degree(S, dir, t)` | Degree elevation along u or v (A5.9 applied to the rows/columns, as in P&T §5.5). |
| `nurbs_surface_decompose(S)` | Bézier decomposition: a list of rows of Bézier patches, each with its own parameter rectangle. |
| `nurbs_surface_project(S, P [, samples=20, maxiter=50])` | Point inversion (P&T eq. 6.6), 2D Newton–Raphson. Returns `[u, v, S(u,v), distance]`. |

### Numeric conversion

| Function | Description |
|---|---|
| `nurbs_float(obj)` | Copy with every datum in double float. The main speed-up for sampling. |
| `nurbs_bfloat(obj)` | Copy in bigfloat (precision `fpprec`). |
| `nurbs_compile()` | Compiles the span search and basis-function kernels (`compile()`); results stay exact for exact input. See `docs/benchmarks.md`. |

---

## Stage 4. Differential geometry

### Curves

| Function | Description |
|---|---|
| `nurbs_curve_curvature(C, u)` | `κ = |C'×C''| / |C'|³`. |
| `nurbs_curve_torsion(C, u)` | `τ = (C'×C'')·C''' / |C'×C''|²`. |
| `nurbs_curve_frenet(C, u)` | Frenet frame `[T, N, B]`. |
| `nurbs_curve_length(C [, nsub])` | Arc length, 5-point Gauss–Legendre with `nsub` sub-intervals per span (default 8). |
| `nurbs_bounding_box(obj)` | Bounding box of the control net (it contains the object). |

### Surfaces

All the local quantities come from one shared helper that makes a single call to `nurbs_surface_derivs(S, u, v, 2)`.

| Function | Description |
|---|---|
| `nurbs_surface_geometry(S, u, v)` | `[point, Su, Sv, n, [E,F,G], [L,M,N], K, H, [k1,k2], singular]`. At a singular point the entries from `n` onwards are `und` and `singular = true`. Use it when several quantities are needed at the same point. |
| `nurbs_fundamental_forms(S, u, v)` | `[[E, F, G], [L, M, N], n]`, or `und` at a singular point. |
| `nurbs_gaussian_curvature(S, u, v)` | `K = (LN − M²)/(EG − F²)`. |
| `nurbs_mean_curvature(S, u, v)` | `H = (EN − 2FM + GL)/(2(EG − F²))`. |
| `nurbs_principal_curvatures(S, u, v)` | `[k1, k2]` with `k1 ≥ k2`. |
| `nurbs_principal_directions(S, u, v)` | Unit principal directions in R³, or `und` at an umbilic or singular point. |
| `nurbs_classify_point(S, u, v [, tol])` | `singular`, `planar`, `umbilic`, `elliptic`, `hyperbolic` or `parabolic`. `tol` overrides `nurbs_tol` locally. |
| `nurbs_is_singular(S, u, v [, tol])` | `true` if `|S_u × S_v| <= tol · max(1, |S_u||S_v|)` (exact test for exact parameters). |
| `nurbs_singular_points(S, nu, nv [, tol])` | Searches a grid; returns `[[u, v, S(u,v)], …]`. |
| `nurbs_curvature_sym(S, i, j, s, t)` | Exact `[K, H]` on patch `(i, j)` as functions of `s, t`. Works with symbolic weights (e.g. `K(s, t; w)`). Costly with irrational weights. |
| `nurbs_surface_area(S [, nsub])` | Area (5×5 Gauss rule, `nsub` = 4 cells per span and direction). |
| `nurbs_total_gaussian_curvature(S [, nsub])` | `∬ K dA`, for checking Gauss–Bonnet (`2πχ`). Singular quadrature nodes contribute 0 (they have zero area element). |
| `nurbs_surface_volume(S [, nsub])` | Enclosed volume of a closed surface (divergence theorem). |

---

## Stage 5. Visualisation and export

| Function | Description |
|---|---|
| `nurbs_plot(obj or [objs] [, options])` | Draws with `draw2d`/`draw3d`, or inline with `wxdraw2d`/`wxdraw3d` in wxMaxima. |
| `nurbs_draw_curve(C [, options])` | List of `draw` objects. Options: `samples=100`, `control=true`, `curve_color=blue`, `control_color=red`. |
| `nurbs_draw_surface(S [, options])` | Options: `nu=30`, `nv=30`, `control=false`, `net_color=red`, `normals=false` or n (an n×n grid of normals), `normal_length=0.3`, `surface_color`. |
| `nurbs_export_plot(obj, "file.ext" [, options])` | `png`, `svg`, `pdf`, `eps` or `jpg` through gnuplot. |
| `nurbs_export_obj(S, file, nu, nv)` | Wavefront OBJ mesh with normals. |
| `nurbs_export_stl(S, file, nu, nv)` | ASCII STL mesh. |
| `nurbs_export_vtk(S, file, nu, nv)` | VTK legacy file with `gaussian_curvature`, `mean_curvature` and normal fields, for ParaView. |
| `nurbs_export_curve_csv(C, file, n)` | Samples `u, x, y[, z]`. |
| `nurbs_export_json(obj, file)` | Exact definition in geomdl's JSON format (`geomdl.exchange.import_json`), 17 significant digits; the `rational` flag reflects the weights. |
| `nurbs_info(obj)` | One-line summary. |

Combining objects with `draw`:

```maxima
wxdraw3d(nurbs_draw_surface(T, nu=40), nurbs_draw_curve(C, control=false));
```

---

## Stage 6. Validation

| File | Purpose |
|---|---|
| `tests/rtest_nurbs.mac` | 95 problems in Maxima's standard rtest format. |
| `tests/run_tests.mac` | Loads the package and runs `batch("rtest_nurbs.mac", test)`. |
| `tests/test_nurbs.mac` | Part A: the rtest suite. Part B: comparison with geomdl reference values (`tests/reference_geomdl.mac`). |
| `tests/compare_geomdl.py` | Regenerates `tests/reference_geomdl.mac` with geomdl. |
| `tests/crossval_export.mac` + `tests/crossval_geomdl.py` | Exports test objects with `nurbs_export_json`, re-evaluates them in geomdl and reports maximum absolute errors for points, derivatives and normals (`tests/crossval/report.md`). |
| `tests/benchmark.mac` + `tests/benchmark_geomdl.py` | Timings (exact, float, bfloat, compiled, geomdl), written to `docs/benchmarks.md`. |

```bash
maxima -b tests/run_tests.mac                 # rtest suite
maxima -b tests/test_nurbs.mac                # rtest + geomdl references
maxima -b tests/crossval_export.mac && python3 tests/crossval_geomdl.py
```
