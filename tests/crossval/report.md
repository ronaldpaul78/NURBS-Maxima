# Cross-validation NURBS-Maxima vs geomdl 5.5.1

Objects exported by Maxima with `nurbs_export_json`, re-imported and re-evaluated by geomdl.
Maximum absolute error.

| case | type | params | point | 1st derivatives | 2nd derivatives | normal | ok |
|---|---|---|---|---|---|---|---|
| c_bezier | curve | 11 | 0.0e+00 | 0.0e+00 | 0.0e+00 | - | yes |
| c_bspline3d | curve | 11 | 8.9e-16 | 1.8e-15 | 2.8e-14 | - | yes |
| c_circle_elevated | curve | 11 | 1.1e-16 | 8.9e-16 | 7.1e-15 | - | yes |
| c_circle | curve | 11 | 1.1e-16 | 1.8e-15 | 7.1e-15 | - | yes |
| c_interp | curve | 11 | 0.0e+00 | 0.0e+00 | 0.0e+00 | - | yes |
| c_rational_elevated | curve | 11 | 4.4e-16 | 3.6e-15 | 1.1e-13 | - | yes |
| c_rational_inserted | curve | 11 | 2.2e-16 | 3.6e-15 | 1.4e-14 | - | yes |
| c_rational | curve | 11 | 8.9e-16 | 3.6e-15 | 1.4e-14 | - | yes |
| c_rational_refined | curve | 11 | 4.4e-16 | 1.8e-15 | 2.3e-13 | - | yes |
| s_bezier | surface | 7 | 2.2e-16 | 4.4e-16 | 2.2e-16 | 2.2e-16 | yes |
| s_cylinder | surface | 7 | 2.2e-16 | 8.9e-16 | 8.3e-15 | 2.2e-16 | yes |
| s_sphere | surface | 7 | 1.1e-16 | 4.4e-16 | 7.1e-15 | 1.1e-16 | yes |
| s_torus_elevated_v | surface | 7 | 5.6e-17 | 1.8e-15 | 1.4e-14 | 1.1e-16 | yes |
| s_torus_inserted_u | surface | 7 | 4.4e-16 | 1.8e-15 | 1.8e-14 | 2.2e-16 | yes |
| s_torus | surface | 7 | 4.4e-16 | 1.8e-15 | 1.8e-14 | 2.2e-16 | yes |
| s_torus_refined_v | surface | 7 | 2.2e-16 | 1.8e-15 | 1.4e-14 | 2.2e-16 | yes |

## Modelling operations (same input in both systems)

| operation | knot vector | homogeneous control points |
|---|---|---|
| knot insertion (A5.1) u=0.4 x2 | 0.0e+00 | 2.2e-16 |
| knot refinement (A5.4) vs repeated insertion | 0.0e+00 | 1.8e-15 |
| degree elevation (A5.9) p=3 -> 5, Bezier form | 0.0e+00 | 1.8e-15 |
| surface knot insertion (A5.3) u=0.3 | 0.0e+00 | 4.4e-16 |
| Bezier decomposition (A5.6), 3 segments | 0.0e+00 | 8.9e-16 |
| Bezier decomposition of the torus, 16 patches | 0.0e+00 | 0.0e+00 |

Note: geomdl.operations.degree_operations applied directly to this non-Bezier curve changes its
shape (max. deviation 3.36e-01 over 51 samples); NURBS-Maxima's A5.9 keeps it unchanged.

RESULT: PASSED
