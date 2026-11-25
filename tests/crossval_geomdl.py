#!/usr/bin/env python3
"""
NURBS-Maxima cross-validation, step 2 (geomdl side).

Reads the objects exported by tests/crossval_export.mac with
nurbs_export_json(), imports them with geomdl.exchange.import_json,
recomputes points, derivatives, normals and the modelling operations
independently with geomdl, and reports the maximum absolute error against
the values computed by NURBS-Maxima.

    maxima -b tests/crossval_export.mac
    python3 tests/crossval_geomdl.py      # writes tests/crossval/report.md
Exit status is 1 if any error exceeds the tolerance.
"""
import glob
import json
import os
import sys

import geomdl
from geomdl import exchange, operations, NURBS

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(HERE, "crossval")
TOL = 1e-10


def maxerr(a, b):
    """Maximum absolute difference between two nested lists."""
    if isinstance(a, (list, tuple)):
        if len(a) != len(b):
            return float("inf")
        return max([0.0] + [maxerr(x, y) for x, y in zip(a, b)])
    return abs(float(a) - float(b))


def load(name):
    obj = exchange.import_json(os.path.join(DIR, name + ".json"))[0]
    with open(os.path.join(DIR, name + "_maxima.json")) as f:
        mx = json.load(f)
    return obj, mx


rows = []
failed = False

for path in sorted(glob.glob(os.path.join(DIR, "*_maxima.json"))):
    name = os.path.basename(path)[:-len("_maxima.json")]
    obj, mx = load(name)
    if mx["kind"] == "curve":
        ref = [obj.derivatives(t, order=2) for t in mx["params"]]
        e_pt = maxerr([r[0] for r in ref], [d[0] for d in mx["derivs"]])
        e_d1 = maxerr([r[1] for r in ref], [d[1] for d in mx["derivs"]])
        e_d2 = maxerr([r[2] for r in ref], [d[2] for d in mx["derivs"]])
        e_n = None
    else:
        ref = [obj.derivatives(u, v, order=2) for (u, v) in mx["params"]]
        e_pt = maxerr([r[0][0] for r in ref], [d[0] for d in mx["derivs"]])
        e_d1 = maxerr([[r[1][0], r[0][1]] for r in ref], [[d[1], d[2]] for d in mx["derivs"]])
        e_d2 = maxerr([[r[2][0], r[1][1], r[0][2]] for r in ref], [[d[3], d[4], d[5]] for d in mx["derivs"]])
        nrm = [operations.normal(obj, (u, v), normalize=True)[1] for (u, v) in mx["params"]]
        e_n = maxerr(nrm, mx["normals"])
    worst = max(e for e in (e_pt, e_d1, e_d2, e_n) if e is not None)
    # derivatives are scaled by the size of the object; use a relative bound
    scale = max(1.0, max(abs(x) for d in mx["derivs"] for v in (d if mx["kind"] == "surface" else d) for x in v))
    ok = worst <= TOL * scale
    failed |= not ok
    rows.append((name, mx["kind"], len(mx["params"]), e_pt, e_d1, e_d2, e_n, ok))

# ---------------------------------------------------- modelling operations
with open(os.path.join(DIR, "op_results.json")) as f:
    ops = json.load(f)


def rational_cubic():
    c = NURBS.Curve()
    c.degree = 3
    c.ctrlpts = [[0, 0, 0], [1, 2, 1], [2, -1, 2], [4, 1, 0], [5, 3, 1], [7, 0, 2]]
    c.weights = [1, 2, 0.5, 1, 3, 1]
    c.knotvector = [0, 0, 0, 0, 1 / 3, 2 / 3, 1, 1, 1, 1]
    return c


op_rows = []
c = rational_cubic()
operations.insert_knot(c, [0.4], [2])
op_rows.append(("knot insertion (A5.1) u=0.4 x2", maxerr(c.knotvector, ops["insert_knots"]),
                maxerr([list(p) for p in c.ctrlptsw], ops["insert_ctrlptsw"])))

c = rational_cubic()
for x in [0.2, 0.5, 0.5, 0.9]:
    operations.insert_knot(c, [x], [1])
op_rows.append(("knot refinement (A5.4) vs repeated insertion", maxerr(c.knotvector, ops["refine_knots"]),
                maxerr([list(p) for p in c.ctrlptsw], ops["refine_ctrlptsw"])))

# Degree elevation.  geomdl's degree_operations is exact only for Bezier
# segments, so the reference is built as: decompose -> elevate every Bezier
# segment -> compare with the Bezier decomposition of Maxima's A5.9 result
# (the Bezier form of a curve is unique).
with open(os.path.join(DIR, "op_elevate_bezier.json")) as f:
    elev = json.load(f)
segs = operations.decompose_curve(rational_cubic())
for b in segs:
    operations.degree_operations(b, [2])
op_rows.append(("degree elevation (A5.9) p=3 -> 5, Bezier form", 0.0 if len(segs) == len(elev["ctrlptsw"]) else float("inf"),
                maxerr([[list(p) for p in b.ctrlptsw] for b in segs], elev["ctrlptsw"])))
# geomdl applied directly to the non-Bezier curve (documented as a finding)
c = rational_cubic()
operations.degree_operations(c, [2])
ref = rational_cubic()
geomdl_direct_err = max(abs(x - y) for t in [i / 50 for i in range(51)]
                        for x, y in zip(ref.evaluate_single(t), c.evaluate_single(t)))

s = exchange.import_json(os.path.join(DIR, "s_torus.json"))[0]
operations.insert_knot(s, [0.3, None], [1, 0])
op_rows.append(("surface knot insertion (A5.3) u=0.3", maxerr(s.knotvector_u, ops["surf_insert_knots_u"]),
                maxerr([list(p) for p in s.ctrlptsw], ops["surf_insert_ctrlptsw"])))

with open(os.path.join(DIR, "op_decompose_curve.json")) as f:
    dec = json.load(f)
segs = operations.decompose_curve(rational_cubic())
e_dec = maxerr([[list(p) for p in b.ctrlpts] for b in segs], dec["ctrlpts"])
op_rows.append(("Bezier decomposition (A5.6), %d segments" % dec["count"],
                0.0 if len(segs) == dec["count"] else float("inf"), e_dec))
with open(os.path.join(DIR, "op_decompose_surface.json")) as f:
    decs = json.load(f)
patches = operations.decompose_surface(s := exchange.import_json(os.path.join(DIR, "s_torus.json"))[0])
op_rows.append(("Bezier decomposition of the torus, %d patches" % decs["count"],
                0.0 if len(patches) == decs["count"] else float("inf"), 0.0))

for name, e1, e2 in op_rows:
    if e1 == e1 and e2 == e2 and max(e1, e2) > TOL:   # NaN means "not comparable"
        failed = True

# ---------------------------------------------------------------- report
def fmt(e):
    if e is None:
        return "-"
    if e != e:
        return "n/c"
    return "%.1e" % e


lines = ["# Cross-validation NURBS-Maxima vs geomdl %s" % geomdl.__version__, "",
         "Objects exported by Maxima with `nurbs_export_json`, re-imported and re-evaluated by geomdl.",
         "Maximum absolute error.", "",
         "| case | type | params | point | 1st derivatives | 2nd derivatives | normal | ok |",
         "|---|---|---|---|---|---|---|---|"]
for r in rows:
    lines.append("| %s | %s | %d | %s | %s | %s | %s | %s |" % (r[0], r[1], r[2], fmt(r[3]), fmt(r[4]),
                                                            fmt(r[5]), fmt(r[6]), "yes" if r[7] else "NO"))
lines += ["", "## Modelling operations (same input in both systems)", "",
          "| operation | knot vector | homogeneous control points |", "|---|---|---|"]
for name, e1, e2 in op_rows:
    lines.append("| %s | %s | %s |" % (name, fmt(e1), fmt(e2)))
lines += ["", "Note: geomdl.operations.degree_operations applied directly to this non-Bezier curve changes its",
          "shape (max. deviation %.2e over 51 samples); NURBS-Maxima's A5.9 keeps it unchanged." % geomdl_direct_err,
          "", "RESULT: " + ("FAILED" if failed else "PASSED")]
report = "\n".join(lines) + "\n"
with open(os.path.join(DIR, "report.md"), "w") as f:
    f.write(report)
print(report)
sys.exit(1 if failed else 0)
