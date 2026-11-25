#!/usr/bin/env python3
"""
Generate reference values with geomdl (NURBS-Python) to validate
NURBS-Maxima.  Writes tests/reference_geomdl.mac, which test_nurbs.mac loads
and compares against the Maxima routines.

Usage (from the NURBS-Maxima folder):
    pip install geomdl==5.5.1      # any 5.x should work
    python3 tests/compare_geomdl.py

Optional: JSON files exported with nurbs_export_json() can be passed as
arguments; they are imported with geomdl and evaluated (interoperability
check).
"""
import math
import os
import sys

import geomdl
from geomdl import NURBS, BSpline, fitting, operations, exchange

S2 = math.sqrt(2.0) / 2.0
HERE = os.path.dirname(os.path.abspath(__file__))


def mx(x):
    """Convert a Python value/list into a Maxima literal."""
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(mx(e) for e in x) + "]"
    return repr(float(x)).replace("e", "e")


# ---------------------------------------------------------------- objects
def curve(p, U, P, W=None):
    c = NURBS.Curve()
    c.degree = p
    c.ctrlpts = P
    c.weights = W if W is not None else [1.0] * len(P)
    c.knotvector = U
    return c


def surface(p, q, U, V, P, W=None):
    s = NURBS.Surface()
    s.degree_u, s.degree_v = p, q
    s.ctrlpts_size_u, s.ctrlpts_size_v = len(P), len(P[0])
    s.ctrlpts = [pt for row in P for pt in row]
    s.weights = [w for row in W for w in row] if W else [1.0] * (len(P) * len(P[0]))
    s.knotvector_u, s.knotvector_v = U, V
    return s


CIRCLE_U = [0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1]
CIRCLE_D = [[1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0], [-1, -1], [0, -1], [1, -1], [1, 0]]
CIRCLE_W = [1, S2, 1, S2, 1, S2, 1, S2, 1]


def circle3d(C, r, e1, e2):
    P = [[C[k] + r * d[0] * e1[k] + r * d[1] * e2[k] for k in range(3)] for d in CIRCLE_D]
    return P, CIRCLE_W[:]


def revolve(prof_P, prof_W, prof_p, prof_U):
    """Surface of revolution about z (same construction as in Maxima)."""
    def rot(X, k):
        x, y, z = X
        k %= 4
        return [[x, y, z], [-y, x, z], [-x, -y, z], [y, -x, z]][k]
    P, W = [], []
    for i in range(9):
        row, wrow = [], []
        for j, X in enumerate(prof_P):
            if i % 2 == 0:
                row.append(rot(X, i // 2))
            else:
                A, B = rot(X, (i - 1) // 2), rot(X, (i + 1) // 2)
                row.append([A[0] + B[0], A[1] + B[1], A[2]])
            wrow.append((1 if i % 2 == 0 else S2) * prof_W[j])
        P.append(row)
        W.append(wrow)
    return surface(2, prof_p, CIRCLE_U, prof_U, P, W), P, W


# Test curves
bez_P = [[0, 0], [1, 2], [3, 3], [4, 0]]
bezier = curve(3, [0, 0, 0, 0, 1, 1, 1, 1], bez_P)

circle = curve(2, CIRCLE_U, [d[:] for d in CIRCLE_D], CIRCLE_W)

bs_P = [[0, 0, 0], [1, 2, 1], [2, -1, 2], [4, 1, 0], [5, 3, 1], [7, 0, 2]]
bs_U = [0, 0, 0, 0, 1/3, 2/3, 1, 1, 1, 1]
bspl = curve(3, bs_U, bs_P)

rat_W = [1, 2, 0.5, 1, 3, 1]
ratc = curve(3, bs_U, bs_P, rat_W)

# Test surfaces
bsurf_P = [[[0, 0, 0], [0, 1, 1], [0, 2, 0], [0, 3, 1]],
           [[1, 0, 1], [1, 1, 2], [1, 2, 1], [1, 3, 0]],
           [[2, 0, 0], [2, 1, -1], [2, 2, 2], [2, 3, 1]],
           [[3, 0, 1], [3, 1, 0], [3, 2, 1], [3, 3, 0]]]
bsurf = surface(3, 3, [0] * 4 + [1] * 4, [0] * 4 + [1] * 4, bsurf_P)

tor_prof_P, tor_prof_W = circle3d([2, 0, 0], 0.5, [1, 0, 0], [0, 0, 1])
torus, _, _ = revolve(tor_prof_P, tor_prof_W, 2, CIRCLE_U)

sph_prof_P = [[0, 0, -1], [1, 0, -1], [1, 0, 0], [1, 0, 1], [0, 0, 1]]
sphere, _, _ = revolve(sph_prof_P, [1, S2, 1, S2, 1], 2, [0, 0, 0, .5, .5, 1, 1, 1])

PARAMS = [0.0, 0.13, 0.25, 0.37, 0.5, 0.61, 0.75, 0.88, 1.0]
UV = [(0.1, 0.2), (0.3, 0.7), (0.5, 0.5), (0.77, 0.41), (0.95, 0.05), (0.0, 0.6)]


def curve_ref(name, c, order=2):
    rows = []
    for u in PARAMS:
        d = c.derivatives(u, order=order)
        rows.append([u] + [list(v) for v in d])
    return f"ref_{name} : {mx(rows)}$\n"


def surf_ref(name, s, order=2):
    rows = []
    for (u, v) in UV:
        d = s.derivatives(u, v, order=order)
        # [u, v, S, Su, Sv, Suu, Suv, Svv]
        rows.append([u, v, d[0][0], d[1][0], d[0][1], d[2][0], d[1][1], d[0][2]])
    return f"ref_{name} : {mx(rows)}$\n"


def normals_ref(name, s):
    rows = []
    for (u, v) in UV:
        n = operations.normal(s, (u, v), normalize=True)
        rows.append([u, v, list(n[1])])
    return f"ref_{name}_normals : {mx(rows)}$\n"


out = ["/* File generated automatically by tests/compare_geomdl.py */\n",
       f"/* geomdl {geomdl.__version__}, Python {sys.version.split()[0]} */\n",
       f"ref_geomdl_version : \"{geomdl.__version__}\"$\n"]
out.append(curve_ref("bezier", bezier))
out.append(curve_ref("circle", circle))
out.append(curve_ref("bspline", bspl))
out.append(curve_ref("rational", ratc))
out.append(surf_ref("bsurf", bsurf))
out.append(surf_ref("torus", torus))
out.append(surf_ref("sphere", sphere))
out.append(normals_ref("torus", torus))
out.append(normals_ref("bsurf", bsurf))

# Global interpolation (chord length + knot averaging)
Q = [[0, 0], [1, 2], [3, 3], [4, 1], [6, 2], [7, 0]]
ic = fitting.interpolate_curve(Q, 3)
out.append(f"ref_interp_Q : {mx(Q)}$\n")
out.append(f"ref_interp_knots : {mx(list(ic.knotvector))}$\n")
out.append(f"ref_interp_ctrlpts : {mx([list(p) for p in ic.ctrlpts])}$\n")

# Knot insertion (geomdl implements algorithm A5.1)
ki = curve(3, bs_U, bs_P, rat_W)
operations.insert_knot(ki, [0.5], [2])
out.append(f"ref_insert_knots : {mx(list(ki.knotvector))}$\n")
out.append(f"ref_insert_ctrlptsw : {mx([list(p) for p in ki.ctrlptsw])}$\n")

# Arc length (geomdl.operations.length_curve)
out.append(f"ref_bspline_length : {mx(operations.length_curve(bspl))}$\n")

path = os.path.join(HERE, "reference_geomdl.mac")
with open(path, "w", encoding="utf-8") as f:
    f.writelines(out)
print("Wrote", path)

# ------------------------------------------- JSON interoperability (optional)
for fn in sys.argv[1:]:
    objs = exchange.import_json(fn)
    o = objs[0]
    if isinstance(o, (NURBS.Curve, BSpline.Curve)):
        print(fn, "-> curve degree", o.degree, " C(0.3) =", o.evaluate_single(0.3))
    else:
        print(fn, "-> surface degrees", o.degree_u, o.degree_v,
              " S(0.3,0.6) =", o.evaluate_single((0.3, 0.6)))
