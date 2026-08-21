#!/usr/bin/env sage
"""Recover one Weierstrass fiber of the resolved CM24 third-q12 pencil.

At the split prime 73 and new-base value 7, Singular's plane-curve
Brill--Noether implementation computes the conductor.  Place 6 is the
canonical simple infinity branch ``xi=-6`` (place 7 lies over the double
branch ``xi=3``).  Its Weierstrass semigroup begins ``0,2,3``.  The corresponding pole-2
and pole-3 functions satisfy a unique seven-term Weierstrass relation.  This
is a bounded finite-field conversion certificate, not yet the generic CM24
Jacobian over the new base.
"""

import argparse
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing
from sage.env import SAGE_SHARE
from sage.interfaces.singular import singular


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--new-base", type=int, default=7)
parser.add_argument("--s-root", type=int, default=33)
parser.add_argument("--j-root", type=int, default=17)
parser.add_argument(
    "--place",
    type=int,
    default=0,
    help="override the automatically selected simple xi=-6 infinity place",
)
parser.add_argument("--adjunction-only", action="store_true")
args = parser.parse_args()

load(str(HERE / "derive_q80_third_q12_cm24_pencil.sage"))

finite = GF(73)
finite_s = finite(args.s_root)
finite_j = finite(args.j_root)
assert finite_s**2 == -6 and finite_j**2 == -3
new_base_value = L(args.new_base)


def reduce_quadratic(value):
    value = K(value)
    return finite(value[0]) + finite_s * finite(value[1])


def reduce_relative(value):
    value = L(value)
    return reduce_quadratic(value[0]) + finite_j * reduce_quadratic(value[1])


finite_plane = PolynomialRing(finite, names=("w", "x"), order="lex")
finite_w, finite_x = finite_plane.gens()
finite_equation = finite_plane(0)
for x_degree, coefficient in enumerate(residual_cubic.list()):
    coefficient = new_old_base(coefficient)
    finite_coefficient = finite_plane(0)
    for w_degree, parameter_coefficient in enumerate(coefficient.list()):
        specialized = new_parameter_ring(parameter_coefficient)(new_base_value)
        finite_coefficient += reduce_relative(specialized) * finite_w**w_degree
    finite_equation += finite_coefficient * finite_x**x_degree

assert finite_equation.degree(finite_w) == 9
assert finite_equation.degree(finite_x) == 3

singular.set_ring(finite_plane._singular_())
singular_equation = finite_equation._singular_()
singular.lib("brnoeth.lib")

# Singular 4.4.1's hnoether.lib uses the local name ``a`` inside extdevelop,
# although ``a`` is also the fixed name of an algebraic-extension parameter.
# A second nonsplit local branch then fails with "object to declare is not a
# name".  Recompile that one upstream procedure with a collision-free local
# name.  The exact replacements make this fail closed if a future Singular
# version changes the procedure.
hnoether_source = Path(SAGE_SHARE) / "singular/LIB/hnoether.lib"
hnoether_text = hnoether_source.read_text()
procedure_start = hnoether_text.index("proc extdevelop (list l, int Exaktheit)")
procedure_end = hnoether_text.index("\nexample\n", procedure_start)
patched_extdevelop = hnoether_text[procedure_start:procedure_end]
renames = {
    "ideal a=hole(lastrow);": "ideal q80row=hole(lastrow);",
    "else { ideal a=lastrow; }": "else { ideal q80row=lastrow; }",
    "a[Q]=delt;": "q80row[Q]=delt;",
    "a[Q+1]=x;": "q80row[Q+1]=x;",
    "lastrow=zurueck(a);": "lastrow=zurueck(q80row);",
    "else { lastrow=a; }": "else { lastrow=q80row; }",
}
for old, new in renames.items():
    assert patched_extdevelop.count(old) == 1
    patched_extdevelop = patched_extdevelop.replace(old, new)
singular.eval("kill extdevelop;")
singular.eval(patched_extdevelop)

singular.eval(
    "printlevel=-1; "
    f"list Q80CURVE=Adj_div({singular_equation.name()}); "
    "def Q80PROJECTIVE=Q80CURVE[1][2]; setring Q80PROJECTIVE;"
)
if args.adjunction_only:
    print(
        f"Q80THIRDADJUNCTION|new_base={args.new_base}|"
        f"genus={singular.eval('Q80CURVE[2][2]')}|"
        f"places={singular.eval('string(Q80CURVE[3])')}|status=PASS_FIBER_ADJUNCTION",
        flush=True,
    )
    quit()
if args.place:
    place = args.place
else:
    # In the degree-one local ring the last two POINTS entries are the two
    # branches over the ordinary projective point at infinity.  The
    # penultimate entry is xi=-6; the final entry is the xi=3 double branch.
    # Convert its local point number back to the global CURVE[3] place index.
    singular.eval(
        "def Q80LOCAL=Q80CURVE[5][1][1]; setring Q80LOCAL; "
        "int Q80LOCALINDEX=size(POINTS)-1; setring Q80PROJECTIVE; "
        "int Q80PLACE=0; int Q80I; "
        "for (Q80I=1;Q80I<=size(Q80CURVE[3]);Q80I=Q80I+1) "
        "{ if (Q80CURVE[3][Q80I][1]==1 && "
        "Q80CURVE[3][Q80I][2]==Q80LOCALINDEX) { Q80PLACE=Q80I; } }"
    )
    place = int(singular.eval("Q80PLACE"))
    assert place > 0
singular.eval(
    f"list Q80WS=Weierstrass({place},3,Q80CURVE); "
    "int Q80S0=Q80WS[1][1]; int Q80S1=Q80WS[1][2]; int Q80S2=Q80WS[1][3]; "
    "poly Q80XN=Q80WS[2][2][1]; poly Q80XD=Q80WS[2][2][2]; "
    "poly Q80YN=Q80WS[2][3][1]; poly Q80YD=Q80WS[2][3][2]; "
    "poly Q80F=CHI;"
)
assert tuple(int(singular.eval(name)) for name in ("Q80S0", "Q80S1", "Q80S2")) == (
    0,
    2,
    3,
)

xn = singular("Q80XN").sage()
xd = singular("Q80XD").sage()
yn = singular("Q80YN").sage()
yd = singular("Q80YD").sage()
projective_equation = singular("Q80F").sage()
projective_ring = xn.parent()
assert all(value.parent() == projective_ring for value in (xd, yn, yd, projective_equation))

# Coefficients multiply Y^2, XY, Y, -X^3, -X^2, -X, -1.
relation_terms = (
    yn**2 * xd**3,
    xn * yn * xd**2 * yd,
    yn * xd**3 * yd,
    -xn**3 * yd**2,
    -xn**2 * xd * yd**2,
    -xn * xd**2 * yd**2,
    -xd**3 * yd**2,
)
remainders = tuple(term.reduce([projective_equation]) for term in relation_terms)
monomials = sorted(set().union(*(remainder.dict().keys() for remainder in remainders)))
relation_matrix = Matrix(
    finite,
    [
        [remainder.dict().get(monomial, finite(0)) for remainder in remainders]
        for monomial in monomials
    ],
)
assert relation_matrix.right_kernel().dimension() == 1
relation = relation_matrix.right_kernel().basis()[0]
is_pinned_fiber = (
    args.new_base == 7 and finite_s == finite(33) and finite_j == finite(17)
    and place == 6
)
if is_pinned_fiber:
    assert tuple(relation) == (1, 12, 43, 21, 51, 61, 33)

# If X'=cX and Y'=cY, where c is the leading cubic coefficient, the
# Weierstrass cubic becomes monic.
cubic_scale = relation[3]
a1 = relation[1]
a2 = relation[4]
a3 = relation[2] * cubic_scale
a4 = relation[5] * cubic_scale
a6 = relation[6] * cubic_scale**2
curve = EllipticCurve(finite, [a1, a2, a3, a4, a6])
assert curve.discriminant() != 0
if is_pinned_fiber:
    assert tuple((a1, a2, a3, a4, a6)) == (12, 51, 27, 40, 26)
    assert curve.discriminant() == 44
    assert curve.j_invariant() == 36
    assert curve.cardinality() == 78
    assert curve.trace_of_frobenius() == -4

status = (
    "PASS_BOUNDED_FIBER"
    if is_pinned_fiber
    else "EXPERIMENTAL_CANONICAL_FIBER"
)

print(
    f"Q80THIRDWEIERSTRASS|prime=73|s={finite_s}|jroot={finite_j}|"
    f"new_base={args.new_base}|"
    f"place={place}|infinity_branch=-6(simple)|semigroup=0,2,3|"
    f"relation={','.join(str(value) for value in relation)}|"
    f"scaled_a1,a2,a3,a4,a6={a1},{a2},{a3},{a4},{a6}|"
    f"Delta={curve.discriminant()}|j={curve.j_invariant()}|"
    f"points={curve.cardinality()}|trace={curve.trace_of_frobenius()}|"
    f"status={status}",
    flush=True,
)
