#!/usr/bin/env python3
"""Exact checks for the full two-pair cubic coordinate model and safe strata.

The general statements are proved in TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md.
Finite monomial/power loops below are regression checks, not full-family
classification. The symbolic basis matrix, quadratic identities and shear
identities are exact certificates for their displayed scopes.
"""
from __future__ import annotations

import hashlib
import json
from math import factorial
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT / "research/foundations/extended-geometry/two_pair_image_cubic_structure_verification.json"
x, y, t = s.symbols("x y t")
q = 2*x*y + t*t


def odd_factorial(n: int) -> int:
    result = 1
    for j in range(1, n+1, 2):
        result *= j
    return result


def sphere(poly):
    result = 0
    for (a, b, c), coeff in s.Poly(s.expand(poly), x, y, t).terms():
        if a == b and c % 2 == 0:
            result += coeff * s.Rational(factorial(a)*odd_factorial(c-1), odd_factorial(2*a+c+1))
    return s.expand(result)


def reduce_sphere(poly):
    result = 0
    for (a, b, c), coeff in s.Poly(s.expand(poly), x, y, t).terms():
        j = min(a, b)
        result += coeff*x**(a-j)*y**(b-j)*t**c*((1-t*t)/2)**j
    return s.expand(result)


def cubic_image_monomial(i, j):
    a, d = (1+t)/2, (1-t)/2
    if i >= j:
        return s.expand(a**j*x**(i-j)*d**(3-i))
    return s.expand(a**i*(y/2)**(j-i)*d**(3-j))


def main():
    assert s.__version__ == "1.14.0", "Use pinned research/requirements.txt"
    basis = [t**j for j in range(4)]
    basis += [v**i*t**j for v in (x, y) for i in range(1, 4) for j in range(4-i)]
    images = [cubic_image_monomial(i, j) for i in range(4) for j in range(4)]
    matrix = s.Matrix([[s.Poly(p, x, y, t).coeff_monomial(b) for p in images] for b in basis])
    determinant = matrix.det()
    assert determinant != 0
    # Every original diagonal monomial contracts to i!*(3-i)!.
    for i in range(4):
        for j in range(4):
            assert sphere(cubic_image_monomial(i, j)) == (s.Rational(factorial(i)*factorial(3-i), 24) if i == j else 0)
    # Independently compare the monomial sphere rule with formal beta integration.
    v = s.symbols("v")
    for a in range(7):
        for b in range(7):
            for c in range(9):
                expected = s.integrate(((1-v*v)/2)**a*v**c, (v, -1, 1))/2 if a == b else 0
                assert sphere(x**a*y**b*t**c) == expected
                assert sphere((q-1)*x**a*y**b*t**c) == 0
    # Parity gives the unique canonical H_3+E_2 model.
    odd = [x**a*y**b*t**(3-a-b) for a in range(4) for b in range(4-a)]
    even = [x**a*y**b*t**(2-a-b) for a in range(3) for b in range(3-a)]
    parity_matrix = s.Matrix([[s.Poly(reduce_sphere(p), x, y, t).coeff_monomial(b) for p in odd+even] for b in basis])
    assert parity_matrix.det() != 0
    # The full six-parameter even family: exact trace identities, in a Cartesian frame.
    u, v, w = s.symbols("u v w")
    aa, bb, cc, dd, ee, ff = s.symbols("a b c d e f")
    a_matrix = s.Matrix([[aa, dd, ee], [dd, bb, ff], [ee, ff, cc]])
    cartesian_e = (s.Matrix([[u, v, w]])*a_matrix*s.Matrix([u, v, w]))[0]
    circular_e = s.expand(cartesian_e.subs({u:(x+y)/s.sqrt(2), v:(x-y)/(s.sqrt(2)*s.I), w:t}))
    trace = s.trace(a_matrix)
    expected_moments = [trace/3, (trace**2+2*s.trace(a_matrix**2))/15,
                        (trace**3+6*trace*s.trace(a_matrix**2)+8*s.trace(a_matrix**3))/105]
    for m, expected in enumerate(expected_moments, 1):
        assert s.expand(sphere(circular_e**m)-expected) == 0
    # Null-conic incidence for the full ten-parameter cubic. This verifies
    # the frame criterion modulo (q-1)*linear, not purity => frame existence.
    frame_s, frame_r = s.symbols("frame_s frame_r")
    h300, h210, h201, h120, h111, h102, h030, h021, h012, h003 = s.symbols(
        "h300 h210 h201 h120 h111 h102 h030 h021 h012 h003")
    cubic_h = (h300*x**3+h210*x*x*y+h201*x*x*t+h120*x*y*y+h111*x*y*t
               +h102*x*t*t+h030*y**3+h021*y*y*t+h012*y*t*t+h003*t**3)
    cone_point = {x:frame_s**2, y:-frame_r**2/2, t:frame_s*frame_r}
    sextic = s.expand(cubic_h.subs(cone_point, simultaneous=True))
    sextic_s, sextic_r = s.diff(sextic, frame_s), s.diff(sextic, frame_r)
    h_gradient = s.Matrix([s.diff(cubic_h,z).subs(cone_point, simultaneous=True) for z in (x,y,t)])
    q_gradient = s.Matrix([s.diff(q,z).subs(cone_point, simultaneous=True) for z in (x,y,t)])
    cross = h_gradient.cross(q_gradient)
    expected_cross = s.Matrix([-2*frame_s*sextic_r, -frame_r*sextic_s,
                               frame_s*sextic_s-frame_r*sextic_r])
    assert all(s.expand(entry)==0 for entry in cross-expected_cross)
    # These identities include both projective endpoints; no affine chart
    # may discard a root at infinity.
    assert s.expand(sextic_s.subs({frame_s:0,frame_r:1})-h021/4) == 0
    assert s.expand(sextic_r.subs({frame_s:0,frame_r:1})+3*h030/4) == 0
    e200, e110, e101, e020, e011, e002 = s.symbols("e200 e110 e101 e020 e011 e002")
    quadratic_e = e200*x*x+e110*x*y+e101*x*t+e020*y*y+e011*y*t+e002*t*t
    quartic = s.expand(quadratic_e.subs(cone_point, simultaneous=True))
    assert s.expand(quartic.subs({frame_s:0,frame_r:1})-e020/4) == 0
    # At the y-axis, repeated-sextic-root incidence kills y^3,y^2*t.
    # The remaining x*y^2 coefficient is removed by (q-1)*h120*y/2.
    canonical = (cubic_h+quadratic_e).subs({h030:0,h021:0,e020:0})
    affine_y = s.expand(canonical-(q-1)*h120*y/2)
    assert s.diff(affine_y,y,2) == 0
    assert reduce_sphere(affine_y-canonical) == 0

    def has_frame(h_phase, e_phase):
        binary_h = s.expand(h_phase.subs(cone_point, simultaneous=True))
        binary_e = s.expand(e_phase.subs(cone_point, simultaneous=True))
        common = s.gcd(s.diff(binary_h,frame_s),s.diff(binary_h,frame_r))
        common = s.gcd(common,binary_e)
        return common == 0 or s.Poly(common,frame_s,frame_r).total_degree() > 0

    frame_controls = {
        "repeated_sextic_root_with_quartic_root": has_frame(x**3,x*x),
        "quartic_rejects_only_repeated_root": has_frame(x**3,y*y),
        "projective_endpoint_at_infinity": has_frame(y**3,y*y),
        "squarefree_sextic_rejects_every_frame": has_frame(x**3+y**3+t**3,s.Integer(0)),
        "identically_zero_sextic_and_quartic": has_frame(q*t,q),
    }
    assert frame_controls == {
        "repeated_sextic_root_with_quartic_root": True,
        "quartic_rejects_only_repeated_root": False,
        "projective_endpoint_at_infinity": True,
        "squarefree_sextic_rejects_every_frame": False,
        "identically_zero_sextic_and_quartic": True,
    }

    # General q-preserving shear and its inverse, as literal polynomial identities.
    a = s.symbols("shear_parameter")
    forward = {x:x, t:t-a*x*x, y:y+a*x*t-a*a*x**3/2}
    inverse = {x:x, t:t+a*x*x, y:y-a*x*t-a*a*x**3/2}
    assert s.expand(q.subs(forward, simultaneous=True)-q) == 0
    assert s.expand(q.subs(inverse, simultaneous=True)-q) == 0
    for z in (x, y, t):
        assert s.expand(forward[z].subs(inverse, simultaneous=True)-z) == 0
    # Infinitesimal functional identity L(x*delta f)=0, checked on monomials.
    # In the proof this follows for every f from SO(3) invariance and delta(x)=0.
    for degree in range(8):
        for i in range(degree+1):
            for j in range(degree-i+1):
                monomial = x**i*y**j*t**(degree-i-j)
                derivative = x*(x*s.diff(monomial,t)-t*s.diff(monomial,y))
                assert sphere(derivative) == 0
    shear = forward[y].subs(a, 1)
    # Positive-degree invariant of the odd component detects semistability.
    assert sphere((q*y-x**3/2)**4) == -s.Rational(4,35)
    sharp_values = {}
    for degree in range(1, 5):
        expected = (-1)**degree*s.Rational(2**(2*degree)*factorial(3*degree)**2, factorial(6*degree+1))
        actual = sphere(y**degree*shear**(3*degree))
        assert actual == expected != 0
        sharp_values[str(degree)] = str(actual)
    for m in range(1, 17):
        assert sphere(shear**m) == 0
    # Deliberately wrong shear coefficient fails the fourth moment.
    assert sphere((y+x*t-x**3)**4) != 0
    artifact = {
        "format": "two-pair-image-cubic-structure-v1",
        "proof": "TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md",
        "coefficient_domain": "Q, symbolic extension by I and sqrt(2) for Cartesian conversion",
        "python": "3.x", "sympy": s.__version__,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "full_family": {"dimension": 16, "basis": [str(z) for z in basis],
                        "image_coordinate_determinant": str(determinant),
                        "parity_coordinate_determinant": str(parity_matrix.det()),
                        "parity_dimensions": [10, 6]},
        "exact_symbolic_certificates": ["invertible 16-dimensional coordinate maps", "three complete quadratic trace identities", "polynomial sphere shear and inverse", "null-conic frame incidence and affine representative identities"],
        "frame_controls": frame_controls,
        "bounded_regressions": {"monomial_exponent_ranges": [[0,6],[0,6],[0,8]],
                                "infinitesimal_monomial_degree_max":7,
                                "shear_pure_moment_max":16,
                                "sharp_shear_multiplier_degree_max":4,
                                "sharp_values":sharp_values},
        "scope": "Coordinate model for all 16 parameters; exact even-parity classification and shear identities as explained in proof. No full mixed-parity classification or universal finite moment cutoff is certified.",
        "outcome": "pass"
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2)+"\n")
    print("PASS full 16-parameter cubic coordinate maps and sphere moment rule")
    print("PASS exact quadratic trace identities, null-conic frame criterion, and q-preserving cubic shear")
    print("PASS bounded regressions and wrong-shear negative control; full mixed family remains open")


if __name__ == "__main__":
    main()
