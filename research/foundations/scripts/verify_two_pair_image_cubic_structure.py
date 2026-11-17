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
    # The full nine-parameter degree-at-most-two sphere family E+ell.
    # M1 removes tr(A). Use Cayley--Hamilton with tr(A)=0 to retain
    # exactly T2,T3 and Bj=b^T A^j b; homogeneous degrees are weighted
    # separately when Gaussian coefficients are turned into sphere moments.
    T2, T3, B0, B1, B2 = s.symbols("T2 T3 B0 B1 B2")
    quadratic_max_moment = 6
    traces = [s.Integer(3),s.Integer(0),T2,T3]
    vector_invariants = [B0,B1,B2]
    for n in range(4,quadratic_max_moment+1):
        traces.append(s.expand(T2*traces[n-2]/2+T3*traces[n-3]/3))
    for n in range(3,quadratic_max_moment+1):
        vector_invariants.append(s.expand(T2*vector_invariants[n-2]/2+T3*vector_invariants[n-3]/3))
    determinant_log = [s.Integer(0)]+[s.Rational(2**(n-1),n)*traces[n]
                                      for n in range(1,quadratic_max_moment+1)]
    determinant_series = [s.Integer(1)]
    for n in range(1,quadratic_max_moment+1):
        determinant_series.append(s.expand(sum(j*determinant_log[j]*determinant_series[n-j]
                                                for j in range(1,n+1))/n))

    def quadratic_series_product(left,right):
        return [s.expand(sum(left[j]*right[n-j]
                            for j in range(max(0,n-len(right)+1),min(n+1,len(left)))))
                for n in range(min(quadratic_max_moment+1,len(left)+len(right)-1))]

    vector_series = [s.Rational(2**n,2)*vector_invariants[n]
                     for n in range(quadratic_max_moment+1)]
    vector_powers = [[s.Integer(1)]]
    for k in range(quadratic_max_moment//2):
        vector_powers.append(quadratic_series_product(vector_powers[-1],vector_series))
    invariant_moments = [s.Integer(1)]
    for m in range(1,quadratic_max_moment+1):
        value = 0
        for k in range(m//2+1):
            j = m-2*k
            coefficient = quadratic_series_product(determinant_series,vector_powers[k])[j]/factorial(k)
            value += factorial(m)*coefficient/s.Integer(odd_factorial(2*j+2*k+1))
        invariant_moments.append(s.expand(value))
    invariant_solutions = {B0:-2*T2/5,B1:-4*T3/21,B2:-4*T2**2/75}
    assert s.expand(invariant_moments[2]-(5*B0+2*T2)/15) == 0
    assert s.expand(invariant_moments[3].subs(B0,invariant_solutions[B0])
                    -2*(21*B1+4*T3)/105) == 0
    assert s.expand(invariant_moments[4].subs({B0:invariant_solutions[B0],B1:invariant_solutions[B1]})
                    -16*(75*B2+4*T2**2)/2625) == 0
    reduced_fifth = -s.Rational(2048,24255)*T2*T3
    reduced_sixth = -s.Rational(8192,165540375)*(147*T2**3+500*T3**2)
    assert s.expand(invariant_moments[5].subs(invariant_solutions)-reduced_fifth) == 0
    assert s.expand(invariant_moments[6].subs(invariant_solutions)-reduced_sixth) == 0
    assert s.expand(invariant_moments[4].subs({B0:invariant_solutions[B0],
                    B1:invariant_solutions[B1],B2:-T2**2/75})) != 0

    # Independent Cartesian monomial evaluation, not the Gaussian
    # generating recurrence or the circular monomial routine.
    def cartesian_sphere(poly):
        result = 0
        for exponents,coefficient in s.Poly(s.expand(poly),u,v,w).terms():
            if any(n%2 for n in exponents):
                continue
            numerator = 1
            for n in exponents:
                numerator *= odd_factorial(n-1)
            result += coefficient*s.Rational(numerator,odd_factorial(sum(exponents)+1))
        return s.expand(result)

    linear_coefficients = s.symbols("linear_u linear_v linear_w")
    generic_linear = sum(coefficient*z for coefficient,z in zip(linear_coefficients,(u,v,w)))
    assert s.expand(cartesian_sphere(cartesian_e+generic_linear)-trace/3) == 0
    traceless_matrix = a_matrix.subs(cc,-aa-bb)
    traceless_e = s.expand(cartesian_e.subs(cc,-aa-bb))
    generic_vector = s.Matrix(linear_coefficients)
    generic_invariants = {
        T2:s.trace(traceless_matrix**2),T3:s.trace(traceless_matrix**3),
        B0:(generic_vector.T*generic_vector)[0],
        B1:(generic_vector.T*traceless_matrix*generic_vector)[0],
        B2:(generic_vector.T*traceless_matrix**2*generic_vector)[0],
    }
    for m in (2,3):
        assert s.expand(cartesian_sphere((traceless_e+generic_linear)**m)
                        -invariant_moments[m].subs(generic_invariants)) == 0
    cartesian_controls = [
        (s.Matrix([[1,2,-1],[2,3,1],[-1,1,-4]]),s.Matrix([2,-1,3])),
        (s.Matrix([[0,0,1],[0,0,s.I],[1,s.I,0]]),s.Matrix([1,s.I,0])),
        (s.zeros(3),s.Matrix([1,s.I,0])),
    ]
    cartesian_coordinates = s.Matrix([u,v,w])
    for control_matrix,control_vector in cartesian_controls:
        control_phase = (cartesian_coordinates.T*control_matrix*cartesian_coordinates
                         +control_vector.T*cartesian_coordinates)[0]
        control_invariants = {
            T2:s.trace(control_matrix**2),T3:s.trace(control_matrix**3),
            B0:(control_vector.T*control_vector)[0],
            B1:(control_vector.T*control_matrix*control_vector)[0],
            B2:(control_vector.T*control_matrix**2*control_vector)[0],
        }
        for m in range(1,quadratic_max_moment+1):
            assert s.expand(cartesian_sphere(control_phase**m)
                            -invariant_moments[m].subs(control_invariants)) == 0
    # Five moments genuinely fail to classify this same nine-dimensional
    # family. The displayed exact complex-coefficient phase first fails at 6.
    sharp_quadratic_phase = u*u-v*v+2*s.I*(s.sqrt(6)*(u+v)+s.sqrt(33)*w)/15
    sharp_quadratic_moments = [cartesian_sphere(sharp_quadratic_phase**m)
                              for m in range(1,quadratic_max_moment+1)]
    sharp_quadratic_value = reduced_sixth.subs({T2:2,T3:0})
    assert sharp_quadratic_moments == [0,0,0,0,0,sharp_quadratic_value]
    assert sharp_quadratic_value != 0


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

    # Graph compactification near the null conic. These are algebraic
    # identities for all sixteen coefficients; the cited tameness and
    # stationary-phase theorems are not computer-verified by this script.
    boundary_u, boundary_z, boundary_c = s.symbols("boundary_u boundary_z boundary_c")
    projective_point = {x:boundary_u**2-boundary_z**2, y:-s.Rational(1,2), t:boundary_u}
    null_point = {x:boundary_u**2, y:-s.Rational(1,2), t:boundary_u}
    affine_point = {z:s.cancel(value/boundary_z) for z,value in projective_point.items()}
    assert s.cancel(q.subs(affine_point,simultaneous=True)-1) == 0
    boundary_f = s.expand(cubic_h.subs(null_point,simultaneous=True))
    boundary_g = s.expand(quadratic_e.subs(null_point,simultaneous=True))
    boundary_equation = s.expand(cubic_h.subs(projective_point,simultaneous=True)
                                 +boundary_z*quadratic_e.subs(projective_point,simultaneous=True)
                                 -boundary_c*boundary_z**3)
    boundary_h = s.cancel((boundary_equation-boundary_f-boundary_z*boundary_g
                           +boundary_c*boundary_z**3)/boundary_z**2)
    assert s.denom(boundary_h) == 1
    s.Poly(boundary_h,boundary_u,boundary_z)
    assert s.expand(boundary_equation-(boundary_f+boundary_z*boundary_g
                    +boundary_z**2*boundary_h-boundary_c*boundary_z**3)) == 0
    assert s.cancel((cubic_h+quadratic_e).subs(affine_point,simultaneous=True)
                    -(boundary_equation+boundary_c*boundary_z**3)/boundary_z**3) == 0
    assert s.expand(s.diff(boundary_equation,boundary_u).subs(boundary_z,0)
                    -s.diff(boundary_f,boundary_u)) == 0
    assert s.expand(s.diff(boundary_equation,boundary_z).subs(boundary_z,0)-boundary_g) == 0
    assert s.diff(boundary_equation,boundary_c) == -boundary_z**3
    boundary_unit = boundary_g+boundary_z*boundary_h-boundary_c*boundary_z**2
    assert s.expand(boundary_equation-boundary_f-boundary_z*boundary_unit) == 0
    # The residue dy wedge dt/(2y) has precisely a second-order boundary pole.
    boundary_jacobian = s.det(s.Matrix([[s.diff(affine_point[y],boundary_u),s.diff(affine_point[y],boundary_z)],
                                      [s.diff(affine_point[t],boundary_u),s.diff(affine_point[t],boundary_z)]]))
    assert s.cancel(boundary_jacobian/(2*affine_point[y])-1/(2*boundary_z**2)) == 0


    # Ordinary boundary nodes in the thirteen-parameter frame family.
    # The displayed identities support the two blow-up chart proofs;
    # their analytic product-pair consequence is a written argument.
    node_u,node_z,node_s,node_v,node_r,node_k=s.symbols("node_u node_z node_s node_v node_r node_k")
    profile_a,profile_b,profile_c=s.symbols("profile_a profile_b profile_c")
    g_exponents=[(i,j) for i in range(4) for j in range(4-i)]
    g_coefficients=s.symbols("g00 g01 g02 g03 g10 g11 g12 g20 g21 g30")
    g_polynomial=sum(coefficient*x**i*t**j for (i,j),coefficient in zip(g_exponents,g_coefficients))
    frame_phase=y*(profile_a+profile_b*t+profile_c*t*t)+g_polynomial
    node_substitution={x:(node_z**2-node_u**2)/(2*node_z),y:1/node_z,t:node_u/node_z}
    assert s.cancel(q.subs(node_substitution,simultaneous=True)-1)==0
    node_equation=s.cancel(node_z**3*frame_phase.subs(node_substitution,simultaneous=True)-node_s*node_z**3)
    s.Poly(node_equation,node_u,node_z)
    node_quadratic=profile_c*node_u**2+profile_b*node_u*node_z+profile_a*node_z**2
    assert all(sum(powers)>=3 for powers,coefficient in s.Poly(s.expand(node_equation-node_quadratic),node_u,node_z).terms())
    z_chart=s.cancel(node_equation.subs(node_u,node_v*node_z)/node_z**2)
    u_chart=s.cancel(node_equation.subs(node_z,node_v*node_u)/node_u**2)
    s.Poly(z_chart,node_v,node_z)
    s.Poly(u_chart,node_v,node_u)
    assert s.expand(z_chart.subs(node_z,0)-(profile_c*node_v**2+profile_b*node_v+profile_a))==0
    assert s.expand(u_chart.subs(node_u,0)-(profile_c+profile_b*node_v+profile_a*node_v**2))==0
    assert s.expand(s.diff(z_chart,node_s)+node_z)==0
    assert s.expand(s.diff(u_chart,node_s)+node_u*node_v**3)==0
    # Repeated profiles retain their shape and have a uniform height-cubic normalization.
    repeated_phase=y*(t-node_r)**2+g_polynomial
    rotated=reduce_sphere(repeated_phase.subs({x:x,t:t+node_k*x,y:y-node_k*t-node_k**2*x/2},simultaneous=True))
    assert s.Poly(rotated,x,y,t).coeff_monomial(t**3)==g_coefficients[3]-2*node_k
    assert s.expand(s.diff(rotated,y)-(t-node_r)**2)==0

    # Complete repeated-profile total-space singularity: one ordinary A1 point.
    repeated_graph=s.expand(node_equation.subs({profile_a:node_r**2,profile_b:-2*node_r,profile_c:1}))
    repeated_strict=s.cancel(repeated_graph.subs(node_u,node_z*(node_v+node_r))/node_z**2)
    assert s.expand(repeated_strict.subs(node_z,0)-node_v**2)==0
    assert s.diff(repeated_strict,node_s)==-node_z
    critical_target=s.expand(s.diff(repeated_strict,node_z).subs({node_v:0,node_z:0})+node_s)
    critical_point={node_v:0,node_z:0,node_s:critical_target}
    assert all(s.expand(s.diff(repeated_strict,w).subs(critical_point))==0 for w in (node_v,node_z,node_s))
    ordinary_hessian=s.hessian(repeated_strict,(node_v,node_z,node_s)).subs(critical_point)
    assert s.expand(ordinary_hessian.det())==-2
    # Volume/residue identities in both section-blow-up charts.
    assert s.cancel((1/(2*node_z))/(1/s.diff(z_chart,node_s)))==-s.Rational(1,2)
    assert s.cancel((1/(2*node_u*node_v**2))/(1/s.diff(u_chart,node_s)))==-node_v/2
    # A point blow-up in a smooth threefold has volume factor lambda^2.
    # The ordinary double point equation has the same factor in all charts.
    odp_lam,odp_a,odp_b,odp_sigma=s.symbols("odp_lam odp_a odp_b odp_sigma")
    centered=s.expand(repeated_strict.subs(node_s,odp_sigma+critical_target))
    ambient_coordinates=(node_v,node_z,odp_sigma)
    for exceptional_coordinate in range(3):
        values=[None,None,None]
        values[exceptional_coordinate]=odp_lam
        remaining=[j for j in range(3) if j!=exceptional_coordinate]
        values[remaining[0]]=odp_lam*odp_a
        values[remaining[1]]=odp_lam*odp_b
        jacobian=s.det(s.Matrix(values).jacobian((odp_lam,odp_a,odp_b)))
        assert s.cancel(jacobian**2/odp_lam**4)==1
        strict_equation=s.cancel(centered.subs(dict(zip(ambient_coordinates,values)),simultaneous=True)/odp_lam**2)
        assert s.denom(strict_equation)==1
        s.Poly(strict_equation,odp_lam,odp_a,odp_b)

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
        "exact_symbolic_certificates": ["invertible 16-dimensional coordinate maps", "three complete quadratic trace identities", "polynomial sphere shear and inverse", "null-conic frame incidence and affine representative identities", "boundary graph Laurent expansion, local jets, and residue form", "full nine-parameter quadratic invariant moments through six", "both boundary blow-up charts and residue multipliers", "universal A1 Hessian determinant -2 and all three crepant point-blow-up charts"],
        "frame_controls": frame_controls,
        "quadratic_sphere_family": {"dimension":9,"exact_pure_prefix":6,
                                    "M5_after_M1_to_M4":str(reduced_fifth),
                                    "M6_after_M1_to_M4":str(reduced_sixth),
                                    "sharp_phase_cartesian":str(sharp_quadratic_phase),
                                    "sharp_moments_1_through_6":[str(value) for value in sharp_quadratic_moments],
                                    "independent_controls":"Cartesian symbolic M1..M3, three exact controls through M6, and altered B2 rejection"},
        "bounded_regressions": {"monomial_exponent_ranges": [[0,6],[0,6],[0,8]],
                                "infinitesimal_monomial_degree_max":7,
                                "shear_pure_moment_max":16,
                                "sharp_shear_multiplier_degree_max":4,
                                "sharp_values":sharp_values},
        "scope": "Coordinate model for all 16 parameters; exact nine-parameter quadratic invariant moments, even-parity trace, frame-incidence, boundary graph, and shear identities as explained in proof. The cited V-filtration theorems and the geometric full-family argument are written proofs, not computer-certified. No universal finite pure-moment cutoff is certified by this script.",
        "outcome": "pass"
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2)+"\n")
    print("PASS full 16-parameter cubic coordinate maps and sphere moment rule")
    print("PASS full nine-parameter quadratic six-moment certificate and sharpness")
    print("PASS exact quadratic trace, null-conic frame, boundary graph, and cubic shear identities")
    print("PASS exact boundary and A1 identities, bounded regressions, and rejection controls")


if __name__ == "__main__":
    main()
