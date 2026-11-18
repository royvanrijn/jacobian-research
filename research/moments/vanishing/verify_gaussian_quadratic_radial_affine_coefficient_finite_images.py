"""Exact controls for the quadratic radial affine-coefficient theorem.

Finite controls cover original Wick identities, the full relative
connection, Gaussian sectors and finite branches. The unrestricted
classification and complete period comparison are the pinned proof.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform

import sympy as sp


s, y, t = sp.symbols("s y t")
alpha, beta, b, c = sp.symbols("alpha beta b c")
H = sp.Function("H")(y)
A = alpha*y+beta
r = 1/t


def zero(value):
    value = sp.cancel(sp.expand(value))
    if value != 0:
        value = sp.simplify(value)
    assert value == 0, value


def source_controls():
    q = 1-t*H-2*t*A*s
    P = A*s*s+H*s+b*y+c

    def DS(value):
        return sp.diff(value, s)+(-1+t*(2*A*s+H))*value

    def DY(value):
        return sp.diff(value, y)+(-y+t*sp.diff(P, y))*value

    reductions = [sp.Integer(1), s]
    for degree in range(2, 7):
        reductions.append(sp.expand(((1-t*H)*reductions[-1]
                           -(degree-1)*reductions[-2])/(2*t*A)))

    def reduce(value):
        result = 0
        for (degree,), coefficient in sp.Poly(sp.expand(value), s).terms():
            result += coefficient*reductions[degree]
        return sp.expand(result)

    def basis(value):
        value = reduce(value)
        constant = value.coeff(s, 0)
        linear = value.coeff(s, 1)
        return (sp.cancel(constant+linear*(1-t*H)/(2*t*A)),
                sp.cancel(-linear/(2*t*A)))

    for degree in range(1, 6):
        zero(reduce(DS(s**degree)))
        zero(reduce(q*s**degree)-degree*reductions[degree-1])
    zero(reduce(DY(q)-(-y+t*b)*q))
    zero(reduce(sp.diff(q,t)+P*q-(b*y+c)*q))

    ay, ky = basis(DY(sp.Integer(1)))
    at, kt = basis(P-c-b*b*t)
    m = sp.sqrt(t)*sp.sqrt(A)
    x = sp.sqrt(t)*(H-r)/(2*sp.sqrt(A))
    G = -y*y/2+b*t*y-b*b*t*t/2
    F = G-x*x
    zero(ay+alpha/(2*A)-sp.diff(F,y))
    zero(at+1/(2*t)-sp.diff(F,t))
    zero(m*ky+sp.diff(x,y))
    zero(m*kt+sp.diff(x,t))

    MY = sp.Matrix([[sp.diff(F,y),-sp.diff(x,y)],[0,sp.diff(G,y)]])
    MT = sp.Matrix([[sp.diff(F,t),-sp.diff(x,t)],[0,sp.diff(G,t)]])
    for entry in sp.diff(MY,t)-sp.diff(MT,y)+MY*MT-MT*MY:
        zero(entry)
    # The ramified cover fixes sqrt(r)=1/sqrt(t); independent principal
    # square-root simplifications would be invalid for complex t.
    zero(-1/(2*sp.sqrt(t)*sp.sqrt(A)*x)+r/(H-r))
    return {"S_exact_relations": 5,
            "retained_boundary_Y_and_t_connections": True,
            "full_triangular_connection": "nabla u=dF*u-dx*e_B; nabla e_B=dG*e_B",
            "flatness": "all four matrix entries exact",
            "endpoint_leading_sign": "-r/(H-r)"}


def wick(value):
    total = 0
    for (ns,ny), coefficient in sp.Poly(sp.expand(value),s,y).terms():
        if ny%2 == 0:
            total += coefficient*sp.factorial(ns)*sp.factorial2(ny-1)
    return sp.expand(total)


def original_controls():
    fixtures = [
        (2*y+1, y**3+y**2-y+2, 1, 2),
        (sp.Integer(3), y**4-y**2+1, -1, 0),
        (y, y**2, 0, 0),
        (y-1, (y-1)**3, 2, -3),
    ]
    markers = [sp.Integer(1),s,y,y*y,s*y,s*s]
    controls = 0
    for aa,hh,bb,cc in fixtures:
        pp = aa*s*s+hh*s+bb*y+cc
        powers = [sp.Integer(1)]
        for _ in range(4):
            powers.append(sp.expand(powers[-1]*pp))

        def coefficient(a0,a1,order):
            value = wick(a0*powers[order])/sp.factorial(order)
            if order:
                value += wick(a1*powers[order-1])/sp.factorial(order-1)
            return value

        for marker in markers:
            for order in range(5):
                # Original Gamma and ordinary-Gaussian Stein controls.
                zero(coefficient(sp.diff(s*marker,s)-s*marker,
                                 (2*aa*s+hh)*s*marker,order))
                zero(coefficient(sp.diff(marker,y)-y*marker,
                                 sp.diff(pp,y)*marker,order))
                controls += 2
        for marker in [1,y,y*y]:
            for order in range(5):
                boundary = wick(marker*(bb*y+cc)**order)/sp.factorial(order)
                zero(coefficient(marker,-(2*aa*s+hh)*marker,order)-boundary)
                controls += 1
    # The correct relative source retains the nonzero endpoint at S=0.
    assert wick(-1) == -1
    return {"phases":len(fixtures),"original_balanced_markers":len(markers),
            "moment_orders":[0,4],"exact_original_controls":controls,
            "common_A_H_roots_and_repeated_H_roots":True,
            "negative_control":"D_S(1) has endpoint -1 at t=0 and is not declared exact"}


def sector_controls():
    x,z = sp.symbols("x z")
    EE = sp.Function("E")(x)
    jp = sp.exp(-x*x)*(sp.I*sp.sqrt(sp.pi)/2-EE)
    jm = sp.exp(-x*x)*(-sp.I*sp.sqrt(sp.pi)/2-EE)
    for jj in [jp,jm]:
        zero((sp.diff(jj,x)+2*x*jj+1).subs(sp.diff(EE,x),sp.exp(x*x)))
    zero(jp-jm-sp.I*sp.sqrt(sp.pi)*sp.exp(-x*x))
    zero(jm.subs(x,-x).subs(sp.Function("E")(-x),-EE)+jp)
    # Exact sector angles are rational multiples of pi.
    plus = (Fraction(-1,6),Fraction(7,6))
    minus = (Fraction(-7,6),Fraction(1,6))
    assert plus[0] < minus[1] and minus[0] < -1 and plus[1] > 1
    assert plus[0]-Fraction(1,2) == Fraction(-2,3)
    assert plus[1]-Fraction(1,2) == Fraction(2,3)
    assert minus[0]+Fraction(1,2) == Fraction(-2,3)
    assert minus[1]+Fraction(1,2) == Fraction(2,3)
    delta = sp.pi/12
    theta = sp.pi/4-delta/2
    zero(sp.cos(2*theta)-sp.sin(delta))
    assert Fraction(2,3)-Fraction(5,24) == Fraction(11,24)
    assert Fraction(1,2)-Fraction(1,24) == Fraction(11,24)
    # The full-sector integral proof gives a uniform bound and O(x^-3)
    # remainder by explicit positive cosines. This check verifies its
    # exact trigonometric identities, not an asymptotic by sample points.
    return {"closed_cover":"[-pi/6,7pi/6] and [-7pi/6,pi/6]",
            "strict_margin":"pi/12 inside both Gaussian asymptotic sectors",
            "positive_quadratic_bound":"sin(pi/12)",
            "positive_linear_bound":"sin(pi/24)",
            "exact_jump":"i*sqrt(pi)*exp(-x^2)",
            "deck_descent":"J_minus(-x)=-J_plus(x), compensating the half-density sign",
            "ODE":"J'=-2xJ-1",
            "unrestricted_bound_proof":"rotated convergent Gaussian integral and integration by parts in written source"}


def critical_and_volume_controls():
    yy,rr,hh = sp.symbols("yy rr hh")
    leading = []
    for degree in range(2,7):
        hp = hh*yy**degree+yy+2
        for ap,label in [(alpha*yy+beta,"affine"),(beta,"constant")]:
            nn = sp.expand(sp.diff(ap,yy)*(hp-rr)**2
                  -2*ap*sp.diff(hp,yy)*(hp-rr)
                  -4*rr*ap**2*yy+4*b*ap**2)
            pol = sp.Poly(nn,yy)
            expected_degree = 2*degree if label=="affine" else 2*degree-1
            expected_lead = ((1-2*degree)*alpha*hh**2 if label=="affine"
                             else -2*degree*beta*hh**2)
            assert pol.degree() == expected_degree
            zero(pol.LC()-expected_lead)
            assert not pol.LC().has(rr)
            leading.append([degree,label,pol.degree()])
    ap = yy
    hp = yy
    bad = sp.Poly(sp.expand((hp-rr)**2-2*ap*(hp-rr)-4*rr*ap**2*yy),yy)
    assert bad.degree()==3 and bad.LC()==-4*rr
    jj = sp.symbols("j",integer=True,nonnegative=True)
    NN = sp.symbols("N",positive=True)
    kk = 2+jj
    endpoint = 1+(kk-1)*NN
    assert endpoint.is_positive
    for dd in [2*kk,4*kk-2]:
        assert (sp.Rational(1,2)+1/dd).is_positive
    return {"leading_coefficient_fixtures":leading,
            "unrestricted_leading_formulas":"(1-2k)*alpha*h_k^2 and -2k*beta*h_k^2, written degree argument",
            "bounded_Y_endpoint_order":"1/2 in Y, or 1/2+v(V) on V^2=A",
            "infinite_Y_endpoint_order":"1+(k-1)N>0 for all k>=2,N>0",
            "absolute_top_degrees":["2k (constant A)","4k-2 (affine A, V^2=A)"],
            "descending_absolute_tail_order":"1/2+1/D>0",
            "negative_control":"k=1 affine A has nonunit leading coefficient -4r"}


def finite_original_marker_controls():
    z,w = sp.symbols("z w")
    remainder = z*z*y+z**3*w
    gamma = 2
    markers = [sp.Integer(1),z,w,y,w*w,w**3,y*w*w]
    order = 6

    def full_wick(value):
        total = 0
        for (nz,nw,ny),coefficient in sp.Poly(sp.expand(value),z,w,y).terms():
            if nz==nw and ny%2==0:
                total += coefficient*sp.factorial(nz)*sp.factorial2(ny-1)
        return sp.expand(total)

    report = []
    for cc,bb in [(0,0),(3,2)]:
        phase = remainder+cc+bb*y
        powers = [sp.Integer(1)]
        for _ in range(order):
            powers.append(sp.expand(powers[-1]*phase))
        exponential = [sp.Integer(1)]
        for degree in range(1,order+1):
            value = cc*exponential[degree-1]
            if degree>=2:
                value += bb*bb*exponential[degree-2]
            exponential.append(value/degree)
        fixture = []
        for marker in markers:
            D = sp.Poly(marker,z,w,y).total_degree()
            normalized = sp.expand(sum(t**n*full_wick(
                (marker*remainder**n).subs(y,y+bb*t))/sp.factorial(n)
                for n in range(D//gamma+1)))
            for degree in range(order+1):
                predicted = sum(exponential[degree-j]*normalized.coeff(t,j)
                                for j in range(degree+1))
                zero(full_wick(marker*powers[degree])/sp.factorial(degree)-predicted)
                if bb==cc==0 and degree>D//gamma:
                    zero(full_wick(marker*powers[degree]))
            if marker==1:
                zero(normalized-1)
            bound = D//gamma if bb==0 else D+2*(D//gamma)
            assert normalized==0 or sp.degree(normalized,t)<=bound
            fixture.append({"marker":str(marker),"normalized_period":str(normalized),
                            "valid_degree_bound":bound})
        report.append({"c":cc,"b":bb,"gamma":gamma,"degree_Y_R":1,
                       "R":"Z^2Y+Z^3W","original_moment_orders":[0,order],
                       "markers":fixture})
    return report


def original_support_controls():
    report = {}
    for degree in [5,6]:
        support = {(i,j) for i in range(degree//2+1)
                   for j in range(degree-2*i+1) if i>0 or j<=1}
        expected = ({(0,0),(0,1)}
                    | {(1,j) for j in range(4 if degree==5 else 5)}
                    | {(2,j) for j in range(2 if degree==5 else 3)}
                    | ({(3,0)} if degree==6 else set()))
        assert support==expected
        report[str(degree)] = [f"S^{i}Y^{j}" for i,j in sorted(support)]
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve()
    repository = here.parents[3]
    files = [here]+[here.with_name(name) for name in [
        "GAUSSIAN_QUADRATIC_RADIAL_AFFINE_COEFFICIENT_FINITE_IMAGES.md",
        "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md",
        "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
        "GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md",
        "GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md",
        "GAUSSIAN_VARIABLE_RADIAL_QUADRATIC_FIBRE_FINITE_IMAGES.md",
        "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md",
        "GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md",
        "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md",
    ]]
    primary_base = repository/"results/gaussian-boundary"
    files += [primary_base/name for name in [
        "analytic_core.tex","seed_period_comparison.tex",
        "polynomial_stationary_balls.tex","sources/laurent_primary_manifest.json",
    ]]
    manifest = json.loads((primary_base/"sources/laurent_primary_manifest.json").read_text(encoding="utf-8"))
    for entry in manifest["sources"]:
        path = primary_base/entry["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest()==entry["sha256"],entry["file"]
        files.append(path)
    result = {
        "status":"passed",
        "checked_on":datetime.now(timezone.utc).date().isoformat(),
        "environment":{"python":platform.python_version(),"sympy":sp.__version__},
        "checks":{"source":source_controls(),"original":original_controls(),
                  "sectors":sector_controls(),"critical_and_volume":critical_and_volume_controls(),
                  "finite_original_markers":finite_original_marker_controls(),
                  "boundary_support":original_support_controls(),
                  "primary_manifest_consistency":{"status":"passed","sources":len(manifest["sources"])}},
        "sha256":{p.relative_to(repository).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        "scope":"Bounded exact original controls and source identities. All multiplicities, complete relative comparison/gluing, scalar polynomial normalization and unrestricted full finite-image classification use the pinned written proof; no finite-prefix extrapolation.",
        "independent_replay":False,"formal_verification":False,"external_review":False,
    }
    output = json.dumps(result,indent=2)+"\n"
    if args.output:
        args.output.write_text(output,encoding="utf-8")
    print(output,end="")


if __name__ == "__main__":
    main()
