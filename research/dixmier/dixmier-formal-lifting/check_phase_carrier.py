#!/usr/bin/env python3
"""QPH1 exact phase, matrix, admission and Riccati controls."""
from collections import defaultdict
from fractions import Fraction as F
from math import comb
from pathlib import Path
import hashlib
import importlib.util
import json
import sympy as s

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE/file)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def bind(name):
    receipt = json.loads((HERE/name).read_text())
    assert receipt["status"] == "PASS"
    for rel, expected in receipt["hashes"].items():
        path = HERE/rel
        if not path.exists():
            path = HERE.parent/"dixmier-audit"/rel
        assert path.is_file() and sha(path) == expected, (name, rel)
    return receipt


def mon(r=0, v=0, a=0, h=0, c=1):
    return {(r, v, a, h): F(c)} if c else {}


def add(*polys):
    out = defaultdict(F)
    for poly in polys:
        for e, c in poly.items():
            out[e] += c
    return {e: c for e, c in out.items() if c}


def scale(poly, c):
    return {e: v*c for e, v in poly.items() if v*c}


def trunc(poly, n):
    return {e: c for e, c in poly.items() if e[3] <= n}


def mul(f, g, n=None):
    out = defaultdict(F)
    for (r, v, a, hh), c in f.items():
        for (rr, vv, aa, hhh), cc in g.items():
            factor = 1
            for k in range(a+1):
                if k:
                    factor *= vv-2*(k-1)
                order = hh+hhh+k
                if n is not None and order > n:
                    continue
                if factor:
                    out[(r+rr, v+vv-2*k, a+aa-k, order)] += (
                        c*cc*comb(a, k)*factor)
    return {e: c for e, c in out.items() if c}


def comm(f, g, n=None):
    return add(mul(f, g, n), scale(mul(g, f, n), -1))


def dr(poly):
    return {(r-1, v, a, hh): c*r
            for (r, v, a, hh), c in poly.items() if r}


def coeff(poly, n):
    return {(r, v, a, 0): c
            for (r, v, a, hh), c in poly.items() if hh == n}


def expr(poly, r, v, a, h):
    return s.Add(*(s.Rational(c.numerator, c.denominator) *
                   r**er*v**ev*a**ea*h**eh
                   for (er, ev, ea, eh), c in poly.items()))


def main():
    bind("painleve_long_controls.json")
    bind("matrix_carrier_controls.json")
    mc = module("phase_matrix_arithmetic", "check_matrix_carrier.py")
    q, r, a, d, h, v, ss = s.symbols("q r a d h v ss")
    mc.x, mc.y, mc.p, mc.z, mc.h = q, r, a, d, h
    mc.positions, mc.momenta = (q, r), (a, d)
    g = 4*a-3*r*a**2
    t0 = a**2/2-r*a**3/4
    c = 3*a**4/4-s.Rational(9, 20)*r*a**5
    c0 = 8*q-g**2
    phase = (-a*ss**2/8+a**2*ss/2-r*a**3*ss/4-2*a**3/3
             +3*r*a**4/4-s.Rational(9, 40)*r**2*a**5)
    assert s.expand(s.diff(phase, a)+(ss-g)**2/8) == 0
    assert s.expand(s.diff(phase, ss)-(t0-a*ss/4)) == 0
    assert s.expand(s.diff(phase, r)-(c-a**3*ss/4)) == 0
    S = [[0, 1], [c0, 2*g]]
    T = [[t0, -a/4], [h-mc.mul(a, c0)/4, t0-a*g/2]]
    C = [[c, -a**3/4],
         [-mc.mul(a**3, c0)/4, c-a**3*g/2]]
    R = mc.diag(r)
    D = mc.add(mc.diag(d), C, -1)
    one = mc.diag(s.Integer(1))
    images = (R, S, T, D)
    relations = 0
    for i in range(4):
        for j in range(i+1, 4):
            expected = -1 if (i, j) == (0, 3) else 1 if (i, j) == (1, 2) else 0
            mc.zero(mc.add(mc.comm(images[i], images[j]),
                           mc.scale(expected*h, one), -1),
                    f"full matrix relation {i},{j}")
            relations += 1
    # One declared higher-order deformation, not a degree search.
    g_alt = g+h*a**3
    integral_g = s.integrate(g_alt, a)
    integral_square = s.integrate(s.expand(g_alt**2), a)
    t_alt = integral_g/4
    z_alt = s.diff(integral_g, r)/4
    c_alt = -s.diff(integral_square, r)/8
    b_alt = 8*q-g_alt**2
    S_alt = [[0, 1], [b_alt, 2*g_alt]]
    T_alt = [[t_alt, -a/4],
             [h-mc.mul(a, b_alt)/4, t_alt-a*g_alt/2]]
    C_alt = [[c_alt, z_alt],
             [mc.mul(z_alt, b_alt), c_alt+2*z_alt*g_alt]]
    D_alt = mc.add(mc.diag(d), C_alt, -1)
    alt_images = (R, S_alt, T_alt, D_alt)
    alt_relations = 0
    for i in range(4):
        for j in range(i+1, 4):
            expected = -1 if (i, j) == (0, 3) else 1 if (i, j) == (1, 2) else 0
            mc.zero(mc.add(mc.comm(alt_images[i], alt_images[j]),
                           mc.scale(expected*h, one), -1),
                    f"higher-order phase relation {i},{j}")
            alt_relations += 1

    xi = s.symbols("xi")
    leadS = s.Matrix([[0, 1], [-9*r**2*xi**4, -6*r*xi**2]])
    leadT = s.Matrix([[-r*xi**3/4, -xi/4],
                      [s.Rational(9, 4)*r**2*xi**5, 5*r*xi**3/4]])
    assert s.expand(leadS.det()-9*r**2*xi**4) == 0
    assert s.expand(leadT.det()-r**2*xi**6/4) == 0

    # Positive split control g=0 after q=v^2/2.
    # Replace q by v and momentum a by a*v: [a*v,v]=h,
    # so the abstract a is represented as v^-1*(a*v).
    V, PP = s.symbols("V PP")
    mc.x, mc.y, mc.p, mc.z = V, r, PP, d
    mc.positions, mc.momenta = (V, r), (PP, d)
    aa = PP/V
    splitS = [[0, 1], [4*V**2, 0]]
    splitT = [[0, -aa/4], [h-mc.mul(aa, V**2), 0]]
    projector = [[s.Rational(1, 2), 1/(4*V)],
                 [V, s.Rational(1, 2)]]
    mc.zero(mc.add(mc.mm(projector, projector), projector, -1),
            "constant-g projector")
    mc.zero(mc.comm(projector, splitS), "constant-g projector S")
    mc.zero(mc.comm(projector, splitT), "constant-g projector T")

    scalarS = aa+2*V
    scalarT = -mc.mul(scalarS, scalarS)/8+V**2/2
    assert s.expand(mc.mul(scalarS, scalarT)-mc.mul(scalarT, scalarS)-h) == 0

    # Rational canonical chart and actual corrected classical fourth image.
    X, Q, P, Z = s.symbols("X Q P Z")
    rr = X*(2-3*X*Q)
    vv = 1/X
    ff = 3*vv-s.Rational(5, 3)*rr*vv**2+s.Rational(1, 3)*rr**2*vv**3
    gg = -s.Rational(5, 12)*vv**4+s.Rational(2, 15)*rr*vv**5
    ac = -X**3*P-s.Rational(2, 3)*X*(1-3*X*Q)*Z+ff
    dc = -Z/(3*X**2)+gg
    def poisson(f, g):
        return s.cancel(s.diff(f, P)*s.diff(g, X)-s.diff(f, X)*s.diff(g, P)
                        +s.diff(f, Z)*s.diff(g, Q)-s.diff(f, Q)*s.diff(g, Z))
    chart = (rr, vv**2/2, ac, dc)
    for i in range(4):
        for j in range(i+1, 4):
            expected = -1 if (i, j) in ((0, 3), (1, 2)) else 0
            assert poisson(chart[i], chart[j]) == expected
    yy = a-v
    beta = 5*v**2-r*v**3-3*a*v
    u = yy/v
    correction = (yy**4*(18*u**2+78*u+125)/20
                  +3*beta*yy**2*(u**3+5*u**2+10*u-5)/10
                  -beta**2*(9*u+2)/6-beta**3/(6*v**2))
    ff = 3*v-s.Rational(5, 3)*r*v**2+s.Rational(1, 3)*r**2*v**3
    gg = -s.Rational(5, 12)*v**4+s.Rational(2, 15)*r*v**5
    diff = (-gg-v**3*(3-r*v)*(a-ff)/2+correction
            -v*a**3/2-a**4/4+s.Rational(3, 10)*r*a**5)
    assert s.cancel(diff) == 0
    assert s.limit(X*2*(4*ac-3*rr*ac**2), X, 0) == -4

    # Four recurrence coefficients; this is a fixed identity control.
    # Normal monomials are r^i*v^j*a^k*h^l with [a,v]=h/v.
    gg = add(mon(a=1, c=4), mon(r=1, a=2, c=-3))
    L = mon(v=1, c=2)
    terms = [L]
    top_samples = []
    eta, chi = s.symbols("eta chi")
    for n in range(1, 5):
        error = add(mul(L, L, n), comm(gg, L, n), mon(v=2, c=-4))
        forcing = coeff(error, n)
        Ln = {(er, ev-1, ea, 0): -cc/4
              for (er, ev, ea, eh), cc in forcing.items()}
        shifted = {(er, ev, ea, n): cc
                   for (er, ev, ea, eh), cc in Ln.items()}
        L = add(L, shifted)
        terms.append(shifted)
        assert not trunc(add(mul(L, L, n), comm(gg, L, n),
                             mon(v=2, c=-4)), n)
        pn = 0
        for (er, ev, ea, eh), cc in Ln.items():
            power_chi = er-ea
            assert power_chi >= 0
            assert ev == 1-3*n+power_chi
            assert ea+power_chi <= n
            pn += s.Rational(cc.numerator, cc.denominator)*eta**ea*chi**power_chi
        lead = s.Poly(pn, eta, chi).coeff_monomial(eta**n)
        top_samples.append(str(lead))
    assert s.expand(expr(coeff(L, 1), r, v, a, h)-(3*r*a-2)/v**2) == 0
    assert top_samples[:3] == ["3", "-45/4", "405/4"]
    SS = add(gg, L)
    TT = add(mon(r=1, a=3, c=F(1, 2)), mon(a=2, c=-F(1, 2)),
             scale(mul(L, mon(a=1)), -F(1, 4)))
    DD = add(mon(a=4, c=F(1, 4)), mon(r=1, a=5, c=-F(3, 10)),
             scale(mul(L, mon(a=3)), F(1, 4)))
    assert not add(comm(SS, TT, 5), mon(h=1, c=-1))
    for target in (SS, TT):
        derivative = {(er, ev, ea, eh+1): cc
                      for (er, ev, ea, eh), cc in dr(target).items()}
        assert not trunc(add(derivative, comm(DD, target, 5)), 5)

    # Complete finite admission seed, before any all-order inference.
    pi = module("phase_admission_arithmetic", "check_painleve_long.py")
    aa = pi.add(pi.mon(b=-1, k=3), pi.mon(a=1, b=-2, k=-F(5, 3)),
                pi.mon(a=2, b=-3, k=F(1, 3)), pi.mon(b=1, c=1, k=-F(1, 3)))
    def substitute(poly):
        out = {}
        for (er, ev, ea, eh), cc in poly.items():
            out = pi.add(out, pi.mul(pi.mon(a=er, b=-ev, d=eh, k=cc),
                                    pi.power(aa, ea)))
        return out
    baseL = add(terms[0], terms[1])
    seedS = pi.add(substitute(gg), substitute(baseL))
    seedT = substitute(add(mon(r=1, a=3, c=F(1, 2)), mon(a=2, c=-F(1, 2)),
                          scale(mul(baseL, mon(a=1)), -F(1, 4))))
    d_minus_d0 = pi.add(pi.mon(b=-4, k=-F(5, 12)),
                        pi.mon(a=1, b=-5, k=F(2, 15)),
                        pi.mon(b=-2, c=1, k=-F(1, 2)),
                        pi.mon(a=1, b=-3, c=1, k=F(1, 6)))
    seedD = pi.add(d_minus_d0, substitute(add(
        mon(a=4, c=F(1, 4)), mon(r=1, a=5, c=-F(3, 10)),
        scale(mul(baseL, mon(a=3)), F(1, 4)))))
    seed = {"S": seedS, "T": seedT, "D_minus_D0": seedD}
    global_seed = {name: pi.globalize(poly) for name, poly in seed.items()}
    retained = json.loads((HERE/"painleve_long_seed.json").read_text())
    for name, poly in global_seed.items():
        old = {tuple(e): F(c) for e, c in retained["global"][name]}
        assert {e: c for e, c in old.items() if e[3] == 0} == {
            e: c for e, c in poly.items() if e[3] == 0}, name+" classical seed"
    seed_receipt = {
        "claim": "QPH1", "scope": "complete base plus h*L1 admission seed",
        "local_exponent_order": ["R", "X", "B", "h"],
        "global_exponent_order": ["X", "Q", "B", "h"],
        "local": {name: pi.rows(poly) for name, poly in seed.items()},
        "global": {name: pi.rows(poly) for name, poly in global_seed.items()}}
    (HERE/"phase_long_seed.json").write_text(json.dumps(seed_receipt, indent=2)+"\n")
    growth_samples = []
    for n in range(2, 5):
        newterm = pi.globalize(substitute(terms[n]))
        leading = {e: c for e, c in newterm.items() if e[3] == n}
        assert max(e[2] for e in leading) == n
        assert max(e[0]+e[1]+3*e[2] for e in leading) == 10*n-1
        growth_samples.append({"n": n, "leading_momentum_order": n,
                               "leading_spatial_degree": 10*n-1})

    names = ["QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md", "check_phase_carrier.py",
             "phase_long_seed.json", "check_matrix_carrier.py",
             "matrix_carrier_controls.json",
             "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
             "ALGEBRAIC_CURVE_STRING_OBSTRUCTION.md",
             "PAINLEVE_GLOBAL_LONG_LIFT.md", "check_painleve_long.py",
             "painleve_long_seed.json", "painleve_long_controls.json"]
    receipt = {
        "claim": "QPH1", "status": "PASS",
        "phase_derivative_identities": 3, "full_matrix_relations": relations,
        "higher_degree_g_plus_h_a3_full_relations": alt_relations,
        "alternative_linear_g_scalar_quantization": "PASS",
        "exact_rational_chart_relations": 6,
        "actual_corrected_classical_fourth_image": "PASS",
        "positive_branch_classical_seed_matches_QPL1": "PASS",
        "negative_branch_trace_X_residue": -4,
        "constant_g_algebraic_split_control": "PASS",
        "weighted_symbol_determinants": ["9*r^2*xi^4", "r^2*xi^6/4"],
        "Riccati_coefficients_checked": 4,
        "Riccati_highest_momentum_coefficients": top_samples,
        "scalar_full_normalized_relations_through_h4": "PASS",
        "finite_global_admission_seed_terms": {
            name: len(poly) for name, poly in global_seed.items()},
        "leading_growth_samples": growth_samples,
        "hashes": {name: sha(HERE/name) for name in names},
        "research_context": [
            "https://msp.org/pjm/2011/254-1/pjm-v254-n1-p08-p.pdf",
            "https://arxiv.org/html/math/0003009v2",
            "https://arxiv.org/html/0712.4124"],
        "coverage_boundary": (
            "Finite exact phase, matrix, symbol, chart, admission seed and selected "
            "Riccati identities only. All-degree freeness, algebraic nonsplitting "
            "including arbitrary nonconstant polynomial g deformations and coefficient "
            "covers, formal diagonalization, global admission "
            "and sharp degree/order growth use written proofs. No finite-left "
            "rank-one Long carrier or matrix endomorphism counterexample is constructed."),
        "independent_replay": False, "formal_verification": False, "external_review": False}
    (HERE/"phase_carrier_controls.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"claim": "QPH1", "status": "PASS",
                      "full_matrix_relations": relations,
                      "finite_admission_seed_terms":
                          receipt["finite_global_admission_seed_terms"],
                      "growth_samples": growth_samples}, indent=2))


if __name__ == "__main__":
    main()
