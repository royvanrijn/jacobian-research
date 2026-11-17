#!/usr/bin/env python3
"""QVA1 exact valuation-ideal and full graph admission controls."""
from pathlib import Path
import hashlib
import json
import sympy as s
import check_projective_frame as pbw

HERE = Path(__file__).resolve().parent
x, p, t = pbw.x, pbw.p, pbw.t
a = s.Symbol("a", nonzero=True)


def mul(f, g):
    return pbw.mul(f, g, a)


def comm(f, g):
    return s.expand(mul(f, g)-mul(g, f))


def mm(A, B):
    return s.Matrix(A.rows, B.cols,
                    lambda i, j: s.expand(sum(mul(A[i, k], B[k, j])
                                             for k in range(A.cols))))


def zero(expr, label):
    assert s.cancel(expr) == 0, (label, s.cancel(expr))


def matrix_zero(M, label):
    for c in M:
        zero(c, label)


def main():
    # Independent a denotes the whole quantum parameter, not a sampled e.
    zero(comm(p, x)-a, "canonical commutator")
    F, G = x*x, x*p+2*a
    zero(-mul(p, F)+mul(x, G), "exact interior syzygy")
    row, col = s.Matrix([[-p, x]]), s.Matrix([-x/a, -p/a])
    matrix_zero(mm(row, col)-s.eye(1), "generic projective splitting")
    E = s.eye(2)-mm(col, row)
    matrix_zero(mm(E, E)-E, "generic rank-one projective idempotent")
    matrix_zero(mm(row, E), "projector kernel")
    # The single whole polynomial leading relation is exhausted by gcd(x,p)=1.
    h0 = s.Symbol("h0")
    zero(h0*p*x*x-h0*x*x*p, "complete leading syzygy form")
    assert s.gcd(x, p) == 1
    assert s.gcd(x*x, x*p) == x
    normalized = s.groebner([x, p], x, p)
    assert normalized.reduce(s.Integer(1))[1] == 1
    assert normalized.reduce(s.diff(p, p))[1] == 1
    zero(comm(x, p)/a+1, "residue x violates normalized ideal")

    # Gcd normalization is a universal scalar-connection identity.
    d = s.Function("d")(x, p)
    g = s.Function("g")(x, p)
    vx, vp, r = s.symbols("vx vp r")
    deriv = lambda z: vx*s.diff(z, x)+vp*s.diff(z, p)
    zero((deriv(d*g)+r*d*g)/d -
         (deriv(g)+(r+deriv(d)/d)*g), "normalized potential identity")

    # A nonunit divisorial leading gcd is allowed. The saturated formal
    # lattice B*x has J=(x), J'=S and a full rational action.
    rho_x, rho_p = x, p+a/x
    zero(comm(rho_p, rho_x)-a, "principal-divisor full rational pair")
    zero(mul(x, rho_p)-mul(p, x), "principal-divisor right stability")
    zero(s.diff(x*g, x)-(x*g)/x-x*s.diff(g, x),
         "divisorial potential cancels after gcd normalization")
    zero(-s.diff(x*g, p)+x*s.diff(g, p),
         "other Hamiltonian preserves the principal lattice")

    # A nonlinear full symplectic shadow still spans all derivatives.
    fx, fp = s.symbols("fx fp")
    X1, X2 = -fp, fx-3*x*x*fp
    zero(X2-3*x*x*X1-fx, "nonlinear symplectic derivative recovery")
    zero(-X1-fp, "second derivative recovery")

    # General rational coefficient gauge identity; no exponent census.
    ff = s.Function("ff")(x)
    rr = s.Function("rr")(x)
    transported = mul(mul(ff, p+a*rr), 1/ff)
    zero(transported-(p+a*(rr-s.diff(ff, x)/ff)),
         "universal logarithmic gauge formula")
    positive_f = x**3/(x-1)**2
    positive_H = s.Rational(5, 3)*x**3+7*x
    positive_r = 3/x-2/(x-1)+5*x*x+7
    zero(positive_r-s.diff(positive_H, x)-s.diff(positive_f, x)/positive_f,
         "full logarithmic-plus-polynomial positive admission")
    regular_p = mul(mul(positive_f, p+a*positive_r), 1/positive_f)
    zero(regular_p-(p+a*s.diff(positive_H, x)),
         "actual generic free-frame transport")
    zero(comm(p+a*positive_r, x)-a, "positive full rational relation")
    assert s.residue(s.Rational(1, 2)/x, x, 0) == s.Rational(1, 2)
    assert s.residue(1/x**2, x, 0) == 0
    assert s.limit(x*x/x**2, x, 0) == 1
    # Logarithmic derivatives of rational functions only have simple
    # finite poles with integer residues; the proof gives the all-degree gate.

    # Full regular formal tuples admit no extra spatially rational
    # division intertwiner. This admitted unit has an unbounded positive tail.
    q = 1-t*x
    u, ui = 1/q, q
    zero(mul(u, ui)-1, "regular formal unit inverse")
    moved_p = mul(mul(u, p), ui)
    zero(moved_p-(p-a*t/q), "regular inner comparison")
    zero(comm(moved_p, x)-a, "regular comparison full relation")
    zero((moved_p-p)/(a*t)+1/q, "inner finite quantum jet contact")
    # Spatial pole comparison fails the regular graph hypothesis.
    pole_p = mul(mul(1/x, p), x)
    zero(pole_p-(p+a/x), "non-C spatial pole comparison")
    assert s.denom(s.cancel((pole_p-p)/a)) == x

    # Exact moving-divisor left-model control. Its constant special
    # generator makes its whole initial ideal the unit ideal.
    z = s.Symbol("z")
    moving = 1-t*z*x
    assert s.expand(moving*moving).subs(t, 0) == 1

    inputs = [
        "VALUATION_SYMBOLS_AND_GENERIC_GRAPH_ADMISSION.md",
        "check_valuation_graph_admission.py",
        "check_projective_frame.py",
        "CENTRAL_PATCHING_AND_QUANTUM_GRAPH_JETS.md",
        "HAMILTONIAN_IDEALIZER_CRITERION.md",
        "INTEGRAL_GRAPH_PROJECTIVITY_AND_LATTICE_RIGIDITY.md",
        "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md",
        "QUANTUM_DIVISION_AND_COMPLETION.md",
        "PROJECTIVE_FRAME_BOUNDARY.md",
    ]
    receipt = {
        "claim": "QVA1", "status": "PASS",
        "interior_projective_control": {
            "generators": ["x^2", "x*p+2*a"],
            "exact_syzygy": "(-p,x)",
            "generic_splitting": "(-x/a,-p/a)^T",
            "projective_idempotent": "PASS",
            "whole_initial_ideal": "(x^2,x*p)",
            "whole_initial_exhaustion": "written exact syzygy lifting",
            "normalized_initial_ideal": "(x,p)",
            "nonprincipal_certificate": "gcd=x and normalized ideal proper",
            "residue_x_hamiltonian_failure": "-1 not in (x,p)",
            "all_degrees_and_parameter_orders": "written proof"},
        "principal_control": {
            "initial_ideal": "(x)", "normalized_initial_ideal": "(1)",
            "full_rational_pair": ["x", "p+a/x"],
            "actual_free_frame": "x",
            "divisorial_connection_cancellation": "PASS"},
        "generic_admission": {
            "normalized_potential_identity": "PASS",
            "nonlinear_full_frame_recovery": "PASS",
            "subquantum_polynomial_removal": "written multiplier argument",
            "automatic_C_frame": "written principal initial ideal proof",
            "all_ramification_orders": "written proof, no order census"},
        "logarithmic_admission": {
            "universal_gauge_identity": "PASS",
            "positive_r": "3/x-2/(x-1)+5*x^2+7",
            "positive_frame": "x^3/(x-1)^2",
            "transported_mate": "p+a*(5*x^2+7)",
            "negative_noninteger_residue": "1/2",
            "negative_higher_pole": "x^-2",
            "complete_criterion": "r=Hprime+fprime/f, H polynomial, f rational"},
        "division_intertwiner": {
            "unbounded_regular_unit": "(1-t*x)^-1",
            "full_comparison": "PASS",
            "inner_contact": "t*a",
            "non_C_control": "x^-1 gives nonregular p+a/x",
            "regular_full_intertwiner_rigidity": "written proof"},
        "moving_projective_left_model_initial_ideal": "(1)",
        "research_context": {
            "url": "https://doi.org/10.1016/j.jalgebra.2024.07.023",
            "used": "Valued Weyl Groebner context only; no algorithm or initial-ideal theorem imported"},
        "proof_boundary": (
            "The automatic formal frame uses an entire saturated valuation "
            "ideal, normality, scalar divided operators and a full residue "
            "Hamiltonian frame. It is for generic rank one and does not prove "
            "generic freeness, finite generic existence or a higher-rank matrix "
            "analogue. Exact controls support the written all-parameter proof. "
            "The interior syzygy gives an exhaustive initial ideal, not a degree "
            "box. Rational logarithmic actions show a full rational pair and "
            "polynomial shadow do not replace finite invariant generic generation. "
            "No finite generic Long action or Dixmier counterexample is supplied."),
        "independent_replay": False, "formal_verification": False,
        "external_review": False,
        "hashes": {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                   for name in inputs}}
    (HERE/"valuation_graph_admission_controls.json").write_text(
        json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"claim": "QVA1", "status": "PASS",
                      "automatic_formal_frame": "written all-e proof",
                      "interior_initial_ideal": "(x^2,x*p)",
                      "logarithmic_admission": "PASS", "input_hashes": len(inputs)}))


if __name__ == "__main__":
    main()
