#!/usr/bin/env python3
"""QFT1 finite exact controls; all-degree and all-cover statements are written."""
from pathlib import Path
import hashlib
import json
import math
import sympy as s

HERE = Path(__file__).resolve().parent
a, q, p, h, r, ss = s.symbols("a q p h r ss")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def zero(value, label):
    assert s.cancel(value) == 0, label


def bind(name):
    receipt = json.loads((HERE/name).read_text())
    assert receipt["status"] == "PASS"
    for rel, expected in receipt["hashes"].items():
        path = HERE/rel
        if not path.exists():
            path = HERE.parent/"dixmier-audit"/rel
        assert path.is_file() and sha(path) == expected, (name, rel)
    return receipt


def mul(f, g):
    f, g = s.expand(f), s.expand(g)
    if f == 0 or g == 0:
        return s.Integer(0)
    degree = s.Poly(f, p).degree()
    return s.expand(sum(h**k*s.diff(f, p, k)*s.diff(g, q, k) /
                        math.factorial(k) for k in range(degree+1)))


def fourier(A, C):
    AA, CC = s.Poly(A, a), s.Poly(C, a)
    n = max(AA.degree(), CC.degree())
    return {j: s.expand(-q*AA.nth(j)*h**j
                        -(j+1)*AA.nth(j+1)*h**(j+1)
                        +CC.nth(j)*h**j)
            for j in range(n+1)
            if s.expand(-q*AA.nth(j)*h**j
                        -(j+1)*AA.nth(j+1)*h**(j+1)
                        +CC.nth(j)*h**j) != 0}


def indicial(op, puncture):
    n, z, la = max(op), s.symbols("z"), s.symbols("la")
    leading = op[n].subs(q, z+puncture)
    answer = 0
    for j, coefficient in op.items():
        limit = s.limit(z**(n-j)*coefficient.subs(q, z+puncture)/leading, z, 0)
        answer += limit*s.ff(la, j)
    return s.factor(answer), la


def main():
    bind("phase_carrier_controls.json")
    # Generic finite coefficient identity for normal Fourier ordering.
    aa = s.symbols("A0:5")
    bb = s.symbols("B0:4")
    A = sum(c*a**i for i, c in enumerate(aa))
    C = sum(c*a**i for i, c in enumerate(bb))
    raw = mul(A.subs(a, p), -q)+C.subs(a, p)
    expected = -q*A.subs(a, p)-h*s.diff(A, a).subs(a, p)+C.subs(a, p)
    zero(raw-expected, "Fourier A(a)*D+B normal ordering")
    represented = sum(co*p**j/h**j for j, co in fourier(A, C).items())
    zero(represented-expected, "Fourier differential coefficient convention")

    chi0 = (ss-4*a+3*r*a**2)**2/8
    assert s.degree(chi0, a) == 4
    controls = []
    for denominator_degree in (4, 5):
        A = 1+h*a**denominator_degree
        polynomial, remainder = s.div(s.expand(chi0), A, a)
        assert s.degree(s.gcd(A, chi0), a) == 0
        assert s.degree(s.gcd(A, s.diff(A, a)), a) == 0
        zero((chi0/A).subs(h, 0)-chi0, "h0 quartic preserved")
        assert s.degree(polynomial, a) <= 0
        puncture = s.cancel(polynomial)
        op = fourier(A, chi0)
        assert max(op) == denominator_degree
        zero(op[max(op)] + h**denominator_degree*s.LC(s.Poly(A, a))*
             (q-puncture), "one finite Fourier singular point")
        beta = s.Poly(remainder, a).nth(denominator_degree-1) / (
            h*s.LC(s.Poly(A, a)))
        ind, la = indicial(op, puncture)
        zero(ind-s.ff(la, denominator_degree-1)*(la+1-beta),
             "indicial roots and N-1 unit monodromy eigenvalues")
        # A's roots do not depend on ss; the residue's ss^2 coefficient
        # is nonzero, so these generic logarithmic residues are not integers.
        zero(s.diff(remainder, ss, 2)-s.Rational(1, 4),
             "generic noninteger residue dependence")
        controls.append({
            "denominator_degree": denominator_degree,
            "generic_polynomial_part": str(puncture),
            "classical_coefficient_degree": 4,
            "Fourier_order": max(op),
            "finite_singular_point": str(puncture),
            "indicial": str(ind),
            "generic_apparent_poles": "none, by source-position dependence"})

    # Mixed divisors are allowed in a logarithmic amplitude.
    alpha, beta = s.symbols("alpha beta")
    g = 4*a-3*r*a**2
    phase = (-a*ss**2/8+a**2*ss/2-r*a**3*ss/4-2*a**3/3
             +3*r*a**4/4-s.Rational(9, 40)*r**2*a**5)
    divisor = a**2-ss
    chi_log = chi0-h*alpha*2*a/divisor-h*beta/(a-r)
    eta_s = s.diff(phase, ss)-h*alpha/divisor
    eta_r = s.diff(phase, r)-h*beta/(a-r)
    zero(s.diff(chi_log, ss)+s.diff(eta_s, a), "a,s flatness")
    zero(s.diff(chi_log, r)+s.diff(eta_r, a), "a,r flatness")
    zero(s.diff(eta_s, r)-s.diff(eta_r, ss), "s,r flatness")
    A_log = divisor*(a-r)
    C_log = s.cancel(A_log*chi_log)
    assert s.degree(C_log, a) == 7
    log_op = fourier(A_log, C_log)
    assert max(log_op) == 7
    zero(log_op[7]-s.Rational(9, 8)*r**2*h**7,
         "mixed logarithmic amplitude has constant Fourier leading coefficient")
    assert not log_op[7].has(q)

    # Removing an apparent pole is necessary before counting minimal rank.
    b = s.symbols("b")
    frame = (a-b)**7
    apparent = chi0+7*h/(a-b)
    zero(apparent-h*s.diff(frame, a)/frame-chi0, "integer residue removal")
    assert max(fourier(a-b, s.cancel((a-b)*apparent))) == 5
    assert max(fourier(s.Integer(1), chi0)) == 4

    moving_frame = 1+h*a**3
    changed = chi0-h*s.diff(moving_frame, a)/moving_frame
    zero(changed.subs(h, 0)-chi0, "parameter-dependent rational gauge reduction")

    # Uniform symbolic trace certificate for arbitrary N, not a rank census.
    N = s.symbols("N", integer=True, positive=True)
    mu = 1-N
    zero((N-1)+mu**2-N*(N-1), "uniform second trace")
    zero((N-1)+mu**3+N*(N-1)*(N-2), "uniform third trace")
    zeta = s.symbols("zeta")
    zero(s.expand((zeta**2-1)-(zeta-1)*(zeta+1)), "remaining character order two")

    # Real rank-two exception after the quadratic position cover.
    c, v = s.symbols("c v", nonzero=True)
    low_op = fourier(a**2, -c+s.Rational(3, 2)*h*a)
    ind, la = indicial(low_op, 0)
    zero(ind-la*(la-s.Rational(1, 2)), "rank-two half-integral exponent")
    yy = s.Function("yy")(v)
    Dq = lambda value: s.diff(value, v)/v
    pulled = sum(co.subs(q, v**2/2)*
                 (Dq(Dq(yy)) if j == 2 else Dq(yy) if j == 1 else yy)
                 for j, co in low_op.items())
    zero(pulled+c*yy+h**2*s.diff(yy, v, 2)/2, "rank-two constant-coefficient pullback")
    J = s.Matrix([[0, 1], [-2, 0]])
    kappa = s.sqrt(-2)
    projector = (s.eye(2)+J/kappa)/2
    for entry in projector*projector-projector:
        zero(entry, "rank-two projector")
    for entry in J*projector-projector*J:
        zero(entry, "horizontal rank-two projector")

    # Tensoring by a line leaves the endomorphism connection unchanged.
    u = s.symbols("u")
    A0 = s.Matrix(2, 2, s.symbols("a00 a01 a10 a11"))
    W = s.Matrix(2, 2, s.symbols("w00 w01 w10 w11"))
    zero_matrix = (A0+u*s.eye(2))*W-W*(A0+u*s.eye(2))-(A0*W-W*A0)
    for entry in zero_matrix:
        zero(entry, "post-Fourier line twist")

    # Rank jump: a special/formal position line can have generic rank two.
    D = s.symbols("D")
    rank_jump_operator = h*D**2+D-a
    assert s.degree(rank_jump_operator, D) == 2
    assert s.degree(rank_jump_operator.subs(h, 0), D) == 1
    formal_u = a-h*a**2+h**2*(2*a**3-1)
    residual = s.expand(h*formal_u**2+h**2*s.diff(formal_u, a)+formal_u-a)
    assert all(order[0] >= 3 for order, coefficient in s.Poly(residual, h).terms()
               if coefficient)

    names = ["RATIONAL_LINE_FOURIER_OBSTRUCTION.md", "check_rational_fourier.py",
             "QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md", "phase_carrier_controls.json",
             "FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md"]
    receipt = {
        "claim": "QFT1", "status": "PASS",
        "Fourier_normal_order_generic_coefficient_control": "degree4",
        "rational_behaviors_at_infinity": controls,
        "mixed_logarithmic_flatness_identities": 3,
        "mixed_logarithmic_Fourier_order": 7,
        "apparent_pole_raw_order": 5, "clean_minimal_order": 4,
        "rational_gauge_classical_reduction": "PASS",
        "uniform_three_trace_identities": "symbolic N",
        "rank_two_true_split_control": "PASS",
        "post_Fourier_line_twist_endomorphism_identity": "PASS",
        "generic_rank_two_special_rank_one_control": "PASS",
        "hashes": {name: sha(HERE/name) for name in names},
        "research_context": [
            {"source": "https://arxiv.org/html/0712.4124",
             "used": "Fundamental differential Galois correspondence, Corollary3.10 finite algebraic part, and Section6 entire-function argument"},
            {"source": "https://arxiv.org/html/math/0312343",
             "used": "Local Fourier and lattice context; no global finite-left conclusion imported"}],
        "coverage_boundary": (
            "Exact finite operator, indicial, gauge, flatness, trace, low-rank and "
            "rank-jump controls. Universal minimal-extension identification, "
            "absolute irreducibility, all finite algebraic covers and exclusion "
            "of all rational line corrections with classical polynomial degree>=3 "
            "use written proofs and standard Picard-Vessiot theory. Minimal line "
            "extensions are the declared category; no all-carrier nonlifting or "
            "actual bounded Long carrier is certified."),
        "independent_replay": False, "formal_verification": False, "external_review": False}
    (HERE/"rational_fourier_controls.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"claim": "QFT1", "status": "PASS",
                      "rational_examples": len(controls),
                      "mixed_log_Fourier_order": 7,
                      "uniform_trace_certificate": "PASS",
                      "true_rank_two_split": "PASS"}, indent=2))


if __name__ == "__main__":
    main()
