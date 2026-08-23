#!/usr/bin/env sage -python
"""
Smoke-test the degree-50 Abel-Jacobi trace route to the missing third q6 section.

Exact Picard input (certified separately):
    S3 = 21 O + 22 A + C - 46 F
with
    A = -P1,
    C = reconstructed -P2,
and q6 fibre
    F6 = O + A - F.

Restricting to a generic q6 fibre (origin old O):
    AJ_q6(C) = S3 + 24 A.

This script works modulo one good prime and one good q6 base value tau.

1. Reconstruct the exact q6 chord pencil modulo p.
2. Restrict the q6 parameter to the already-explicit Hensel section C=-P2.
   Its old-base map must have degree 50.
3. For each irreducible factor of T_C(u)-tau, evaluate C on the corresponding
   finite extension and transport that point through the SAME certified
   binary-quartic -> child-Jacobian map.
4. Sum Frobenius conjugates.  This gives the standard-child group trace of the
   50 points without adjoining/symbolically adding 50 roots.
5. The two quartic signs over the pole branch u=u0 are old O and A=-P1 in
   some order.  Therefore the two orientation candidates are

       Trace(C) - 25*Q0 - 24*QA
       Trace(C) - 25*QA - 24*Q0.

No claim is made yet which orientation is S3.  The next interpolation stage
will resolve that by the global section profile / lattice checks.

This is intentionally only a smoke test.  It does not interpolate a global
section or write a generated-results certificate.

Run:
  sage -python ~/Downloads/probe_h92_q6_third_p2_trace_modp.sage
"""

import argparse
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ
)


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--tau", type=int, help="force one q6 base value instead of scanning")
parser.add_argument("--scan-limit", type=int, default=50)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1FILE = GEN / "elkies-k3-h92-p1-lift.json"
P2FILE = GEN / "elkies-k3-h92-p2-hensel-100003-p1024.json"
RRFILE = GEN / "elkies-k3-h92-q6-global-rr.json"
CHILDFILE = GEN / "elkies-k3-h92-q6-child-jacobian.json"

for path in (CORE, ANCHOR, H92, P1FILE, P2FILE, RRFILE, CHILDFILE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
squarefree_binary_quartic = scope["squarefree_binary_quartic"]
binary_quartic_invariants = scope["binary_quartic_invariants"]
transport_binary_quartic_point_to_jacobian = (
    scope["transport_binary_quartic_point_to_jacobian"]
)

p = ZZ(args.prime)
if not p.is_prime() or p in (2, 3):
    raise ValueError("prime must be good and different from 2,3")
F = GF(p)


def modq(value):
    value = QQ(value)
    if ZZ(value.denominator()) % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {value}")
    return F(ZZ(value.numerator())) / F(ZZ(value.denominator()))


p1data = json.loads(P1FILE.read_text())
p2data = json.loads(P2FILE.read_text())
rr = json.loads(RRFILE.read_text())
child = json.loads(CHILDFILE.read_text())
assert p1data["status"] == "PASS_EXACT_H92_P1"
assert p2data["schema"] == "elkies-k3.h92-p2-hensel-lift.v1"
assert p2data["complete"]
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"

R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()


def poly(values):
    return R([modq(value) for value in values])


def rf_poly(values):
    return K(poly(values))


# P1 is stored directly in entrance u=1/t coordinates.
xp = K(poly(p1data["x_entrance_base"]["numerator_coefficients"])) / K(
    poly(p1data["x_entrance_base"]["denominator_coefficients"])
)
yp = K(poly(p1data["y_entrance_base"]["numerator_coefficients"])) / K(
    poly(p1data["y_entrance_base"]["denominator_coefficients"])
)

# The Hensel artifact is in t.  Reciprocal evaluation puts the reconstructed
# coordinate (whose pinned H3 sign is -P2) into entrance u=1/t.
def reciprocal(values):
    answer = K(0)
    for index, value in enumerate(values):
        answer += K(modq(value)) / K(u**index)
    return answer


z2 = reciprocal(p2data["Z"])
xC = reciprocal(p2data["X"]) / z2**2
yC = reciprocal(p2data["Y"]) / z2**3

anchor = SourceFileLoader("q6_p2_trace_anchor", str(ANCHOR)).load_module()
r0, s0 = anchor.EXPECTED_H92
unused_ring, formulas = anchor.parse_h92(H92)
A1q, Aq, B1q, Bq, B2q = tuple(QQ(value(r0, s0)) for value in formulas)
A1, A, B1, B, B2 = map(modq, (A1q, Aq, B1q, Bq, B2q))
old_a = K(A1) / u**3 + K(A) / u**4
old_b = K(B1) / u**5 + K(B) / u**6 + K(B2) / u**7

assert yp**2 == xp**3 + old_a*xp + old_b
assert yC**2 == xC**3 + old_a*xC + old_b

h = poly(p1data["structured_denominator"]["Z4_coefficients"])


def coefficient_pair(entry):
    Ap = poly(entry["A_coefficients_low_to_high"])
    Bp = poly(entry["B_coefficients_low_to_high"])
    return K(Ap) / K(h**2), K(Bp) / K(h)


(a0, b0), (a1, b1) = tuple(
    coefficient_pair(entry) for entry in rr["kernel"]["sections"]
)

# Restriction of the q6 parameter to C=-P2.
mC = (yC - yp) / (xC - xp)
TC = (a1 + b1*mC) / (a0 + b0*mC)
TC_degree = max(TC.numerator().degree(), TC.denominator().degree())
assert TC_degree == 50, TC_degree

print(
    f"Q6P2TRACE|prime={p}|C_parameter_degree={TC_degree}|stage=setup|status=PASS",
    flush=True,
)

Tchild = child["minimal_short_weierstrass"]


def child_coeff_at(key, tau):
    coeffs = Tchild[key]
    return sum(modq(value) * tau**index for index, value in enumerate(coeffs))


def eval_poly_at(poly_value, alpha):
    parent = alpha.parent()
    answer = parent(0)
    for coefficient in reversed(poly_value.list()):
        answer = answer*alpha + parent(coefficient)
    return answer


def eval_rf_at(value, alpha):
    top = eval_poly_at(R(value.numerator()), alpha)
    bottom = eval_poly_at(R(value.denominator()), alpha)
    if not bottom:
        raise ZeroDivisionError("rational-function denominator vanished")
    return top / bottom


def to_base(value):
    if value.parent() is F:
        return F(value)
    if value**p != value:
        raise ArithmeticError("Frobenius trace coordinate did not descend")
    polynomial = value.polynomial()
    if polynomial.degree() > 0:
        raise ArithmeticError("fixed extension element was not constant")
    return F(0) if polynomial.is_zero() else F(polynomial[0])


def fourth_sixth_unit(std_a, std_b, target_a, target_b):
    if not std_a or not std_b or not target_a or not target_b:
        raise ArithmeticError("degenerate Jacobian coefficient at specialization")
    ZR = PolynomialRing(F, "z")
    z = ZR.gen()
    candidates = [
        root for root, multiplicity in (z**4 - target_a/std_a).roots()
        if multiplicity == 1
    ]
    candidates = [
        root for root in candidates
        if std_b * root**6 == target_b
    ]
    if not candidates:
        raise ArithmeticError("no fourth/sixth-power minimizing unit")
    return candidates[0]


def point_from_old_section(
    alpha, quartic, square_factor, minimizing_unit,
    q6_m, child_a, child_b, extension
):
    xpa = eval_rf_at(xp, alpha)
    xca = eval_rf_at(xC, alpha)
    yca = eval_rf_at(yC, alpha)
    mca = eval_rf_at(mC, alpha)
    mqa = eval_rf_at(q6_m, alpha)
    if mca != mqa:
        raise ArithmeticError("C point does not lie on specialized q6 fibre")
    raw_sqrt = 2*xca + xpa - mca**2
    sf = eval_rf_at(square_factor, alpha)
    if not sf:
        raise ArithmeticError("square-factor vanished at C point")
    w = raw_sqrt / sf
    quartic_ext = quartic.change_ring(extension)
    if w**2 != quartic_ext(alpha):
        # The opposite convention for the quadratic discriminant differs only
        # by sign; squaring must still agree, so this is a genuine failure.
        raise ArithmeticError("recovered quartic square root failed")
    mapped = transport_binary_quartic_point_to_jacobian(
        quartic_ext, alpha, extension(1), w, extension(minimizing_unit)
    )
    if mapped["child_a"] != extension(child_a):
        raise ArithmeticError("transported child A mismatch")
    if mapped["child_b"] != extension(child_b):
        raise ArithmeticError("transported child B mismatch")
    E = EllipticCurve(extension, [0, 0, 0, extension(child_a), extension(child_b)])
    Q = E(mapped["child_x"], mapped["child_y"])
    return E, Q


def trace_factor_point(factor, factor_index, quartic, square_factor,
                       minimizing_unit, q6_m, child_a, child_b):
    degree = factor.degree()
    if degree == 1:
        alpha = -factor[0] / factor[1]
        E, Q = point_from_old_section(
            F(alpha), quartic, square_factor, minimizing_unit,
            q6_m, child_a, child_b, F
        )
        return EllipticCurve(F, [0, 0, 0, child_a, child_b])(Q)

    name = f"a{factor_index}_{degree}"
    extension = GF(p**degree, name=name, modulus=factor.monic())
    alpha = extension.gen()
    Eext, Q = point_from_old_section(
        alpha, quartic, square_factor, minimizing_unit,
        q6_m, child_a, child_b, extension
    )
    total = Eext(0)
    current = Q
    for index in range(degree):
        total += current
        if index + 1 != degree:
            if current.is_zero():
                current = Eext(0)
            else:
                xx, yy = current.xy()
                current = Eext(xx**p, yy**p)

    Ebase = EllipticCurve(F, [0, 0, 0, child_a, child_b])
    if total.is_zero():
        return Ebase(0)
    xx, yy = total.xy()
    return Ebase(to_base(xx), to_base(yy))


def quartic_pole_points(tau, quartic, minimizing_unit, child_a, child_b):
    # At both old O and A=-P1 the chord has a pole, so T=b1/b0.
    B0 = poly(rr["kernel"]["sections"][0]["B_coefficients_low_to_high"])
    B1p = poly(rr["kernel"]["sections"][1]["B_coefficients_low_to_high"])
    line = B1p - tau*B0
    quotient, remainder = line.quo_rem(u**3)
    if remainder or quotient.degree() != 1:
        raise ArithmeticError("q6 pole branch did not reduce to one linear u-root")
    u0 = -quotient[0] / quotient[1]
    value = quartic(u0)
    roots = value.sqrt(all=True)
    if len(roots) != 2:
        raise ArithmeticError("quartic pole value did not have two signs")

    Ebase = EllipticCurve(F, [0, 0, 0, child_a, child_b])
    points = []
    for w in roots:
        mapped = transport_binary_quartic_point_to_jacobian(
            quartic, F(u0), F(1), F(w), F(minimizing_unit)
        )
        assert mapped["child_a"] == child_a
        assert mapped["child_b"] == child_b
        points.append(Ebase(mapped["child_x"], mapped["child_y"]))
    if points[0] == points[1]:
        raise ArithmeticError("the two pole signs collapsed")
    return F(u0), tuple(points)


tau_candidates = [args.tau] if args.tau is not None else range(1, args.scan_limit + 1)
last_error = None

for tau_integer in tau_candidates:
    try:
        tau = F(tau_integer)

        child_a = child_coeff_at("A_coefficients_low_to_high", tau)
        child_b = child_coeff_at("B_coefficients_low_to_high", tau)
        Ebase = EllipticCurve(F, [0, 0, 0, child_a, child_b])
        if not Ebase.discriminant():
            raise ArithmeticError("singular child specialization")

        # The degree-50 old-base fibre of C over this q6 base value.
        H = R(TC.numerator() - tau*TC.denominator())
        if H.degree() != 50:
            raise ArithmeticError(f"C fibre degree dropped to {H.degree()}")
        if H.gcd(H.derivative()).degree() != 0:
            raise ArithmeticError("C fibre polynomial is not squarefree")
        H = H.monic()

        # Reconstruct the specialized q6 quartic from the certified pencil.
        q6_m = (a1 - tau*a0) / (tau*b0 - b1)
        radicand = (
            q6_m**4
            - 6*xp*q6_m**2
            + 8*yp*q6_m
            - 3*xp**2
            - 4*old_a
        )
        quartic, square_factor = squarefree_binary_quartic(radicand, R)
        if quartic.degree() != 4:
            raise ArithmeticError(f"specialized quartic degree {quartic.degree()}")

        invariant_i, invariant_j = binary_quartic_invariants(quartic)
        std_a = -27*invariant_i
        std_b = -27*invariant_j
        minimizing_unit = fourth_sixth_unit(
            std_a, std_b, child_a, child_b
        )

        factorization = H.factor()
        factor_degrees = []
        traceC = Ebase(0)
        factor_index = 0
        for factor, multiplicity in factorization:
            if multiplicity != 1:
                raise ArithmeticError("non-squarefree factorization multiplicity")
            factor_index += 1
            factor_degrees.append(factor.degree())
            traceC += trace_factor_point(
                factor, factor_index, quartic, square_factor,
                minimizing_unit, q6_m, child_a, child_b
            )
        if sum(factor_degrees) != 50:
            raise ArithmeticError("factor degrees do not sum to 50")

        u0, pole_points = quartic_pole_points(
            tau, quartic, minimizing_unit, child_a, child_b
        )
        Qplus, Qminus = pole_points

        # One orientation is old O / A, the other is A / old O.
        candidate_1 = traceC - 25*Qplus - 24*Qminus
        candidate_2 = traceC - 25*Qminus - 24*Qplus
        if candidate_1 == candidate_2:
            raise ArithmeticError("orientation candidates unexpectedly coincide")

        def point_text(point):
            if point.is_zero():
                return "O"
            xx, yy = point.xy()
            return f"{int(xx)},{int(yy)}"

        print(
            "Q6P2TRACE_SPECIALIZATION|"
            f"prime={p}|tau={int(tau)}|u0={int(u0)}|"
            f"factor_degrees={','.join(map(str, factor_degrees))}|"
            f"traceC={point_text(traceC)}|"
            f"candidate1={point_text(candidate_1)}|"
            f"candidate2={point_text(candidate_2)}|status=PASS",
            flush=True,
        )
        print(
            "Q6P2TRACE_RESULT|"
            "degree50_trace=PASS|quartic_transport=PASS|frobenius_descent=PASS|"
            "orientation=UNRESOLVED_TWO_SIGNS|status=PASS_SMOKE_TEST",
            flush=True,
        )
        break

    except Exception as error:
        last_error = error
        if args.tau is not None:
            raise
        print(
            f"Q6P2TRACE_SKIP|tau={tau_integer}|"
            f"reason={type(error).__name__}:{error}",
            flush=True,
        )
else:
    raise RuntimeError(f"no good specialization found; last error: {last_error}")
