#!/usr/bin/env sage -python
"""Recover the rational I5/I6 supports by simple-root p-adic lifting.

The degree-22 discriminant has roots of multiplicity 5 and 6.  Its fourth
derivative has the I5 support as a simple root, and its fifth derivative has
the I6 support as a simple root.  Lifting those two roots avoids a rational
polynomial gcd on million-bit coefficients.  Candidate rational supports are
accepted only after one exact QQ division by their prescribed powers.
"""

import json
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CHECKPOINT = ROOT / "artifacts/local/elkies-k3/q24-2a5-q4o230-jacobian-checkpoint.json"
OUTPUT = ROOT / "artifacts/local/elkies-k3/q24-2a5-q4o230-repeated-supports.json"
CANDIDATES = ROOT / "artifacts/local/elkies-k3/q24-2a5-q4o230-support-candidates.json"
P = 103
MAX_EXPONENT = 1048576
started = time.monotonic()


def log(stage, **fields):
    suffix = "|".join(f"{key}={value}" for key, value in fields.items())
    print(
        f"A5A5Q4O230SUPPORTS|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{suffix}" if suffix else ""),
        flush=True,
    )


data = json.loads(CHECKPOINT.read_text())
RU = PolynomialRing(QQ, "U")
U = RU.gen()
Delta = RU([QQ(value) for value in data["Delta_coefficients_low_to_high"]])
assert Delta.degree() == 22

Fp = GF(P)
Rp = PolynomialRing(Fp, "u")
delta_p = Rp([Fp(value) for value in Delta.list()])
squarefree_p = list(delta_p.squarefree_decomposition())
profile_p = sorted((factor.degree(), int(multiplicity))
                   for factor, multiplicity in squarefree_p)
assert profile_p == [(1, 5), (1, 6), (11, 1)]
roots_mod_p = {
    int(multiplicity): int(-factor[0] / factor[1])
    for factor, multiplicity in squarefree_p
    if factor.degree() == 1
}
assert set(roots_mod_p) == {5, 6}
log("MODULAR", prime=P, profile=profile_p, root5=roots_mod_p[5], root6=roots_mod_p[6])


def coefficients_of_derivative(poly, order):
    coefficients = poly.list()
    for _ in range(order):
        coefficients = [(index + 1) * coefficients[index + 1]
                        for index in range(len(coefficients) - 1)]
    return coefficients


def evaluate_mod(coefficients, value, modulus):
    answer = ZZ(0)
    for coefficient in reversed(coefficients):
        numerator = ZZ(coefficient.numerator()) % modulus
        denominator = ZZ(coefficient.denominator()) % modulus
        answer = (answer * value + numerator * denominator.inverse_mod(modulus)) % modulus
    return answer


def lift_simple_root(multiplicity, root):
    # For a root of exact multiplicity m, Delta^(m-1) has a simple root.
    h = coefficients_of_derivative(Delta, multiplicity - 1)
    dh = [(index + 1) * h[index + 1] for index in range(len(h) - 1)]
    exponent = 1
    residue = ZZ(root)
    previous_candidate = None
    while exponent < MAX_EXPONENT:
        exponent *= 2
        modulus = ZZ(P) ** exponent
        f_value = evaluate_mod(h, residue, modulus)
        derivative_value = evaluate_mod(dh, residue, modulus)
        assert derivative_value % P
        residue = (residue - f_value * derivative_value.inverse_mod(modulus)) % modulus
        assert evaluate_mod(h, residue, modulus) == 0
        try:
            candidate = residue.rational_reconstruction(modulus)
        except ArithmeticError:
            log("LIFT", multiplicity=multiplicity, exponent=exponent,
                modulus_bits=modulus.nbits(), candidate="none")
            previous_candidate = None
            continue
        bits = max(abs(candidate.numerator()).nbits(), candidate.denominator().nbits())
        log("LIFT", multiplicity=multiplicity, exponent=exponent,
            modulus_bits=modulus.nbits(), candidate_bits=bits)
        if candidate == previous_candidate:
            return candidate, exponent
        previous_candidate = candidate
    raise RuntimeError(f"support I{multiplicity} did not reconstruct")


if CANDIDATES.exists():
    saved = json.loads(CANDIDATES.read_text())
    supports = {5: QQ(saved["I5"]), 6: QQ(saved["I6"])}
    log("CANDIDATES_RESUME", output=CANDIDATES)
else:
    supports = {}
    for multiplicity in (5, 6):
        support, exponent = lift_simple_root(multiplicity, roots_mod_p[multiplicity])
        supports[multiplicity] = support
        log("RECONSTRUCTED", multiplicity=multiplicity,
            numerator_bits=abs(support.numerator()).nbits(),
            denominator_bits=support.denominator().nbits(), exponent=exponent)
    CANDIDATES.write_text(json.dumps({
        "I5": str(supports[5]),
        "I6": str(supports[6]),
    }, sort_keys=True) + "\n")
    log("CANDIDATES", output=CANDIDATES)

assert supports[5] != supports[6]

# One exact division certifies both lower derivative vanishings at once.  The
# good-prime residual is squarefree and nonzero at both supports, proving the
# orders are exactly 5 and 6 (rather than merely at least those values).
repeated = (U - supports[5])**5 * (U - supports[6])**6
residual, remainder = Delta.quo_rem(repeated)
assert not remainder and residual.degree() == 11
repeated_p = (Rp.gen() - Fp(roots_mod_p[5]))**5 * (Rp.gen() - Fp(roots_mod_p[6]))**6
residual_p, remainder_p = delta_p.quo_rem(repeated_p)
assert not remainder_p and residual_p.degree() == 11
assert residual_p.gcd(residual_p.derivative()) == 1
assert residual_p(Fp(roots_mod_p[5])) and residual_p(Fp(roots_mod_p[6]))
assert Rp([Fp(value) for value in residual.list()]).monic() == residual_p.monic()
log("EXACT_FACTOR", repeated_degree=11, residual_degree=11,
    residual_mod_p_squarefree=True)

payload = {
    "status": "PASS_EXACT_Q24_2A5_Q4O230_REPEATED_SUPPORTS",
    "method": (
        "simple-root Hensel lifting of Delta^(m-1), rational reconstruction, "
        "one exact QQ factor division, and a good-prime residual squarefree gate"
    ),
    "prime": P,
    "modular_squarefree_profile": [
        [int(degree), int(multiplicity)] for degree, multiplicity in profile_p
    ],
    "supports": {
        f"I{multiplicity}": {
            "U": str(supports[multiplicity]),
            "numerator_bits": int(abs(supports[multiplicity].numerator()).nbits()),
            "denominator_bits": int(supports[multiplicity].denominator().nbits()),
            "exact_discriminant_order": multiplicity,
        }
        for multiplicity in (5, 6)
    },
    "residual_nodal_factor_coefficients_low_to_high": [
        str(value) for value in residual.monic().list()
    ],
    "large_rational_polynomial_gcd_required": False,
    "elapsed_seconds": round(time.monotonic() - started, 6),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
log("DONE", status=payload["status"], output=OUTPUT)
