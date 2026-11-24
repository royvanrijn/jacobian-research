#!/usr/bin/env python3
"""Exact finite shifted-resultant screens; no unrestricted positivity claim."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit("FAIL: optimized Python execution is not permitted")

HERE = Path(__file__).resolve().parent
NOTE = HERE / "SFC_SHIFTED_RESULTANT_POSITIVITY_AUDIT.md"
NOTE_PIN = 'c4360e066478e6651b8dbd684f8fd67865a36d698508ecc0648201b58bc44195'
SYMPY_VERSION = "1.14.0"
SOURCE_PINS = {
    "SFC_DEAMPLIFICATION_BRIDGES.md": "df9f19e4418c2abdf47c0f2696a60637dcfb761f93a9b91457aad4b5c671e6b0",
    "SFC_INITIAL_RESULTANT_AND_SHIFTED_WINDOWS.md": "9df960f24475ad06f6447739dc30c8d25f9a00b629f10b9f32c88a00b2631751",
    "UNRESTRICTED_SFC_MECHANISMS.tex": "f8b3dfc5947601fb002aee0ba942b0fc20b7ef46ca8e0d4f69b5a3e5e1361785"
}
EXPECTED = {
    "1": {
        "degree": 3,
        "residual_degree": 0,
        "residual_sha256": "47d62cf50f4414a1b6669e711cbaf7bb5a21a6196ba5fa67bbd8e2b43a977b48",
        "residual_terms": 1,
        "resultant_sha256": "a32b5281e96eebc51ec7418f919ee56dd7f9a9298ce8d1530dca02520be579a7",
        "schur_terms": 1,
        "terms": 8
    },
    "2": {
        "degree": 9,
        "residual_degree": 5,
        "residual_sha256": "04f8d0fe7eddf6c4b60b34c6801d769b806f1a0914467cc1a6c169041da4950a",
        "residual_terms": 15,
        "resultant_sha256": "07cc04849e2d2a951f0fb55345fac8aac43a63affb74bfc18fb2cbdc3ace6f6a",
        "schur_terms": 9,
        "terms": 43
    },
    "3": {
        "degree": 18,
        "residual_degree": 13,
        "residual_sha256": "610e90cad112531274b2ef638cca4341c01917ed4d14fa56722af42f52693141",
        "residual_terms": 75,
        "resultant_sha256": "f406e2321ab6c814e1d700cfd5a6090cd8ba35185b4eecb2dfba833aa5240254",
        "schur_terms": 41,
        "terms": 148
    },
    "4": {
        "degree": 30,
        "residual_degree": 24,
        "residual_sha256": "287a813a3c6a3812849e1c24ea668fc1673fd96374f21f306e21c6f2962d6cb3",
        "residual_terms": 235,
        "resultant_sha256": "bcca81e18e275ebdbf8df24c8363200da6e2c456c0ee6919925722c679f9912b",
        "schur_terms": 124,
        "terms": 386
    },
    "5": {
        "degree": 45,
        "residual_degree": 38,
        "residual_sha256": "0490ac6e8e3f26aced35d6c47a7cf7d62b4aa7d5886e14676f79d306bf9010b0",
        "residual_terms": 570,
        "resultant_sha256": "d68c67218fb76566e9d037ca4d203cdd25686a768a57b9dbd9f20dbb7b776e2b",
        "schur_terms": 295,
        "terms": 841
    }
}


class VerificationError(RuntimeError):
    pass


def need(condition, message):
    if not bool(condition):
        raise VerificationError(message)


def poly_digest(poly):
    data = [[int(i), int(j), int(c.p), int(c.q)]
            for (i, j), c in poly.terms()]
    return sha256(json.dumps(data, separators=(",", ":")).encode()).hexdigest()


def bareiss(matrix):
    """Independent integer determinant, including exact-division guards."""
    a = [list(row) for row in matrix]
    n = len(a)
    need(all(len(row) == n for row in a), "nonsquare determinant input")
    if not n:
        return 1
    previous = 1
    sign = 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k]), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = pivot * a[i][j] - a[i][k] * a[k][j]
                need(value % previous == 0, "inexact Bareiss division")
                a[i][j] = value // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def integer_form(m, first, second):
    """Independent binomial construction, coefficients in ascending order."""
    return [comb(m * first + k, k)
            * comb(m * second + m - k, m - k)
            for k in range(m + 1)]


def sylvester(first, second):
    f, g = first[::-1], second[::-1]
    m, n = len(f) - 1, len(g) - 1
    need(f[0] and g[0], "unexpected degree drop in numeric control")
    rows = []
    for shift in range(n):
        rows.append([0] * shift + f + [0] * (n - 1 - shift))
    for shift in range(m):
        rows.append([0] * shift + g + [0] * (m - 1 - shift))
    return bareiss(rows)


def check_resultant(poly, n, r, u):
    expected = EXPECTED[str(n)]
    need(poly.total_degree() == expected["degree"], f"n={n}: degree mismatch")
    need(len(poly.terms()) == expected["terms"], f"n={n}: term-count mismatch")
    need(all(bool(c > 0) for c in poly.coeffs()), f"n={n}: nonpositive coefficient")
    need(poly.eval({r: 0, u: 0}) == 1, f"n={n}: orientation/constant mismatch")
    need(poly_digest(poly) == expected["resultant_sha256"],
         f"n={n}: full coefficient digest mismatch")


def must_reject(function, message):
    try:
        function()
    except VerificationError:
        return
    raise VerificationError(message)


def symbolic_screens(S):
    r, u, t = S.symbols("r u t")

    def form(m):
        return S.Poly(sum(
            S.rf(m * r + 1, k) * S.rf(m * u + 1, m - k) * t ** k
            / (S.factorial(k) * S.factorial(m - k))
            for k in range(m + 1)), t, r, u).as_expr()

    records = []
    results, residuals = {}, {}
    for n in range(1, 6):
        resultant = S.Poly(S.resultant(form(n), form(n + 1), t), r, u)
        check_resultant(resultant, n, r, u)
        factor = S.Poly((r + 1) * (u + 1) * (r + u + 1) ** n, r, u)
        residual, remainder = resultant.div(factor)
        need(remainder.is_zero, f"n={n}: common factor did not divide exactly")
        expected = EXPECTED[str(n)]
        need(residual.total_degree() == expected["residual_degree"],
             f"n={n}: residual degree mismatch")
        need(len(residual.terms()) == expected["residual_terms"],
             f"n={n}: residual term-count mismatch")
        need(all(bool(c > 0) for c in residual.coeffs()),
             f"n={n}: nonpositive residual coefficient")
        need(poly_digest(residual) == expected["residual_sha256"],
             f"n={n}: full residual digest mismatch")

        reconstructed = {}
        schur_count = 0
        for degree in range(residual.total_degree() + 1):
            previous = S.Integer(0)
            for k in range(degree // 2 + 1):
                current = residual.coeff_monomial(r ** (degree - k) * u ** k)
                coefficient = current - previous
                need(coefficient >= 0, f"n={n}: negative Schur coefficient")
                if coefficient:
                    schur_count += 1
                    for j in range(k, degree - k + 1):
                        monomial = (degree - j, j)
                        reconstructed[monomial] = reconstructed.get(monomial, 0) + coefficient
                previous = current
        need(S.Poly.from_dict(reconstructed, (r, u)) == residual,
             f"n={n}: Schur reconstruction mismatch")
        need(schur_count == expected["schur_terms"], f"n={n}: Schur term count")

        controls = []
        for first, second in ((0, 0), (1, 2), (2, 1)):
            direct = sylvester(integer_form(n, first, second),
                               integer_form(n + 1, first, second))
            symbolic = resultant.eval({r: first, u: second})
            need(symbolic == direct, f"n={n}: independent determinant mismatch")
            controls.append({"exponents": [first, second], "resultant": str(direct)})

        records.append({"n": n, **expected, "negative_ordinary_coefficients": 0,
                        "negative_schur_coefficients": 0,
                        "independent_sylvester_controls": controls})
        results[n], residuals[n] = resultant, residual
        print(f"PASS binary n={n}: {len(resultant.terms())} positive terms; "
              f"{schur_count} positive Schur terms", flush=True)

    need(results[1].as_expr() == S.expand((r + 1) * (u + 1) * (r + u + 1)),
         "initial-resultant calibration mismatch")
    need(S.Poly(form(3).subs({r: -S.Rational(1, 3), u: -S.Rational(1, 3)}), t).is_zero,
         "negative-exponent shifted control: g3 is not zero")
    need(results[2].eval({r: -S.Rational(1, 3), u: -S.Rational(1, 3)}) == 0,
         "negative-exponent shifted control: resultant is nonzero")
    need(results[1].eval({r: -S.Rational(1, 3), u: -S.Rational(1, 3)}) == S.Rational(4, 27),
         "negative-exponent initial control mismatch")

    q = 4 * residuals[2].as_expr()
    restricted = S.Poly(q.subs(u, 0), r)
    need(restricted.all_coeffs() == [17, 37, 22, 4], "cubic restriction mismatch")
    a, b, c, d = map(int, restricted.all_coeffs())
    discriminant = b*b*c*c - 4*a*c*c*c - 4*b*b*b*d - 27*a*a*d*d + 18*a*b*c*d
    need(discriminant == -428 and S.discriminant(restricted.as_expr(), r) == -428,
         "nonreal-root negative control failed")
    need(S.gcd_list(S.Poly(q, u).all_coeffs()) == 1, "vertical-factor gcd mismatch")
    need(S.gcd_list(S.Poly(q, r).all_coeffs()) == 1, "horizontal-factor gcd mismatch")
    top = sum(coefficient * r ** i * u ** j
              for (i, j), coefficient in S.Poly(q, r, u).terms() if i + j == 5)
    need(S.expand(top - 81*r*r*u*u*(r+u)) == 0, "top homogeneous part mismatch")
    need(S.expand(q.subs(u, r) - (2*r+1)*(3*r+1)**2*(3*r+2)**2) == 0,
         "diagonal calibration mismatch")
    must_reject(lambda: check_resultant(results[2] + S.Poly(r*u, r, u), 2, r, u),
                "altered coefficient was accepted")
    must_reject(lambda: check_resultant(-results[2], 2, r, u),
                "wrong resultant orientation was accepted")
    return records, {
        "cubic_coefficients_descending": [17, 37, 22, 4],
        "cubic_discriminant": discriminant,
        "vertical_and_horizontal_coefficient_gcds": [1, 1],
        "negative_exponent_initial_resultant": "4/27",
        "negative_exponent_shifted_resultant": "0",
        "altered_coefficient_rejected": True,
        "wrong_orientation_rejected": True,
    }


def shape_contiguity_control(S):
    x, y = S.symbols("x y")
    p2 = x*x + x*y + y*y
    p3 = x**3 + x*x*y + x*y*y + y**3
    transform = lambda p: S.expand(p + x*S.diff(p, x))
    q2, q3 = transform(p2), transform(p3)
    need(q2 == 3*x*x + 2*x*y + y*y, "shape-two quadratic mismatch")
    need(q3 == 4*x**3 + 3*x*x*y + 2*x*y*y + y**3, "shape-two cubic mismatch")
    original = int(S.resultant(p2.subs(y, 1), p3.subs(y, 1), x))
    shifted = int(S.resultant(q2.subs(y, 1), q3.subs(y, 1), x))
    need((original, shifted) == (1, 16), "Gamma-contiguity resultants mismatch")
    inverse = S.Poly.from_dict({monomial: coefficient/(monomial[0]+1)
                               for monomial, coefficient in S.Poly(x*q2, x, y).terms()}, (x, y))
    witness = lambda p: S.Poly(p, x, y).coeff_monomial(x*x*y) - S.Poly(p, x, y).coeff_monomial(x*y*y)
    need(all(witness(p) == 0 for p in (x*p2, y*p2, p3)),
         "claimed nonmembership functional does not annihilate the ideal component")
    need(witness(inverse.as_expr()) == S.Rational(1, 6), "nonmembership witness mismatch")
    return {"base_resultant": original, "shape_shifted_resultant": shifted,
            "inverse_transport_witness": "1/6", "both_endpoint_pairs_regular": True,
            "regularity_counterexample_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--note", type=Path, default=NOTE)
    parser.add_argument("--source-dir", type=Path, default=HERE)
    parser.add_argument("--output", type=Path,
                        default=HERE / "sfc_shifted_resultant_positivity_verification.json")
    args = parser.parse_args()
    need(args.note.is_file(), "missing pinned note")
    need(sha256(args.note.read_bytes()).hexdigest() == NOTE_PIN, "note SHA-256 mismatch")
    need(set(EXPECTED) == {str(n) for n in range(1, 6)}, "incomplete expected screen data")
    need(len(SOURCE_PINS) == 3, "incomplete source pins")
    for name, digest in SOURCE_PINS.items():
        path = args.source_dir / name
        need(path.is_file(), f"missing pinned source: {name}")
        need(sha256(path.read_bytes()).hexdigest() == digest, f"source SHA-256 mismatch: {name}")
    import sympy as S
    need(S.__version__ == SYMPY_VERSION, "SymPy version mismatch")
    screens, negative_controls = symbolic_screens(S)
    shape_control = shape_contiguity_control(S)
    receipt = {
        "status": "PASS",
        "note_sha256": NOTE_PIN,
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": SOURCE_PINS,
        "sympy_version": S.__version__,
        "scope": "Complete exact binary symbolic screens for n=1,...,5; independent fixed-parameter determinants and exact obstruction controls.",
        "claims_not_made": ["all-order resultant positivity", "all-support resultant positivity",
                             "unrestricted SFC", "formal verification", "external review"],
        "binary_screens": screens,
        "independent_numeric_sylvester_determinants": 15,
        "negative_controls": negative_controls,
        "positive_shape_contiguity_control": shape_control,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    temporary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print("PASS: five complete symbolic binary screens, 15 independent Sylvester determinants, "
          "full Schur reconstruction and exact negative controls.")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"FAIL: {error}", file=sys.stderr)
        raise SystemExit(1)
