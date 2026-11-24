#!/usr/bin/env python3
"""Exact support for one primary-source coefficient correction.

This checks an unrestricted symbolic coefficient identity and explicit Weyl
operators.  It does not refute the starred Dixmier theorem or inspect the
unseen publisher proof.  The commuting witness deliberately does not satisfy
the proposition's nonzero-commutator hypothesis.
"""

from __future__ import annotations

if not __debug__:
    raise RuntimeError("This verifier requires assertions; do not use -O or -OO.")

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform

import sympy as sp


X, Y = sp.symbols("X Y")


def multiply(left, right):
    """Normal symbols X^i Y^j, with [Y,X]=1."""
    left, right = sp.expand(left), sp.expand(right)
    if left == 0 or right == 0:
        return sp.Integer(0)
    bound = min(sp.degree(left, Y), sp.degree(right, X))
    return sp.expand(sum(
        sp.diff(left, Y, k)*sp.diff(right, X, k)/sp.factorial(k)
        for k in range(int(bound) + 1)
    ))


def commutator(left, right):
    return sp.expand(multiply(left, right) - multiply(right, left))


def beta(value):
    """Coefficient-fixing antiinvolution beta(X)=X, beta(Y)=-Y."""
    value = sp.expand(value)
    if value == 0:
        return value
    reflected = value.subs(Y, -Y)
    bound = min(sp.degree(value, X), sp.degree(value, Y))
    return sp.expand(sum(
        sp.diff(reflected, X, k, Y, k)/sp.factorial(k)
        for k in range(int(bound) + 1)
    ))


def verify():
    assert commutator(Y, X) == 1
    assert beta(X) == X and beta(Y) == -Y
    assert beta(X*Y) == -X*Y - 1

    # Formal jet symbols are algebraically independent.  Thus expansion of
    # this polynomial verifies the stated coefficient identity, rather than
    # sampling particular polynomial degrees or exponent pairs.
    r, s = sp.symbols("r s")
    g, gp, gpp, h, hp, hpp = sp.symbols("g gp gpp h hp hpp")
    odd = (
        s*(s - 1)*gpp*h/2 - r*(r - 1)*hpp*g/2
        + (s - 1)*gp*(s*hp/2) - r*(s*hpp/2)*g
        + s*(r*gpp/2)*h - (r - 1)*hp*(r*gp/2)
    )
    derivative_of_top = s*gpp*h + (s - r)*gp*hp - r*g*hpp
    corrected = (r + s - 1)*derivative_of_top/2
    assert sp.expand(odd - corrected) == 0
    # The inspected display has a different prefactor and an extra term.
    inspected_display = (r + s - 1)*derivative_of_top + 2*r*g*hpp
    assert sp.expand(odd - inspected_display) != 0

    euler = X*Y + sp.Rational(1, 2)
    p0 = multiply(euler, euler)
    p1 = multiply(p0, euler)
    assert p0 == X**2*Y**2 + 2*X*Y + sp.Rational(1, 4)
    assert p1 == (X**3*Y**3 + sp.Rational(9, 2)*X**2*Y**2
                  + sp.Rational(13, 4)*X*Y + sp.Rational(1, 8))
    assert beta(euler) == -euler
    assert beta(p0) == p0 and beta(p1) == -p1
    assert commutator(p0, p1) == 0
    assert commutator(p0, p1) != sp.Rational(1, 2)
    assert sp.degree(p0, X) == 2 and sp.degree(p1, X) == 3
    assert sp.Poly(p0, X).coeff_monomial(X**2) == Y**2
    assert sp.Poly(p1, X).coeff_monomial(X**3) == Y**3
    assert sp.Poly(p0, X).coeff_monomial(X) == 2*Y
    assert sp.Poly(p1, X).coeff_monomial(X**2) == sp.Rational(9, 2)*Y**2

    gg, hh = Y**2, Y**3
    top = 3*sp.diff(gg, Y)*hh - 2*sp.diff(hh, Y)*gg
    assert top == 0 and sp.diff(top, Y) == 0
    alleged_extra = 2*2*gg*sp.diff(hh, Y, 2)
    assert alleged_extra == 24*Y**3 and alleged_extra != 0
    witness_substitution = {
        r: 2, s: 3, g: gg, gp: sp.diff(gg, Y),
        gpp: sp.diff(gg, Y, 2), h: hh, hp: sp.diff(hh, Y),
        hpp: sp.diff(hh, Y, 2),
    }
    assert sp.expand(odd.subs(witness_substitution)) == 0
    assert sp.expand(inspected_display.subs(witness_substitution)) == alleged_extra

    # An actual nonzero-commutator pair, using the source's convention.
    positive0 = -X/6 + Y**6 + 2*Y**2 + 3
    positive1 = 3*Y
    assert beta(positive0) == positive0 and beta(positive1) == -positive1
    assert commutator(positive0, positive1) == sp.Rational(1, 2)
    assert sp.degree(positive0, X) == 1 and sp.degree(positive1, X) == 0

    counts = {
        "unrestricted_symbolic_coefficient_identities": 1,
        "exact_commuting_witnesses": 1,
        "genuine_nonzero_commutator_controls": 1,
        "negative_controls": 2,
    }
    details = {
        "normal_order": "X^i Y^j with [Y,X]=1",
        "correct_odd_coefficient": "(r+s-1)/2 * (s*g'*h-r*h'*g)'",
        "witness": {
            "H": str(euler), "P0": str(p0), "P1": str(p1),
            "r": 2, "s": 3, "g": "Y^2", "h": "Y^3",
            "actual_commutator": "0",
            "actual_next_odd_coefficient": "0",
            "alleged_remaining_coefficient": "24*Y^3",
            "satisfies_nonzero_commutator_hypothesis": False,
        },
        "positive_control": {
            "P0": str(positive0), "P1": str(positive1),
            "commutator": "1/2",
            "satisfies_all_displayed_proposition_hypotheses": True,
        },
    }
    return counts, details


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    directory = Path("research/dixmier/dixmier-formal-lifting")
    stem = "starred_dixmier_primary_coefficient_audit"
    ledger_path = directory / (stem + "_source_ledger.json")
    ledger = json.loads((root / ledger_path).read_text())
    for source in ledger["sources"]:
        assert digest(root / source["path"]) == source["sha256"], source["path"]
    counts, details = verify()
    proof = directory / "STARRED_DIXMIER_PRIMARY_COEFFICIENT_AUDIT.md"
    checker = directory / Path(__file__).name
    receipt = {
        "schema": "starred_dixmier_primary_coefficient_audit_verification_v1",
        "status": "PASS",
        "arithmetic": "Exact SymPy polynomial identities and Weyl normal ordering",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "counts": counts,
        "details": details,
        "source_pins_checked": len(ledger["sources"]),
        "packet_pins": {str(path): digest(root / path)
                        for path in (proof, checker, ledger_path)},
        "scope": "Correction of one displayed intermediate coefficient "
                 "identity in arXiv:1401.5141v1, with exact witnesses. "
                 "The separate read-only foundation audit is recorded in "
                 "the ledger and is not performed by this checker.",
        "not_claimed": [
            "A counterexample to Proposition 2.6 or the starred Dixmier theorem",
            "An assertion about the uninspected publisher proof",
            "Impossibility of proving the same conclusion by another argument",
            "A new extension beyond the retained power-of-two classification",
            "Formal verification or external human review",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "counts": counts,
                      "receipt": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
