#!/usr/bin/env python3
"""Exact controls for ELLIPTIC_EXCEPTION_ARITHMETIC.md (SymPy)."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path
import platform

import sympy as sp


def squarefree_part(value: int) -> int:
    sign = -1 if value < 0 else 1
    value = abs(value)
    result = sign
    for prime, exponent in sp.factorint(value).items():
        if exponent % 2:
            result *= int(prime)
    return result


def primitive_quartic_residues(d: int, a: int, b: int, p: int, k: int) -> dict:
    """The two projective charts contain every primitive pair modulo p^k."""
    modulus = p**k
    assert b % d == 0
    squares = {n*n % modulus for n in range(modulus)}
    pairs = [(m, 1) for m in range(modulus)]
    pairs += [(1, e) for e in range(0, modulus, p)]
    values = {(d*m**4+a*m*m*e*e+(b//d)*e**4) % modulus for m, e in pairs}
    return {"soluble": bool(values & squares), "quartic_residues": sorted(values),
            "square_residues": sorted(squares), "primitive_projective_charts": len(pairs)}


def run() -> dict:
    checks = []

    def equal(name, left, right=0):
        assert sp.cancel(sp.expand(left-right)) == 0, name
        checks.append(name)

    v, W, x, X, y, Z, m, n = sp.symbols("v W x X y Z m n")
    pi, ell = sp.Rational(128, 25), sp.Rational(9, 10)
    r, s0 = -sp.Rational(12, 5), sp.Rational(1816347, 819200)
    s = s0-v
    t = sp.Rational(10, 9)*v*v-sp.Rational(10, 9)*s0*v+sp.Rational(125, 1024)
    h = W**3+r*W*W+s*W+t
    q = W*W+ell*W+v
    A, C = sp.Rational(3, 2)/pi**2, sp.Rational(9, 16)/pi
    equal("specialized_linear_cubic_coefficient", s, A-ell*r-v)
    equal("specialized_constant_cubic_coefficient", t, (C-v*s)/ell)
    beta, delta = t+ell*s+v*r, v*t
    inverse = W**5-sp.Rational(3, 2)*W**4+A*W**3+beta*W*W+C*W+delta
    equal("complete_fixed_map_factorization", inverse, q*h)

    L = 819200*v-243483
    vL, v0 = sp.Rational(243483, 819200), sp.Rational(81, 400)
    Q = 20480000*v*v-110592000*v+143467857
    J = (400*v-81)*Q
    K = sp.Integer(41231686041600000000)
    Delta = -L*L*Q/412316860416000000
    D = sp.Rational(81, 100)-4*v
    equal("specialized_cubic_discriminant", sp.discriminant(h, W), Delta)
    equal("specialized_quadratic_discriminant", sp.discriminant(q, W), D)
    equal("matching_cover_normalization", D*Delta, L*L*J/K)
    a, b, prime0 = -10656, 27091730, 2709173
    equal("translated_integral_elliptic_model",
          ((X**3-10758382*X+6600759744).subs(X, x-3552)), x**3+a*x*x+b*x)
    equal("quadratic_discriminant_squareclass_coordinate",
          D.subs(v, sp.Rational(3, 6400)*(x+432)), -sp.Rational(3, 1600)*x)
    assert b == 10*prime0 and all(prime0 % d for d in range(2, isqrt(prime0)+1))
    checks.append("exact_prime_and_squarefree_coefficient_factorization")
    assert a < 0 and b > 0
    checks.append("real_negativity_for_x_negative")
    line_polynomial = sp.Poly(x**3+a*x*x+b*x-(m*x+n)**2, x)
    equal("chord_squareclass_product_identity", -line_polynomial.TC(), n*n)
    equal("chord_through_two_torsion_product_identity",
          sp.Poly(sp.cancel(line_polynomial.as_expr().subs(n, 0)/x), x).TC(), b)

    px, py = sp.Rational(199200, 49), sp.Rational(59459616, 343)
    shifted_px = px+3552
    equal("nontorsion_witness_on_translated_curve", py*py, shifted_px**3+a*shifted_px**2+b*shifted_px)
    equal("nontorsion_witness_squareclass_two", shifted_px, 2*sp.Rational(432, 7)**2)
    positive_divisors = sorted([1, 2, 5, 10, prime0, 2*prime0, 5*prime0, 10*prime0])
    known = {1, 2, 5*prime0, 10*prime0}
    assert {squarefree_part(d*e) for d in known for e in known} == known
    assert {squarefree_part(5*d) for d in known} == set(positive_divisors)-known
    checks.append("full_positive_squareclass_partition_by_known_subgroup")
    obstruction = primitive_quartic_residues(5, a, b, 2, 4)
    assert obstruction["quartic_residues"] == [5, 10, 15] and not obstruction["soluble"]
    checks.append("primitive_class_five_quartic_obstruction_mod_sixteen")
    field_classes = sorted(squarefree_part(-3*d) for d in known)
    assert field_classes == [-81275190, -40637595, -6, -3]
    checks.append("exact_four_quadratic_field_squareclasses")

    # Dual-isogeny descent: all signed squarefree classes, with every
    # primitive projective residue checked at both moduli.
    ad, bd = -2*a, a*a-4*b
    assert (ad, bd) == (21312, 5183416)
    assert bd == 2**3*7**3*1889
    checks.append("dual_isogeny_coefficients_and_factorization")
    dual_candidates = sorted(sign*d for sign in [-1, 1]
                             for d in [1, 2, 7, 14, 1889, 3778, 13223, 26446])
    dual16 = {d: primitive_quartic_residues(d, ad, bd, 2, 4) for d in dual_candidates}
    excluded16 = {d for d, result in dual16.items() if not result["soluble"]}
    assert excluded16 == {2, 7, 3778, 13223, -1, -14, -1889, -26446}
    checks.append("complete_dual_descent_primitive_mod_sixteen_exclusions")
    remaining = set(dual_candidates)-excluded16
    dual25 = {d: primitive_quartic_residues(d, ad, bd, 5, 2) for d in remaining}
    excluded25 = {d for d, result in dual25.items() if not result["soluble"]}
    assert excluded25 == {-2, -7, -3778, -13223}
    dual_image = remaining-excluded25
    assert dual_image == {1, 14, 1889, 26446}
    checks.append("complete_dual_descent_primitive_mod_twenty_five_exclusions")
    qx, qy = sp.Rational(3584, 25), sp.Rational(4301248, 125)
    equal("dual_descent_class_fourteen_witness_on_curve", qy*qy, qx**3+ad*qx*qx+bd*qx)
    equal("dual_descent_class_fourteen_witness_squareclass", qx, 14*sp.Rational(16, 5)**2)
    assert {squarefree_part(d*e) for d in [1, 14] for e in [1, bd]} == dual_image
    checks.append("dual_descent_all_four_survivors_realized")
    assert len(known)*len(dual_image) == 2**(2+2)
    checks.append("imported_two_isogeny_rank_formula_gives_rank_two")
    model_discriminant = 16*b*b*(a*a-4*b)
    reduction_counts = {}
    for p in [3, 11]:
        assert model_discriminant % p != 0
        reduction_counts[p] = 1+sum((yy*yy-xx**3-a*xx*xx-b*xx) % p == 0
                                    for xx in range(p) for yy in range(p))
    assert reduction_counts == {3: 4, 11: 10}
    checks.append("exact_good_reduction_point_counts_three_and_eleven")
    assert gcd(*reduction_counts.values()) == 2
    checks.append("torsion_reduction_bound_equals_known_two_torsion")

    # The displayed triple-root expansion establishes index three, rather
    # than merely inferring it from a double discriminant zero.
    z, w = sp.symbols("z w")
    expansion = w**3-w*z+sp.Rational(10, 9)*z*z-sp.Rational(42649, 16384)*z
    equal("exact_triple_root_local_expansion", h.subs({v: vL+z, W: sp.Rational(4, 5)+w}), expansion)
    assert sp.Rational(42649, 16384) != 0
    checks.append("triple_root_nonzero_transverse_derivative")
    assert sp.discriminant(Q, v) == 477703618560000 != 0
    checks.append("two_distinct_simple_quadratic_discriminant_branches")
    assert Q.subs(v, v0) == 121912785 != 0
    assert Q.subs(v, vL) == sp.Rational(3683347832025, 32768) != 0
    assert vL != v0
    checks.append("all_elliptic_and_cubic_branch_values_separated")
    equal("triple_branch_cubic_discriminant_order_two", sp.cancel(Delta/L**2), -Q/412316860416000000)
    coefficient_valuations = [-sp.Poly(sp.Poly(h, W).nth(j), v).degree() for j in range(4)]
    assert coefficient_valuations == [-2, -1, 0, 0]
    assert all(sp.Rational(2*j, 3)-2 <= value for j, value in enumerate(coefficient_valuations))
    assert gcd(2, 3) == 1
    checks.append("infinity_newton_edge_coprime_degree_three")
    assert sp.Poly(J, v).degree() == 3 and sp.discriminant(J, v) != 0
    checks.append("elliptic_cover_exactly_four_geometric_branch_values")
    ramification_over_elliptic = [3, 3, 3]
    ramification_degree = sum(e-1 for e in ramification_over_elliptic)
    assert ramification_degree == 6 and 3*(2*1-2)+ramification_degree == 2*4-2
    checks.append("cubic_root_cover_genus_four_riemann_hurwitz")
    target_b = -32*pi*pi/9*(beta+sp.Rational(5, 4)/pi)
    assert sp.Poly(target_b, v).degree() == 2
    checks.append("finite_target_map_fibres_from_quadratic_coordinate")
    assert sp.Poly(sp.resultant(q, h, W), v).degree() == 4
    checks.append("finite_deleted_quadratic_cubic_collision_parameters")

    here = Path(__file__).resolve().parent
    sources = []
    for name, statement in [
        ("dujella-number-theory-2021-section15-5.pdf",
         "pp. 507--509, eqs. 15.29--15.31: two-isogeny squareclass maps, primitive quartics, and rank formula"),
        ("dujella-number-theory-2021-section15-3.pdf",
         "Proposition 15.7, p. 490: rational torsion injects into odd good reduction"),
        ("faltings-1983-endlichkeitssaetze.pdf",
         "Satz 7, p. 365: genus at least two implies finitely many points over every fixed number field"),
    ]:
        path = here/"sources"/name
        assert path.read_bytes().startswith(b"%PDF-")
        sources.append({"path": "sources/"+name, "sha256": sha256(path.read_bytes()).hexdigest(),
                        "statement": statement, "proof_imported_not_replayed": True})
    return {
        "status": "pass", "python": platform.python_version(), "sympy": sp.__version__,
        "command": ".venv/bin/python research/arithmetic/verify_elliptic_exception_arithmetic.py",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256((here/"ELLIPTIC_EXCEPTION_ARITHMETIC.md").read_bytes()).hexdigest(),
        "exact_controls": checks,
        "quadratic_field_radicands": field_classes,
        "class_five_mod_sixteen": obstruction,
        "dual_descent": {"all_candidates": dual_candidates,
                         "excluded_mod_sixteen": sorted(excluded16),
                         "excluded_mod_twenty_five_after_first_test": sorted(excluded25),
                         "exact_image": sorted(dual_image), "rank": 2},
        "torsion": {"finite_field_counts": reduction_counts, "exact_order": 2},
        "cubic_root_cover": {"degree_over_elliptic": 3,
                             "geometric_ramification_indices": ramification_over_elliptic,
                             "ramification_degree": ramification_degree, "genus": 4},
        "external_primary_sources": sources,
        "scope": ("Exact identities and complete primitive residue controls support the written "
                  "squareclass, valuation and ramification proofs. The classical rank formula "
                  "and Faltings finiteness theorem are explicit imports. No Mordell--Weil "
                  "basis, saturation, infinite Hasse family, uniform field multiplicity "
                  "bound, or quantitative distinct-field asymptotic is asserted."),
        "formal_verification": False, "independent_whole_proof_replay": False,
        "external_review": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = run()
    output = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        Path(__file__).resolve().with_name("elliptic_exception_arithmetic_checks.json").write_text(
            output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
