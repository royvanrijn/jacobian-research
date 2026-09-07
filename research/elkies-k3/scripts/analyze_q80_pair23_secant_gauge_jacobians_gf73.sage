#!/usr/bin/env sage
"""Classify unit-corrected Jacobians of pair23 secant gauges over GF(73)."""

import hashlib
import json
from pathlib import Path

from sage.all import GF, FunctionField, PolynomialRing, gcd


ROOT = Path(__file__).resolve().parents[2]
SOURCE = (
    ROOT / "artifacts/local/"
    "q80-pair23-secant-reference-gauge-search-gf73.json"
)
source_bytes = SOURCE.read_bytes()
source_hash = hashlib.sha256(source_bytes).hexdigest()
source = json.loads(source_bytes)
assert source["schema"] == "q80-pair23-secant-reference-gauge-search-gf73-v1"
hits = source["genus_one_hits"]
assert len(hits) == 6

finite = GF(73, impl="modn")
base = FunctionField(finite, "s")
s = base.gen()
cover_ring = PolynomialRing(base, "tau")
tau = cover_ring.gen()


def kodaira_data(ord_a, ord_b, ord_delta):
    if ord_a == 0 or ord_b == 0:
        n = int(ord_delta)
        return (f"I{n}", n-1, n*(n-1), n, n)
    if ord_delta == 2:
        return ("II", 0, 0, 1, 2)
    if ord_delta == 3:
        return ("III", 1, 2, 2, 3)
    if ord_delta == 4:
        return ("IV", 2, 6, 3, 4)
    if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
        return ("I0*", 4, 24, 4, 6)
    if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
        n = int(ord_delta-6)
        rank = n+4
        return (f"I{n}*", rank, 2*rank*(rank-1), 4, n+6)
    if ord_delta == 8:
        return ("IV*", 6, 72, 3, 8)
    if ord_delta == 9:
        return ("III*", 7, 126, 2, 9)
    if ord_delta == 10:
        return ("II*", 8, 240, 1, 10)
    raise ArithmeticError((ord_a, ord_b, ord_delta))


rows = []
for index, hit in enumerate(hits):
    cover = cover_ring(sum(
        base(coefficient)*tau**tau_degree*s**s_degree
        for tau_degree, s_degree, coefficient
        in hit["integral_double_cover_terms_T_U_coefficient"]
    ))
    factorization_object = cover.factor()
    factorization = tuple(factorization_object)
    twist = base(factorization_object.unit())
    odd_part = cover_ring(1)
    for factor, exponent in factorization:
        if int(exponent) % 2:
            odd_part *= factor
    quartic = odd_part.monic()
    assert quartic.degree() in (3, 4)
    coefficients = list(quartic.list())+[base(0)]*5
    e, d, c, b, a = coefficients[:5]
    invariant_I = 12*a*e-3*b*d+c**2
    invariant_J = (
        72*a*c*e+9*b*c*d-27*a*d**2-27*b**2*e-2*c**3
    )
    jacobian_A = twist**2*(-27*invariant_I)
    jacobian_B = twist**3*(-27*invariant_J)
    delta_core = twist**6*(4*invariant_I**3-invariant_J**2)
    assert jacobian_A.denominator() == 1
    assert jacobian_B.denominator() == 1
    assert delta_core.denominator() == 1
    A = jacobian_A.numerator()
    B = jacobian_B.numerator()
    Delta = delta_core.numerator()

    finite_scalings = []
    for factor, _ in gcd(A, B).factor():
        scale_order = min(A.valuation(factor)//4, B.valuation(factor)//6)
        if scale_order <= 0:
            continue
        A //= factor**(4*scale_order)
        B //= factor**(6*scale_order)
        Delta //= factor**(12*scale_order)
        finite_scalings.append((str(factor), int(scale_order)))
    assert A.degree() <= 8 and B.degree() <= 12 and Delta.degree() <= 24

    finite_signature = []
    root_rank = 0
    root_count = 0
    root_determinant = 1
    euler_number = 0
    for factor, exponent in Delta.factor():
        data = kodaira_data(
            int(A.valuation(factor)), int(B.valuation(factor)), int(exponent)
        )
        kind, rank, count, determinant, euler = data
        degree = int(factor.degree())
        finite_signature.append((
            str(factor), degree, int(A.valuation(factor)),
            int(B.valuation(factor)), int(exponent), kind,
        ))
        root_rank += degree*rank
        root_count += degree*count
        root_determinant *= determinant**degree
        euler_number += degree*euler
    infinity_orders = (8-A.degree(), 12-B.degree(), 24-Delta.degree())
    infinity_kind, rank, count, determinant, euler = kodaira_data(
        *infinity_orders
    )
    root_rank += rank
    root_count += count
    root_determinant *= determinant
    euler_number += euler
    assert euler_number == 24

    record = {
        "index": index,
        "left_keys": hit["left_keys"],
        "right_keys": hit["right_keys"],
        "twist": str(twist),
        "finite_scalings": [list(row) for row in finite_scalings],
        "degrees_A_B_Delta": [A.degree(), B.degree(), Delta.degree()],
        "A_coefficients_low_to_high": list(map(int, A.list())),
        "B_coefficients_low_to_high": list(map(int, B.list())),
        "Delta_coefficients_low_to_high": list(map(int, Delta.list())),
        "finite_signature": [list(row) for row in finite_signature],
        "infinity_orders": list(map(int, infinity_orders)),
        "infinity_fiber": infinity_kind,
        "root_data": [root_rank, root_count, root_determinant],
        "geometric_CM24_MW_rank": 18-root_rank,
        "euler_number": euler_number,
        "matches_deforming_target": (
            (root_rank, root_count, root_determinant) == (16, 66, 2048)
        ),
    }
    rows.append(record)
    print(
        "Q80PAIR23SECANTJAC|"
        f"index={index}|left={tuple(hit['left_keys'])}|right={tuple(hit['right_keys'])}|"
        f"scalings={tuple(finite_scalings)}|degrees={(A.degree(), B.degree(), Delta.degree())}|"
        f"finite_signature={tuple(finite_signature)}|infinity={infinity_kind}|"
        f"root_data={(root_rank, root_count, root_determinant)}|"
        f"MW={18-root_rank}|target={int(record['matches_deforming_target'])}|status=PASS",
        flush=True,
    )

targets = [row for row in rows if row["matches_deforming_target"]]
artifact = {
    "schema": "q80-pair23-secant-gauge-jacobian-audit-gf73-v1",
    "status": "exact_unit_corrected_k3_jacobian_classification",
    "prime": 73,
    "source_artifact": str(SOURCE.relative_to(ROOT)),
    "source_sha256": source_hash,
    "rows": rows,
    "target_hits": targets,
    "rank_claim": None,
    "reproduce": (
        "sage elkies-k3/scripts/"
        "analyze_q80_pair23_secant_gauge_jacobians_gf73.sage"
    ),
}
output = (
    ROOT / "artifacts/local/"
    "q80-pair23-secant-gauge-jacobian-audit-gf73.json"
)
output.write_text(json.dumps(artifact, indent=2, sort_keys=True, default=int)+"\n")
print(
    "Q80PAIR23SECANTJAC|"
    f"classified={len(rows)}|target_hits={len(targets)}|artifact={output}|"
    "status=PASS_EXACT_CLASSIFICATION",
    flush=True,
)
