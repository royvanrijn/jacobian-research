#!/usr/bin/env sage -python
"""Bounded, fail-closed fixed-resolvent quartic-orbit experiment for E302.

The input cubic displayed by the descent equation defines a very nonmaximal
order.  This program first recomputes the maximal order and its oriented
binary-cubic form.  Only that form is used in 4*det(x*A-y*B).

The search is deliberately a small normalized affine slice, not a claim that
the slice is a fundamental domain.  In particular, rational GL(3) changes of
one seed are not reported as new quartic objects: they preserve its rational
quartic algebra.  A field can advance only after all of the arithmetic gates
in CURVE302_FIXED_RESOLVENT_ORBIT_PROTOCOL.json pass.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import itertools
import json
from math import gcd
from pathlib import Path
import subprocess
import sys

from sage.all import AA, GF, QQ, ZZ, EllipticCurve, PolynomialRing, matrix, pari, prime_range, vector
from sage.version import version as SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
PROTOCOL = Path(__file__).with_name("CURVE302_FIXED_RESOLVENT_ORBIT_PROTOCOL.json")
MAXIMUM_FORM = ART / "rank_jump_curve302_maximal_norm_form_v1.json"
STRICT = ART / "rank_jump_curve302_strict_constructor_arithmetic_v1.json"
FILTRATION = ART / "curve302_recovered_quotient_local_filtration_v1.json"
CURVE = CAS / "icarm_curve302.py"
HIDDEN_SOLVER = CAS / "benchmark_curve302_strict_covers_hidden.sage"
OUTPUT = ART / "curve302_fixed_resolvent_orbits_v1.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-fixed-resolvent-orbits-v1"

sys.path.insert(0, str(CAS))
import icarm_curve302 as curve  # noqa: E402


def read(path: Path):
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def put_new(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def string_matrix(value):
    return [[str(value[i, j]) for j in range(value.ncols())] for i in range(value.nrows())]


def cubic_discriminant(coefficients):
    a, b, c, d = map(ZZ, coefficients)
    return b * b * c * c - 4 * a * c**3 - 4 * b**3 * d - 27 * a * a * d * d + 18 * a * b * c * d


def maximal_cubic_preflight():
    """Independently recover the normal maximal-order binary cubic form."""
    maximum, strict = read(MAXIMUM_FORM), read(STRICT)
    assert maximum["status"] == strict["status"] == "PASS"
    record = next(case for case in maximum["cases"] if case["source"].endswith("rank_jump_curve302_strict_constructor_arithmetic_v1.json"))
    R = PolynomialRing(QQ, "z")
    raw = R([QQ(value) for value in strict["cubic_ascending"]])
    pari.allocatemem(64_000_000, read(PROTOCOL)["bounds"]["pari_stack_bytes"], silent=True)
    nf = pari.nfinit([pari(raw), strict["S_finite"]])
    assert str(nf.disc()) == strict["field_discriminant"]
    assert ZZ(nf[3]) > 1
    zk = list(nf.nf_get_zk())
    assert zk[0] == 1
    pair = list(map(ZZ, pari.nfalgtobasis(nf, zk[1] * zk[2])))
    w = zk[1] - pair[2]
    t = zk[2] - pair[1]

    def coords(value):
        row = list(map(ZZ, pari.nfalgtobasis(nf, value)))
        return [row[0] + row[1] * pair[2] + row[2] * pair[1], row[1], row[2]]

    table = [[coords(left * right) for right in (1, w, t)] for left in (1, w, t)]
    a, b, c, d = -table[1][1][2], table[1][1][1], -table[2][2][2], table[2][2][1]
    form = [a, b, c, d]
    assert table[1][1] == [-a * c, b, -a]
    assert table[1][2] == [-a * d, 0, 0]
    assert table[2][2] == [-b * d, d, -c]
    assert [str(value) for value in form] == record["binary_cubic_descending"]
    assert cubic_discriminant(form) == ZZ(nf.disc())
    assert str(ZZ(nf[3])) == record["defining_order_index"]
    return {
        "status": "PASS",
        "raw_defining_cubic_ascending": strict["cubic_ascending"],
        "raw_order_index_in_maximal_order": str(nf[3]),
        "field_discriminant": str(nf.disc()),
        "maximal_order_basis_GP": [str(value) for value in zk],
        "normal_basis_GP": [str(value) for value in (1, w, t)],
        "normal_basis_multiplication_table": [[[str(x) for x in row] for row in block] for block in table],
        "binary_cubic_descending": [str(value) for value in form],
        "binary_cubic_discriminant": str(cubic_discriminant(form)),
        "artifact_crosscheck": relative(MAXIMUM_FORM),
        "assertion": "The determinant fibre below is attached to this maximal cubic ring, never to the nonmaximal displayed cubic polynomial.",
    }, raw, nf


def matrix_pair_for_affine_slice(form, px, py):
    """Return A,B whose fixed invariant is the supplied binary cubic form."""
    a, b, c, d = map(ZZ, form)
    qx, qy, nx, ny = -a, 2 * px - b, 2 * py - c, -d
    A = matrix(QQ, [[0, QQ(1) / 2, 0], [QQ(1) / 2, nx, px], [0, px, qx]])
    B = matrix(QQ, [[0, 0, -QQ(1) / 2], [0, -ny, -py], [-QQ(1) / 2, -py, -qy]])
    return A, B


def pair_coefficients(A, B):
    """Integral quadratic-form coefficients in the convention stated in protocol."""
    order = ((0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2))
    return [ZZ(A[i, j] if i == j else 2 * A[i, j]) for i, j in order] + [ZZ(B[i, j] if i == j else 2 * B[i, j]) for i, j in order]


def determinant_coefficients(A, B):
    R = PolynomialRing(QQ, names=("x", "y"))
    x, y = R.gens()
    D = 4 * (x * A - y * B).det()
    return [ZZ(D.monomial_coefficient(x ** (3 - i) * y**i)) for i in range(4)]


def primitive(values):
    content = 0
    for value in values:
        content = gcd(content, abs(int(value)))
    return content == 1


def permutation_sign(values):
    assert sorted(values) == [1, 2, 3]
    return -1 if sum(values[i] > values[j] for i in range(3) for j in range(i + 1, 3)) % 2 else 1


def quartic_multiplication(coefficients):
    """Direct HCL-III lambda/c invariant reconstruction; no external ring code."""
    assert len(coefficients) == 12
    aa, bb = [QQ(value) for value in coefficients[:6]], [QQ(value) for value in coefficients[6:]]
    pos = {(1, 1): 0, (2, 2): 1, (3, 3): 2, (1, 2): 3, (2, 1): 3,
           (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5}

    def lam(i, j, k, ell):
        return aa[pos[i, j]] * bb[pos[k, ell]] - bb[pos[i, j]] * aa[pos[k, ell]]

    def ci(i):
        return {1: lam(2, 3, 1, 1), 2: -lam(1, 3, 2, 2), 3: lam(1, 2, 3, 3)}[i]

    memo = {}

    def cijk(i, j, k):
        key = (i, j, k)
        if key in memo:
            return memo[key]
        original = key
        if k == 0:
            if i > j:
                i, j = j, i
            rem = [1, 2, 3]
            rem.remove(j)
            if i != j:
                rem.remove(i)
            k = rem[0]
            answer = sum(cijk(j, k, r) * cijk(r, i, k) - cijk(i, j, r) * cijk(r, k, k) for r in range(1, 4))
        else:
            if j == k and i != j:
                j, i = i, k
            if i == j:
                if j == k:
                    rem = [1, 2, 3]
                    rem.remove(i)
                    j, k = rem
                    answer = permutation_sign([i, j, k]) * lam(i, k, i, j) + ci(i)
                else:
                    j = k
                    rem = [1, 2, 3]
                    rem.remove(i)
                    rem.remove(j)
                    k = rem[0]
                    answer = permutation_sign([i, j, k]) * lam(i, i, i, k)
            elif i == k:
                rem = [1, 2, 3]
                rem.remove(i)
                rem.remove(j)
                k = rem[0]
                answer = permutation_sign([i, j, k]) * lam(i, k, j, j) / 2 + ci(j) / 2
            else:
                answer = permutation_sign([i, j, k]) * lam(j, j, i, i)
        memo[original] = answer
        return answer

    table = []
    for i in range(4):
        row = []
        for j in range(4):
            if i == 0:
                row.append(vector(QQ, [ZZ(k == j) for k in range(4)]))
            elif j == 0:
                row.append(vector(QQ, [ZZ(k == i) for k in range(4)]))
            else:
                row.append(vector(QQ, [cijk(i, j, k) for k in range(4)]))
        table.append(row)
    assert all(value in ZZ for row in table for product in row for value in product)
    table = [[vector(ZZ, product) for product in row] for row in table]
    for i, j, k in itertools.product(range(4), repeat=3):
        left = sum((table[r][k] * table[i][j][r] for r in range(4)), vector(ZZ, 4))
        right = sum((table[i][r] * table[j][k][r] for r in range(4)), vector(ZZ, 4))
        assert left == right
    left_multiplication = []
    for i in range(4):
        left_multiplication.append(matrix(ZZ, 4, 4, lambda row, col: table[i][col][row]))
    return table, left_multiplication


def short_vectors(bound):
    """One sign representative of primitive nonzero vectors in Z^3."""
    vectors = []
    for row in itertools.product(range(-bound, bound + 1), repeat=3):
        if not any(row) or gcd(gcd(abs(row[0]), abs(row[1])), abs(row[2])) != 1:
            continue
        first = next(value for value in row if value)
        if first < 0:
            continue
        vectors.append(row)
    return sorted(vectors, key=lambda row: (sum(abs(value) for value in row), row))


def element_matrix(row, multiplication):
    return sum((ZZ(row[i]) * multiplication[i + 1] for i in range(3)), matrix.zero(ZZ, 4))


def algebra_probe(multiplication, bound):
    """Either exhibit an exact zero divisor or an exact primitive field generator."""
    zero_divisor = None
    for row in short_vectors(bound):
        M = element_matrix(row, multiplication)
        if M.det() == 0:
            zero_divisor = row
            break
        q = M.charpoly()
        if q.is_irreducible():
            return {"status": "FIELD_GENERATOR", "coordinates": list(row), "polynomial_ascending": [str(value) for value in q.list()]}
    if zero_divisor is not None:
        return {"status": "ZERO_DIVISOR", "coordinates": list(zero_divisor), "norm": "0"}
    return {"status": "UNKNOWN_NO_SHORT_WITNESS", "bound": bound}


def factor_pattern(poly, p):
    if poly.discriminant() % p == 0:
        return None
    reduced = poly.change_ring(GF(p))
    return sorted(int(factor.degree()) for factor, multiplicity in reduced.factor() for _ in range(multiplicity))


def structural_field_certificate(poly, field_discriminant, prime_bound):
    """Field/maximality/S4/real gates, all exact once a field generator is found."""
    if not poly.is_irreducible():
        return {"status": "REJECTED_NOT_A_FIELD"}
    nf = pari.nfinit(pari(poly))
    if ZZ(nf.disc()) != ZZ(field_discriminant):
        return {"status": "REJECTED_NONMAXIMAL_QUARTIC_ORDER", "field_discriminant": str(nf.disc()), "expected": str(field_discriminant)}
    witnesses = {}
    for p in prime_range(2, prime_bound + 1):
        pattern = factor_pattern(poly, p)
        if pattern is None:
            continue
        if pattern == [4] and "four_cycle" not in witnesses:
            witnesses["four_cycle"] = {"prime": int(p), "factor_degrees": pattern}
        if pattern == [1, 3] and "three_cycle" not in witnesses:
            witnesses["three_cycle"] = {"prime": int(p), "factor_degrees": pattern}
        if len(witnesses) == 2:
            break
    if len(witnesses) != 2:
        return {"status": "UNKNOWN_GALOIS_WITNESSES_NOT_FOUND", "witnesses": witnesses}
    real_count = int(pari.polsturm(pari(poly)))
    if real_count != 4:
        return {"status": "REJECTED_NOT_COMPLETELY_REAL", "real_root_count": real_count, "galois_witnesses": witnesses}
    return {"status": "PASS", "field_discriminant": str(nf.disc()), "quartic_order_maximal": True,
            "galois_group": "S4", "galois_witnesses": witnesses, "real_root_count": real_count,
            "proof": "An irreducible quartic gives transitivity; the two unramified factorisation patterns provide a 4-cycle and a 3-cycle, forcing S4."}


def multiply_mod_cubic(left, right, coefficients):
    a0, a1, a2, leading = coefficients
    assert leading == 1
    product = [QQ(0)] * 5
    for i, u in enumerate(left):
        for j, v in enumerate(right):
            product[i + j] += u * v
    return [product[0] - a0 * product[3] + a0 * a2 * product[4],
            product[1] - a1 * product[3] + (-a0 + a1 * a2) * product[4],
            product[2] - a2 * product[3] + (-a1 + a2**2) * product[4]]


def point_beta(point, f, nf):
    x, y = map(QQ, point.xy())
    X, Y = 4 * x, 8 * y + 4 * x + 4
    assert Y * Y == f(X)
    den = ZZ(X.denominator()).sqrt()
    assert den in ZZ and den * den == X.denominator()
    beta = [QQ(X * den**2), -QQ(den**2), QQ(0)]
    theta = pari.Mod("z", pari(f))
    assert pari.nfeltnorm(nf, sum(pari(value) * theta**i for i, value in enumerate(beta))) == (Y * den**3) ** 2
    return beta


def m24_strict_character_basis(f, nf):
    """Reconstruct the three existing M24 strict squareclasses in raw coordinates."""
    data = read(FILTRATION)
    words = data["local_filtration_mod_2"]["M24_strict_kernel_public_words"]
    assert len(words) == 3
    E = EllipticCurve(QQ, [QQ(value) for value in curve.GENERAL_WEIERSTRASS_COEFFICIENTS])
    betas = [point_beta(E(point), f, nf) for point in curve.POINTS]
    coefficients = [QQ(value) for value in f.list()]
    answer = []
    for word in words:
        alpha = [QQ(1), QQ(0), QQ(0)]
        for bit, beta in zip(word, betas):
            if bit:
                alpha = multiply_mod_cubic(alpha, beta, coefficients)
        answer.append(alpha)
    return answer


def squareclass_character_certificate(alpha, f, nf, prime_bound):
    """Exact split-good-prime square characters against the three M24 directions."""
    basis = m24_strict_character_basis(f, nf)
    rows = basis + [alpha]
    chars = [[] for _ in rows]
    blocks = []
    for p in prime_range(3, prime_bound + 1):
        if f.discriminant() % p == 0:
            continue
        roots = [QQ(root) for root in f.change_ring(GF(p)).roots(multiplicities=False)]
        if len(roots) != 3:
            continue
        try:
            values = [[int(sum(value * root**i for i, value in enumerate(row)) % p) for root in roots] for row in rows]
        except (TypeError, ZeroDivisionError):
            continue
        if any(value == 0 for row in values for value in row):
            continue
        for target, value_row in zip(chars, values):
            target.extend(int(pow(value, (p - 1) // 2, p) == p - 1) for value in value_row)
        old_rank, rank = matrix(GF(2), chars[:3]).rank(), matrix(GF(2), chars).rank()
        blocks.append({"prime": int(p), "roots": [int(root) for root in roots], "M24_rank": int(old_rank), "with_candidate_rank": int(rank)})
        if old_rank == 3 and rank == 4:
            return {"status": "PASS", "character_blocks": blocks, "rank_M24": 3, "rank_with_candidate": 4}
    return {"status": "UNKNOWN_NO_BOUNDED_INDEPENDENCE_WITNESS", "character_blocks": blocks}


def strict_class_certificate(alpha, raw, nf, prime_bound):
    """Exact ideal-square and S-local square tests for an extracted alpha."""
    S = [int(p) for p in read(STRICT)["S_finite"]]
    theta = pari.Mod("z", pari(raw))
    value = sum(pari(QQ(coefficient)) * theta**i for i, coefficient in enumerate(alpha))
    norm = QQ(pari.nfeltnorm(nf, value))
    if not (norm > 0 and norm.sqrt() in QQ):
        return {"status": "REJECTED_NORM_NOT_POSITIVE_SQUARE", "norm": str(norm)}
    factor = pari.idealfactor(nf, value)
    half = pari.idealhnf(nf, 1)
    entries = []
    for col in range(factor.ncols()):
        prime, exponent = factor[0, col], int(factor[1, col])
        if exponent % 2:
            return {"status": "REJECTED_PRINCIPAL_IDEAL_NOT_A_SQUARE", "norm": str(norm), "odd_exponent": exponent}
        half = pari.idealmul(nf, half, pari.idealpow(nf, prime, exponent // 2))
        entries.append({"prime_hnf": str(pari.idealhnf(nf, prime)), "exponent": exponent})
    half = pari.idealhnf(nf, half)
    assert pari.idealpow(nf, half, 2) == pari.idealhnf(nf, value)
    basis_denominator = pari.denominator(pari.nfalgtobasis(nf, value))
    local_value = value * basis_denominator**2
    local = []
    for p in S:
        for index, prime in enumerate(pari.idealprimedec(nf, p)):
            square = bool(pari.nfislocalpower(nf, prime, local_value, 2))
            local.append({"rational_prime": p, "index": index, "prime_hnf": str(pari.idealhnf(nf, prime)), "square": square})
    if not all(row["square"] for row in local):
        return {"status": "REJECTED_NOT_S_SPLIT", "norm": str(norm), "local_squares": local}
    roots = raw.roots(AA, multiplicities=False)
    real_values = [sum(QQ(coefficient) * root**i for i, coefficient in enumerate(alpha)) for root in roots]
    if not all(value > 0 for value in real_values):
        return {"status": "REJECTED_NOT_REAL_SPLIT", "norm": str(norm), "real_signs": [int(value > 0) - int(value < 0) for value in real_values]}
    independent = squareclass_character_certificate(alpha, raw, nf, prime_bound)
    return {"status": "PASS" if independent["status"] == "PASS" else "UNKNOWN_INDEPENDENCE", "alpha_raw_power_basis": [str(value) for value in alpha],
            "norm": str(norm), "norm_square_root": str(norm.sqrt()), "half_ideal_hnf": string_matrix(half),
            "ideal_factorization": entries, "local_square_certificates": local,
            "real_signs": [1, 1, 1], "independence": independent}


def classical_resolvent(poly):
    """The standard cubic resolvent of a monic quartic polynomial."""
    R = poly.parent()
    u = R.gen()
    a, b, c, d = poly[3], poly[2], poly[1], poly[0]
    return u**3 - b * u**2 + (a * c - 4 * d) * u + (4 * b * d - a * a * d - c * c)


def extract_squareclass(poly, raw):
    """Try the norm-square quartic chart; failure is explicit rather than guessed."""
    resolvent = classical_resolvent(poly)
    if not resolvent.is_irreducible():
        return {"status": "REJECTED_REDUCIBLE_CLASSICAL_RESOLVENT"}
    Kc = resolvent.number_field("u")
    Kr = raw.number_field("z")
    outcome = Kc.is_isomorphic(Kr, isomorphism_maps=True)
    if not outcome or outcome[0] is not True:
        return {"status": "UNKNOWN_CUBIC_FIELD_IDENTIFICATION"}
    to_raw = outcome[1]
    # In the norm-square quartic chart, a resolvent root is alpha-Tr(alpha)/2.
    # The trace of alpha is the negative quartic x^3 coefficient divided by two.
    trace_alpha = -QQ(poly[3]) / 2
    alpha = to_raw(Kc.gen()) + trace_alpha
    coefficients = [QQ(alpha.polynomial()[i]) for i in range(3)]
    return {"status": "PASS", "alpha_raw_power_basis": coefficients,
            "classical_resolvent_ascending": [str(value) for value in resolvent.list()],
            "identification": "Exact Sage number-field isomorphism from the classical resolvent field to the raw cubic field."}


def integral_quadrics(alpha, raw):
    """Build the two standard 2-cover quadrics, preserving the existing cover convention."""
    R = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    u, v, w, z = R.gens()
    product = multiply_mod_cubic(alpha, multiply_mod_cubic([u, v, w], [u, v, w], [QQ(x) for x in raw.list()]), [QQ(x) for x in raw.list()])
    def primitive_integral(poly):
        scale = 1
        for coefficient in poly.coefficients():
            scale = scale.lcm(coefficient.denominator())
        integral = poly * scale
        content = 0
        for coefficient in integral.coefficients():
            content = gcd(content, abs(int(coefficient)))
        assert content
        return integral / content
    return [str(primitive_integral(product[1] + z**2)), str(primitive_integral(product[2]))]


def handoff_to_sealed_cover(candidate, raw, protocol):
    """This is called only after strict+independence PASS; no pre-acceptance PointsQI."""
    alpha = [QQ(value) for value in candidate["strict"]["alpha_raw_power_basis"]]
    token = "fixed-resolvent-" + sha256(json.dumps(candidate, sort_keys=True).encode()).hexdigest()[:16]
    sealed = WORK / (token + ".json")
    job = {"schema": "elliptic-curves.sealed-two-quadric-job.v1", "token": token,
           "integral_quadrics": integral_quadrics(alpha, raw),
           "solver": {"engine": "Magma GenusOneModel Minimise/Reduce followed by PointsQI", "pointsqi_height": 10000000,
                      "remote_request_timeout_seconds": 50}}
    put_new(sealed, job)
    # The existing isolated worker receives only this sealed file.  The parent
    # records its output for the follow-on exact E302 map/recheck.
    subprocess.run(["sage", "-python", str(HIDDEN_SOLVER), "--solver", "--sealed", str(sealed)], cwd=ROOT, check=False, timeout=58)
    result = sealed.with_name(sealed.stem + "-result.json")
    return {"status": "DISPATCHED", "sealed_job": relative(sealed), "sealed_job_sha256": digest(sealed),
            "solver_result": relative(result) if result.exists() else None,
            "solver_result_sha256": digest(result) if result.exists() else None}


def compute():
    protocol = read(PROTOCOL)
    preflight, raw, nf = maximal_cubic_preflight()
    assert preflight["status"] == "PASS"
    form = [ZZ(value) for value in preflight["binary_cubic_descending"]]
    field_discriminant = ZZ(preflight["field_discriminant"])
    bound = protocol["bounds"]["affine_slice_parameter_bound"]
    candidates, accepted = [], []
    for px, py in itertools.product(range(-bound, bound + 1), repeat=2):
        A, B = matrix_pair_for_affine_slice(form, px, py)
        coefficients = pair_coefficients(A, B)
        assert primitive(coefficients)
        determinant = determinant_coefficients(A, B)
        assert determinant == form
        table, multiplication = quartic_multiplication(coefficients)
        probe = algebra_probe(multiplication, protocol["bounds"]["quartic_generator_l1_bound"])
        row = {"slice_parameters": {"p_x": px, "p_y": py}, "pair_coefficients": [str(value) for value in coefficients],
               "determinant_form_descending": [str(value) for value in determinant], "primitive": True,
               "quartic_ring_associative_integral": True, "algebra_probe": probe}
        if probe["status"] != "FIELD_GENERATOR":
            row["status"] = "REJECTED_ZERO_DIVISOR" if probe["status"] == "ZERO_DIVISOR" else "UNKNOWN_NO_FIELD_WITNESS"
            candidates.append(row)
            continue
        Q = PolynomialRing(QQ, "X")
        polynomial = Q([QQ(value) for value in probe["polynomial_ascending"]])
        structural = structural_field_certificate(polynomial, field_discriminant, protocol["bounds"]["finite_field_witness_prime_bound"])
        row["quartic_generator_polynomial_ascending"] = probe["polynomial_ascending"]
        row["structural"] = structural
        if structural["status"] != "PASS":
            row["status"] = structural["status"]
            candidates.append(row)
            continue
        extraction = extract_squareclass(polynomial, raw)
        row["squareclass_extraction"] = {key: ([str(value) for value in value] if key == "alpha_raw_power_basis" and isinstance(value, list) else value) for key, value in extraction.items()}
        if extraction["status"] != "PASS":
            row["status"] = extraction["status"]
            candidates.append(row)
            continue
        strict = strict_class_certificate(extraction["alpha_raw_power_basis"], raw, nf, protocol["bounds"]["finite_field_witness_prime_bound"])
        row["strict"] = strict
        if strict["status"] != "PASS":
            row["status"] = strict["status"]
            candidates.append(row)
            continue
        row["status"] = "ACCEPTED_STRICT_INDEPENDENT"
        accepted.append(row)
        candidates.append(row)
        # Protocol stops at the first independent strict character.  The cover
        # handoff is intentionally the only branch allowed to call PointsQI.
        row["cover_handoff"] = handoff_to_sealed_cover(row, raw, protocol)
        break
    status = "PASS_NEW_STRICT_CHARACTER_HANDED_OFF" if accepted else "BOUNDED_SLICE_NO_ACCEPTED_CHARACTER"
    return {"schema": "elliptic-curves.curve302-fixed-resolvent-orbits.v1", "status": status,
            "software": {"sage": SAGE_VERSION, "pari": str(pari.version())}, "protocol": relative(PROTOCOL), "protocol_sha256": digest(PROTOCOL),
            "bindings": {relative(path): digest(path) for path in (PROTOCOL, MAXIMUM_FORM, STRICT, FILTRATION, CURVE, HIDDEN_SOLVER, Path(__file__))},
            "maximal_cubic_preflight": preflight,
            "bounded_search": {"normalization": protocol["normalization"], "parameter_box": [-bound, bound],
                                "raw_pairs_enumerated": len(candidates), "raw_pair_candidates": candidates,
                                "accepted_count": len(accepted), "completeness": "NOT_CLAIMED"},
            "next_action": "SEALED_2COVER_POINTSQI" if accepted else "STOP_BOUNDED_SLICE_NO_CLASS_CREATED",
            "boundary": protocol["boundary"]}


def capture():
    WORK.mkdir(parents=True, exist_ok=True)
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    put_new(OUTPUT, compute())
    print(read(OUTPUT)["status"], flush=True)


def check():
    stored = read(OUTPUT)
    assert stored["bindings"] == {relative(path): digest(path) for path in (PROTOCOL, MAXIMUM_FORM, STRICT, FILTRATION, CURVE, HIDDEN_SOLVER, Path(__file__))}
    # A check intentionally recomputes the bounded arithmetic but never dispatches
    # a solver unless the recorded output itself contains a genuine accepted hit.
    if stored["status"] != "PASS_NEW_STRICT_CHARACTER_HANDED_OFF":
        replay = compute()
        assert replay == stored
    print("PASS fixed-resolvent maximal-order preflight and bounded slice replay", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("capture", "check"))
    args = parser.parse_args()
    if args.mode == "capture":
        capture()
    else:
        check()
