#!/usr/bin/env sage-python
"""Exact, witness-free preprocessing for the eight strict E302 four-covers.

``--capture`` is the public preparation boundary.  It makes the fixed eight
strict-coset inputs, attaches only the public discriminant-prime support, and
launches one fresh ``--worker`` process for each cover.  A worker reads just a
sealed pair of quadrics and that support: it does not import the curve, point
table, Kummer words, alpha, a control point, or a cover witness.

For each worker the pipeline is:

* verify the integral ternary conic and its determinant decomposition;
* ask PARI ``qfsolve``/``qfparam`` for a conic parametrization using that
  exact decomposition;
* derive the binary quartic; and
* apply PARI ``hyperellred``, ``hyperellminimalmodel``, ``hyperellred``.

Every parametrization and every binary-quartic transport is checked as a
polynomial identity over QQ.  A PARI factorization-interface rejection is a
method status, not an assertion about rational points or solubility.
"""

from __future__ import annotations

import argparse
import importlib.machinery
import json
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import time

from sage.all import QQ, ZZ, PolynomialRing, matrix, vector
from sage.libs.pari.all import pari
from sage.version import version as SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
SOURCE = CAS / "benchmark_curve302_strict_covers_hidden.sage"
FACTOR_PROOFS = ART / "record_prime_factor_proofs_20260904.json"
OUTPUT = ART / "curve302_strict_cover_preprocessing_witness_free_v1.json"
WORK = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-strict-cover-preprocessing-witness-free-v1"

POLICY = {
    "fixed_primary_tokens": [f"case-{index:02d}" for index in range(8, 16)],
    "conic_algorithm": "PARI qfsolve followed by qfparam with the supplied exact determinant decomposition",
    "quartic_algorithms": ["PARI hyperellred", "PARI hyperellminimalmodel", "PARI hyperellred"],
    "worker_timeout_seconds": 120,
    "selection": "all eight strict-coset controls; no witness-dependent ordering, tuning, or stopping",
}


def read(path: Path):
    return json.loads(path.read_text())


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def put_new(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def source_module():
    return importlib.machinery.SourceFileLoader("curve302_hidden_preprocess_source", str(SOURCE)).load_module()


def coefficient_strings(poly):
    return [str(poly[index]) for index in range(5)]


def binary_poly(values):
    ring = PolynomialRing(QQ, "x")
    return ring([QQ(value) for value in values])


def matrix_strings(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def matrix_from_strings(value):
    return matrix(QQ, [[QQ(entry) for entry in row] for row in value])


def polynomial_bits(poly):
    nonzero = [coefficient for coefficient in poly.coefficients() if coefficient]
    assert nonzero
    return max(max(abs(coefficient.numerator()).nbits(), coefficient.denominator().nbits())
               for coefficient in nonzero)


def public_factor_support():
    """Read the existing primality certificate; do not factor in a worker."""
    record = read(FACTOR_PROOFS)
    assert record["method"] == "PARI isprime: unconditional primality test"
    proved = set(record["proved_primes"])
    factors = record["records"]["302"]["factorizations"]["DISCRIMINANT_FACTORIZATION"]
    support = []
    for prime, exponent in factors:
        assert prime in proved and int(exponent) > 0
        if prime not in support:
            support.append(prime)
    assert support
    return support


def parse_quadrics(quadrics):
    ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))
    u, v, w, z = ring.gens()
    first, second = [ring(text) for text in quadrics]
    assert first.total_degree() == second.total_degree() == 2
    assert first.monomial_coefficient(z * z) == 1
    assert first == first(u, v, w, 0) + z * z
    assert second == second(u, v, w, 0)
    xyz = (u, v, w)
    gram = matrix(QQ, 3, 3,
                  lambda row, column: second.monomial_coefficient(xyz[row] * xyz[column]) /
                  (1 if row == column else 2))
    assert all(entry.denominator() == 1 for entry in gram.list())
    gram = matrix(ZZ, gram)
    gram_form = sum(gram[row, column] * xyz[row] * xyz[column]
                    for row in range(3) for column in range(3))
    assert gram_form == second
    assert gram.det()
    return ring, first, second, gram


def determinant_certificate(gram, support):
    """Exact prime-support times square decomposition of abs(det(gram)).

    The support primes have their primality proved by FACTOR_PROOFS before the
    sealed job is made.  The final base is intentionally *not* called prime:
    only its square contribution is certified here.
    """
    absolute = abs(ZZ(gram.det()))
    remaining = absolute
    valuations = []
    for text in support:
        prime = ZZ(text)
        assert prime > 1
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if exponent:
            valuations.append([str(prime), exponent])
    root = remaining.sqrt()
    assert root * root == remaining
    reconstructed = root * root
    for text, exponent in valuations:
        reconstructed *= ZZ(text) ** exponent
    assert reconstructed == absolute
    return {
        "absolute_determinant": str(absolute),
        "support_valuations": valuations,
        "square_cofactor_base": str(root),
        "square_cofactor_exponent": 2,
        "exact_product_verified": True,
    }


def factor_matrix(certificate):
    rows = [[ZZ(prime), ZZ(exponent)] for prime, exponent in certificate["support_valuations"]]
    root = ZZ(certificate["square_cofactor_base"])
    if root != 1:
        rows.append([root, ZZ(certificate["square_cofactor_exponent"])])
    assert rows
    return pari.matrix(len(rows), 2, [entry for row in rows for entry in row])


def homogeneous_pullback(quartic, change):
    ring = quartic.parent()
    x = ring.gen()
    s, t = change * vector([x, 1])
    return ring(sum(quartic[index] * s**index * t**(4 - index) for index in range(5)))


def transport_ratio(source, target, change):
    """Certify source(s,t) = square * target(x) exactly over QQ."""
    assert source.degree() == target.degree() == 4 and change.det()
    pull = homogeneous_pullback(source, change)
    index = next(index for index in range(5) if target[index])
    ratio = QQ(pull[index] / target[index])
    assert ratio and ratio.is_square() and pull == ratio * target
    return ratio


def conic_and_quartic(first, second, gram, certificate):
    """Run and verify qfsolve/qfparam, then return the exact raw quartic."""
    solution = pari.qfsolve([pari(gram), factor_matrix(certificate)])
    if solution.type() != "t_COL":
        return None, {"qfsolve_return": str(solution)}
    point = vector(ZZ, solution)
    assert point * gram * point == 0
    parametrization = matrix(QQ, pari.qfparam(pari(gram), solution, 1))
    assert parametrization.dimensions() == (3, 3) and parametrization.det()
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    gamma = parametrization * vector([x * x, x, 1])
    old_ring = first.parent()
    u, v, w, _ = old_ring.gens()
    conic_identity = ring(second(gamma[0], gamma[1], gamma[2], 0))
    assert not conic_identity
    raw = ring(-first(gamma[0], gamma[1], gamma[2], 0))
    assert raw.degree() == 4
    assert ring(first(gamma[0], gamma[1], gamma[2], 0) + raw) == 0
    return {
        "qfsolve_point": [str(entry) for entry in point],
        "qfsolve_point_verified": True,
        "parametrization": matrix_strings(parametrization),
        "conic_parameter_identity_verified": True,
        "raw_quartic": coefficient_strings(raw),
        "raw_quartic_coefficient_bits": polynomial_bits(raw),
        "raw_cover_identity_verified": True,
    }, None


def minimize_and_reduce(raw):
    """Apply the fixed exact PARI minimization/reduction sequence."""
    ring = raw.parent()
    current = raw
    composite = matrix.identity(QQ, 2)
    steps = []
    for function_name in ("hyperellred", "hyperellminimalmodel", "hyperellred"):
        denominator_square = current.denominator() ** 2
        source = ring(current * denominator_square)
        assert all(coefficient.denominator() == 1 for coefficient in source.coefficients())
        function = pari(f"(q)->{{my(m);my(z={function_name}(q,&m));[z,m]}}")
        z, transformation = function(source)
        target = ring(z[0]) + ring(z[1]) ** 2 / 4
        change = matrix(QQ, transformation[1])
        ratio = transport_ratio(source, target, change)
        steps.append({
            "algorithm": function_name,
            "input_quartic": coefficient_strings(source),
            "input_denominator_square": str(denominator_square),
            "output_quartic": coefficient_strings(target),
            "parameter_change": matrix_strings(change),
            "pullback_ratio": str(ratio),
            "pullback_ratio_is_rational_square": True,
            "exact_pullback_identity_verified": True,
            "output_coefficient_bits": polynomial_bits(target),
        })
        current = target
        composite *= change
    final_ratio = transport_ratio(raw, current, composite)
    return {
        "steps": steps,
        "final_quartic": coefficient_strings(current),
        "final_quartic_coefficient_bits": polynomial_bits(current),
        "composite_parameter_change": matrix_strings(composite),
        "raw_to_final_pullback_ratio": str(final_ratio),
        "raw_to_final_ratio_is_rational_square": True,
        "raw_to_final_exact_pullback_identity_verified": True,
    }


def verify_success(sealed, result):
    """Replay every stored algebraic identity without solving or searching."""
    ring, first, second, gram = parse_quadrics(sealed["integral_quadrics"])
    certificate = determinant_certificate(gram, sealed["determinant_factor_support"])
    assert result["determinant_certificate"] == certificate
    data = result["preprocessing"]
    point = vector(ZZ, [ZZ(entry) for entry in data["qfsolve_point"]])
    assert point * gram * point == 0
    parametrization = matrix_from_strings(data["parametrization"])
    assert parametrization.dimensions() == (3, 3) and parametrization.det()
    binary = PolynomialRing(QQ, "x")
    x = binary.gen()
    gamma = parametrization * vector([x * x, x, 1])
    assert binary(second(gamma[0], gamma[1], gamma[2], 0)) == 0
    raw = binary_poly(data["raw_quartic"])
    assert raw == binary(-first(gamma[0], gamma[1], gamma[2], 0))
    current = raw
    composite = matrix.identity(QQ, 2)
    for step in data["reduction"]["steps"]:
        source = binary_poly(step["input_quartic"])
        assert source == binary(current * current.denominator() ** 2)
        change = matrix_from_strings(step["parameter_change"])
        target = binary_poly(step["output_quartic"])
        assert transport_ratio(source, target, change) == QQ(step["pullback_ratio"])
        current = target
        composite *= change
    assert current == binary_poly(data["reduction"]["final_quartic"])
    assert composite == matrix_from_strings(data["reduction"]["composite_parameter_change"])
    assert transport_ratio(raw, current, composite) == QQ(data["reduction"]["raw_to_final_pullback_ratio"])
    return True


def worker(sealed_path: Path) -> None:
    """Worker arithmetic boundary: no curve-302 source is imported here."""
    sealed = read(sealed_path)
    assert set(sealed) == {"schema", "token", "integral_quadrics", "determinant_factor_support"}
    assert sealed["schema"] == "elliptic-curves.sealed-strict-cover-preprocess-job.v1"
    started = time.monotonic()
    result = {
        "schema": "elliptic-curves.sealed-strict-cover-preprocess-result.v1",
        "token": sealed["token"],
        "sealed_input_sha256": digest(sealed_path),
        "engine": {"sage": SAGE_VERSION, "pari": str(pari("version()"))},
        "status": "WORKER_ERROR",
    }
    try:
        _, first, second, gram = parse_quadrics(sealed["integral_quadrics"])
        certificate = determinant_certificate(gram, sealed["determinant_factor_support"])
        result["determinant_certificate"] = certificate
        try:
            preprocessing, no_parameter = conic_and_quartic(first, second, gram, certificate)
        except Exception as error:
            # A factor-matrix interface failure is not a local/global conic claim.
            result["status"] = "CONIC_SOLVER_REJECTED_DETERMINANT_ENCODING"
            result["conic_solver_exception_type"] = type(error).__name__
        else:
            if preprocessing is None:
                result["status"] = "CONIC_SOLVER_NO_PARAMETRIZATION"
                result["conic_solver_return"] = no_parameter
            else:
                preprocessing["reduction"] = minimize_and_reduce(binary_poly(preprocessing["raw_quartic"]))
                result["preprocessing"] = preprocessing
                result["status"] = "PREPROCESSED_WITH_EXACT_TRANSPORTS"
                assert verify_success(sealed, result)
    except Exception as error:
        result["status"] = "WORKER_ERROR"
        result["exception_type"] = type(error).__name__
    result["wall_seconds"] = time.monotonic() - started
    put_new(sealed_path.with_name(sealed_path.stem + "-result.json"), result)
    print(sealed["token"], result["status"], flush=True)


def sealed_case(case, support):
    return {
        "schema": "elliptic-curves.sealed-strict-cover-preprocess-job.v1",
        "token": case["token"],
        "integral_quadrics": case["integral_quadrics"],
        "determinant_factor_support": support,
    }


def capture():
    source = source_module()
    panel = source.build_panels()
    cases = [case for case in panel["cases"] if case["token"] in POLICY["fixed_primary_tokens"]]
    assert [case["token"] for case in cases] == POLICY["fixed_primary_tokens"]
    support = public_factor_support()
    WORK.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in cases:
        sealed_path = WORK / f"{case['token']}.json"
        result_path = sealed_path.with_name(sealed_path.stem + "-result.json")
        expected = sealed_case(case, support)
        if sealed_path.exists():
            assert read(sealed_path) == expected
        else:
            put_new(sealed_path, expected)
        if not result_path.exists():
            try:
                subprocess.run(["sage", "-python", str(Path(__file__)), "--worker", "--sealed", str(sealed_path)],
                               check=False, timeout=POLICY["worker_timeout_seconds"])
            except subprocess.TimeoutExpired:
                put_new(result_path, {
                    "schema": "elliptic-curves.sealed-strict-cover-preprocess-result.v1",
                    "token": case["token"], "sealed_input_sha256": digest(sealed_path),
                    "status": "WORKER_TIMEOUT",
                })
        result = read(result_path)
        assert result["token"] == case["token"] and result["sealed_input_sha256"] == digest(sealed_path)
        if result["status"] == "PREPROCESSED_WITH_EXACT_TRANSPORTS":
            assert verify_success(expected, result)
        rows.append({
            "token": case["token"],
            "sealed_input": str(sealed_path.relative_to(ROOT)),
            "sealed_input_sha256": digest(sealed_path),
            "result": str(result_path.relative_to(ROOT)),
            "result_sha256": digest(result_path),
            "worker_status": result["status"],
            "raw_quartic_coefficient_bits": result.get("preprocessing", {}).get("raw_quartic_coefficient_bits"),
            "final_quartic_coefficient_bits": result.get("preprocessing", {}).get("reduction", {}).get("final_quartic_coefficient_bits"),
        })
    completed = [row for row in rows if row["worker_status"] == "PREPROCESSED_WITH_EXACT_TRANSPORTS"]
    output = {
        "schema": "elliptic-curves.curve302-strict-cover-preprocessing-witness-free.v1",
        "status": "PREPROCESSING_COMPLETED",
        "protocol": "fixed eight known-soluble strict-coset controls; a worker sees only its quadrics and public determinant-prime support",
        "policy": POLICY,
        "bindings": {str(path.relative_to(ROOT)): digest(path) for path in (SOURCE, FACTOR_PROOFS, Path(__file__))},
        "factor_certificate_source": str(FACTOR_PROOFS.relative_to(ROOT)),
        "cases": rows,
        "exactly_preprocessed_count": len(completed),
        "boundary": "A preprocessing or conic-interface failure is not a claim about a cover's rational points, solubility, rank, or Sha.",
    }
    put_new(OUTPUT, output)
    print(output["status"], len(completed), flush=True)


def check():
    stored = read(OUTPUT)
    bindings = {str(path.relative_to(ROOT)): digest(path) for path in (SOURCE, FACTOR_PROOFS, Path(__file__))}
    assert stored["bindings"] == bindings
    assert stored["policy"] == POLICY and len(stored["cases"]) == 8
    assert public_factor_support()
    for row in stored["cases"]:
        sealed_path = ROOT / row["sealed_input"]
        result_path = ROOT / row["result"]
        assert digest(sealed_path) == row["sealed_input_sha256"]
        assert digest(result_path) == row["result_sha256"]
        result = read(result_path)
        if result["status"] == "PREPROCESSED_WITH_EXACT_TRANSPORTS":
            assert verify_success(read(sealed_path), result)
    print("PASS exact replay of witness-free strict-cover preprocessing", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--sealed", type=Path)
    args = parser.parse_args()
    assert sum((args.capture, args.check, args.worker)) == 1
    if args.worker:
        assert args.sealed is not None
        worker(args.sealed)
    elif args.capture:
        capture()
    else:
        check()


if __name__ == "__main__":
    main()
