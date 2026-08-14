#!/usr/bin/env python3
"""Exhaust the 80 independent Fermigier bidegree-(2,1) pencils exactly.

For every pair in the rank-10-by-rank-8 exceptional quotient product, this
script studies

    x(T) = (a(c)T+b(c)+k(T-T_22)(T-T_20))/(cT+1),

where the endpoint conditions determine ``a(c),b(c)``.  The exact algebra
exhibits five universal rational discriminant components.  Their squareclass
kernels are factored over the appropriate rational function field.  After
those factors are stripped, irreducibility of the degree-32 residual over QQ
is certified by an irreducible good reduction modulo a recorded prime.

A characteristic-zero bivariate factorization is attempted only if the
pinned modular-prime list supplies no irreducibility witness.  Thus this is a
bounded structural classification, not a rational-point or score sweep.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import platform
import re
import shutil
import subprocess
from typing import Any, Iterable, Sequence

import sympy as sp

from analyze_fermigier_bidegree21_pilot import component_record, pencil_polynomial
from analyze_fermigier_exceptional_transport import (
    E22_INDEPENDENT_EXCEPTIONAL_LABELS,
    E22_T,
    E22_U,
    RANK20_T,
    RANK20_U,
    e22_exceptional_quotient,
    rank20_exceptional_quotient,
    rational_text,
)


Q = Fraction
OUTPUT_RELATIVE = Path(
    "artifacts/generated-results/elliptic_fermigier_bidegree21_all80.json"
)
MODULAR_PRIMES = (101, 103, 107, 109, 113, 127, 131)
FAMILY_COEFFICIENTS = (
    18103855887324900,
    102302344648,
    -879500,
    -257843832010380,
    -650709150,
    1860,
    1195214262641,
    1718550,
    -2,
    -2051321790,
    -1860,
    1149050,
)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def sha256_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update((line + "\n").encode())
    return digest.hexdigest()


def primitive_linear_coefficients(
    left_x: Fraction, right_x: Fraction
) -> tuple[int, int]:
    """Primitive coefficients A,B for A*c+B*k=0, the cancellation locus."""

    first = left_x - right_x
    second = RANK20_T - E22_T
    denominator_lcm = sp.ilcm(first.denominator, second.denominator)
    a = first.numerator * (denominator_lcm // first.denominator)
    b = second.numerator * (denominator_lcm // second.denominator)
    divisor = sp.igcd(abs(a), abs(b))
    a //= divisor
    b //= divisor
    if a < 0 or (a == 0 and b < 0):
        a, b = -a, -b
    return int(a), int(b)


def independent_pairs(
    root: Path,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    tuple[tuple[str, Fraction, str, Fraction], ...],
]:
    e22_record, e22_points = e22_exceptional_quotient()
    rank20_record, rank20_points = rank20_exceptional_quotient(root)
    allowed = set(E22_INDEPENDENT_EXCEPTIONAL_LABELS)
    e22_independent = tuple(row for row in e22_points if row[0] in allowed)
    if tuple(label for label, _ in e22_independent) != E22_INDEPENDENT_EXCEPTIONAL_LABELS:
        raise AssertionError("the independent E22 exceptional population changed")
    if len(rank20_points) != 8:
        raise AssertionError("the independent rank-20 exceptional population changed")
    pairs = tuple(
        (left_label, left_x, right_label, right_x)
        for left_label, left_x in e22_independent
        for right_label, right_x in rank20_points
    )
    if len(pairs) != 80:
        raise AssertionError("the bidegree-(2,1) population is no longer 10*8")
    return e22_record, rank20_record, pairs


def exact_component_classification(
    pair: tuple[str, Fraction, str, Fraction],
) -> dict[str, Any]:
    """Classify all five rational discriminant components for one pencil."""

    left_label, left_x, right_label, right_x = pair
    polynomial, symbols = pencil_polynomial(left_x, right_x)
    T, c, k = symbols["T"], symbols["c"], symbols["k"]
    generic = sp.Poly(polynomial.as_expr(), T, domain=sp.QQ.frac_field(c, k))
    leading_scalar = sp.cancel(generic.LC() / ((c - k) ** 2 * (c + k) ** 2))
    if (
        generic.degree() != 10
        or leading_scalar == 0
        or leading_scalar.free_symbols
        or sp.factor(generic.LC())
        != sp.factor(leading_scalar * (c - k) ** 2 * (c + k) ** 2)
    ):
        raise AssertionError(f"generic degree/leading coefficient changed for {left_label}x{right_label}")

    first = sp.Rational(E22_T.numerator, E22_T.denominator)
    second = sp.Rational(RANK20_T.numerator, RANK20_T.denominator)
    slope = -sp.Rational(
        (left_x - right_x).numerator, (left_x - right_x).denominator
    ) / sp.Rational(
        (RANK20_T - E22_T).numerator,
        (RANK20_T - E22_T).denominator,
    )
    cancellation_a, cancellation_b = primitive_linear_coefficients(left_x, right_x)
    cancellation_label = f"{cancellation_a}*c+{cancellation_b}*k=0"

    numerator = symbols["numerator"]
    denominator = symbols["denominator"]
    for c_value, anchor in ((-1 / first, first), (-1 / second, second)):
        if denominator.subs({c: c_value, T: anchor}) != 0:
            raise AssertionError("anchor-pole denominator identity failed")
        if sp.cancel(numerator.subs({c: c_value, T: anchor})) != 0:
            raise AssertionError("anchor-pole numerator identity failed")
    cancellation_numerator = sp.Poly(
        sp.expand(numerator.subs(k, slope * c)),
        T,
        domain=sp.QQ.frac_field(c),
    )
    cancellation_denominator = sp.Poly(
        denominator, T, domain=sp.QQ.frac_field(c)
    )
    if not sp.rem(cancellation_numerator, cancellation_denominator).is_zero:
        raise AssertionError("the declared cancellation component does not cancel")

    components = [
        {
            "type": "degree_drop_plus",
            **component_record(
                polynomial,
                symbols,
                label="k-c=0",
                substitutions={k: c},
                validity="valid degree-drop locus",
                interpretation=(
                    "The T^10 and T^9 coefficients vanish; the remaining exact "
                    "squareclass is computed over QQ(c)."
                ),
            ),
        },
        {
            "type": "degree_drop_minus",
            **component_record(
                polynomial,
                symbols,
                label="k+c=0",
                substitutions={k: -c},
                validity="valid degree-drop locus",
                interpretation=(
                    "The T^10 and T^9 coefficients vanish; the remaining exact "
                    "squareclass is computed over QQ(c)."
                ),
            ),
        },
        {
            "type": "numerator_denominator_cancellation",
            **component_record(
                polynomial,
                symbols,
                label=cancellation_label,
                substitutions={k: slope * c},
                validity="valid but representation is non-reduced",
                interpretation=(
                    "The numerator is divisible by c*T+1; cancellation returns "
                    "the already-classified affine transport."
                ),
            ),
        },
        {
            "type": "rank20_anchor_pole",
            **component_record(
                polynomial,
                symbols,
                label="28917*c+10=0",
                substitutions={c: -1 / second},
                validity="invalid at rank20 anchor",
                interpretation=(
                    "Numerator and denominator both vanish at the rank-20 anchor."
                ),
            ),
        },
        {
            "type": "E22_anchor_pole",
            **component_record(
                polynomial,
                symbols,
                label="39508*c+39=0",
                substitutions={c: -1 / first},
                validity="invalid at E22 anchor",
                interpretation="Numerator and denominator both vanish at the E22 anchor.",
            ),
        },
    ]
    low_genus = [
        item
        for item in components
        if item["validity"].startswith("valid") and item["genus_at_most_one"]
    ]
    return {
        "left": left_label,
        "left_x": rational_text(left_x),
        "right": right_label,
        "right_x": rational_text(right_x),
        "generic_T_degree": generic.degree(),
        "cancellation_factor": {
            "primitive_coefficients_c_k": [cancellation_a, cancellation_b],
            "equation": cancellation_label,
            "derived_from": "(left_x-right_x)*c+(T_rank20-T_E22)*k=0",
        },
        "exact_identity_checks": {
            "generic_leading_coefficient": (
                f"{rational_text(leading_scalar)}*(c-k)^2*(c+k)^2"
            ),
            "anchor_endpoint_numerator_denominator_vanishing": True,
            "cancellation_divisibility_by_cT_plus_1": True,
        },
        "rational_components": components,
        "valid_genus_at_most_one_components": low_genus,
    }


def _fraction_mod(value: Fraction, prime: int) -> int:
    return value.numerator * pow(value.denominator, -1, prime) % prime


def good_reduction_prime(
    left_x: Fraction, right_x: Fraction, prime: int
) -> tuple[bool, str]:
    values = (E22_T, RANK20_T, left_x, right_x)
    if prime <= 2 or any(value.denominator % prime == 0 for value in values):
        return False, "input denominator or characteristic two"
    first = _fraction_mod(E22_T, prime)
    second = _fraction_mod(RANK20_T, prime)
    left = _fraction_mod(left_x, prime)
    right = _fraction_mod(right_x, prime)
    difference_t = (second - first) % prime
    difference_x = (left - right) % prime
    if first == 0 or second == 0 or difference_t == 0:
        return False, "anchor or anchor difference degenerates"
    if (difference_x - difference_t) % prime == 0:
        return False, "cancellation and k+c components collide"
    if (difference_x + difference_t) % prime == 0:
        return False, "cancellation and k-c components collide"
    return True, "good reduction"


def _singular_integer(output: str, marker: str) -> int:
    match = re.search(rf"^{re.escape(marker)}\n(-?\d+)$", output, re.MULTILINE)
    if match is None:
        raise AssertionError(f"missing Singular marker {marker}")
    return int(match.group(1))


def _modular_program(
    left_x: Fraction, right_x: Fraction, prime: int
) -> str:
    first = _fraction_mod(E22_T, prime)
    second = _fraction_mod(RANK20_T, prime)
    left = _fraction_mod(left_x, prime)
    right = _fraction_mod(right_x, prime)
    difference_x = (left - right) % prime
    difference_t = (second - first) % prime
    constants = [value % prime for value in FAMILY_COEFFICIENTS]
    return f"""
ring r={prime},(c,k,T),dp;
number t1={first}; number t2={second};
number x1={left}; number x2={right};
poly a=(x1*(c*t1+1)-x2*(c*t2+1))/(t1-t2);
poly b=x1*(c*t1+1)-a*t1;
poly D=c*T+1;
poly N=a*T+b+k*(T-t1)*(T-t2);
poly f0={constants[0]}+{constants[1]}*T^2+{constants[2]}*T^4+T^6;
poly f1={constants[3]}+{constants[4]}*T^2+{constants[5]}*T^4;
poly f2={constants[6]}+{constants[7]}*T^2+{constants[8]}*T^4;
poly f3={constants[9]}+{constants[10]}*T^2;
poly f4={constants[11]}+T^2;
poly P=f0*D^4+f1*N*D^3+f2*N^2*D^2+f3*N^3*D+f4*N^4;
poly rr=resultant(P,diff(P,T),T);
poly leading_shape=(c-k)^2*(c+k)^2;
poly disc=rr/leading_shape;
if (rr-leading_shape*disc!=0) {{ "@@ERROR_LEADING_DIVISION"; }}
poly known=(k-c)^2*(k+c)^2*(c*t1+1)^12*(c*t2+1)^12
           *({difference_x}*c+{difference_t}*k)^12;
poly residual=disc/known;
if (disc-known*residual!=0) {{ "@@ERROR_KNOWN_DIVISION"; }}
list L=factorize(residual);
int nonconstant=0;
int idx=0;
for (idx=1;idx<=size(L[1]);idx++) {{
  if (deg(L[1][idx])>0) {{ nonconstant=nonconstant+1; }}
}}
intvec wc=1,0,0;
intvec wk=0,1,0;
intvec wt=0,0,1;
"@@RESULTANT_TOTAL_DEGREE"; deg(rr);
"@@DISCRIMINANT_TOTAL_DEGREE"; deg(disc);
"@@RESIDUAL_TOTAL_DEGREE"; deg(residual);
"@@RESIDUAL_C_DEGREE"; deg(residual,wc);
"@@RESIDUAL_K_DEGREE"; deg(residual,wk);
"@@RESIDUAL_T_DEGREE"; deg(residual,wt);
"@@RESIDUAL_TERMS"; size(residual);
"@@NONCONSTANT_FACTOR_COUNT"; nonconstant;
for (idx=1;idx<=size(L[1]);idx++) {{
  if (deg(L[1][idx])>0) {{
    "@@FACTOR_DEGREE"; deg(L[1][idx]);
    "@@FACTOR_EXPONENT"; L[2][idx];
  }}
}}
"@@RESIDUAL_BEGIN"; residual; "@@RESIDUAL_END";
"""


def modular_attempt(
    left_x: Fraction,
    right_x: Fraction,
    prime: int,
    *,
    timeout: int,
) -> dict[str, Any]:
    good, reason = good_reduction_prime(left_x, right_x, prime)
    if not good:
        return {"prime": prime, "status": "skipped", "reason": reason}
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for modular factorization")
    try:
        completed = subprocess.run(
            [executable, "-q"],
            input=_modular_program(left_x, right_x, prime),
            text=True,
            capture_output=True,
            timeout=timeout,
            check=True,
        )
    except subprocess.TimeoutExpired:
        return {"prime": prime, "status": "timeout", "timeout_seconds": timeout}
    if completed.stderr.strip():
        raise RuntimeError(f"Singular wrote to stderr: {completed.stderr.strip()}")
    output = completed.stdout
    if "@@ERROR_" in output or "   ?" in output:
        raise AssertionError(f"Singular modular calculation failed:\n{output}")
    degrees = [int(value) for value in re.findall(r"^@@FACTOR_DEGREE\n(-?\d+)$", output, re.MULTILINE)]
    exponents = [int(value) for value in re.findall(r"^@@FACTOR_EXPONENT\n(-?\d+)$", output, re.MULTILINE)]
    if len(degrees) != len(exponents):
        raise AssertionError("Singular factor degree/exponent output is inconsistent")
    residual_match = re.search(
        r"^@@RESIDUAL_BEGIN\n(.*?)\n@@RESIDUAL_END$",
        output,
        re.MULTILINE | re.DOTALL,
    )
    if residual_match is None:
        raise AssertionError("missing modular residual polynomial")
    residual_text = residual_match.group(1).strip()
    record = {
        "prime": prime,
        "status": "factored",
        "good_reduction_checks": reason,
        "resultant_total_degree": _singular_integer(output, "@@RESULTANT_TOTAL_DEGREE"),
        "discriminant_total_degree": _singular_integer(output, "@@DISCRIMINANT_TOTAL_DEGREE"),
        "residual_total_degree": _singular_integer(output, "@@RESIDUAL_TOTAL_DEGREE"),
        "residual_degree_in_c": _singular_integer(output, "@@RESIDUAL_C_DEGREE"),
        "residual_degree_in_k": _singular_integer(output, "@@RESIDUAL_K_DEGREE"),
        "residual_degree_in_T": _singular_integer(output, "@@RESIDUAL_T_DEGREE"),
        "residual_term_count": _singular_integer(output, "@@RESIDUAL_TERMS"),
        "nonconstant_factor_count": _singular_integer(output, "@@NONCONSTANT_FACTOR_COUNT"),
        "factor_signature": [[degree, exponent] for degree, exponent in zip(degrees, exponents, strict=True)],
        "residual_polynomial_sha256": sha256_text(residual_text + "\n"),
    }
    if record["resultant_total_degree"] > 76:
        raise AssertionError("the resultant exceeded the exact degree bound 76")
    if record["discriminant_total_degree"] > 72:
        raise AssertionError("the discriminant exceeded the exact degree bound 72")
    if record["residual_total_degree"] > 32 or record["residual_degree_in_T"] != 0:
        raise AssertionError("the stripped residual exceeded its exact degree bound")
    record["irreducible_degree_32_witness"] = (
        record["residual_total_degree"] == 32
        and record["nonconstant_factor_count"] == 1
        and record["factor_signature"] == [[32, 1]]
    )
    return record


def _qq_factorization_fallback(
    left_x: Fraction,
    right_x: Fraction,
    *,
    timeout: int,
) -> dict[str, Any]:
    """Run the expensive QQ factorization only after modular witnesses fail."""

    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for characteristic-zero fallback")
    cancellation_a, cancellation_b = primitive_linear_coefficients(left_x, right_x)
    program = f"""
ring r=0,(c,k,T),dp;
number t1={rational_text(E22_T)}; number t2={rational_text(RANK20_T)};
number x1={rational_text(left_x)}; number x2={rational_text(right_x)};
poly a=(x1*(c*t1+1)-x2*(c*t2+1))/(t1-t2);
poly b=x1*(c*t1+1)-a*t1;
poly D=c*T+1;
poly N=a*T+b+k*(T-t1)*(T-t2);
poly f0=18103855887324900+102302344648*T^2-879500*T^4+T^6;
poly f1=-257843832010380-650709150*T^2+1860*T^4;
poly f2=1195214262641+1718550*T^2-2*T^4;
poly f3=-2051321790-1860*T^2;
poly f4=1149050+T^2;
poly P=f0*D^4+f1*N*D^3+f2*N^2*D^2+f3*N^3*D+f4*N^4;
poly rr=resultant(P,diff(P,T),T);
poly leading_shape=(c-k)^2*(c+k)^2;
poly disc=rr/leading_shape;
poly known=(k-c)^2*(k+c)^2*(28917*c+10)^12*(39508*c+39)^12
           *({cancellation_a}*c+{cancellation_b}*k)^12;
poly residual=disc/known;
if (rr-leading_shape*disc!=0 || disc-known*residual!=0) {{ "@@ERROR_DIVISION"; }}
list L=factorize(residual);
int idx=0;
for (idx=1;idx<=size(L[1]);idx++) {{
  if (deg(L[1][idx])>0) {{
    "@@FACTOR_DEGREE"; deg(L[1][idx]);
    "@@FACTOR_EXPONENT"; L[2][idx];
    "@@FACTOR_BEGIN"; L[1][idx]; "@@FACTOR_END";
  }}
}}
"""
    try:
        completed = subprocess.run(
            [executable, "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=True,
        )
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "timeout_seconds": timeout}
    if completed.stderr.strip():
        raise RuntimeError(f"Singular wrote to stderr: {completed.stderr.strip()}")
    if "@@ERROR_" in completed.stdout or "   ?" in completed.stdout:
        raise AssertionError(f"Singular QQ fallback failed:\n{completed.stdout}")
    degrees = [int(value) for value in re.findall(r"^@@FACTOR_DEGREE\n(-?\d+)$", completed.stdout, re.MULTILINE)]
    exponents = [int(value) for value in re.findall(r"^@@FACTOR_EXPONENT\n(-?\d+)$", completed.stdout, re.MULTILINE)]
    factors = re.findall(
        r"^@@FACTOR_BEGIN\n(.*?)\n@@FACTOR_END$",
        completed.stdout,
        re.MULTILINE | re.DOTALL,
    )
    return {
        "status": "complete",
        "factor_signature": [[degree, exponent] for degree, exponent in zip(degrees, exponents, strict=True)],
        "factor_sha256": [sha256_text(value.strip() + "\n") for value in factors],
        "irreducible_degree_32": len(degrees) == 1 and degrees == [32] and exponents == [1],
        "rational_linear_residual_factor_count": sum(degree == 1 for degree in degrees),
    }


def modular_or_qq_pair(
    work: tuple[tuple[str, Fraction, str, Fraction], tuple[int, ...], int, int]
) -> dict[str, Any]:
    pair, primes, modular_timeout, qq_timeout = work
    left_label, left_x, right_label, right_x = pair
    attempts = []
    witness = None
    for prime in primes:
        attempt = modular_attempt(left_x, right_x, prime, timeout=modular_timeout)
        attempts.append(attempt)
        if attempt.get("irreducible_degree_32_witness"):
            witness = attempt
            break
    fallback = None
    if witness is None:
        fallback = _qq_factorization_fallback(left_x, right_x, timeout=qq_timeout)
    classified = witness is not None or bool(
        fallback and fallback.get("irreducible_degree_32")
    )
    return {
        "left": left_label,
        "right": right_label,
        "modular_attempts": attempts,
        "irreducibility_witness": witness,
        "qq_factorization_fallback": fallback,
        "residual_irreducible_over_QQ": classified,
        "proof": (
            "Gauss lemma: a primitive characteristic-zero factorization would "
            "reduce nontrivially at this good prime because total degree 32 is "
            "preserved; the recorded residual reduction is irreducible."
            if witness is not None
            else (
                "direct characteristic-zero factorization"
                if classified
                else "unresolved within the pinned modular and QQ caps"
            )
        ),
    }


def run(
    root: Path,
    *,
    workers: int = 4,
    modular_timeout: int = 60,
    qq_timeout: int = 600,
) -> dict[str, Any]:
    if workers <= 0:
        raise ValueError("workers must be positive")
    e22_record, rank20_record, pairs = independent_pairs(root)
    population_lines = (
        f"{left}|{rational_text(left_x)}|{right}|{rational_text(right_x)}"
        for left, left_x, right, right_x in pairs
    )
    population_sha256 = sha256_lines(population_lines)

    exact_records = []
    breakthrough = None
    for pair in pairs:
        record = exact_component_classification(pair)
        exact_records.append(record)
        if record["valid_genus_at_most_one_components"]:
            breakthrough = {
                "pair": [record["left"], record["right"]],
                "components": record["valid_genus_at_most_one_components"],
            }
            break

    modular_records: list[dict[str, Any]] = []
    if breakthrough is None:
        work = [
            (pair, MODULAR_PRIMES, modular_timeout, qq_timeout)
            for pair in pairs
        ]
        if workers == 1:
            modular_records = [modular_or_qq_pair(item) for item in work]
        else:
            start_method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"
            with mp.get_context(start_method).Pool(workers) as pool:
                modular_records = pool.map(modular_or_qq_pair, work)
        modular_records.sort(key=lambda item: (item["left"], item["right"]))

    exact_histogram = Counter(
        (
            item["type"],
            item["validity"],
            item["T_degree_after_specialization"],
            tuple(map(tuple, item["factor_signature_over_parameter_function_field"])),
            item["squareclass_kernel_degree"],
            item["squareclass_genus"],
        )
        for record in exact_records
        for item in record["rational_components"]
    )
    unresolved = [
        [record["left"], record["right"]]
        for record in modular_records
        if not record["residual_irreducible_over_QQ"]
    ]
    qq_fallback_pairs = [
        [record["left"], record["right"]]
        for record in modular_records
        if record["qq_factorization_fallback"] is not None
    ]
    witness_prime_histogram = Counter(
        record["irreducibility_witness"]["prime"]
        for record in modular_records
        if record["irreducibility_witness"] is not None
    )
    component_manifest_sha256 = sha256_lines(
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in exact_records
    )
    modular_manifest_sha256 = sha256_lines(
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in modular_records
    )

    completed = (
        breakthrough is None
        and len(exact_records) == 80
        and len(modular_records) == 80
        and not unresolved
    )
    artifact: dict[str, Any] = {
        "schema_version": "elliptic-curves.fermigier-bidegree21-all80.v1",
        "status": (
            "complete exact bounded structural classification"
            if completed
            else (
                "stopped on a valid genus-at-most-one component"
                if breakthrough is not None
                else "bounded run complete with unresolved residuals"
            )
        ),
        "claim_level": (
            "exact computation; no new section or specialization found"
            if completed
            else "partial exact computation; inspect unresolved or breakthrough records"
        ),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "family_id": "fermigier-mestre-v1",
        "anchors": {
            "E22": {
                "canonical_parameter": f"u={rational_text(E22_U)}",
                "alias": f"literal_shift_T={rational_text(E22_T)}",
                "certified_rank_lower_bound": 22,
                "exceptional_quotient_rank_lower_bound": e22_record["exceptional_quotient_rank_lower_bound"],
                "independent_labels": list(E22_INDEPENDENT_EXCEPTIONAL_LABELS),
            },
            "rank20": {
                "canonical_parameter": f"u={rational_text(RANK20_U)}",
                "alias": f"literal_shift_T={rational_text(RANK20_T)}",
                "certified_rank_lower_bound": 20,
                "exceptional_quotient_rank_lower_bound": rank20_record["exceptional_quotient_rank_lower_bound"],
                "independent_labels": [f"R20E{index}" for index in range(1, 9)],
                "source_artifact": rank20_record["source_artifact"],
                "source_artifact_sha256": rank20_record["source_artifact_sha256"],
            },
        },
        "population": {
            "definition": "10 independent E22 quotient representatives times 8 independent rank20 quotient representatives",
            "expected_pair_count": 80,
            "pair_count": len(pairs),
            "manifest_sha256": population_sha256,
            "pairs": [
                {
                    "left": left,
                    "left_x": rational_text(left_x),
                    "right": right,
                    "right_x": rational_text(right_x),
                }
                for left, left_x, right, right_x in pairs
            ],
        },
        "universal_discriminant_factor_proof": {
            "cleared_polynomial_parameter_degree_bound": 4,
            "degree_10_discriminant_coefficient_weight": 18,
            "discriminant_total_degree_upper_bound": 72,
            "resultant_total_degree_upper_bound": 76,
            "known_factor_total_degree": 40,
            "residual_total_degree_upper_bound": 32,
            "known_factors": [
                {"factor": "k-c", "multiplicity": 2, "reason": "two-root degree drop at infinity"},
                {"factor": "k+c", "multiplicity": 2, "reason": "two-root degree drop at infinity"},
                {"factor": "c*T_E22+1", "multiplicity": 12, "reason": "four roots collide linearly at the E22 anchor pole"},
                {"factor": "c*T_rank20+1", "multiplicity": 12, "reason": "four roots collide linearly at the rank20 anchor pole"},
                {"factor": "(left_x-right_x)*c+(T_rank20-T_E22)*k", "multiplicity": 12, "reason": "four roots follow the cancelled denominator linearly"},
            ],
            "exact_checks_per_pair": (
                "The script verifies the leading coefficient, both endpoint "
                "numerator/denominator identities, exact numerator divisibility "
                "on the cancellation locus, and all five function-field squareclasses."
            ),
            "multiplicity_exactness": (
                "The local weighted identities give the displayed lower bounds. "
                "At every recorded good prime, exact division leaves an irreducible "
                "degree-32 quotient; together with the degree-72 upper bound this "
                "proves the multiplicities are exact and the list exhaustive."
            ),
        },
        "rational_component_classification": {
            "completed_pair_count": len(exact_records),
            "histogram": [
                {
                    "type": key[0],
                    "validity": key[1],
                    "T_degree": key[2],
                    "factor_signature": [list(item) for item in key[3]],
                    "squareclass_kernel_degree": key[4],
                    "squareclass_genus": key[5],
                    "count": count,
                }
                for key, count in sorted(exact_histogram.items())
            ],
            "manifest_sha256": component_manifest_sha256,
            "records": exact_records,
        },
        "residual_irreducibility": {
            "method": (
                "Try the pinned good primes in order.  Irreducibility of the "
                "degree-preserving reduction certifies irreducibility over QQ by "
                "Gauss lemma.  Use full QQ factorization only if all modular attempts fail."
            ),
            "pinned_prime_list": list(MODULAR_PRIMES),
            "completed_pair_count": len(modular_records),
            "irreducible_pair_count": sum(record["residual_irreducible_over_QQ"] for record in modular_records),
            "witness_prime_histogram": [
                {"prime": prime, "count": count}
                for prime, count in sorted(witness_prime_histogram.items())
            ],
            "qq_fallback_pairs": qq_fallback_pairs,
            "unresolved_pairs": unresolved,
            "manifest_sha256": modular_manifest_sha256,
            "records": modular_records,
        },
        "breakthrough": breakthrough,
        "scope": {
            "all_80_independent_pairs_classified": completed,
            "finite_denominator_chart": "d=1",
            "infinity_chart": (
                "not repeated here; the complete k=0 projective Mobius chart was "
                "classified separately, while this run exhausts the declared finite "
                "bidegree-(2,1) pencils"
            ),
            "not_claimed": (
                "No rational-point search is made on the irreducible degree-32 "
                "discriminant curves, nor on intersections of their components."
            ),
        },
        "outcome": {
            "rational_discriminant_components_per_pair": 5,
            "valid_rational_components_per_pair": 3,
            "valid_genus_at_most_one_components": sum(
                len(record["valid_genus_at_most_one_components"])
                for record in exact_records
            ),
            "new_base_changes": 0,
            "new_sections": 0,
            "new_specializations": 0,
            "target_met": False,
        },
        "software": {
            "python": platform.python_version(),
            "sympy": sp.__version__,
            "singular": subprocess.run(
                [shutil.which("Singular") or "Singular", "--version"],
                text=True,
                capture_output=True,
                check=True,
            ).stdout.splitlines()[0],
        },
        "bounds": {
            "workers": workers,
            "modular_timeout_seconds_per_attempt": modular_timeout,
            "qq_timeout_seconds_per_fallback_pair": qq_timeout,
            "modular_prime_count": len(MODULAR_PRIMES),
        },
        "reproducing_command": (
            "PYTHONPATH=elliptic-curves:elliptic-curves/cas .venv/bin/python "
            "elliptic-curves/cas/analyze_fermigier_bidegree21_all80.py --workers 4"
        ),
    }
    stable = dict(artifact)
    stable.pop("generated_at_utc")
    artifact["result_sha256"] = sha256_text(
        json.dumps(stable, sort_keys=True, separators=(",", ":"))
    )
    return artifact


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--modular-timeout", type=int, default=60)
    parser.add_argument("--qq-timeout", type=int, default=600)
    parser.add_argument("--output", type=Path, default=root / OUTPUT_RELATIVE)
    args = parser.parse_args()
    artifact = run(
        root,
        workers=args.workers,
        modular_timeout=args.modular_timeout,
        qq_timeout=args.qq_timeout,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["outcome"], sort_keys=True))
    print(f"result_sha256={artifact['result_sha256']}")


if __name__ == "__main__":
    main()
