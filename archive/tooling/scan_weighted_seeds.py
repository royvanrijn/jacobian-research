#!/usr/bin/env python3
"""Scan canonical and small deformed weighted seeds for uniform geometry.

The scan is exploratory.  It proves the recorded one-variable identities for
each candidate, but finite-field histograms are diagnostics rather than a
characteristic-zero monodromy theorem.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.weighted import (  # noqa: E402
    WeightedSeedModel,
    canonical_seed,
    deformation_basis,
    w,
)


def rational_text(value):
    return str(sp.factor(value))


def coefficients_mod(poly, prime):
    coefficients = {}
    for (exponent,), coefficient in sp.Poly(poly, w).terms():
        numerator, denominator = map(int, sp.fraction(coefficient))
        if denominator % prime == 0:
            return None
        reduced = numerator * pow(denominator, -1, prime) % prime
        if reduced:
            coefficients[exponent] = reduced
    return coefficients


def evaluate_mod(coefficients, value, prime):
    return sum(coefficient*pow(value, exponent, prime)
               for exponent, coefficient in coefficients.items()) % prime


def simple_root_histogram(model, prime):
    primitive = coefficients_mod(model.primitive, prime)
    seed = coefficients_mod(model.seed, prime)
    if not primitive or not seed:
        return None
    if max(primitive) != model.fiber_degree:
        return None
    if max(seed) != model.seed_degree:
        return None
    histogram = Counter()
    for slope in range(prime):
        for intercept in range(prime):
            simple = 0
            for root in range(prime):
                value = (evaluate_mod(primitive, root, prime)
                         - slope*root + intercept) % prime
                derivative = (evaluate_mod(seed, root, prime) - slope) % prime
                if value == 0 and derivative != 0:
                    simple += 1
            histogram[simple] += 1
    return dict(sorted(histogram.items()))


def candidate_record(label, degree, parameters, seed, primes):
    try:
        model = WeightedSeedModel(seed)
    except AssertionError:
        return None
    m0, m1, extra = model.zero_profile()
    extra_degree = sp.Poly(extra, w).degree()
    critical = sp.diff(model.seed, w)
    critical_squarefree = sp.gcd(
        sp.Poly(critical, w), sp.Poly(sp.diff(critical, w), w)
    ).degree() == 0
    rational_extra_roots = {
        str(root): multiplicity
        for root, multiplicity in sp.roots(extra, w).items()
        if root.is_Rational
    }
    histograms = {}
    for prime in primes:
        histogram = simple_root_histogram(model, prime)
        histograms[str(prime)] = histogram
    distinguished_only = extra_degree == 0
    score = (
        100*int(distinguished_only)
        + 20*int(critical_squarefree)
        - 5*extra_degree
        - len(rational_extra_roots)
    )
    return {
        "label": label,
        "requested_seed_degree": degree,
        "parameters": parameters,
        "seed": rational_text(model.seed),
        "primitive": rational_text(model.primitive),
        "seed_degree": model.seed_degree,
        "fiber_degree": model.fiber_degree,
        "kappa": rational_text(model.kappa),
        "a": rational_text(model.a),
        "zero_multiplicity": m0,
        "one_multiplicity": m1,
        "extra_factor": rational_text(extra),
        "extra_degree": extra_degree,
        "distinguished_zeros_only": distinguished_only,
        "rational_extra_roots": rational_extra_roots,
        "critical_squarefree": critical_squarefree,
        "simple_root_histograms": histograms,
        "score": score,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-degree", type=int, default=6)
    parser.add_argument("--coefficient-bound", type=int, default=2)
    parser.add_argument("--primes", default="5,7,11")
    parser.add_argument("--output", default="artifacts/generated-results/weighted_seed_scan.json")
    args = parser.parse_args()
    assert 2 <= args.max_degree <= 8
    primes = tuple(int(value) for value in args.primes.split(",") if value)
    values = range(-args.coefficient_bound, args.coefficient_bound+1)

    records = []
    seen = set()
    base = 2*w - 3*w**2

    def add(label, degree, parameters, seed):
        key = sp.srepr(sp.expand(seed))
        if key in seen:
            return
        seen.add(key)
        record = candidate_record(label, degree, parameters, seed, primes)
        if record is not None:
            records.append(record)

    for degree in range(2, args.max_degree+1):
        add(f"canonical_H_{degree}", degree, {}, canonical_seed(degree))
        if degree == 2:
            continue
        indices = tuple(range(1, degree-1))
        for coefficients in itertools.product(values, repeat=len(indices)):
            if coefficients[-1] == 0:
                continue
            seed = base + sum(
                coefficient*deformation_basis(index)
                for index, coefficient in zip(indices, coefficients)
            )
            parameters = {f"theta_{index}": coefficient
                          for index, coefficient in zip(indices, coefficients)}
            add("deformation", degree, parameters, seed)

    records.sort(key=lambda record: (
        -record["score"], record["fiber_degree"], record["seed"]
    ))
    summary = {
        "candidate_count": len(records),
        "canonical_count": sum(record["label"].startswith("canonical") for record in records),
        "distinguished_only_count": sum(record["distinguished_zeros_only"] for record in records),
        "squarefree_critical_count": sum(record["critical_squarefree"] for record in records),
        "by_seed_degree": dict(sorted(Counter(
            record["seed_degree"] for record in records
        ).items())),
    }
    output = {
        "scope": {
            "max_seed_degree": args.max_degree,
            "coefficient_bound": args.coefficient_bound,
            "primes": primes,
            "warning": "Finite-field histograms are diagnostics, not a monodromy proof.",
        },
        "summary": summary,
        "top_candidates": records[:30],
        "records": records,
    }
    destination = ROOT / args.output
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("top candidates:")
    for record in records[:12]:
        print(
            record["label"],
            "seed_degree", record["seed_degree"],
            "fiber_degree", record["fiber_degree"],
            "extra_degree", record["extra_degree"],
            "critical_squarefree", record["critical_squarefree"],
            "seed", record["seed"],
        )
    print("wrote", destination)


if __name__ == "__main__":
    main()
