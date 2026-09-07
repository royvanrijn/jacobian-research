
from sage.all import *
from pathlib import Path
import json
import math
import sys
import time

sys.path.insert(0, "elliptic-curves/cas")

from anonymous_rank_candidate import (
    GENERAL_WEIERSTRASS_COEFFICIENTS,
    DISCRIMINANT_FACTORS,
    E29_COEFFICIENTS,
)

OUT = Path("artifacts/local/elliptic-curves")
OUT.mkdir(parents=True, exist_ok=True)


def integer_ainvs(E):
    return [ZZ(v) for v in E.ainvs()]


def verify_discriminant_factorization(delta):
    product_value = ZZ(1)

    for p, e in DISCRIMINANT_FACTORS:
        p = ZZ(p)
        e = ZZ(e)

        assert p.is_prime(), f"factor is not prime: {p}"
        product_value *= p ** e

    assert product_value == abs(ZZ(delta)), (
        "supplied factorization does not equal discriminant"
    )

    return product_value


def radical(factors):
    result = ZZ(1)

    for p, _ in factors:
        result *= ZZ(p)

    return result


def bad_prime_set(factors):
    return {ZZ(p) for p, _ in factors}


def trace_features(E, bounds, known_bad=None):
    """
    Compute several complementary Mestre--Nagao fingerprints.

    legacy:
        The same basic quantity used by the existing repo search:
            sum (2-a_p)/(p+1-a_p) log p

    S0:
        1/log(B) * sum_{p<=B, good} a_p log(p)/p

    S5:
        sum_good log((p+1-a_p)/p)
        + sum_bad log(1.5*(p-1)/p)

    We retain the full vector over many B rather than reducing it to one
    scalar.  This is important for held-out testing and avoids optimizing
    against a single cutoff.
    """

    max_bound = max(bounds)
    bad = set(known_bad or [])

    primes = list(prime_range(2, max_bound + 1))

    legacy = RR(0)
    s0_raw = RR(0)
    s5 = RR(0)

    result = {}
    next_bound_index = 0

    t0 = time.monotonic()

    for p in primes:
        p = ZZ(p)
        ap = ZZ(E.ap(p))

        legacy += RR(2 - ap) / RR(p + 1 - ap) * log(RR(p))

        if p not in bad:
            s0_raw += RR(ap) * log(RR(p)) / RR(p)
            s5 += log(RR(p + 1 - ap) / RR(p))
        else:
            s5 += log(RR(1.5) * RR(p - 1) / RR(p))

        while (
            next_bound_index < len(bounds)
            and p >= bounds[next_bound_index]
        ):
            B = bounds[next_bound_index]

            result[B] = {
                "legacy": legacy,
                "S0": s0_raw / log(RR(B)),
                "S5": s5,
            }

            next_bound_index += 1

    # prime_range may end just below B.
    while next_bound_index < len(bounds):
        B = bounds[next_bound_index]

        result[B] = {
            "legacy": legacy,
            "S0": s0_raw / log(RR(B)),
            "S5": s5,
        }

        next_bound_index += 1

    return result, time.monotonic() - t0


def print_features(label, features):
    for B, values in features.items():
        print(
            "ANONRANK"
            f"|stage=trace"
            f"|curve={label}"
            f"|B={B}"
            f"|legacy={values['legacy']:.12f}"
            f"|S0={values['S0']:.12f}"
            f"|S5={values['S5']:.12f}",
            flush=True,
        )


print("ANONRANK|stage=start", flush=True)

E = EllipticCurve(QQ, GENERAL_WEIERSTRASS_COEFFICIENTS)
E29 = EllipticCurve(QQ, E29_COEFFICIENTS)

print(
    "ANONRANK|stage=input"
    f"|ainvs={integer_ainvs(E)}",
    flush=True,
)


# ------------------------------------------------------------
# Minimal-model check
# ------------------------------------------------------------

Emin = E.global_minimal_model()

minimal_same = integer_ainvs(Emin) == integer_ainvs(E)

print(
    "ANONRANK|stage=minimal"
    f"|same={int(minimal_same)}"
    f"|minimal_ainvs={integer_ainvs(Emin)}",
    flush=True,
)


# Work with the actual minimal model from here onward.
E = Emin

delta = ZZ(E.discriminant())
c4 = ZZ(E.c4())
c6 = ZZ(E.c6())

print(
    "ANONRANK|stage=invariants"
    f"|c4={c4}"
    f"|c6={c6}"
    f"|delta={delta}"
    f"|log_abs_delta={RR(log(abs(delta))):.15f}",
    flush=True,
)


# ------------------------------------------------------------
# Exact supplied discriminant factorization
# ------------------------------------------------------------

verify_discriminant_factorization(delta)

bad = bad_prime_set(DISCRIMINANT_FACTORS)

print(
    "ANONRANK|stage=disc_factor"
    f"|prime_count={len(DISCRIMINANT_FACTORS)}"
    f"|factors="
    + ",".join(f"{p}^{e}" for p, e in DISCRIMINANT_FACTORS),
    flush=True,
)


# ------------------------------------------------------------
# Semistability / conductor
# ------------------------------------------------------------

g = gcd(abs(c4), abs(delta))

print(
    "ANONRANK|stage=semistable"
    f"|gcd_c4_delta={g}",
    flush=True,
)

assert g == 1, (
    "gcd(c4,Delta) != 1; cannot infer all bad fibers are multiplicative"
)

N_from_delta = radical(DISCRIMINANT_FACTORS)

print(
    "ANONRANK|stage=conductor_from_semistable"
    f"|N={N_from_delta}"
    f"|logN={RR(log(N_from_delta)):.15f}",
    flush=True,
)


# Ask Sage independently as well.
try:
    sage_N = ZZ(E.conductor())

    print(
        "ANONRANK|stage=conductor_sage"
        f"|N={sage_N}"
        f"|matches={int(sage_N == N_from_delta)}",
        flush=True,
    )

    assert sage_N == N_from_delta

except Exception as exc:
    print(
        "ANONRANK|stage=conductor_sage"
        f"|status=ERROR"
        f"|error={type(exc).__name__}:{exc}",
        flush=True,
    )


# ------------------------------------------------------------
# Root number
# ------------------------------------------------------------

try:
    root_number = ZZ(E.root_number())

    print(
        "ANONRANK|stage=root_number"
        f"|value={root_number}",
        flush=True,
    )

except Exception as exc:
    root_number = None

    print(
        "ANONRANK|stage=root_number"
        f"|status=ERROR"
        f"|error={type(exc).__name__}:{exc}",
        flush=True,
    )


# ------------------------------------------------------------
# Torsion
# ------------------------------------------------------------

try:
    torsion = E.torsion_subgroup()

    print(
        "ANONRANK|stage=torsion"
        f"|order={torsion.order()}"
        f"|invariants={torsion.invariants()}",
        flush=True,
    )

except Exception as exc:
    torsion = None

    print(
        "ANONRANK|stage=torsion"
        f"|status=ERROR"
        f"|error={type(exc).__name__}:{exc}",
        flush=True,
    )


# ------------------------------------------------------------
# Multi-window trace fingerprints
# ------------------------------------------------------------

bounds = [1000, 2000, 5000, 10000, 20000]

anonymous_features, anonymous_seconds = trace_features(
    E,
    bounds,
    known_bad=bad,
)

print_features("anonymous", anonymous_features)

# For E29 let Sage decide good/bad by Delta divisibility.
e29_delta = ZZ(E29.discriminant())
e29_bad = {
    ZZ(p)
    for p, _ in factor(abs(e29_delta))
}

e29_features, e29_seconds = trace_features(
    E29,
    bounds,
    known_bad=e29_bad,
)

print_features("E29", e29_features)

for B in bounds:
    a = anonymous_features[B]
    b = e29_features[B]

    print(
        "ANONRANK"
        f"|stage=trace_compare"
        f"|B={B}"
        f"|legacy_delta={a['legacy'] - b['legacy']:.12f}"
        f"|S0_delta={a['S0'] - b['S0']:.12f}"
        f"|S5_delta={a['S5'] - b['S5']:.12f}",
        flush=True,
    )


# ------------------------------------------------------------
# Machine-readable checkpoint
# ------------------------------------------------------------

result = {
    "curve": {
        "ainvs": [str(v) for v in integer_ainvs(E)],
        "minimal_same_as_input": bool(minimal_same),
        "c4": str(c4),
        "c6": str(c6),
        "discriminant": str(delta),
        "log_abs_discriminant": float(log(abs(delta))),
        "semistable": g == 1,
        "conductor": str(N_from_delta),
        "log_conductor": float(log(N_from_delta)),
        "root_number": None if root_number is None else int(root_number),
    },
    "trace_bounds": bounds,
    "trace_features": {
        str(B): {
            key: float(value)
            for key, value in anonymous_features[B].items()
        }
        for B in bounds
    },
    "e29_features": {
        str(B): {
            key: float(value)
            for key, value in e29_features[B].items()
        }
        for B in bounds
    },
}

output = OUT / "anonymous-high-rank-candidate-profile.json"
output.write_text(json.dumps(result, indent=2) + "\n")

print(
    "ANONRANK|stage=done"
    f"|status=PASS"
    f"|out={output}"
    f"|anonymous_trace_seconds={anonymous_seconds:.3f}"
    f"|e29_trace_seconds={e29_seconds:.3f}",
    flush=True,
)
