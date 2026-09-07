
from sage.all import *
from pathlib import Path
import json
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

PROTOCOL = "ANONRANK2"


def ainvs_int(E):
    return [ZZ(x) for x in E.ainvs()]


def factor_product(factors):
    x = ZZ(1)
    for p, e in factors:
        p = ZZ(p)
        e = ZZ(e)
        assert p.is_prime(), f"not prime: {p}"
        x *= p**e
    return x


def radical(factors):
    x = ZZ(1)
    for p, _ in factors:
        x *= ZZ(p)
    return x


def trace_features(E, bounds, bad):
    """
    Produce several partial-sum fingerprints.

    legacy_all:
        historical repo-style sum over all primes.

    legacy_good:
        same expression restricted to good reduction.

    S0:
        current repo convention:
          (1/log B) sum_good a_p log(p)/p

        Hence strongly NEGATIVE is the high-rank direction.

    S5:
        secondary Mestre-Nagao-style feature retained from v1.

    IMPORTANT:
        snapshots happen BEFORE processing the first prime > B,
        so the B cutoff is exact.
    """

    bounds = sorted(bounds)
    maxB = max(bounds)

    legacy_all = RR(0)
    legacy_good = RR(0)
    s0_raw = RR(0)
    s5 = RR(0)

    results = {}
    next_bound = 0

    def snapshot(B):
        results[B] = {
            "legacy_all": RR(legacy_all),
            "legacy_good": RR(legacy_good),
            "S0": RR(s0_raw / log(RR(B))),
            "S5": RR(s5),
        }

    for p in prime_range(2, maxB + 1):
        p = ZZ(p)

        # Snapshot all bounds strictly below this prime.
        while next_bound < len(bounds) and bounds[next_bound] < p:
            snapshot(bounds[next_bound])
            next_bound += 1

        ap = ZZ(E.ap(p))

        term = (
            RR(2 - ap)
            / RR(p + 1 - ap)
            * log(RR(p))
        )

        legacy_all += term

        if p not in bad:
            legacy_good += term
            s0_raw += RR(ap) * log(RR(p)) / RR(p)
            s5 += log(RR(p + 1 - ap) / RR(p))
        else:
            # Preserve the v1 convention separately from legacy_good.
            s5 += log(RR(1.5) * RR(p - 1) / RR(p))

    while next_bound < len(bounds):
        snapshot(bounds[next_bound])
        next_bound += 1

    return results


def print_features(label, features):
    for B in sorted(features):
        f = features[B]
        print(
            f"{PROTOCOL}|stage=trace"
            f"|curve={label}"
            f"|B={B}"
            f"|legacy_all={f['legacy_all']:.12f}"
            f"|legacy_good={f['legacy_good']:.12f}"
            f"|S0={f['S0']:.12f}"
            f"|S5={f['S5']:.12f}",
            flush=True,
        )


print(f"{PROTOCOL}|stage=start", flush=True)

E0 = EllipticCurve(QQ, GENERAL_WEIERSTRASS_COEFFICIENTS)
E29 = EllipticCurve(QQ, E29_COEFFICIENTS)

E = E0.global_minimal_model()

minimal_same = ainvs_int(E0) == ainvs_int(E)

print(
    f"{PROTOCOL}|stage=minimal"
    f"|same={int(minimal_same)}"
    f"|ainvs={ainvs_int(E)}",
    flush=True,
)

delta = ZZ(E.discriminant())
c4 = ZZ(E.c4())
c6 = ZZ(E.c6())

assert factor_product(DISCRIMINANT_FACTORS) == abs(delta)

bad = {ZZ(p) for p, _ in DISCRIMINANT_FACTORS}

print(
    f"{PROTOCOL}|stage=invariants"
    f"|c4={c4}"
    f"|c6={c6}"
    f"|delta={delta}"
    f"|logDelta={RR(log(abs(delta))):.15f}",
    flush=True,
)

g = gcd(abs(c4), abs(delta))
assert g == 1

N_expected = radical(DISCRIMINANT_FACTORS)
N_sage = ZZ(E.conductor())

assert N_expected == N_sage

print(
    f"{PROTOCOL}|stage=conductor"
    f"|semistable=1"
    f"|N={N_sage}"
    f"|logN={RR(log(N_sage)):.15f}",
    flush=True,
)

root_number = ZZ(E.root_number())

print(
    f"{PROTOCOL}|stage=root_number"
    f"|value={root_number}",
    flush=True,
)

torsion = E.torsion_subgroup()

print(
    f"{PROTOCOL}|stage=torsion"
    f"|order={torsion.order()}"
    f"|invariants={torsion.invariants()}",
    flush=True,
)

# ------------------------------------------------------------
# E29 exact comparison data
# ------------------------------------------------------------

E29 = E29.global_minimal_model()

delta29 = ZZ(E29.discriminant())
c429 = ZZ(E29.c4())

# E29 is also semistable. Factoring Delta gives its bad-prime support.
fac29 = factor(abs(delta29))
bad29 = {ZZ(p) for p, _ in fac29}

assert gcd(abs(c429), abs(delta29)) == 1

N29 = ZZ(E29.conductor())

print(
    f"{PROTOCOL}|stage=e29"
    f"|N={N29}"
    f"|logN={RR(log(N29)):.15f}"
    f"|root_number={E29.root_number()}",
    flush=True,
)

ratio = QQ(N29, N_sage)

print(
    f"{PROTOCOL}|stage=conductor_compare"
    f"|ratio={ratio}"
    f"|ratio_RR={RR(ratio):.15f}"
    f"|log_ratio={RR(log(RR(ratio))):.15f}"
    f"|ratio_over_79={RR(ratio / 79):.15f}",
    flush=True,
)

# ------------------------------------------------------------
# Trace fingerprints
# ------------------------------------------------------------

bounds = [
    250,
    500,
    1000,
    2000,
    5000,
    10000,
    20000,
    50000,
    100000,
]

t0 = time.monotonic()

anon_features = trace_features(E, bounds, bad)

print(
    f"{PROTOCOL}|stage=anon_trace_done"
    f"|seconds={time.monotonic()-t0:.3f}",
    flush=True,
)

t1 = time.monotonic()

e29_features = trace_features(E29, bounds, bad29)

print(
    f"{PROTOCOL}|stage=e29_trace_done"
    f"|seconds={time.monotonic()-t1:.3f}",
    flush=True,
)

print_features("anonymous", anon_features)
print_features("E29", e29_features)

# ------------------------------------------------------------
# Cross-curve deltas
# ------------------------------------------------------------

for B in bounds:
    a = anon_features[B]
    b = e29_features[B]

    print(
        f"{PROTOCOL}|stage=compare"
        f"|B={B}"
        f"|legacy_good_delta={a['legacy_good'] - b['legacy_good']:.12f}"
        f"|legacy_all_delta={a['legacy_all'] - b['legacy_all']:.12f}"
        f"|S0_delta={a['S0'] - b['S0']:.12f}"
        f"|S5_delta={a['S5'] - b['S5']:.12f}",
        flush=True,
    )

# ------------------------------------------------------------
# Persistence
# ------------------------------------------------------------

data = {
    "ainvs": [str(x) for x in ainvs_int(E)],
    "minimal_same": bool(minimal_same),
    "discriminant": str(delta),
    "conductor": str(N_sage),
    "log_conductor": float(log(N_sage)),
    "root_number": int(root_number),
    "torsion_order": int(torsion.order()),
    "e29_conductor": str(N29),
    "conductor_ratio": str(ratio),
    "bounds": bounds,
    "anonymous": {
        str(B): {
            k: float(v)
            for k, v in anon_features[B].items()
        }
        for B in bounds
    },
    "E29": {
        str(B): {
            k: float(v)
            for k, v in e29_features[B].items()
        }
        for B in bounds
    },
}

out = OUT / "anonymous-high-rank-candidate-profile-v2.json"
out.write_text(json.dumps(data, indent=2) + "\n")

print(
    f"{PROTOCOL}|stage=done|status=PASS|out={out}",
    flush=True,
)
