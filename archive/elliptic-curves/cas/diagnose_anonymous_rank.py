#!/usr/bin/env python3
"""
Long-range rank diagnostics for the anonymous elliptic curve.

This performs TWO distinct calculations.

1. Delta=2 Bober-calibrated explicit formula.

   We use the same difference method already used elsewhere in this repo:

       bound(candidate)
         = 21.70
         + (log N_candidate - log N_E20)/(2*pi*Delta)
         - (PS_candidate - PS_E20)/(pi*Delta)

   with Delta=2 and prime support

       floor(exp(2*pi*Delta)) = 286751.

   The value is an analytic-rank upper diagnostic conditional on GRH.
   It is NOT an unconditional algebraic rank bound.

2. Long-range Mestre/Nagao fingerprints for anonymous vs E29.

   These can optionally continue through 10^7.

No point search or descent is performed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import mpmath

sys.path.insert(0, str(Path(__file__).resolve().parent))

from anonymous_rank_candidate import (
    GENERAL_WEIERSTRASS_COEFFICIENTS,
)

from elkies_klagsbrun_rank29 import (
    GENERAL_WEIERSTRASS_COEFFICIENTS as E29_MODEL,
)


Q = Fraction

PROTOCOL = "ANONDIAG"

DELTA = 2
EF_LIMIT = 286_751

# Published Bober calibration used by the existing repository diagnostic.
BOBER_E20_BOUND = mpmath.mpf("21.70")

BOBER_E20_MODEL = (
    Q(1),
    Q(0),
    Q(0),
    Q(-431092980766333677958362095891166),
    Q(5156283555366643659035652799871176909391533088196),
)

DEFAULT_BOUNDS = (
    10_000,
    20_000,
    50_000,
    100_000,
    200_000,
    500_000,
    1_000_000,
)


def gp_rat(value):
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def gp_model(coefficients):
    return "[" + ",".join(
        gp_rat(value)
        for value in coefficients
    ) + "]"


def build_gp_program(bounds):
    max_bound = max(max(bounds), EF_LIMIT)

    bounds_gp = ",".join(str(x) for x in bounds)

    anon = gp_model(GENERAL_WEIERSTRASS_COEFFICIENTS)
    e29 = gp_model(E29_MODEL)
    e20 = gp_model(BOBER_E20_MODEL)

    return f"""
default(realprecision,80);

D={DELTA};
EFLIM={EF_LIMIT};
MAXB={max_bound};
BDS=[{bounds_gp}];

EA=ellminimalmodel(ellinit({anon}));
E29=ellminimalmodel(ellinit({e29}));
E20=ellminimalmodel(ellinit({e20}));

NA=ellglobalred(EA)[1];
N29=ellglobalred(E29)[1];
N20=ellglobalred(E20)[1];

print("META|anonymous|",NA,"|",log(NA));
print("META|E29|",N29,"|",log(N29));
print("META|E20|",N20,"|",log(N20));

efterm(E,p,a,D)=
{{
    my(K,ans,sk,s0,s1,k,lp);

    lp=log(p);
    K=floor(2*Pi*D/lp);
    ans=0;

    if(E.disc % p == 0,
        sk=1;
        for(k=1,K,
            sk*=a;
            ans += lp*sk/p^k*(1-k*lp/(2*Pi*D));
        );
    ,
        s0=2;
        s1=a;

        for(k=1,K,
            if(k==1,
                sk=s1;
            ,
                sk=a*s1-p*s0;
                s0=s1;
                s1=sk;
            );

            ans += lp*sk/p^k*(1-k*lp/(2*Pi*D));
        );
    );

    return(ans);
}};

LAall=0;
LAgood=0;
S0A=0;
S5A=0;

L29all=0;
L29good=0;
S029=0;
S529=0;

PSA=0;
PS29=0;
PS20=0;

bi=1;
nb=#BDS;

forprime(p=2,MAXB,

    while(bi<=nb && p>BDS[bi],
        bb=BDS[bi];

        print(
            "TRACE|anonymous|",bb,"|",
            LAall,"|",
            LAgood,"|",
            S0A/log(bb),"|",
            S5A
        );

        print(
            "TRACE|E29|",bb,"|",
            L29all,"|",
            L29good,"|",
            S029/log(bb),"|",
            S529
        );

        bi++;
    );

    aa=ellap(EA,p);
    a29=ellap(E29,p);

    ta=(2-aa)/(p+1-aa)*log(p);
    t29=(2-a29)/(p+1-a29)*log(p);

    LAall += ta;
    L29all += t29;

    if(EA.disc % p != 0,
        LAgood += ta;
        S0A += aa*log(p)/p;
        S5A += log((p+1-aa)/p);
    ,
        S5A += log((3/2)*(p-1)/p);
    );

    if(E29.disc % p != 0,
        L29good += t29;
        S029 += a29*log(p)/p;
        S529 += log((p+1-a29)/p);
    ,
        S529 += log((3/2)*(p-1)/p);
    );

    if(p<=EFLIM,
        a20=ellap(E20,p);

        PSA += efterm(EA,p,aa,D);
        PS29 += efterm(E29,p,a29,D);
        PS20 += efterm(E20,p,a20,D);
    );
);

while(bi<=nb,
    bb=BDS[bi];

    print(
        "TRACE|anonymous|",bb,"|",
        LAall,"|",
        LAgood,"|",
        S0A/log(bb),"|",
        S5A
    );

    print(
        "TRACE|E29|",bb,"|",
        L29all,"|",
        L29good,"|",
        S029/log(bb),"|",
        S529
    );

    bi++;
);

print(
    "EF|anonymous|",
    NA,"|",
    log(NA),"|",
    PSA
);

print(
    "EF|E29|",
    N29,"|",
    log(N29),"|",
    PS29
);

print(
    "EF|E20|",
    N20,"|",
    log(N20),"|",
    PS20
);

quit;
"""


def parse_output(text):
    meta = {}
    traces = {
        "anonymous": {},
        "E29": {},
    }
    ef = {}

    for raw in text.splitlines():
        line = raw.strip()

        if line.startswith("META|"):
            _, label, conductor, logn = line.split("|")
            meta[label] = {
                "conductor": conductor,
                "logN": logn,
            }

        elif line.startswith("TRACE|"):
            (
                _,
                label,
                bound,
                legacy_all,
                legacy_good,
                s0,
                s5,
            ) = line.split("|")

            traces[label][int(bound)] = {
                "legacy_all": legacy_all,
                "legacy_good": legacy_good,
                "S0": s0,
                "S5": s5,
            }

        elif line.startswith("EF|"):
            _, label, conductor, logn, prime_sum = line.split("|")

            ef[label] = {
                "conductor": conductor,
                "logN": logn,
                "prime_sum": prime_sum,
            }

    return meta, traces, ef


def calibrated_bound(ef, label):
    mpmath.mp.dps = 80

    ref = ef["E20"]
    candidate = ef[label]

    log_term = (
        mpmath.mpf(candidate["logN"])
        - mpmath.mpf(ref["logN"])
    ) / (
        2 * mpmath.pi * DELTA
    )

    prime_term = -(
        mpmath.mpf(candidate["prime_sum"])
        - mpmath.mpf(ref["prime_sum"])
    ) / (
        mpmath.pi * DELTA
    )

    difference = log_term + prime_term

    upper = (
        BOBER_E20_BOUND
        + difference
    )

    return {
        "log_term": log_term,
        "prime_term": prime_term,
        "difference": difference,
        "upper": upper,
    }


def parity_ceiling(value, parity):
    """
    Largest nonnegative integer <= value
    having parity 0 (even) or 1 (odd).
    """
    n = math.floor(float(value))

    while n >= 0 and n % 2 != parity:
        n -= 1

    return n


def nstr(x, digits=20):
    return mpmath.nstr(x, digits)


def run_gp(bounds, timeout, stack_bytes):
    program = build_gp_program(bounds)

    started = time.monotonic()

    proc = subprocess.run(
        [
            "gp",
            "-q",
            "-s",
            str(stack_bytes),
        ],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout,
    )

    elapsed = time.monotonic() - started

    if proc.returncode != 0:
        raise RuntimeError(
            "PARI failed:\n"
            + proc.stderr[-4000:]
        )

    if "***" in proc.stderr:
        raise RuntimeError(
            "PARI error:\n"
            + proc.stderr[-4000:]
        )

    return proc.stdout, elapsed


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--max-bound",
        type=int,
        default=1_000_000,
        help=(
            "maximum long-range trace cutoff; "
            "try 10000000 for the deep pass"
        ),
    )

    ap.add_argument(
        "--timeout",
        type=float,
        default=7200,
    )

    ap.add_argument(
        "--stack-bytes",
        type=int,
        default=1_000_000_000,
    )

    ap.add_argument(
        "--output",
        type=Path,
        default=Path(
            "artifacts/local/elliptic-curves/"
            "anonymous-rank-diagnostic.json"
        ),
    )

    args = ap.parse_args()

    base_bounds = list(DEFAULT_BOUNDS)

    if args.max_bound > base_bounds[-1]:
        x = 2_000_000

        while x < args.max_bound:
            base_bounds.append(x)

            # 2m -> 5m -> 10m -> 20m ...
            if str(x).startswith("2"):
                x = x // 2 * 5
            else:
                x *= 2

        base_bounds.append(args.max_bound)

    bounds = sorted(
        set(
            b
            for b in base_bounds
            if b <= args.max_bound
        )
    )

    # Delta=2 explicit formula needs this much support regardless
    # of requested trace maximum.
    if args.max_bound < EF_LIMIT:
        bounds.append(EF_LIMIT)

    bounds = sorted(set(bounds))

    print(
        f"{PROTOCOL}|stage=start"
        f"|trace_max={args.max_bound}"
        f"|ef_delta={DELTA}"
        f"|ef_limit={EF_LIMIT}"
        f"|bounds={bounds}",
        flush=True,
    )

    stdout, elapsed = run_gp(
        bounds,
        args.timeout,
        args.stack_bytes,
    )

    meta, traces, ef = parse_output(stdout)

    required = {
        "anonymous",
        "E29",
        "E20",
    }

    if set(ef) != required:
        print(stdout)
        raise RuntimeError(
            f"missing EF records: got {set(ef)}"
        )

    anon_bound = calibrated_bound(
        ef,
        "anonymous",
    )

    e29_bound = calibrated_bound(
        ef,
        "E29",
    )

    print(
        f"{PROTOCOL}|stage=explicit_formula"
        f"|curve=anonymous"
        f"|Delta={DELTA}"
        f"|prime_limit={EF_LIMIT}"
        f"|log_term={nstr(anon_bound['log_term'])}"
        f"|prime_term={nstr(anon_bound['prime_term'])}"
        f"|upper={nstr(anon_bound['upper'],30)}",
        flush=True,
    )

    print(
        f"{PROTOCOL}|stage=explicit_formula"
        f"|curve=E29"
        f"|Delta={DELTA}"
        f"|prime_limit={EF_LIMIT}"
        f"|log_term={nstr(e29_bound['log_term'])}"
        f"|prime_term={nstr(e29_bound['prime_term'])}"
        f"|upper={nstr(e29_bound['upper'],30)}",
        flush=True,
    )

    print(
        f"{PROTOCOL}|stage=explicit_compare"
        f"|anonymous_minus_E29="
        f"{nstr(anon_bound['upper'] - e29_bound['upper'],30)}",
        flush=True,
    )

    # Analytic functional-equation parity:
    #
    # anonymous root number +1 -> even analytic order
    # E29 root number -1       -> odd analytic order
    anon_parity_max = parity_ceiling(
        anon_bound["upper"],
        0,
    )

    e29_parity_max = parity_ceiling(
        e29_bound["upper"],
        1,
    )

    print(
        f"{PROTOCOL}|stage=grh_parity"
        f"|curve=anonymous"
        f"|root_number=+1"
        f"|analytic_rank_parity=even"
        f"|largest_even_below_bound={anon_parity_max}",
        flush=True,
    )

    print(
        f"{PROTOCOL}|stage=grh_parity"
        f"|curve=E29"
        f"|root_number=-1"
        f"|analytic_rank_parity=odd"
        f"|largest_odd_below_bound={e29_parity_max}",
        flush=True,
    )

    for B in sorted(
        set(traces["anonymous"])
        & set(traces["E29"])
    ):
        # Don't clutter output with EF-only 286751 unless user
        # explicitly requested that trace cutoff.
        if B == EF_LIMIT and args.max_bound < EF_LIMIT:
            continue

        a = traces["anonymous"][B]
        e = traces["E29"][B]

        la = mpmath.mpf(a["legacy_good"])
        le = mpmath.mpf(e["legacy_good"])

        s0a = mpmath.mpf(a["S0"])
        s0e = mpmath.mpf(e["S0"])

        s5a = mpmath.mpf(a["S5"])
        s5e = mpmath.mpf(e["S5"])

        print(
            f"{PROTOCOL}|stage=trace"
            f"|B={B}"
            f"|anon_legacy_good={nstr(la,16)}"
            f"|e29_legacy_good={nstr(le,16)}"
            f"|legacy_delta={nstr(la-le,16)}"
            f"|anon_S0={nstr(s0a,16)}"
            f"|e29_S0={nstr(s0e,16)}"
            f"|S0_delta={nstr(s0a-s0e,16)}"
            f"|anon_S5={nstr(s5a,16)}"
            f"|e29_S5={nstr(s5e,16)}"
            f"|S5_delta={nstr(s5a-s5e,16)}",
            flush=True,
        )

    artifact = {
        "status": "complete",
        "method": {
            "explicit_formula_delta": DELTA,
            "explicit_formula_prime_limit": EF_LIMIT,
            "bober_E20_reference_upper": str(
                BOBER_E20_BOUND
            ),
            "explicit_formula_conditional_on_GRH": True,
            "algebraic_interpretation_requires_BSD": True,
            "trace_bounds": bounds,
        },
        "meta": meta,
        "explicit_formula": {
            "anonymous": {
                "raw": ef["anonymous"],
                "log_term": nstr(
                    anon_bound["log_term"],
                    70,
                ),
                "prime_term": nstr(
                    anon_bound["prime_term"],
                    70,
                ),
                "calibrated_upper": nstr(
                    anon_bound["upper"],
                    70,
                ),
                "root_number": 1,
                "analytic_parity": "even",
                "largest_parity_compatible_integer_below_bound":
                    anon_parity_max,
            },
            "E29": {
                "raw": ef["E29"],
                "log_term": nstr(
                    e29_bound["log_term"],
                    70,
                ),
                "prime_term": nstr(
                    e29_bound["prime_term"],
                    70,
                ),
                "calibrated_upper": nstr(
                    e29_bound["upper"],
                    70,
                ),
                "root_number": -1,
                "analytic_parity": "odd",
                "largest_parity_compatible_integer_below_bound":
                    e29_parity_max,
            },
            "E20_reference": ef["E20"],
        },
        "traces": traces,
        "elapsed_seconds": elapsed,
    }

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.output.write_text(
        json.dumps(
            artifact,
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )

    print(
        f"{PROTOCOL}|stage=done"
        f"|status=PASS"
        f"|seconds={elapsed:.3f}"
        f"|out={args.output}",
        flush=True,
    )

    print(
        f"{PROTOCOL}|note="
        "explicit_formula_is_GRH_conditional_analytic_bound;"
        "do_not_claim_algebraic_rank_upper_without_BSD",
        flush=True,
    )


if __name__ == "__main__":
    main()
