#!/usr/bin/env python3
"""Exact two-remainder gcd controls for fixed homogeneous SFC supports.

Build the existing polynomial chart over ZZ[L,H], retaining termwise
factor products. No height cutoff or numerical root approximation is used.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import platform
from pathlib import Path
import time

import flint
from flint import fmpz, fmpz_poly, nmod_mpoly_ctx

from verify_sfc_homogeneous_trinomial_factored_sign import (
    LEFT, RIGHT, RING, coefficient, expand_factors,
)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def primitive(poly):
    require(bool(poly), "Zero polynomial needs a separate common-component analysis.")
    return poly // poly.content()


def exact_divide(poly, divisor):
    quotient, remainder = divmod(poly, divisor)
    require(not remainder, "Claimed factor does not divide the exact polynomial.")
    return quotient


def factored_sum(terms):
    """Each term is (sign, [(linear-factor multiset, quotient), ...])."""
    products = []
    for sign, factors in terms:
        common = Counter()
        value = RING.constant(sign)
        for part, quotient in factors:
            common += part
            value *= quotient
        products.append((common, value))
    common = products[0][0].copy()
    for factors, _ in products[1:]:
        common &= factors
    value = sum(expand_factors(factors - common) * poly for factors, poly in products)
    return common, primitive(value)


def remainder_pair(middle, span):
    aa, bb, cc = [coefficient(2, k, middle, span) for k in range(3)]
    dd, ee, ff, gg = [coefficient(3, k, middle, span) for k in range(4)]
    ss = factored_sum(((1, (cc, cc, dd)), (1, (aa, bb, gg)), (-1, (aa, cc, ff))))
    tt = factored_sum(((1, (cc, cc, ee)), (-1, (bb, cc, ff)),
                       (1, (bb, bb, gg)), (-1, (aa, cc, gg))))
    return ss, tt


def summary(poly):
    terms = poly.to_dict()
    return {"degree": int(max(map(sum, terms))),
            "bidegree": [int(x) for x in poly.degrees()], "terms": len(terms)}


def digest(poly):
    # Stable numeric serialization, independent of FLINT pretty-printing.
    result = hashlib.sha256()
    for powers, value in sorted(poly.to_dict().items()):
        result.update(f"{powers[0]},{powers[1]}:{value}\n".encode("ascii"))
    return result.hexdigest()


def positive_coefficients(poly):
    terms = poly.to_dict()
    return terms.get((0, 0), 0) > 0 and min(terms.values()) >= 0


def coprimality_certificate(ss, tt, prime):
    require(fmpz(prime).is_prime(), "Coprimality reduction requires a prime field.")
    modring = nmod_mpoly_ctx.get(("L", "H"), prime)
    sm = modring.from_dict({k: int(v) % prime for k, v in ss.to_dict().items()})
    tm = modring.from_dict({k: int(v) % prime for k, v in tt.to_dict().items()})
    require(summary(ss)["degree"] == summary(sm)["degree"] and
            summary(tt)["degree"] == summary(tm)["degree"], "Bad reduction loses total degree.")
    require(sm.gcd(tm).is_constant(), "No modular coprimality certificate at this prime.")


def factor_blocks(factors):
    grouped = defaultdict(list)
    for (axis, step, shift), multiplicity in sorted(factors.items()):
        grouped[axis, step, multiplicity].append(shift)
    blocks = []
    for (axis, step, multiplicity), shifts in sorted(grouped.items()):
        stride = 2 if step == 2 else 1
        first = last = shifts[0]
        for shift in shifts[1:] + [None]:
            if shift is not None and shift == last + stride:
                last = shift
                continue
            blocks.append({"axis": "LH"[axis], "step": step,
                           "first_shift": first, "last_shift": last,
                           "shift_stride": stride, "multiplicity": multiplicity})
            first = last = shift
    return blocks


def nonnegative_integer_roots(poly, axis):
    require(bool(poly), "The eliminant is zero.")
    values = poly.to_dict()
    require(all(k[1-axis] == 0 for k in values), "The eliminant is not univariate.")
    coefficients = [0] * (max(k[axis] for k in values) + 1)
    for powers, value in values.items():
        coefficients[powers[axis]] = int(value)
    _, factors = fmpz_poly(coefficients).factor()
    roots = set()
    for factor, _ in factors:
        if factor.degree() == 1:
            numerator, denominator = -int(factor[0]), int(factor[1])
            if numerator % denominator == 0 and numerator // denominator >= 0:
                roots.add(numerator // denominator)
    return sorted(roots)


def integer_elimination(ss, tt):
    """Both projections are supersets; test candidates on the same fibre."""
    require(ss.degrees()[1] + tt.degrees()[1] > 0 and
            ss.degrees()[0] + tt.degrees()[0] > 0, "Use a univariate branch for this pair.")
    rl, rh = ss.resultant(tt, "H"), ss.resultant(tt, "L")
    left, right = nonnegative_integer_roots(rl, 0), nonnegative_integer_roots(rh, 1)
    candidates = [(l, h) for l in left for h in right]
    common = [(l, h) for l, h in candidates
              if ss.subs({"L": l, "H": h}) == 0 and tt.subs({"L": l, "H": h}) == 0]
    return {"L_eliminant": {**summary(rl), "sha256": digest(rl)},
            "H_eliminant": {**summary(rh), "sha256": digest(rh)},
            "nonnegative_integer_L_roots": left, "nonnegative_integer_H_roots": right,
            "candidate_pairs": candidates, "common_pairs": common}


def verify_controls():
    from verify_sfc_trinomial_reduced_sieve import projected_coefficients

    checked = 0
    trinomial_eliminations = []
    for span in range(2, 5):
        for middle in range(1, span):
            (sf, ss), (tf, tt) = remainder_pair(middle, span)
            aa, bb, cc = projected_coefficients(2, middle, span)
            dd, ee, ff, gg = projected_coefficients(3, middle, span)
            direct_s = cc*cc*dd+aa*bb*gg-aa*cc*ff
            direct_t = cc*cc*ee-bb*cc*ff+(bb*bb-aa*cc)*gg
            require(primitive(expand_factors(sf)*ss) == primitive(direct_s),
                    "S differs from the independent direct moment expansion.")
            require(primitive(expand_factors(tf)*tt) ==
                    primitive(direct_t),
                    "T differs from the independent direct moment expansion.")
            checked += 1
            if span == 3:
                common = direct_s.gcd(direct_t)
                result = integer_elimination(primitive(exact_divide(direct_s, common)),
                                             primitive(exact_divide(direct_t, common)))
                original_pairs = [(l, h) for l, h in result["candidate_pairs"]
                                  if direct_s.subs({"L": l, "H": h}) == 0 and
                                  direct_t.subs({"L": l, "H": h}) == 0]
                require(original_pairs == result["common_pairs"] == [],
                        "Small-support candidates were not excluded in the original remainders.")
                trinomial_eliminations.append({"middle": middle, "span": span, **result,
                                              "original_common_pairs": original_pairs})

    # Leading coefficients can vanish together without a common affine zero.
    spurious = integer_elimination(LEFT*RIGHT+1, LEFT*RIGHT+2)
    require(spurious["candidate_pairs"] == [(0, 0)] and not spurious["common_pairs"],
            "A degree-drop projection was mistaken for a common affine root.")
    actual = integer_elimination(LEFT+RIGHT-3, LEFT*RIGHT-2)
    require(actual["common_pairs"] == [(1, 2), (2, 1)],
            "The two projections did not retain genuine nonnegative common roots.")
    try:
        coprimality_certificate((LEFT+RIGHT+1)*(LEFT+2),
                               (LEFT+RIGHT+1)*(RIGHT+3), 1000000007)
    except ArithmeticError:
        pass
    else:
        raise ArithmeticError("A corrupted gcd omission was not rejected.")
    return {"direct_moment_profiles": checked,
            "trinomial_eliminations": trinomial_eliminations,
            "degree_drop_control": spurious, "common_root_control": actual,
            "corrupted_gcd_control": "omitting a common linear factor is rejected"}


def run_profile(middle, span, prime=1000000007, eliminate=False):
    started = time.monotonic()
    print(f"({middle},{span}): building exact factored remainders", flush=True)
    (sf, ss), (tf, tt) = remainder_pair(middle, span)
    common = sf & tf
    ss = primitive(expand_factors(sf - common) * ss)
    tt = primitive(expand_factors(tf - common) * tt)
    total = LEFT + RIGHT + span
    ss = primitive(exact_divide(ss, total**3))
    tt = primitive(exact_divide(tt, total**3))
    print(f"({middle},{span}): after termwise factors and N^3: "
          f"{summary(ss)}, {summary(tt)}", flush=True)
    modring = nmod_mpoly_ctx.get(("L", "H"), prime)
    sm = modring.from_dict({k: int(v) % prime for k, v in ss.to_dict().items()})
    tm = modring.from_dict({k: int(v) % prime for k, v in tt.to_dict().items()})
    modular_gcd = sm.gcd(tm)
    print(f"({middle},{span}): modular residual gcd: {summary(modular_gcd)}", flush=True)
    if modular_gcd.is_constant():
        common_extra = RING.constant(1)
    else:
        print(f"({middle},{span}): computing characteristic-zero residual gcd", flush=True)
        common_extra = primitive(ss.gcd(tt))
    print(f"({middle},{span}): residual gcd: {summary(common_extra)}", flush=True)
    sr, tr = primitive(exact_divide(ss, common_extra)), primitive(exact_divide(tt, common_extra))
    coprimality_certificate(sr, tr, prime)
    scalar, factors = common_extra.factor()
    harmless = True
    extra = []
    extra_axes = Counter()
    for factor, multiplicity in factors:
        coefficients = factor.to_dict()
        positive = (max(map(sum, coefficients)) == 1 and
                    coefficients.get((0, 0), 0) > 0 and min(coefficients.values()) >= 0)
        harmless &= positive
        if positive and ((1, 0) not in coefficients or (0, 1) not in coefficients):
            axis = 0 if (1, 0) in coefficients else 1
            powers = (1, 0) if axis == 0 else (0, 1)
            extra_axes[axis, int(coefficients[powers]), int(coefficients[0, 0])] += int(multiplicity)
        else:
            extra.append({"factor": str(factor), "multiplicity": int(multiplicity),
                          "certified_positive_linear": positive})
    if span == 81 and middle in (1, 2, 40, 80):
        expected = Counter({(0, 1, j): 3 + 4*(j <= middle) for j in range(1, 82)})
        expected.update((0, 2, j) for j in range(1, middle+1, 2))
        expected.update((1, 1, j) for j in range(1, max(81-middle, 40)+1))
        expected.update({(1, 1, j): 2 for j in range(1, 82-middle)})
        expected.update((0, 1, j) for j in range(middle+1, 41))
        require(not extra and common + extra_axes == expected,
                "The computed gcd differs from the displayed span-81 formula.")
    record = {"middle": middle, "span": span,
              "termwise_common_factor_blocks": factor_blocks(common),
              "total_degree_factor_multiplicity_removed": 3,
              "additional_axis_factor_blocks": factor_blocks(extra_axes),
              "additional_common_factors": extra,
              "gcd_certified_positive_on_nonnegative_quadrant": harmless,
              "reduced_S": {**summary(sr), "sha256": digest(sr)},
              "reduced_T": {**summary(tr), "sha256": digest(tr)},
              "coprimality_prime": prime,
              "degree_preserving_modular_gcd": "1",
              "reduced_S_nonnegative_coefficients_positive_constant": positive_coefficients(sr),
              "elapsed_seconds": round(time.monotonic() - started, 3)}
    if eliminate:
        require(harmless, "Retain the unresolved common-factor zero set before elimination.")
        record["integer_elimination"] = integer_elimination(sr, tr)
    print(f"({middle},{span}): finished in {record['elapsed_seconds']}s; "
          f"harmless={harmless}", flush=True)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--span", type=int, default=81)
    parser.add_argument("--gaps", type=int, nargs="+", default=[1, 2, 40, 80])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--eliminate", action="store_true",
                        help="Compute both characteristic-zero eliminants; use small controls first.")
    args = parser.parse_args()
    require(args.span >= 2 and all(0 < a < args.span for a in args.gaps), "Invalid support.")
    controls = verify_controls()
    records = [run_profile(a, args.span, eliminate=args.eliminate) for a in args.gaps]
    result = {"status": "PASS", "python": platform.python_version(),
              "python_flint": flint.__version__,
              "scope": "exact gcd factorization and coprimality at the listed fixed supports only",
              "controls": controls,
              "profiles": records,
              "assurance": "internal computation; no independent review or formal verification"}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
