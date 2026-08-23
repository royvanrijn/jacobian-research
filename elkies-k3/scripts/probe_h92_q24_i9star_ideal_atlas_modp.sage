#!/usr/bin/env sage -python
'''Pull I24=(u^8,xP*y+yP*x) through the recorded modular I9* chart atlas.

This is a diagnostic companion to probe_h92_q24_d12_reduced16_modp.sage.
It reports the common exceptional order of the two generators in every chart
written by derive_h92_q24_i9star_resolution_modp.sage, and factors the
strict-transform slice on each new exceptional coordinate.

It does not yet glue irreducible factors across chart overlaps or identify
them with the deterministic D13 simple-root basis.
'''

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
    ]
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate repository")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--resolution", type=Path)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
p = ZZ(args.prime)
F = GF(p)

q8_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8 = next(
    (
        path
        for path in q8_candidates
        if path.exists()
        and json.loads(path.read_text()).get("status")
        == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
    ),
    None,
)
if Q8 is None:
    raise SystemExit("No passing q8 D13 child artifact")

RESOLUTION = (
    args.resolution.resolve()
    if args.resolution
    else LOCAL / f"q24-i9star-resolution-mod-{p}.json"
)
if not RESOLUTION.exists():
    raise SystemExit(
        f"Missing {RESOLUTION}\nRun:\n"
        f"  sage -python elkies-k3/scripts/derive_h92_q24_i9star_resolution_modp.sage --prime {p}"
    )


def red_q(value):
    value = QQ(value)
    d = ZZ(value.denominator())
    if d % p == 0:
        raise ZeroDivisionError(f"bad denominator at p={p}")
    return F(ZZ(value.numerator())) / F(d)


q8 = json.loads(Q8.read_text())
resolution = json.loads(RESOLUTION.read_text())
child = q8["child"]

BU = PolynomialRing(F, "U")
U = BU.gen()
Acurve = BU([red_q(v) for v in child["minimal_A_coefficients_low_to_high"]])
Bcurve = BU([red_q(v) for v in child["minimal_B_coefficients_low_to_high"]])

EXACT = LOCAL / "q8-q24-horizontal-section-qq.json"
MOD = LOCAL / f"q24-degree46-direct-global-mod-{p}.json"
source = None
if EXACT.exists():
    exact = json.loads(EXACT.read_text())
    if exact.get("status") == "PASS_EXACT_Q24_HORIZONTAL_SECTION":
        sec = exact["section"]
        Z = BU([red_q(v) for v in sec["Z_coefficients_low_to_high"]])
        X = BU([red_q(v) for v in sec["X_coefficients_low_to_high"]])
        Y = BU([red_q(v) for v in sec["Y_coefficients_low_to_high"]])
        source = str(EXACT.relative_to(ROOT))
if source is None:
    if not MOD.exists():
        raise SystemExit(f"Missing q24 section: {EXACT} and {MOD}")
    modular = json.loads(MOD.read_text())
    sec = modular["section_mod_p"]
    Z = BU([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
    X = BU([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
    Y = BU([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
    source = str(MOD.relative_to(ROOT))

assert Y**2 == X**3 + Acurve * X * Z**4 + Bcurve * Z**6

i9 = next(item for item in child["finite_fibres"] if item["kodaira"] == "I9*")
QU = PolynomialRing(QQ, "U")
fQ = QU(str(i9["factor"]))
f = BU([red_q(c) for c in fQ.list()])
alpha = -f[0] / f[1]

T = PolynomialRing(F, "t")
t = T.gen()
KT = T.fraction_field()


def shift(poly):
    return T(BU(poly)(alpha + t))


xloc = KT(shift(X)) / KT(shift(Z) ** 2)
yloc = KT(shift(Y)) / KT(shift(Z) ** 3)
assert xloc.valuation() == 0
assert F(xloc(0)) != 0

S = PolynomialRing(F, names=("u", "x", "y"), order="degrevlex")
u, x, y = S.gens()
KS = S.fraction_field()


def eval_univariate_rf(value, expression):
    value = KT(value)
    n = S.zero()
    d = S.zero()
    for exponent, coefficient in enumerate(T(value.numerator()).list()):
        n += S(coefficient) * expression**exponent
    for exponent, coefficient in enumerate(T(value.denominator()).list()):
        d += S(coefficient) * expression**exponent
    return KS(n) / KS(d)


def poly_var_valuation(poly, index):
    poly = S(poly)
    if not poly:
        return 10**9
    return min(exp[index] for exp, coefficient in poly.dict().items() if coefficient)


def rf_var_valuation(value, index):
    value = KS(value)
    return poly_var_valuation(value.numerator(), index) - poly_var_valuation(
        value.denominator(), index
    )


def set_exceptional_zero(poly, index):
    substitutions = [u, x, y]
    substitutions[index] = F.zero()
    return S(poly)(*substitutions)


def normalized_restriction(value, exceptional, index, order):
    value = KS(value) / KS(exceptional**order)
    n0 = set_exceptional_zero(value.numerator(), index)
    d0 = set_exceptional_zero(value.denominator(), index)
    return {
        "numerator_nonzero": bool(n0),
        "denominator_unit_generically": bool(d0),
        "numerator_restriction": str(n0),
        "denominator_restriction": str(d0),
    }


def exceptional_slice(strict, index):
    restricted = set_exceptional_zero(S(strict), index)
    if not restricted:
        return {"restriction": "0", "factors": [], "bad": True}
    try:
        factors = [
            {"factor": str(factor), "multiplicity": int(multiplicity)}
            for factor, multiplicity in restricted.factor()
        ]
        bad = len(factors) != 1 or any(v["multiplicity"] != 1 for v in factors)
    except Exception as exc:
        factors = [{"error": f"{type(exc).__name__}: {exc}"}]
        bad = True
    return {"restriction": str(restricted), "factors": factors, "bad": bool(bad)}


records = []
chart_count = 0
principal_count = 0
bad_slice_count = 0

for center in resolution.get("centers", []):
    for chart in center.get("charts", []):
        chart_count += 1
        kind = chart["chart"]
        index = {"u": 0, "x": 1, "y": 2}[kind]
        exceptional = (u, x, y)[index]
        u0, x0, y0 = tuple(S(v) for v in chart["origin_map"])

        xPpull = eval_univariate_rf(xloc, u0)
        yPpull = eval_univariate_rf(yloc, u0)
        g1 = KS(u0**8)
        g2 = xPpull * KS(y0) + yPpull * KS(x0)

        v1 = rf_var_valuation(g1, index)
        v2 = rf_var_valuation(g2, index)
        common = min(v1, v2)
        r1 = normalized_restriction(g1, exceptional, index, common)
        r2 = normalized_restriction(g2, exceptional, index, common)
        principal = (
            r1["denominator_unit_generically"]
            and r2["denominator_unit_generically"]
            and (r1["numerator_nonzero"] or r2["numerator_nonzero"])
        )
        principal_count += int(principal)

        slice_data = exceptional_slice(chart["strict_transform"], index)
        bad_slice_count += int(slice_data["bad"])
        records.append(
            {
                "center": center["label"],
                "depth": int(center["depth"]),
                "chart": kind,
                "origin_map": list(chart["origin_map"]),
                "valuation_u8": int(v1),
                "valuation_L": int(v2),
                "common_exceptional_order": int(common),
                "normalized_u8": r1,
                "normalized_L": r2,
                "generically_principal": bool(principal),
                "exceptional_slice": slice_data,
            }
        )
        print(
            "Q24IDEALCHART|"
            f"center={center['label']}|depth={center['depth']}|chart={kind}|"
            f"v_u8={v1}|v_L={v2}|common={common}|"
            f"principal={int(principal)}|slice_bad={int(slice_data['bad'])}|"
            "status=PASS_DIAGNOSTIC",
            flush=True,
        )

if not chart_count:
    raise ArithmeticError("resolution artifact has no chart records")
all_principal = principal_count == chart_count

payload = {
    "schema": "elkies-k3.h3-q24-i9star-ideal-atlas-modp.v1",
    "status": (
        "PASS_Q24_IDEAL_GENERICALLY_PRINCIPAL_ON_RECORDED_CHARTS"
        if all_principal
        else "Q24_IDEAL_NOT_PRINCIPAL_ON_ALL_RECORDED_CHARTS"
    ),
    "proof_boundary": (
        "This checks recorded charts only. It does not glue reduced exceptional "
        "factors or match them to the deterministic D13 root basis."
    ),
    "prime": int(p),
    "inputs": {
        "q8_child": str(Q8.relative_to(ROOT)),
        "q24_section": source,
        "resolution": str(RESOLUTION.relative_to(ROOT)),
    },
    "ideal": ["u^8", "xP*y+yP*x"],
    "summary": {
        "charts": int(chart_count),
        "generically_principal": int(principal_count),
        "all_generically_principal": bool(all_principal),
        "bad_exceptional_slices": int(bad_slice_count),
    },
    "charts": records,
}
OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-i9star-ideal-atlas-mod-{p}.json"
)
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24IDEALATLAS_RESULT|"
    f"prime={p}|charts={chart_count}|principal={principal_count}|"
    f"bad_slices={bad_slice_count}|status={payload['status']}",
    flush=True,
)
