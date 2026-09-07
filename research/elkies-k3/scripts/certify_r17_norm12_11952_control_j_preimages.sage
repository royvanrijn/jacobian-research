#!/usr/bin/env sage -python
"""Certify rational j-preimages of the four published R17 control curves."""

import argparse
import hashlib
import json
import shlex
import sys
from pathlib import Path

from sage.all import PolynomialRing, QQ, ZZ


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ALTERNATE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
DEFAULT_PUBLISHED = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-control-j-preimages-v1.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--alternate", type=Path, default=DEFAULT_ALTERNATE)
parser.add_argument("--published", type=Path, default=DEFAULT_PUBLISHED)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()
for name in ("alternate", "published", "output"):
    setattr(args, name, getattr(args, name).resolve())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


alternate = json.loads(args.alternate.read_text())
published = json.loads(args.published.read_text())
if alternate.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
    raise ValueError("alternate direct equation is not certified")
if published.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
    raise ValueError("published R17 model is not certified")

Ru = PolynomialRing(QQ, "u")
Rt = PolynomialRing(QQ, "t")
alt_model = alternate["weierstrass_model"]
A_alt = Ru([QQ(value) for value in alt_model["A_coefficients_low_to_high"]])
B_alt = Ru([QQ(value) for value in alt_model["B_coefficients_low_to_high"]])
A_old = Rt([QQ(value) for value in published["A_coefficients_low_to_high"]])
B_old = Rt([QQ(value) for value in published["B_coefficients_low_to_high"]])

records = []
for label, parameter in published["published_high_rank_fibre_parameters"].items():
    t0 = QQ(parameter)
    a_old = A_old(t0)
    b_old = B_old(t0)
    old_discriminant_core = 4 * a_old**3 + 27 * b_old**2
    if not old_discriminant_core:
        raise ArithmeticError(f"published control {parameter} is singular")
    # Equality of 6912*A^3/(4*A^3+27*B^2) cross-multiplies to
    # A_alt^3*B_old^2-A_old^3*B_alt^2=0.
    preimage = A_alt**3 * b_old**2 - a_old**3 * B_alt**2
    if preimage.degree() != 24:
        raise ArithmeticError("unexpected projective j-preimage degree")
    factorization = preimage.factor()
    degree_exponents = sorted(
        [[int(factor.degree()), int(exponent)] for factor, exponent in factorization]
    )
    rational_roots = preimage.roots(QQ)
    if rational_roots or any(degree == 1 for degree, unused in degree_exponents):
        raise ArithmeticError("a published control unexpectedly has a rational alternate parameter")
    records.append(
        {
            "rank_label": label,
            "published_parameter": parameter,
            "j_preimage_degree": 24,
            "factor_degree_exponents_over_QQ": degree_exponents,
            "finite_rational_roots": [],
            "root_at_infinity": False,
            "rational_alternate_fibre": False,
        }
    )

payload = {
    "schema": "elkies-k3.r17-norm12-11952-control-j-preimages.v1",
    "status": "PASS_NO_RATIONAL_ALTERNATE_J_PREIMAGE_FOR_PUBLISHED_RANK25_28_CONTROLS",
    "method": {
        "equation": "A_alt(u)^3*B_old(t0)^2-A_old(t0)^3*B_alt(u)^2=0",
        "projective_degree": 24,
        "finite_gate": "no linear factor over QQ",
        "infinity_gate": "the degree remains 24",
    },
    "controls": records,
    "conclusion": (
        "None of the four published rank-25--28 control curves occurs, even up to "
        "geometric isomorphism or quadratic twist, as a rational fibre of the direct "
        "alternate-Q80 family."
    ),
    "inputs": {
        "alternate": {"path": str(args.alternate.relative_to(ROOT)), "sha256": sha256(args.alternate)},
        "published": {"path": str(args.published.relative_to(ROOT)), "sha256": sha256(args.published)},
    },
    "worker": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
    "claim_boundary": {
        "proved": [
            "no rational parameter of the direct alternate-Q80 family has the same j-invariant as any published rank-25--28 control",
            "the four old control curves therefore cannot be used as literal rational-fibre specialization controls in the alternate chart",
        ],
        "not_proved": [
            "absence of algebraic alternate parameters of degree above one",
            "any rank statistic for rational fibres of the alternate family",
            "any comparison of low-degree multisection visibility after a different calibration set is chosen",
        ],
    },
    "reproduce": shlex.join(
        [
            str(Path(sys.executable)),
            str(Path(__file__).resolve().relative_to(ROOT)),
            "--output",
            str(args.output.relative_to(ROOT)),
        ]
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"control-j-preimage artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "R17NORM12O11952CONTROLS|controls=4|degree=24|rational_preimages=0|"
    "status=PASS_NO_RATIONAL_ALTERNATE_J_PREIMAGE_FOR_PUBLISHED_RANK25_28_CONTROLS"
)
