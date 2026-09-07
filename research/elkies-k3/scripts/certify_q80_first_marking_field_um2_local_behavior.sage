#!/usr/bin/env sage
"""Compare the u=-2 first-marking field with the p=19 and p=61 q12 fields."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--cover", type=Path,
    default=RESULTS / "q80-slope-8-87-first-marked-cover-qq.json",
)
parser.add_argument(
    "--local-parameter", type=Path,
    default=RESULTS / "q80-cm24-slope-8-87-qq-local-parameter.json",
)
parser.add_argument(
    "--p19-pencil", type=Path,
    default=RESULTS / "q80-third-q12-um2-p19-resolved-pencil.json",
)
parser.add_argument(
    "--p61-pencil", type=Path,
    default=RESULTS / "q80-third-q12-um2-p61-resolved-pencil.json",
)
parser.add_argument(
    "--output", type=Path,
    default=RESULTS / "q80-first-marking-field-um2-local-behavior.json",
)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


cover_payload = json.loads(args.cover.read_text())
local_payload = json.loads(args.local_parameter.read_text())
if cover_payload.get("status") != "PASS_EXACT_FIRST_MARKED_COVER":
    raise ValueError("uncertified first-marking cover")
if local_payload.get("status") != "PASS_EXACT_CM24_CENTERED_PARAMETER":
    raise ValueError("uncertified local parameter")

ring = PolynomialRing(QQ, "t")
t = ring.gen()
cover = ring(cover_payload["squarefree_cover"]["polynomial"])
u_cm24 = QQ(local_payload["cm24_global_parameter"])
u_value = QQ(-2)
t_value = u_value - u_cm24
radicand = QQ(cover(t_value))
if not radicand or radicand.is_square():
    raise ArithmeticError("u=-2 did not define a nonsplit quadratic number field")


def decimal_hash(integer):
    return hashlib.sha256(str(integer).encode()).hexdigest()


def local_record(prime, pencil_path):
    payload = json.loads(pencil_path.read_text())
    if payload["specialization"] != {
        "prime": prime,
        "u": "-2",
        "extension_modulus": payload["specialization"]["extension_modulus"],
    }:
        raise ValueError(f"p={prime}: unexpected pencil specialization")
    finite = GF(prime)
    finite_ring = PolynomialRing(finite, "r")
    modulus = finite_ring(payload["specialization"]["extension_modulus"])
    if modulus.degree() != 2 or not modulus.is_irreducible():
        raise ArithmeticError(f"p={prime}: q12 coefficient field is not quadratic")
    target_discriminant = modulus[1] ** 2 - 4 * modulus[0] * modulus[2]
    denominator = finite(radicand.denominator())
    if not denominator:
        raise ZeroDivisionError(f"p={prime}: first-marking radicand is not integral")
    reduced = finite(radicand.numerator()) / denominator
    if not reduced:
        raise ArithmeticError(f"p={prime}: first-marking field is ramified")
    split = bool(reduced.is_square())
    ratio = reduced / target_discriminant
    scales = [value for value in finite if value**2 == ratio]
    return {
        "prime": prime,
        "first_marking_radicand_mod_p": int(reduced),
        "first_marking_field_behavior": "split" if split else "inert",
        "target_q12_modulus": str(modulus),
        "target_q12_discriminant": int(target_discriminant),
        "target_q12_discriminant_is_square": bool(target_discriminant.is_square()),
        "radicand_to_target_discriminant_ratio": int(ratio),
        "anti_invariant_scale_roots": [int(value) for value in scales],
        "residue_fields_match": not split and bool(scales),
    }


records = {
    "19": local_record(19, args.p19_pencil),
    "61": local_record(61, args.p61_pencil),
}
if records["19"]["first_marking_radicand_mod_p"] != 16:
    raise ArithmeticError("unexpected p=19 first-marking squareclass")
if records["19"]["residue_fields_match"]:
    raise ArithmeticError("p=19 unexpectedly matched the first-marking field")
if records["61"]["first_marking_radicand_mod_p"] != 2:
    raise ArithmeticError("unexpected p=61 first-marking squareclass")
if not records["61"]["residue_fields_match"]:
    raise ArithmeticError("p=61 did not match the first-marking residue field")

output = {
    "schema": "elkies-k3.q80-first-marking-field-um2-local-behavior.v1",
    "status": "PASS_EXACT_FIRST_MARKING_FIELD_LOCAL_COMPARISON",
    "inputs": {
        "cover": {"path": str(args.cover.relative_to(ROOT)), "sha256": sha256(args.cover)},
        "local_parameter": {
            "path": str(args.local_parameter.relative_to(ROOT)),
            "sha256": sha256(args.local_parameter),
        },
        "p19_pencil": {"path": str(args.p19_pencil.relative_to(ROOT)), "sha256": sha256(args.p19_pencil)},
        "p61_pencil": {"path": str(args.p61_pencil.relative_to(ROOT)), "sha256": sha256(args.p61_pencil)},
    },
    "specialization": {
        "global_parameter": "u=-2",
        "local_parameter": "t=u-u_CM24",
        "u_CM24": str(u_cm24),
        "cover_degree": int(cover.degree()),
        "radicand_is_nonzero_nonsquare_over_QQ": True,
        "radicand_numerator_bits": int(abs(radicand.numerator()).nbits()),
        "radicand_denominator_bits": int(radicand.denominator().nbits()),
        "radicand_numerator_decimal_sha256": decimal_hash(radicand.numerator()),
        "radicand_denominator_decimal_sha256": decimal_hash(radicand.denominator()),
    },
    "local_behavior": records,
    "conclusion": (
        "The exact first-marking quadratic number field at u=-2 splits at 19, "
        "so it cannot be the source of the quadratic p=19 horizontal field. "
        "It is inert at 61 and has the same residue field as the p=61 horizontal."
    ),
    "claim_boundary": {
        "proved": [
            "the exact first-marking cover specializes to a quadratic number field at u=-2",
            "that number field splits at p=19 while the q12 horizontal needs GF(19^2)",
            "that number field is inert at p=61 and its residue field is isomorphic to the q12 field",
        ],
        "not_proved": [
            "the p=61 horizontal is the reduction of a section over the first-marking number field",
            "any characteristic-zero third-q12 horizontal or child equation",
        ],
    },
    "reproduce": "sage -python elkies-k3/scripts/certify_q80_first_marking_field_um2_local_behavior.sage",
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if args.check:
    if not args.output.exists() or args.output.read_text() != serialized:
        raise SystemExit(f"local-behavior artifact is stale: {args.output}")
else:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
print(
    "Q80FIRSTMARKINGFIELD|u=-2|p19=split|p61=inert|"
    "status=PASS_EXACT_FIRST_MARKING_FIELD_LOCAL_COMPARISON"
)
