#!/usr/bin/env python3
"""Certify F_1 for record curve 356 or 385 from a class-group quotient.

Only the 2-primary class quotient and the classes of the bad prime ideals are
used.  Fundamental units, a full BNF certification, Simon local images, and
higher descent are deliberately outside this calculation.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
PUBLIC = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
PRESSURE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-kummer-classgroup-pressure-v1.json"
)
PUBLIC_STATUS = "PASS_PINNED_PUBLIC_POINT_PROJECTION_FOR_69_RECOGNIZED_FIBRES"
PRESSURE_STATUS = "PROVED_KUMMER_FORCED_CUBIC_CLASS_GROUP_2RANK_LOWER_BOUNDS"
BNF_SCHEMA = "elliptic-curves.elkies-2026-record-pari219-bnf.v2"
BNF_STATUS = "completed_certified_class_quotient_upper"
SCHEMA = "elliptic-curves.record-complete-f1-from-class-quotient.v1"
PROTOCOL = "RECORDF1CLQ-v1"
GENERIC_RANK = 17
KNOWN_RANK = 29
KNOWN_RESIDUAL_DIMENSION = KNOWN_RANK - GENERIC_RANK
A356 = 24391876744717707263532695900840552395172973498186560300
B356 = 46943906433780620456844832699051340439698711588743845207309557656274241785479710000
A385 = 5309239946790504992658629933056863415942952216170388559928487467
B385 = 148662610051436076306955509240772635466805071470938222844950640163632678492682210346622193359526
TARGETS = {
    356: {
        "reduced_cubic": f"x^3 - x^2 - {A356}*x - {B356}",
        "pressure_cubic": f"z^3 + 4*z^2 - {16 * A356}*z + {64 * B356}",
        "field_isomorphism": "z = -4*x",
    },
    385: {
        "reduced_cubic": f"x^3 - {A385}*x - {B385}",
        "pressure_cubic": f"z^3 - 3*z^2 + {3 - A385}*z + {A385 + B385 - 1}",
        "field_isomorphism": "z = 1 - x",
    },
}


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def resolve_recorded_path(value: str, *, relative_to: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = relative_to / path
    return path.resolve()


def f2_rank(rows: list[list[int]], width: int) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        if len(row) != width:
            raise ArithmeticError("binary class-coordinate row has the wrong width")
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def gp_quote(path: Path) -> str:
    return str(path).replace("\\", "\\\\").replace('"', '\\"')


def class_coordinate_program(
    *, checkpoint: Path, polynomial: str, bad_primes: list[int], stack_bytes: int
) -> str:
    primes = ",".join(str(value) for value in bad_primes)
    return f'''default(parisizemax,{stack_bytes});
principal_famat(nf,t) =
{{
  my(I=idealhnf(nf,1),s=matsize(t));
  for(i=1,s[1],
    I=idealmul(nf,I,idealpow(nf,idealhnf(nf,t[i,1]),t[i,2]))
  );
  return(I)
}};
main() =
{{
addprimes([{primes}]);
f={polynomial};
iferr(b=read("{gp_quote(checkpoint)}"),E,
  print("{PROTOCOL}|stage=input|status=ERROR|message=",E);quit(2));
if(type(b)!="t_VEC",print("{PROTOCOL}|stage=input|status=ERROR|message=non_bnf");quit(3));
if(b.pol!=f,print("{PROTOCOL}|stage=input|status=ERROR|message=polynomial_mismatch");quit(4));
iferr(c=bnfcertify(b,1),E,
  print("{PROTOCOL}|stage=class_quotient|status=ERROR|message=",E);quit(5));
if(!c,print("{PROTOCOL}|stage=class_quotient|status=ERROR|message=certify_zero");quit(6));
print("{PROTOCOL}|stage=class_quotient|status=PASS|cyc=",b.cyc,
      "|signature=",b.nf.sign);
S=[{primes}];
for(ip=1,#S,
  p=S[ip];
  dec=idealprimedec(b,p);
  print("{PROTOCOL}|stage=rational_prime|status=PASS|p=",p,"|count=",#dec);
  for(j=1,#dec,
    P=dec[j];
    iferr(v=bnfisprincipal(b,P,5),E,
      print("{PROTOCOL}|stage=prime_ideal|status=ERROR|p=",p,
            "|index=",j,"|message=",E);quit(7));
    if(type(v)!="t_VEC" || #v!=2 || type(v[2])!="t_MAT",
      print("{PROTOCOL}|stage=prime_ideal|status=ERROR|p=",p,
            "|index=",j,"|message=missing_compact_identity");quit(8));
    iferr(J=idealmul(b,principal_famat(b,v[2]),idealfactorback(b,b.gen,v[1])),E,
      print("{PROTOCOL}|stage=prime_ideal|status=ERROR|p=",p,
            "|index=",j,"|message=",E);quit(9));
    if(idealhnf(b,J)!=idealhnf(b,P),
      print("{PROTOCOL}|stage=prime_ideal|status=ERROR|p=",p,
            "|index=",j,"|message=identity_mismatch");quit(10));
    print("{PROTOCOL}|stage=prime_ideal|status=PASS|p=",p,
          "|index=",j,"|e=",P.e,"|f=",P.f,"|coordinates=",v[1]);
  )
);
print("{PROTOCOL}|stage=complete|status=PASS");
return(1)
}};
if(!main(),quit(11));
quit(0)
'''


CLASS_RE = re.compile(
    rf"^{re.escape(PROTOCOL)}\|stage=class_quotient\|status=PASS"
    r"\|cyc=\[(?P<cyc>[^]]*)\]\|signature=\[(?P<signature>[^]]*)\]$"
)
PRIME_RE = re.compile(
    rf"^{re.escape(PROTOCOL)}\|stage=rational_prime\|status=PASS"
    r"\|p=(?P<p>\d+)\|count=(?P<count>\d+)$"
)
IDEAL_RE = re.compile(
    rf"^{re.escape(PROTOCOL)}\|stage=prime_ideal\|status=PASS"
    r"\|p=(?P<p>\d+)\|index=(?P<index>\d+)\|e=(?P<e>\d+)"
    r"\|f=(?P<f>\d+)\|coordinates=\[(?P<coordinates>[^]]*)\]~?$"
)


def integer_list(text: str) -> list[int]:
    if not text.strip():
        return []
    return [int(value.strip()) for value in text.split(",")]


def parse_gp_log(text: str, bad_primes: list[int]) -> dict[str, Any]:
    if f"{PROTOCOL}|stage=complete|status=PASS" not in text:
        raise ArithmeticError("class-coordinate replay has no terminal PASS marker")
    if f"{PROTOCOL}|stage=" in text and "|status=ERROR" in text:
        raise ArithmeticError("class-coordinate replay contains an error marker")
    class_records = [CLASS_RE.match(line) for line in text.splitlines()]
    class_records = [match for match in class_records if match]
    if len(class_records) != 1:
        raise ArithmeticError("class-coordinate replay has no unique class record")
    class_match = class_records[0]
    assert class_match is not None
    invariants = integer_list(class_match.group("cyc"))
    signature = integer_list(class_match.group("signature"))

    prime_counts: dict[int, int] = {}
    ideals = []
    for line in text.splitlines():
        match = PRIME_RE.match(line)
        if match:
            prime_counts[int(match.group("p"))] = int(match.group("count"))
            continue
        match = IDEAL_RE.match(line)
        if match:
            coordinates = integer_list(match.group("coordinates"))
            if len(coordinates) != len(invariants):
                raise ArithmeticError("a bad-prime ideal coordinate has the wrong length")
            ideals.append(
                {
                    "rational_prime": int(match.group("p")),
                    "index_one_based": int(match.group("index")),
                    "ramification_index": int(match.group("e")),
                    "residue_degree": int(match.group("f")),
                    "class_coordinates": coordinates,
                    "compact_principal_identity_verified": True,
                }
            )
    if list(prime_counts) != bad_primes:
        raise ArithmeticError("bad rational primes are missing or misordered in replay")
    expected_pairs = [
        (prime, index)
        for prime in bad_primes
        for index in range(1, prime_counts[prime] + 1)
    ]
    observed_pairs = [
        (record["rational_prime"], record["index_one_based"])
        for record in ideals
    ]
    if observed_pairs != expected_pairs:
        raise ArithmeticError("bad-prime ideal records are incomplete or misordered")
    return {
        "class_group_invariants": invariants,
        "signature": signature,
        "rational_prime_decomposition_counts": prime_counts,
        "prime_ideals": ideals,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve-id", type=int, choices=tuple(TARGETS), required=True)
    parser.add_argument("--gp", type=Path, required=True)
    parser.add_argument("--class-quotient-metadata", type=Path, required=True)
    parser.add_argument("--public-certificate", type=Path, default=PUBLIC)
    parser.add_argument("--pressure-certificate", type=Path, default=PRESSURE)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stack-bytes", type=int, default=4_000_000_000)
    parser.add_argument("--timeout-seconds", type=float, default=1800.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.stack_bytes <= 0 or args.timeout_seconds <= 0:
        parser.error("resource bounds must be positive")
    for path in (args.gp, args.class_quotient_metadata, args.public_certificate, args.pressure_certificate):
        if not path.is_file():
            raise FileNotFoundError(path)
    for path in (args.log, args.output):
        if path.exists() and not args.overwrite:
            raise FileExistsError(path)

    metadata_path = args.class_quotient_metadata.resolve()
    public_path = args.public_certificate.resolve()
    pressure_path = args.pressure_certificate.resolve()
    metadata = json.loads(metadata_path.read_text())
    public = json.loads(public_path.read_text())
    pressure = json.loads(pressure_path.read_text())
    if metadata.get("schema") != BNF_SCHEMA or metadata.get("status") != BNF_STATUS:
        raise ArithmeticError("class-quotient metadata is not certified complete")
    curve_id = args.curve_id
    target = TARGETS[curve_id]
    if int(metadata.get("curve_id", -1)) != curve_id:
        raise ArithmeticError("class-quotient metadata belongs to another curve")
    bnf_input = metadata.get("input", {})
    if bnf_input.get("mode") != "class-quotient-upper" or int(bnf_input.get("bnfcertify_flag", -1)) != 1:
        raise ArithmeticError("metadata does not certify a one-sided class quotient")
    if metadata.get("checkpoint_scope") != "certified_class_group_quotient_upper_only":
        raise ArithmeticError("metadata checkpoint has the wrong certification scope")
    backend = metadata.get("backend", {})
    if backend.get("binary_sha256") != digest(args.gp.resolve()):
        raise ArithmeticError("GP binary differs from the class-quotient backend")
    if public.get("status") != PUBLIC_STATUS or pressure.get("status") != PRESSURE_STATUS:
        raise ArithmeticError("a canonical input certificate is not passing")
    public_record = next(row for row in public["records"] if int(row["id"]) == curve_id)
    pressure_record = next(
        row for row in pressure["curves"] if int(row["curve_id"]) == curve_id
    )
    bad_primes = [int(value) for value in public_record["bad_primes"]]
    if [int(value) for value in bnf_input.get("factor_hint_primes", [])] != bad_primes:
        raise ArithmeticError("class-quotient factor hints differ from the exact bad set")
    if len(public_record.get("points", [])) != KNOWN_RANK:
        raise ArithmeticError("public certificate does not contain MW29")
    if (
        int(pressure_record["point_count"]) != KNOWN_RANK
        or int(pressure_record["residual_gain_over_mw17"]) != KNOWN_RESIDUAL_DIMENSION
        or int(pressure_record["rational_two_torsion_dimension"]) != 0
    ):
        raise ArithmeticError("pressure certificate has incompatible MW/F1 gates")

    checkpoint_value = metadata.get("checkpoint")
    checkpoint_hash = metadata.get("checkpoint_sha256")
    if not checkpoint_value or not checkpoint_hash:
        raise ArithmeticError("class-quotient metadata has no retained checkpoint")
    checkpoint = resolve_recorded_path(
        str(checkpoint_value), relative_to=ROOT
    )
    if not checkpoint.is_file() or digest(checkpoint) != checkpoint_hash:
        raise ArithmeticError("class-quotient checkpoint hash mismatch")
    polynomial = str(bnf_input["reduced_cubic"])
    if polynomial != target["reduced_cubic"]:
        raise ArithmeticError("unexpected reduced cubic for the requested record")
    if str(pressure_record["two_division_cubic"]) != target["pressure_cubic"]:
        raise ArithmeticError("the pressure cubic does not match the exact field transform")

    program = class_coordinate_program(
        checkpoint=checkpoint,
        polynomial=polynomial,
        bad_primes=bad_primes,
        stack_bytes=args.stack_bytes,
    )
    completed = subprocess.run(
        [str(args.gp.resolve()), "-q", "-f", "-s", str(args.stack_bytes)],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=args.timeout_seconds,
        check=False,
    )
    args.log.parent.mkdir(parents=True, exist_ok=True)
    args.log.write_text(completed.stdout)
    if completed.returncode != 0:
        raise ArithmeticError(
            f"class-coordinate GP replay failed with return code {completed.returncode}"
        )
    replay = parse_gp_log(completed.stdout, bad_primes)

    invariants = replay["class_group_invariants"]
    even_positions = [index for index, value in enumerate(invariants) if value % 2 == 0]
    computed_class_two_rank = len(even_positions)
    recorded_upper = metadata.get("measurement", {}).get(
        "computed_class_group_mod2_dimension"
    )
    if recorded_upper is None or int(recorded_upper) != computed_class_two_rank:
        raise ArithmeticError("replayed and recorded class-group 2-ranks differ")
    proved_class_lower = int(pressure_record["proved_class_group_2rank_lower_bound"])
    if computed_class_two_rank < proved_class_lower:
        raise ArithmeticError("class-group upper bound contradicts the known lower bound")
    class_mod2_exact = computed_class_two_rank == proved_class_lower

    class_rows = [
        [record["class_coordinates"][index] & 1 for index in even_positions]
        for record in replay["prime_ideals"]
    ]
    computed_bad_class_rank = f2_rank(class_rows, computed_class_two_rank)
    exact_bad_class_rank = computed_bad_class_rank if class_mod2_exact else None
    bad_prime_ideal_count = len(replay["prime_ideals"])
    if bad_prime_ideal_count != int(pressure_record["bad_prime_ideal_count"]):
        raise ArithmeticError("bad-prime ideal count differs from the pressure certificate")
    pressure_decomposition = [
        (
            int(record["rational_prime"]),
            int(record["prime_ideal_index_one_based"]),
            int(record["ramification_index"]),
            int(record["residue_degree"]),
        )
        for record in pressure_record["bad_prime_ideal_columns"]
    ]
    replayed_decomposition = [
        (
            record["rational_prime"],
            record["index_one_based"],
            record["ramification_index"],
            record["residue_degree"],
        )
        for record in replay["prime_ideals"]
    ]
    if replayed_decomposition != pressure_decomposition:
        raise ArithmeticError("bad-prime decomposition differs under the field transform")
    signature = replay["signature"]
    if signature != [int(value) for value in pressure_record["field_signature"]]:
        raise ArithmeticError("replayed field signature differs from the pressure certificate")
    unit_squareclass_dimension = sum(signature)
    rational_norm_target_dimension = len(bad_primes) + 1
    if class_mod2_exact:
        k_s_2_dimension_upper = (
            unit_squareclass_dimension
            + computed_class_two_rank
            + bad_prime_ideal_count
            - computed_bad_class_rank
        )
    else:
        # A one-sided presentation can overestimate class rank and the rank
        # of the bad-ideal span can collapse in the true quotient.  Hence no
        # subtraction is allowed until upper and lower class ranks coincide.
        k_s_2_dimension_upper = (
            unit_squareclass_dimension
            + computed_class_two_rank
            + bad_prime_ideal_count
        )
    norm_square_envelope_upper = (
        k_s_2_dimension_upper - rational_norm_target_dimension
    )
    if norm_square_envelope_upper < KNOWN_RANK:
        raise ArithmeticError("global norm envelope contradicts certified MW29")
    total_selmer_upper = norm_square_envelope_upper
    parity_sharpened_total_selmer_upper = total_selmer_upper - (
        (total_selmer_upper - int(pressure_record["proved_total_two_selmer_dimension_mod_2"])) & 1
    )
    relative_upper = parity_sharpened_total_selmer_upper - KNOWN_RANK
    if relative_upper < 0:
        raise ArithmeticError("parity-sharpened bound contradicts certified MW29")
    complete_f1 = relative_upper == 0
    status = (
        "CERTIFIED_COMPLETE_F1_DIMENSION_12_AND_EXACT_RANK_29"
        if complete_f1
        else "CERTIFIED_F1_DIMENSION_UPPER_BOUND"
    )

    result = {
        "schema": SCHEMA,
        "protocol": PROTOCOL,
        "status": status,
        "curve_id": curve_id,
        "inputs": {
            "class_quotient_metadata": str(metadata_path),
            "class_quotient_metadata_sha256": digest(metadata_path),
            "class_quotient_checkpoint": str(checkpoint),
            "class_quotient_checkpoint_sha256": digest(checkpoint),
            "public_certificate": str(public_path),
            "public_certificate_sha256": digest(public_path),
            "pressure_certificate": str(pressure_path),
            "pressure_certificate_sha256": digest(pressure_path),
            "gp_binary": str(args.gp.resolve()),
            "gp_binary_sha256": digest(args.gp.resolve()),
            "class_coordinate_log": str(args.log.resolve()),
            "class_coordinate_log_sha256": digest(args.log.resolve()),
        },
        "class_group_mod_2": {
            "field_isomorphism_to_pressure_certificate": target[
                "field_isomorphism"
            ],
            "computed_quotient_invariants": invariants,
            "certified_dimension_upper_bound": computed_class_two_rank,
            "certified_dimension_lower_bound": proved_class_lower,
            "dimension_exact": class_mod2_exact,
            "exact_dimension": computed_class_two_rank if class_mod2_exact else None,
        },
        "bad_prime_class_span": {
            "bad_rational_prime_count": len(bad_primes),
            "bad_prime_ideal_count": bad_prime_ideal_count,
            "rank_in_computed_class_quotient_mod_2": computed_bad_class_rank,
            "rank_in_true_class_group_mod_2": exact_bad_class_rank,
            "prime_ideals": replay["prime_ideals"],
        },
        "dimension_argument": {
            "unit_squareclass_dimension": unit_squareclass_dimension,
            "k_s_2_dimension_upper_bound": k_s_2_dimension_upper,
            "rational_norm_target_dimension": rational_norm_target_dimension,
            "norm_map_surjective": True,
            "norm_surjectivity_witnesses": ["-1", *[str(value) for value in bad_primes]],
            "norm_square_global_envelope_dimension_upper_bound": norm_square_envelope_upper,
            "total_two_selmer_dimension_lower_bound": KNOWN_RANK,
            "total_two_selmer_dimension_upper_bound": total_selmer_upper,
            "total_two_selmer_dimension_parity": int(
                pressure_record["proved_total_two_selmer_dimension_mod_2"]
            ),
            "parity_sharpened_total_two_selmer_dimension_upper_bound": (
                parity_sharpened_total_selmer_upper
            ),
            "mw29_relative_two_selmer_dimension_upper_bound": relative_upper,
            "formula": (
                "dim K(S,2) = dim(O_K^*/O_K^{*2}) + dim Cl(K)[2] + "
                "#S_K - rank(<S_K> in Cl(K)/2); norms of -1 and the bad "
                "rational primes make K(S,2) -> Q(S_Q,2) surjective"
            ),
        },
        "f1": {
            "definition": (
                f"Sel_2(E_{curve_id}/Q) modulo the specialized generic MW17 image"
            ),
            "known_subspace_dimension": KNOWN_RESIDUAL_DIMENSION,
            "known_subspace_labels": [f"P{index}" for index in range(18, 30)],
            "dimension": KNOWN_RESIDUAL_DIMENSION if complete_f1 else None,
            "dimension_upper_bound": (
                parity_sharpened_total_selmer_upper - GENERIC_RANK
            ),
            "basis": (
                [f"P{index} mod <P1,...,P17>" for index in range(18, 30)]
                if complete_f1
                else None
            ),
            "equals_known_w": complete_f1,
        },
        "rank_conclusion": {
            "exact_rank": KNOWN_RANK if complete_f1 else None,
            "rank_lower_bound": KNOWN_RANK,
            "rank_upper_bound": parity_sharpened_total_selmer_upper,
        },
        "certification": {
            "class_group_is_true_quotient_of_computed_presentation": True,
            "class_coordinate_principal_ideal_identities_verified": True,
            "class_mod2_identified_by_matching_upper_and_lower_dimensions": class_mod2_exact,
            "bad_prime_class_span_exact": class_mod2_exact,
            "unit_dimension_from_dirichlet_signature": True,
            "norm_surjectivity_exact": True,
            "known_mw29_mod2_independence_certified_by_inputs": True,
            "selmer_parity_certified_by_input": True,
            "grh_assumption": None,
            "local_conditions_needed": False if complete_f1 else None,
            "higher_descent_used": False,
        },
        "claim_boundary": (
            "The class quotient is used only as an unconditional mod-2 upper "
            "bound. Its bad-prime span is treated as exact only when that upper "
            "dimension equals the independent Kummer lower bound. The 2-Selmer "
            "group is contained in the resulting norm-square envelope; no local "
            "or higher-descent claim is needed when that envelope already equals MW29."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".tmp")
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print(
        f"{PROTOCOL}|status={status}|class2_upper={computed_class_two_rank}"
        f"|bad_class_rank={computed_bad_class_rank}"
        f"|norm_envelope_upper={norm_square_envelope_upper}"
        f"|relative_upper={relative_upper}|output={args.output}"
    )


if __name__ == "__main__":
    main()
