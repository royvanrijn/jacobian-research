#!/usr/bin/env sage
"""Run post-discovery closure after quotienting all 29 known Kummer classes.

This is the certified-BNF backend for the quotient-native descent.  It does
not construct a full Selmer basis or enumerate coverings.  Its order is:

1. construct the global S-squareclass norm envelope from a certified BNF;
2. embed all known rational points in that envelope and prove their exact
   squareclass coordinates;
3. quotient those coordinates immediately;
4. add exact local conditions monotonically, checkpointing after each place;
5. stop as soon as the residual upper bound is zero.

The BNF is only one provider for the global envelope.  The emitted ambient
manifest is also the interchange format for a future F2-only class-relation
provider, so the local quotient machinery does not depend on how the global
upper bound was certified.

This backend assumes that all twelve directions beyond MW17 are already
known.  It is not a prospective MW17-relative replay and cannot calibrate or
authorize the Selmer candidate gate.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
import time
from typing import Any

from sage.all import QQ, pari, prime_range
from sage.env import SAGE_EXTCODE
from sage.misc.randstate import set_random_seed


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
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
SCHEMA = "elliptic-curves.mw29-relative-2selmer-from-bnf.v1"
PROTOCOL = "MW29REL2BNF-v1"
HECKE_PROTOCOL = "ELKIES356HECKEGLOBAL2-v1"
HECKE_PROVIDER = CAS / "run_elkies_2026_curve356_hecke_s_squareclasses.jl"
HECKE_REDUCED_CUBIC = (
    "x^3 - x^2 - "
    "24391876744717707263532695900840552395172973498186560300*x - "
    "46943906433780620456844832699051340439698711588743845207309557656274241785479710000"
)
HECKE_BASIS_RE = re.compile(
    rf"^{re.escape(HECKE_PROTOCOL)}\|stage=norm_basis\|status=PASS"
    r"\|index=(?P<index>\d+)\|coefficients=(?P<coefficients>[^|]+)"
    r"\|norm=(?P<norm>.+)$"
)
HECKE_COMPLETE_RE = re.compile(
    rf"^{re.escape(HECKE_PROTOCOL)}\|stage=complete\|status=PASS"
    r"\|global_s_squareclass_dimension=(?P<sclass>\d+)"
    r"\|norm_square_envelope_dimension=(?P<norm>\d+)"
    r"\|class_group_two_rank=(?P<class2>\d+)$"
)

sys.path.insert(0, str(CAS))
from build_mw29_relative_2selmer_matrix import build_certificate  # noqa: E402
from residual_selmer_quotient import (  # noqa: E402
    f2_dot,
    f2_nullspace_basis,
    f2_rank_rows,
)
from run_elkies_2026_relative_2selmer_checkpointed import (  # noqa: E402
    SIMON_GP_FUNCTION,
)
from run_fermigier_rank20_auxiliary_fingerprints import (  # noqa: E402
    prime_local_rows,
)


RELATIVE_GP_FUNCTIONS = r'''
ell2global_norm_envelope(ell,bnf,K,curve_theta) =
{
  my(A,B,C,polrel,polprime,badideal,badprimes,S,LS2,normspace);
  if(#ell < 13,ell=ellinit(ell));
  if(ell.a1 != 0 || ell.a3 != 0,error("ell2global_norm_envelope: nonzero a1/a3"));
  A=ell.a2; B=ell.a4; C=ell.a6;
  polrel=Pol([1,A,B,C]); polprime=polrel';
  polprime=subst(polprime,variable(polprime),curve_theta);
  badideal=abs(K*idealadd(bnf,polprime,bnf.index));
  S=bnfpSelmer(bnf,badideal,2); LS2=S[1]; S=S[2];
  normspace=kernorm(LS2,vector(#S,i,S[i].p),2);
  badprimes=factorint(badideal[1,1]*2*numerator(ell.disc))[,1];
  return([LS2,lift(normspace),badprimes]);
};

/* MW29_RELATIVE_GP_DEFINITION_SPLIT */

ell2allowed_at_place(ell,bnf,K,curve_theta,LS2,p) =
{
  my(A,B,C,polrel,theta_embeddings,real_place,signs,pp,prank,locimage,LS2image,allowed);
  if(#ell < 13,ell=ellinit(ell));
  A=ell.a2; B=ell.a4; C=ell.a6; polrel=Pol([1,A,B,C]);
  if(p==-1,
    if(bnf.r1==3,
      theta_embeddings=nfeltembed(bnf.nf,curve_theta);real_place=1;
      for(i=2,#theta_embeddings,if(theta_embeddings[i]<theta_embeddings[real_place],real_place=i));
      signs=vector(#LS2,i,nfeltsign(bnf.nf,LS2[i],real_place)<0);
      allowed=matker(Mat(signs*Mod(1,2)))*Mod(1,2);
      return([lift(allowed),1,1,matrank(Mat(signs*Mod(1,2)))])
    , return([matid(#LS2),0,0,0]))
  );
  pp=ppinit(bnf.nf,p); prank=#pp-(p!=2);
  locimage=elllocalimage_mapped(bnf.nf,pp,K,polrel,curve_theta);
  LS2image=LS2localimage(bnf.nf,LS2,pp);
  locimage=matintersect(LS2image,locimage);
  allowed=concat(matker(LS2image),matinverseimage(LS2image,locimage)*Mod(1,2));
  allowed=matimage(allowed*Mod(1,2));
  return([lift(allowed),prank,#locimage,#LS2image]);
};

/* MW29_RELATIVE_GP_DEFINITION_SPLIT */

ell2allowed_at_place_nf(ell,nf,K,curve_theta,LS2,p) =
{
  my(A,B,C,polrel,theta_embeddings,real_place,signs,pp,prank,locimage,LS2image,allowed);
  if(#ell < 13,ell=ellinit(ell));
  A=ell.a2; B=ell.a4; C=ell.a6; polrel=Pol([1,A,B,C]);
  if(p==-1,
    if(nf.r1==3,
      theta_embeddings=nfeltembed(nf,curve_theta);real_place=1;
      for(i=2,#theta_embeddings,if(theta_embeddings[i]<theta_embeddings[real_place],real_place=i));
      signs=vector(#LS2,i,nfeltsign(nf,LS2[i],real_place)<0);
      allowed=matker(Mat(signs*Mod(1,2)))*Mod(1,2);
      return([lift(allowed),1,1,matrank(Mat(signs*Mod(1,2)))])
    , return([matid(#LS2),0,0,0]))
  );
  pp=ppinit(nf,p); prank=#pp-(p!=2);
  locimage=elllocalimage_mapped(nf,pp,K,polrel,curve_theta);
  LS2image=LS2localimage(nf,LS2,pp);
  locimage=matintersect(LS2image,locimage);
  allowed=concat(matker(LS2image),matinverseimage(LS2image,locimage)*Mod(1,2));
  allowed=matimage(allowed*Mod(1,2));
  return([lift(allowed),prank,#locimage,#LS2image]);
};
'''


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def binary_columns(matrix: Any) -> list[list[int]]:
    return [[int(entry) & 1 for entry in column] for column in matrix]


def solve_rows(vectors: list[list[int]], target: list[int]) -> list[int] | None:
    """Solve a GF(2) row combination, retaining exact generator indices."""

    basis: dict[int, tuple[int, int]] = {}
    for index, vector in enumerate(vectors):
        value = sum((int(bit) & 1) << column for column, bit in enumerate(vector))
        combination = 1 << index
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot][0]
                combination ^= basis[pivot][1]
            else:
                basis[pivot] = (value, combination)
                break
    value = sum((int(bit) & 1) << column for column, bit in enumerate(target))
    combination = 0
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            return None
        value ^= basis[pivot][0]
        combination ^= basis[pivot][1]
    return [(combination >> index) & 1 for index in range(len(vectors))]


def allowed_in_norm_coordinates(
    allowed_columns_in_s: list[list[int]], norm_columns_in_s: list[list[int]]
) -> list[list[int]]:
    """Pull an allowed S-squareclass subspace back to norm coordinates."""

    s_width = len(norm_columns_in_s[0]) if norm_columns_in_s else (
        len(allowed_columns_in_s[0]) if allowed_columns_in_s else 0
    )
    allowed_annihilator = f2_nullspace_basis(allowed_columns_in_s, s_width)
    equations_on_norm = [
        [f2_dot(functional, column) for column in norm_columns_in_s]
        for functional in allowed_annihilator
    ]
    return f2_nullspace_basis(equations_on_norm, len(norm_columns_in_s))


def transform_point(meta: dict[str, Any], values: list[str]) -> Any:
    point = pari([QQ(values[0]), QQ(values[1])])
    for change in meta["point_change_sequence"]:
        point = pari.ellchangepoint(
            point, pari([QQ(value) for value in change])
        )
    return point


def point_coordinates(
    *,
    nf: Any,
    field_polynomial: Any,
    curve_theta: Any,
    ls2: Any,
    norm_columns: list[list[int]],
    point_values: list[list[str]],
    meta: dict[str, Any],
    auxiliary_prime_bound: int,
) -> tuple[list[list[int]], list[list[int]], list[dict[str, int]]]:
    ls2_alphas = [alpha for alpha in ls2]
    points = [transform_point(meta, values) for values in point_values]
    point_alphas = [pari.Mod(point[0], field_polynomial) - curve_theta for point in points]
    all_alphas = ls2_alphas + point_alphas
    signatures: list[list[int]] = [[] for _ in all_alphas]
    ls2_coordinates: list[list[int] | None] = [None] * len(points)
    prime_records: list[dict[str, int]] = []
    bad_discriminant = int(pari.poldisc(field_polynomial))
    nfissquare = pari("nfissquare")

    for q_value in prime_range(3, auxiliary_prime_bound + 1):
        q = int(q_value)
        if bad_discriminant % q == 0:
            continue
        local, places = prime_local_rows(pari, nf, all_alphas, q)
        for index, row in enumerate(local):
            signatures[index].extend(int(bit) & 1 for bit in row)
        ls2_signatures = signatures[: len(ls2_alphas)]
        signature_rank = f2_rank_rows(ls2_signatures, len(signatures[0]))
        prime_records.append(
            {"prime": q, "place_count": len(places), "ls2_signature_rank": signature_rank}
        )
        if signature_rank < len(ls2_alphas):
            continue
        for point_index, alpha in enumerate(point_alphas):
            if ls2_coordinates[point_index] is not None:
                continue
            candidate = solve_rows(
                ls2_signatures, signatures[len(ls2_alphas) + point_index]
            )
            if candidate is None:
                continue
            ratio = alpha
            for bit, basis_alpha in zip(candidate, ls2_alphas):
                if bit:
                    ratio /= basis_alpha
            if bool(nfissquare(nf, ratio)):
                ls2_coordinates[point_index] = candidate
        if all(row is not None for row in ls2_coordinates):
            break

    if any(row is None for row in ls2_coordinates):
        missing = [index + 1 for index, row in enumerate(ls2_coordinates) if row is None]
        raise ArithmeticError(
            f"auxiliary fingerprints did not embed point indices {missing}"
        )
    exact_ls2 = [list(row) for row in ls2_coordinates if row is not None]
    norm_coordinates = []
    for index, row in enumerate(exact_ls2):
        coordinates = solve_rows(norm_columns, row)
        if coordinates is None:
            raise ArithmeticError(f"known point P{index + 1} misses the norm envelope")
        norm_coordinates.append(coordinates)
    if f2_rank_rows(norm_coordinates, len(norm_columns)) != len(norm_coordinates):
        raise ArithmeticError("known Kummer images are not independent in the norm envelope")
    return exact_ls2, norm_coordinates, prime_records


def prioritized_places(curve_id: int, bad_primes: list[int]) -> list[int]:
    """Cheap prior; the emitted matrix supplies the exact posterior order."""

    finite = list(dict.fromkeys(int(prime) for prime in bad_primes))
    # These prefixes are the exact minimum distinguishing cuts for the known
    # rigid-invisible ten-direction presentations.  They are scheduling priors
    # only: the unknown MW29 quotient can have different local behavior.
    prior = (
        [2, 3, 13, 23, 751]
        if curve_id == 356
        else [13, 29, 47, 89]
    )
    preferred = [prime for prime in prior if prime in finite]
    ordered = [*preferred, *(prime for prime in finite if prime not in preferred)]
    if curve_id == 385 and 2 in ordered:
        # Defer the conspicuous I_16 fibre until after the odd-place prior.
        ordered.remove(2)
        ordered.append(2)
    return [*preferred, -1, *(prime for prime in ordered if prime not in preferred)]


def parse_hecke_global_log(path: Path) -> dict[str, Any]:
    """Accept only the terminal output of the pinned unconditional provider."""

    text = path.read_text(errors="replace")
    required = (
        f"{HECKE_PROTOCOL}|stage=input|status=PASS|hecke_version=0.37.6"
        "|random_seed=20260904|curve=356|grh=false",
        f"{HECKE_PROTOCOL}|stage=maximal_order|status=PASS",
        f"{HECKE_PROTOCOL}|stage=class_group|status=PASS|grh=false",
        f"{HECKE_PROTOCOL}|stage=unit_group|status=PASS|grh=false",
    )
    missing = [marker for marker in required if marker not in text]
    if missing:
        raise ArithmeticError(
            "Hecke global-provider log lacks unconditional completion markers: "
            + ", ".join(missing)
        )
    basis = []
    complete = None
    for line in text.splitlines():
        match = HECKE_BASIS_RE.match(line)
        if match:
            coefficients = match.group("coefficients").split(",")
            if len(coefficients) != 3:
                raise ArithmeticError("Hecke norm-basis row is not cubic")
            basis.append(
                {
                    "index": int(match.group("index")),
                    "coefficients": coefficients,
                    "norm": match.group("norm"),
                }
            )
            continue
        match = HECKE_COMPLETE_RE.match(line)
        if match:
            complete = {
                "global_s_squareclass_dimension": int(match.group("sclass")),
                "norm_square_envelope_dimension": int(match.group("norm")),
                "class_group_two_rank": int(match.group("class2")),
            }
    if complete is None:
        raise ArithmeticError("Hecke global-provider log has no terminal PASS record")
    norm_dimension = complete["norm_square_envelope_dimension"]
    if [row["index"] for row in basis] != list(range(1, norm_dimension + 1)):
        raise ArithmeticError("Hecke norm-basis rows are incomplete or misordered")
    complete["basis"] = basis
    return complete


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve-id", type=int, choices=(356, 385), required=True)
    provider = parser.add_mutually_exclusive_group(required=True)
    provider.add_argument("--bnf-metadata", type=Path)
    provider.add_argument(
        "--hecke-global-log",
        type=Path,
        help="completed output of the pinned unconditional curve-356 Hecke provider",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--public-certificate", type=Path, default=PUBLIC)
    parser.add_argument("--classgroup-pressure-certificate", type=Path, default=PRESSURE)
    parser.add_argument("--auxiliary-prime-bound", type=int, default=10000)
    parser.add_argument("--pari-stack-bytes", type=int, default=4_000_000_000)
    parser.add_argument("--random-seed", type=int, default=20260904)
    parser.add_argument("--maximum-cut-size", type=int, default=8)
    parser.add_argument("--maximum-cut-subsets", type=int, default=1_000_000)
    parser.add_argument("--stop-after-places", type=int)
    parser.add_argument(
        "--place-order",
        help="comma-separated finite primes and/or infinity; omitted places follow",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if min(
        args.auxiliary_prime_bound,
        args.pari_stack_bytes,
        args.maximum_cut_subsets,
    ) <= 0:
        parser.error("resource and search bounds must be positive")
    if args.maximum_cut_size < 0:
        parser.error("--maximum-cut-size must be nonnegative")
    if args.stop_after_places is not None and args.stop_after_places <= 0:
        parser.error("--stop-after-places must be positive")
    if args.resume and args.overwrite:
        parser.error("--resume and --overwrite are mutually exclusive")
    return args


def main() -> None:
    args = parse_args()
    if args.output.exists() and not (args.resume or args.overwrite):
        raise FileExistsError(args.output)
    public_path = args.public_certificate.resolve()
    pressure_path = args.classgroup_pressure_certificate.resolve()
    provider_kind = "bnf" if args.bnf_metadata is not None else "hecke"
    provider_path = (
        args.bnf_metadata if args.bnf_metadata is not None else args.hecke_global_log
    ).resolve()
    public = json.loads(public_path.read_text())
    pressure = json.loads(pressure_path.read_text())
    meta = json.loads(provider_path.read_text()) if provider_kind == "bnf" else None
    if public.get("status") != PUBLIC_STATUS:
        raise ArithmeticError("public point certificate is not passing")
    if pressure.get("status") != PRESSURE_STATUS:
        raise ArithmeticError("class-group pressure/parity certificate is not passing")
    record = next(
        row for row in public["records"] if int(row["id"]) == args.curve_id
    )
    pressure_record = next(
        row for row in pressure["curves"] if int(row["curve_id"]) == args.curve_id
    )
    point_values = [[str(value) for value in point] for point in record["points"]]
    if len(point_values) != 29:
        raise ArithmeticError("record fibre does not supply exactly 29 points")
    if (
        int(pressure_record["point_count"]) != 29
        or int(pressure_record["rational_two_torsion_dimension"]) != 0
    ):
        raise ArithmeticError("pressure certificate has incompatible MW29 gates")
    residual_parity = (
        int(pressure_record["proved_total_two_selmer_dimension_mod_2"]) - 29
    ) % 2
    expected_model = [str(value) for value in record["ainvs"]]
    if provider_kind == "bnf":
        assert meta is not None
        if meta.get("global_minimal_model") != expected_model:
            raise ArithmeticError("BNF metadata belongs to a different curve model")
        checkpoint_path = Path(meta["bnf_checkpoint"])
        if not checkpoint_path.is_absolute():
            checkpoint_path = (ROOT / checkpoint_path).resolve()
        if digest(checkpoint_path) != meta["bnf_checkpoint_sha256"]:
            raise ArithmeticError("BNF binary hash mismatch")
        transformed_meta = meta
    else:
        if args.curve_id != 356:
            raise ArithmeticError("the pinned Hecke global provider is only for curve 356")
        hecke_global = parse_hecke_global_log(provider_path)
        transformed_meta = {
            "global_minimal_model": expected_model,
            "transformed_model": expected_model,
            "point_change_sequence": [],
            "field_cubic": HECKE_REDUCED_CUBIC,
            "curve_theta_in_field": f"Mod(-x, {HECKE_REDUCED_CUBIC})",
        }
        checkpoint_path = None

    pari.allocatemem(args.pari_stack_bytes)
    set_random_seed(args.random_seed)
    if provider_kind == "bnf":
        assert checkpoint_path is not None and meta is not None
        from research_runtime.pari_context import certified_bnf_checkpoint
        bnf = certified_bnf_checkpoint(checkpoint_path, digest(checkpoint_path), discover=True)
        if str(bnf.nf_get_pol()) != meta["field_cubic"]:
            raise ArithmeticError("BNF field polynomial mismatch")
        nf = bnf
    else:
        pari.addprimes([int(value) for value in record["bad_primes"]])
        from research_runtime.pari_context import prepared_nf
        nf = prepared_nf(pari(HECKE_REDUCED_CUBIC), record["bad_primes"])
        bnf = None
    simon = Path(SAGE_EXTCODE) / "pari" / "simon"
    for name in ("ellQ.gp", "ell.gp", "qfsolve.gp", "resultant3.gp"):
        pari.read(simon / name)
    pari("DEBUGLEVEL_ell=0;LIMBIGPRIME=0;LIM1=0;LIM3=0;LIMTRIV=0;")
    local_definition = SIMON_GP_FUNCTION.split(
        "/* ELKIES_R17_GP_DEFINITION_SPLIT */"
    )[0]
    pari(local_definition)
    for definition in RELATIVE_GP_FUNCTIONS.split(
        "/* MW29_RELATIVE_GP_DEFINITION_SPLIT */"
    ):
        pari(definition)

    curve = pari.ellinit([QQ(value) for value in transformed_meta["transformed_model"]])
    curve_theta = pari(transformed_meta["curve_theta_in_field"])
    started = time.monotonic()
    if provider_kind == "bnf":
        global_data = pari("ell2global_norm_envelope")(curve, bnf, 1, curve_theta)
        ls2, normspace, bad_primes_raw = global_data
        norm_columns = binary_columns(normspace)
        ls2_dimension = len(ls2)
        norm_dimension = len(norm_columns)
    else:
        assert hecke_global is not None
        field_polynomial = pari(HECKE_REDUCED_CUBIC)
        alphas = []
        for row in hecke_global["basis"]:
            c0, c1, c2 = (QQ(value) for value in row["coefficients"])
            alphas.append(pari.Mod(c0 + c1 * pari("x") + c2 * pari("x") ** 2, field_polynomial))
        ls2 = pari(alphas)
        norm_dimension = len(alphas)
        ls2_dimension = norm_dimension
        norm_columns = [
            [int(row == column) for row in range(norm_dimension)]
            for column in range(norm_dimension)
        ]
        bad_primes_raw = [int(value) for value in record["bad_primes"]]
    print(
        f"{PROTOCOL}|curve={args.curve_id}|stage=global_envelope|status=complete"
        f"|sclass={ls2_dimension}|norm={norm_dimension}",
        flush=True,
    )

    old = None
    if args.resume and args.output.exists():
        candidate = json.loads(args.output.read_text())
        compatible = (
            candidate.get("schema") == SCHEMA
            and candidate.get("curve_id") == args.curve_id
            and candidate.get("inputs", {}).get("global_provider_sha256")
            == digest(provider_path)
            and candidate.get("inputs", {}).get("public_certificate_sha256") == digest(public_path)
            and candidate.get("inputs", {}).get(
                "classgroup_pressure_certificate_sha256"
            ) == digest(pressure_path)
        )
        if not compatible:
            raise ArithmeticError("existing resume checkpoint has incompatible inputs")
        old = candidate

    if old is not None and old.get("known_mw", {}).get(
        "global_square_verification_completed"
    ) is True:
        ls2_point_rows = old["known_mw"]["rows_in_global_s_squareclasses"]
        known_norm_rows = old["known_mw"]["rows_in_norm_envelope"]
        auxiliary_primes = old["known_mw"]["auxiliary_primes"]
    else:
        ls2_point_rows, known_norm_rows, auxiliary_primes = point_coordinates(
            nf=nf,
            field_polynomial=pari(transformed_meta["field_cubic"]),
            curve_theta=curve_theta,
            ls2=ls2,
            norm_columns=norm_columns,
            point_values=point_values,
            meta=transformed_meta,
            auxiliary_prime_bound=args.auxiliary_prime_bound,
        )
    if len(known_norm_rows) != 29 or f2_rank_rows(known_norm_rows, norm_dimension) != 29:
        raise ArithmeticError("MW29 quotient gate lost a known Kummer direction")
    print(
        f"{PROTOCOL}|curve={args.curve_id}|stage=mw_quotient|status=complete"
        f"|known=29|ambient_quotient={norm_dimension - 29}",
        flush=True,
    )

    default_order = prioritized_places(
        args.curve_id, [int(value) for value in bad_primes_raw]
    )
    if args.place_order:
        requested = [
            -1 if token.strip().lower() == "infinity" else int(token)
            for token in args.place_order.split(",")
            if token.strip()
        ]
        if len(set(requested)) != len(requested):
            raise ArithmeticError("--place-order contains duplicates")
        invalid = sorted(set(requested) - set(default_order))
        if invalid:
            raise ArithmeticError(f"--place-order contains irrelevant places {invalid}")
        place_order = [*requested, *(place for place in default_order if place not in requested)]
    else:
        place_order = default_order

    completed_places: list[dict[str, Any]] = []
    if old is not None:
        old_places = old.get("completed_local_places", [])
        if isinstance(old_places, list):
            completed_places = [
                place for place in old_places
                if isinstance(place, dict) and int(place.get("place_integer", -2)) in place_order
            ]
    completed_by_integer = {
        int(place["place_integer"]): place for place in completed_places
    }

    provider_inputs = {
        "global_provider_kind": provider_kind,
        "global_provider": str(provider_path),
        "global_provider_sha256": digest(provider_path),
    }
    if provider_kind == "bnf":
        assert checkpoint_path is not None
        provider_inputs.update(
            {
                "bnf_checkpoint": str(checkpoint_path),
                "bnf_checkpoint_sha256": digest(checkpoint_path),
            }
        )
    else:
        provider_inputs.update(
            {
                "hecke_provider_script": str(HECKE_PROVIDER),
                "hecke_provider_script_sha256": digest(HECKE_PROVIDER),
            }
        )

    def document(status: str, certificate: dict[str, Any]) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "protocol": PROTOCOL,
            "status": status,
            "curve_id": args.curve_id,
            "operational_role": "post-discovery MW29-relative closure only",
            "prospective_mw17_candidate_gate_evidence": False,
            "inputs": {
                "public_certificate": str(public_path),
                "public_certificate_sha256": digest(public_path),
                "classgroup_pressure_certificate": str(pressure_path),
                "classgroup_pressure_certificate_sha256": digest(pressure_path),
                **provider_inputs,
            },
            "global_envelope": {
                "provider": (
                    "certified PARI BNF + Simon bnfpSelmer/kernorm"
                    if provider_kind == "bnf"
                    else (
                        "Hecke 0.37.6 class/unit groups with GRH=false; "
                        "Hecke pselmer_group norm kernel"
                    )
                ),
                "global_s_squareclass_dimension": ls2_dimension,
                "norm_square_envelope_dimension": norm_dimension,
                "ls2_generators": [str(value) for value in ls2],
                "norm_basis_columns_in_global_s_squareclasses": norm_columns,
            },
            "known_mw": {
                "labels": [f"P{index}" for index in range(1, 30)],
                "rows_in_global_s_squareclasses": ls2_point_rows,
                "rows_in_norm_envelope": known_norm_rows,
                "rank": 29,
                "certified_residual_selmer_dimension_parity": residual_parity,
                "auxiliary_primes": auxiliary_primes,
                "global_square_verification_completed": True,
            },
            "requested_place_order": [
                "infinity" if place == -1 else str(place) for place in place_order
            ],
            "completed_local_places": [
                completed_by_integer[place]
                for place in place_order
                if place in completed_by_integer
            ],
            "ambient_manifest": certificate["ambient_manifest"],
            "relative_certificate": certificate["relative_certificate"],
            "elapsed_seconds": time.monotonic() - started,
            "claim_boundary": (
                "The zero-excess status is an unconditional relative Selmer upper "
                "bound from the certified global envelope, exact local subset, and "
                "certified parity sharpening when used. "
                "A nonzero kernel is exact only after every relevant place is complete. "
                "Because P18..P29 are quotient inputs, this result cannot calibrate a "
                "prospective MW17-only candidate gate."
            ),
        }

    def rebuild_certificate(all_places_complete: bool) -> dict[str, Any]:
        places = [completed_by_integer[place] for place in place_order if place in completed_by_integer]
        ambient_manifest = {
            "schema": "elliptic-curves.mw29-relative-2selmer-ambient.v1",
            "case_id": f"record-r29-{args.curve_id}",
            "known_mw_target_rank": 29,
            "ambient_norm_square_dimension": norm_dimension,
            "residual_selmer_dimension_parity": residual_parity,
            "known_mw_rows": [
                {"label": f"P{index + 1}", "row": row}
                for index, row in enumerate(known_norm_rows)
            ],
            "places": [
                {
                    "place": place["place"],
                    "allowed_subspace_basis": place[
                        "allowed_subspace_basis_in_norm_envelope"
                    ],
                }
                for place in places
            ],
            "certification": {
                "method": (
                    (
                        "PARI bnfcertify; Simon bnfpSelmer norm envelope"
                        if provider_kind == "bnf"
                        else (
                            "Hecke 0.37.6 class_group/unit_group with GRH=false; "
                            "Hecke pselmer_group norm kernel; PARI nfcertify"
                        )
                    )
                    + "; exact nfissquare Kummer coordinates; "
                    "full-dimensional Simon local images"
                ),
                "hypothesis": None,
                "global_ambient_upper_envelope_certified": True,
                "global_ambient_exact": True,
                "norm_condition_incorporated": True,
                "known_mw_kummer_coordinates_certified": True,
                "supplied_local_conditions_certified": True,
                "supplied_subspaces_are_necessary_selmer_conditions": True,
                "all_required_local_conditions_complete": all_places_complete,
                "residual_dimension_parity_certified": True,
            },
        }
        relative = dict(
            build_certificate(
                ambient_manifest,
                maximum_cut_size=args.maximum_cut_size,
                maximum_cut_subsets=args.maximum_cut_subsets,
            )
        )
        return {"ambient_manifest": ambient_manifest, "relative_certificate": relative}

    processed_this_run = 0
    initial = rebuild_certificate(len(completed_by_integer) == len(place_order))
    if initial["relative_certificate"]["status"].startswith(
        "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO"
    ):
        final = document("CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO", initial)
        atomic_json(args.output, final)
        print(f"{PROTOCOL}|curve={args.curve_id}|stage=complete|status={final['status']}", flush=True)
        return

    for place in place_order:
        if place in completed_by_integer:
            continue
        if args.stop_after_places is not None and processed_this_run >= args.stop_after_places:
            break
        print(
            f"{PROTOCOL}|curve={args.curve_id}|stage=local|place="
            f"{'infinity' if place == -1 else place}|status=start",
            flush=True,
        )
        local_raw = (
            pari("ell2allowed_at_place")(
                curve, bnf, 1, curve_theta, ls2, place
            )
            if provider_kind == "bnf"
            else pari("ell2allowed_at_place_nf")(
                curve, nf, 1, curve_theta, ls2, place
            )
        )
        allowed_s = binary_columns(local_raw[0])
        allowed_norm = allowed_in_norm_coordinates(allowed_s, norm_columns)
        completed_by_integer[place] = {
            "place": "infinity" if place == -1 else str(place),
            "place_integer": place,
            "allowed_subspace_basis_in_norm_envelope": allowed_norm,
            "allowed_subspace_dimension_in_norm_envelope": len(allowed_norm),
            "allowed_subspace_basis_columns_in_global_s_squareclasses": allowed_s,
            "ambient_local_kummer_dimension": int(local_raw[1]),
            "computed_local_kummer_image_dimension": int(local_raw[2]),
            "localized_global_s_squareclass_image_dimension": int(local_raw[3]),
        }
        processed_this_run += 1
        complete_coverage = len(completed_by_integer) == len(place_order)
        certificate = rebuild_certificate(complete_coverage)
        raw_kernel = certificate["relative_certificate"]["relative_local_matrix"][
            "unexplained_selmer_excess_kernel_dimension"
        ]
        residual_upper = certificate["relative_certificate"]["relative_selmer_bound"][
            "parity_sharpened_upper_bound"
        ]
        status = (
            "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO"
            if residual_upper == 0
            else (
                "COMPLETE_RELATIVE_2SELMER_KERNEL"
                if complete_coverage
                else "INCREMENTAL_RELATIVE_2SELMER_UPPER_BOUND"
            )
        )
        atomic_json(args.output, document(status, certificate))
        print(
            f"{PROTOCOL}|curve={args.curve_id}|stage=local|place="
            f"{'infinity' if place == -1 else place}|status=complete"
            f"|raw_kernel={raw_kernel}|residual_upper={residual_upper}",
            flush=True,
        )
        if residual_upper == 0:
            break

    complete_coverage = len(completed_by_integer) == len(place_order)
    certificate = rebuild_certificate(complete_coverage)
    raw_kernel = certificate["relative_certificate"]["relative_local_matrix"][
        "unexplained_selmer_excess_kernel_dimension"
    ]
    residual_upper = certificate["relative_certificate"]["relative_selmer_bound"][
        "parity_sharpened_upper_bound"
    ]
    status = (
        "CERTIFIED_RELATIVE_2SELMER_QUOTIENT_ZERO"
        if residual_upper == 0
        else (
            "COMPLETE_RELATIVE_2SELMER_KERNEL"
            if complete_coverage
            else "INCREMENTAL_RELATIVE_2SELMER_UPPER_BOUND"
        )
    )
    final = document(status, certificate)
    atomic_json(args.output, final)
    print(
        f"{PROTOCOL}|curve={args.curve_id}|stage=complete|status={status}"
        f"|raw_kernel={raw_kernel}|residual_upper={residual_upper}"
        f"|places={len(completed_by_integer)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
