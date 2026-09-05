#!/usr/bin/env python3
"""Checkpointed open-source R17 2-descent via PARI and Simon's GP routines.

This route separates the expensive cubic BNF from the Selmer computation.
A fully certified BNF is persisted with PARI ``writebin`` and may be reused by
later local-condition, Kummer-embedding, and cover-construction workers.  It
therefore avoids repeating monolithic ``ellrankinit`` after the global field
stage has once completed.

The descent worker is blind: it receives only the curve and certified BNF.
The generic MW17 points are introduced by a later worker to align a quotient
basis.  A third worker then enumerates quotient representatives and
constructs/searches their explicit intersections of quadrics.  Held-out
exceptional control points are introduced only after that search has ended.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]

from run_elkies_2026_relative_2selmer_open import (  # noqa: E402
    GENERIC_RANK,
    factor_hints,
    f2_linear_combination_coefficients,
    f2_rank,
    file_sha256,
    load_authoritative_cases,
    query_open_source_versions,
    selected_manifest_cases,
    supervise_source,
)


INPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-input.v1"
OUTPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-checkpointed-run.v3"
PROTOCOL = "ELKIESR17CHECKREL2"
DEFAULT_SAGE_PYTHON = Path(
    "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python"
)


BNF_WORKER = r'''from __future__ import annotations
from hashlib import sha256
import json
import os
from pathlib import Path
import time

from sage.all import EllipticCurve, QQ, pari
from sage.version import version as sage_version

payload = json.loads(Path(INPUT_PATH).read_text())
pari.allocatemem(int(payload["pari_stack_bytes"]))
pari.default("debug", int(payload["pari_debug"]))
hints = [int(value) for value in payload.get("factor_hint_primes", [])]
if hints:
    pari.addprimes(hints)

def stage(name, status, **fields):
    suffix = "".join(f"|{key}={value}" for key, value in sorted(fields.items()))
    print(f"ELKIESR17CHECKREL2|case={payload['case_id']}|stage={name}|status={status}{suffix}", flush=True)

def transformed_model(curve):
    steps = []
    if curve[0] != 0 or curve[2] != 0:
        change = pari([1, 0, -curve[0]/2, -curve[2]/2])
        curve = pari.ellchangecurve(curve, change)
        steps.append([str(value) for value in change])
    while True:
        denominator = int(pari.denominator(pari([curve[1], curve[3], curve[4]])))
        if denominator == 1:
            break
        factorization = pari.factor(denominator)
        radical = 1
        for index in range(len(factorization[0])):
            radical *= int(factorization[0][index])
        change = pari([QQ(1)/radical, 0, 0, 0])
        curve = pari.ellchangecurve(curve, change)
        steps.append([str(value) for value in change])
    return curve, steps

curve = pari(EllipticCurve(QQ, [QQ(value) for value in payload["global_minimal_model"]]))
curve, changes = transformed_model(curve)
transformed = [str(curve[index]) for index in range(5)]
if transformed[0] != "0" or transformed[2] != "0":
    raise ArithmeticError("failed to remove a1 and a3")
curve_cubic = pari(f"x^3+({curve[1]})*x^2+({curve[3]})*x+({curve[4]})")
reduction = pari.polredbest(curve_cubic, 1)
cubic, curve_theta = reduction[0], reduction[1]
field_generator_candidates = pari.nfisisom(cubic, curve_cubic)
curve_field_theta = pari.Mod("x", curve_cubic)
field_generator_in_curve_field = None
for candidate in field_generator_candidates:
    if pari.subst(pari.lift(curve_theta), "x", candidate) == curve_field_theta:
        field_generator_in_curve_field = candidate
        break
if field_generator_in_curve_field is None:
    raise ArithmeticError("could not invert the reduced-field isomorphism")
stage("polredbest", "complete", reduced_polynomial=cubic)

stage("nfinit", "start", factor_hints=len(hints))
started = time.monotonic()
nf = pari.nfinit([cubic, hints]) if hints else pari.nfinit(cubic)
nf_seconds = time.monotonic() - started
stage("nfinit", "complete", seconds=f"{nf_seconds:.6f}")
stage("nfcertify", "start")
started = time.monotonic()
obstructions = list(pari.nfcertify(nf))
if obstructions:
    raise ArithmeticError(f"maximal-order certification failed: {obstructions}")
nf_certify_seconds = time.monotonic() - started
stage("nfcertify", "complete", seconds=f"{nf_certify_seconds:.6f}")

stage("bnfinit", "start", flag=payload["bnf_flag"], tech=":".join(str(v) for v in payload["bnf_tech"]))
started = time.monotonic()
bnf = pari.bnfinit(nf, int(payload["bnf_flag"]), payload["bnf_tech"])
bnf_seconds = time.monotonic() - started
stage("bnfinit", "complete", seconds=f"{bnf_seconds:.6f}")
stage("bnfcertify", "start", flag=0)
started = time.monotonic()
if not bool(pari.bnfcertify(bnf)):
    raise ArithmeticError("full BNF certification failed")
bnf_certify_seconds = time.monotonic() - started
stage("bnfcertify", "complete", seconds=f"{bnf_certify_seconds:.6f}")

checkpoint = Path(payload["bnf_checkpoint"])
checkpoint.parent.mkdir(parents=True, exist_ok=True)
temporary = checkpoint.with_suffix(checkpoint.suffix + ".tmp")
temporary.unlink(missing_ok=True)
pari.writebin(str(temporary), bnf)
os.replace(temporary, checkpoint)
checkpoint_hash = sha256(checkpoint.read_bytes()).hexdigest()
result = {
    "schema": "elliptic-curves.elkies-2026-relative-2selmer-bnf-checkpoint.v2",
    "case_id": payload["case_id"],
    "global_minimal_model": payload["global_minimal_model"],
    "transformed_model": transformed,
    "point_change_sequence": changes,
    "cubic": str(curve_cubic),
    "field_cubic": str(cubic),
    "curve_theta_in_field": str(curve_theta),
    "field_generator_in_curve_field": str(field_generator_in_curve_field),
    "cubic_coefficients_ascending": [str(curve[4]), str(curve[3]), str(curve[1]), "1"],
    "field_discriminant": str(nf[2]),
    "defining_order_index": str(nf[3]),
    "field_signature": [int(value) for value in nf.nf_get_sign()],
    "class_number": str(bnf.bnf_get_no()),
    "class_group_cyclic_invariants": [str(value) for value in bnf.bnf_get_cyc()],
    "fundamental_unit_count": (
        0 if str(bnf.bnf_get_fu()) == "0" else len(bnf.bnf_get_fu())
    ),
    "bnf_flag": int(payload["bnf_flag"]),
    "nfcertify_completed": True,
    "bnfcertify_completed": True,
    "factor_hint_primes": [str(value) for value in hints],
    "bnf_checkpoint": str(checkpoint),
    "bnf_checkpoint_sha256": checkpoint_hash,
    "sage_version": str(sage_version),
    "pari_version": str(pari.version()),
    "timings": {
        "nfinit_seconds": nf_seconds,
        "nfcertify_seconds": nf_certify_seconds,
        "bnfinit_seconds": bnf_seconds,
        "bnfcertify_seconds": bnf_certify_seconds,
    },
}
Path(OUTPUT_PATH).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
stage("bnf_checkpoint", "complete", output=OUTPUT_PATH, checkpoint=checkpoint)
'''


SIMON_GP_FUNCTION = r'''
elllocalimage_mapped(nf,pp,K,polrel,curve_theta) =
{
  my(p,prank,rac,pts,bound,essai,mrank,r,xx,delta,ph,delta2,local_prec,ival);
  p=pp[1][1][1]; prank=#pp-(p!=2);
  rac=polrootsmodpn(K*polrel,p);
  pts=matrix(0,0); bound=p+6; essai=0; mrank=0;
  while(mrank<prank,
    essai++;
    if(essai%16==0,pts=matimage(pts);bound*=p);
    r=random(#rac)+1; local_prec=random(rac[r][2]+3)-2;
    xx=rac[r][1]+p^local_prec*random(bound);
    delta=K*(xx-curve_theta);
    if(!psquare(K*subst(polrel,variable(polrel),xx),p),next);
    ph=[];
    for(i=1,#pp,
      ph=concat(ph,[ival=idealval(nf,delta,pp[i][1])]);
      delta2=delta/pp[i][2]^ival;
      if(p==2,
        ph=concat(ph,ideallog(nf,delta2,pp[i][3])~)
      , ph=concat(ph,[1-nfpsquareodd(nf,delta2,pp[i][4])])));
    pts=concat(pts,ph~*Mod(1,2));
    mrank=matrank(pts*Mod(1,2)));
  return(matimage(pts));
};

/* ELKIES_R17_GP_DEFINITION_SPLIT */

ell2selmer_basis_gen(ell,bnf,K,curve_theta) =
{
  my(A,B,C,polrel,polprime,badideal,badprimes,descentprimes,S,LS2,normspace,selmer,theta_embeddings,real_place,signs,p,pp,prank,locimage,LS2image,allowed,localspaces,localplaces,localmeta,audit,omitted,standalone,j,selmerinnorm);
  if(#ell < 13,ell=ellinit(ell));
  if(ell.a1 != 0 || ell.a3 != 0,error("ell2selmer_basis_gen: nonzero a1/a3"));
  A=ell.a2; B=ell.a4; C=ell.a6;
  polrel=Pol([1,A,B,C]); polprime=polrel';
  polprime=subst(polprime,variable(polprime),curve_theta);
  badideal=abs(K*idealadd(bnf,polprime,bnf.index));
  S=bnfpSelmer(bnf,badideal,2); LS2=S[1]; S=S[2];
  normspace=kernorm(LS2,vector(#S,i,S[i].p),2);
  selmer=normspace;
  localspaces=List();localplaces=List();localmeta=List();
  if(bnf.r1==3,
    theta_embeddings=nfeltembed(bnf.nf,curve_theta);real_place=1;
    for(i=2,#theta_embeddings,if(theta_embeddings[i]<theta_embeddings[real_place],real_place=i));
    signs=vector(#LS2,i,nfeltsign(bnf.nf,LS2[i],real_place)<0);
    allowed=matker(Mat(signs*Mod(1,2)))*Mod(1,2);
    listput(localspaces,allowed);listput(localplaces,-1);listput(localmeta,[1,1,matrank(Mat(signs*Mod(1,2)))]);
    selmer=matintersect(selmer,allowed)*Mod(1,2)
  ,
    allowed=matid(#LS2)*Mod(1,2);
    listput(localspaces,allowed);listput(localplaces,-1);listput(localmeta,[0,0,0])
  );
  /* Detector v2 records every bad rational prime, even when Simon's
     smaller S-support makes that local condition redundant. */
  descentprimes=factorint(badideal[1,1]*2)[,1];
  badprimes=factorint(badideal[1,1]*2*numerator(ell.disc))[,1];
  for(i=1,#badprimes,
    p=badprimes[i];
    print("ELKIESR17CHECKREL2|stage=local_condition|status=start|prime=",p);
    pp=ppinit(bnf.nf,p);
    prank=#pp-(p!=2);
    locimage=elllocalimage_mapped(bnf.nf,pp,K,polrel,curve_theta);
    LS2image=LS2localimage(bnf.nf,LS2,pp);
    locimage=matintersect(LS2image,locimage);
    allowed=concat(matker(LS2image),matinverseimage(LS2image,locimage)*Mod(1,2));
    allowed=matimage(allowed*Mod(1,2));
    listput(localspaces,allowed);listput(localplaces,p);listput(localmeta,[prank,#locimage,#LS2image]);
    selmer=matintersect(selmer,allowed);
    selmer=matimage(selmer*Mod(1,2));
    print("ELKIESR17CHECKREL2|stage=local_condition|status=complete|prime=",p,"|dimension=",#selmer);
  );
  audit=List();
  for(i=1,#localspaces,
    omitted=normspace;
    for(j=1,#localspaces,if(i!=j,omitted=matintersect(omitted,localspaces[j])*Mod(1,2)));
    standalone=matintersect(normspace,localspaces[i])*Mod(1,2);
    listput(audit,[localplaces[i],#localspaces[i],#standalone,#omitted,#normspace-#omitted]);
  );
  selmerinnorm=matinverseimage(normspace,selmer)*Mod(1,2);
  return([LS2,lift(selmer),badprimes,Vec(audit),#LS2,#normspace,#normspace-#selmer,Vec(localspaces),Vec(localplaces),lift(normspace),Vec(localmeta),lift(selmerinnorm),descentprimes]);
};
'''


SELMER_WORKER = r'''from __future__ import annotations
import json
from pathlib import Path
import sys
import time

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, pari
from sage.env import SAGE_EXTCODE
from sage.misc.randstate import set_random_seed

payload = json.loads(Path(INPUT_PATH).read_text())
sys.path.insert(0, payload["cas_directory"])
from build_bnf_free_two_covers import cover_for, multiply_mod_cubic
meta = json.loads(Path(payload["bnf_metadata"]).read_text())
pari.allocatemem(int(payload["pari_stack_bytes"]))
pari.default("debug", int(payload["pari_debug"]))
set_random_seed(int(payload["random_seed"]))

def stage(name, status, **fields):
    suffix = "".join(f"|{key}={value}" for key, value in sorted(fields.items()))
    print(f"ELKIESR17CHECKREL2|case={payload['case_id']}|stage={name}|status={status}{suffix}", flush=True)

checkpoint = Path(meta["bnf_checkpoint"])
if __import__('hashlib').sha256(checkpoint.read_bytes()).hexdigest() != meta["bnf_checkpoint_sha256"]:
    raise ArithmeticError("BNF checkpoint hash mismatch")
bnf = pari.read(str(checkpoint))
if not bool(pari.bnfcertify(bnf)):
    raise ArithmeticError("reloaded BNF failed certification")
if str(bnf.nf_get_pol()) != meta["field_cubic"]:
    raise ArithmeticError("BNF defining polynomial mismatch")

simon = Path(SAGE_EXTCODE) / "pari" / "simon"
for name in ("ellQ.gp", "ell.gp", "qfsolve.gp", "resultant3.gp"):
    pari.read(simon / name)
pari(f"DEBUGLEVEL_ell={int(payload['simon_verbose'])}; LIMBIGPRIME=0; LIM1=0; LIM3=0; LIMTRIV=0;")
for definition in SIMON_FUNCTION.split("/* ELKIES_R17_GP_DEFINITION_SPLIT */"):
    pari(definition)
curve = pari.ellinit([QQ(value) for value in meta["transformed_model"]])
global_curve = EllipticCurve(QQ, [QQ(value) for value in meta["global_minimal_model"]])
elliptic_bad_local_data = {}
for data in global_curve.local_data():
    prime = int(data.prime().gens()[0])
    elliptic_bad_local_data[str(prime)] = {
        "kodaira_symbol": str(data.kodaira_symbol()),
        "tamagawa_number": int(data.tamagawa_number()),
        "conductor_exponent": int(data.conductor_valuation()),
        "minimal_discriminant_valuation": int(data.discriminant_valuation()),
    }
elliptic_bad_primes = sorted(int(value) for value in elliptic_bad_local_data)
curve_theta = pari(meta["curve_theta_in_field"])
stage("selmer_basis", "start")
started = time.monotonic()
raw = pari("ell2selmer_basis_gen")(curve, bnf, 1, curve_theta)
elapsed = time.monotonic() - started
LS2, matrix, bad_primes = raw[0], raw[1], raw[2]
local_audit = raw[3]
global_squareclass_dimension = int(raw[4])
norm_kernel_dimension = int(raw[5])
full_local_condition_matrix_rank = int(raw[6])
local_allowed_subspaces = raw[7]
local_places = raw[8]
normspace = raw[9]
local_metadata = raw[10]
selmer_in_normspace = raw[11]
descent_support_primes = [int(value) for value in raw[12]]
finite_local_condition_primes = [int(value) for value in bad_primes]
auxiliary_descent_primes = sorted(
    set(descent_support_primes) - set(elliptic_bad_primes) - {2}
)
dimension = len(matrix)
stage("selmer_basis", "complete", seconds=f"{elapsed:.6f}", dimension=dimension)

def binary_basis_columns(value):
    return [[int(entry) for entry in column] for column in value]

if (
    len(local_audit) != len(local_allowed_subspaces)
    or len(local_audit) != len(local_places)
    or len(local_audit) != len(local_metadata)
):
    raise ArithmeticError("local-condition audit vectors have inconsistent lengths")

curve_f = pari(meta["cubic"])
field_f = pari(meta["field_cubic"])
field_generator_in_curve_field = pari(meta["field_generator_in_curve_field"])
nf = bnf
cubic_coefficients = [QQ(value) for value in meta["cubic_coefficients_ascending"]]
ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))

def square_root_q(value):
    value = QQ(value)
    if value < 0 or not ZZ(value.numerator()).is_square() or not ZZ(value.denominator()).is_square():
        return None
    return QQ(ZZ(value.numerator()).sqrt()) / QQ(ZZ(value.denominator()).sqrt())

def search_cover(alpha, bound):
    for u in range(-bound, bound + 1):
      for v in range(-bound, bound + 1):
       for w in range(-bound, bound + 1):
        if u == v == w == 0:
            continue
        constant, theta, theta2 = multiply_mod_cubic(
            alpha,
            multiply_mod_cubic([QQ(u), QQ(v), QQ(w)], [QQ(u), QQ(v), QQ(w)], cubic_coefficients),
            cubic_coefficients,
        )
        if theta2 != 0:
            continue
        z = square_root_q(-theta)
        if z is None or z == 0:
            continue
        x_value = constant / (z*z)
        y_value = square_root_q(
            cubic_coefficients[0] + cubic_coefficients[1]*x_value
            + cubic_coefficients[2]*x_value*x_value + x_value*x_value*x_value
        )
        if y_value is not None:
            return [str(QQ(u)), str(QQ(v)), str(QQ(w)), str(z)], [str(x_value), str(y_value)]
    return None, None

basis = []
for column in range(dimension):
    alpha = pari.Mod(1, field_f)
    exponent = matrix[column]
    for row in range(len(LS2)):
        if int(exponent[row]) & 1:
            alpha *= LS2[row]
    field_lifted = pari.lift(alpha)
    field_coefficients = [QQ(str(field_lifted.polcoef(degree))) for degree in range(3)]
    curve_alpha = pari.subst(field_lifted, "x", field_generator_in_curve_field)
    curve_lifted = pari.lift(curve_alpha)
    coefficients = [QQ(str(curve_lifted.polcoef(degree))) for degree in range(3)]
    norm = QQ(str(pari.nfeltnorm(nf, alpha)))
    if square_root_q(norm) is None:
        raise ArithmeticError("a Selmer basis class has nonsquare norm")
    search_started = time.monotonic()
    witness, image = search_cover(coefficients, int(payload["raw_basis_search_bound"]))
    basis.append({
        "basis_index": column + 1,
        "alpha_coefficients": [str(value) for value in coefficients],
        "field_alpha_coefficients": [str(value) for value in field_coefficients],
        "norm": str(norm),
        "norm_square_root": str(square_root_q(norm)),
        "cover": cover_for(coefficients, cubic_coefficients, ring),
        "blind_search": {
            "generic_points_supplied": False,
            "exceptional_points_supplied": False,
            "coefficient_bound": int(payload["raw_basis_search_bound"]),
            "status": "point_found" if witness else "no_point_within_bound",
            "cover_point": witness,
            "transformed_elliptic_point": image,
            "seconds": time.monotonic() - search_started,
        },
    })

result = {
    "schema": "elliptic-curves.elkies-2026-relative-2selmer-simon-basis.v3",
    "case_id": payload["case_id"],
    "global_minimal_model": meta["global_minimal_model"],
    "transformed_model": meta["transformed_model"],
    "point_change_sequence": meta["point_change_sequence"],
    "cubic": meta["cubic"],
    "field_cubic": meta["field_cubic"],
    "curve_theta_in_field": meta["curve_theta_in_field"],
    "cubic_coefficients_ascending": meta["cubic_coefficients_ascending"],
    "two_selmer_dimension": dimension,
    "bad_rational_primes": [str(value) for value in elliptic_bad_primes],
    "descent_support_rational_primes": [str(value) for value in descent_support_primes],
    "auxiliary_descent_primes": [str(value) for value in auxiliary_descent_primes],
    "finite_local_condition_primes": [
        str(value) for value in finite_local_condition_primes
    ],
    "local_condition_matrix": {
        "ambient": "the global norm-square subspace of the S-squareclass group",
        "global_s_squareclass_dimension": global_squareclass_dimension,
        "global_norm_square_subspace_dimension": norm_kernel_dimension,
        "global_norm_square_subspace_basis_columns_in_s_squareclasses": binary_basis_columns(normspace),
        "selmer_basis_columns_in_global_norm_square_subspace": binary_basis_columns(selmer_in_normspace),
        "global_norm_condition_rank": global_squareclass_dimension - norm_kernel_dimension,
        "full_finite_and_archimedean_matrix_rank_on_norm_subspace": full_local_condition_matrix_rank,
        "rank_nullity_verified": (
            full_local_condition_matrix_rank + dimension == norm_kernel_dimension
        ),
        "places": [
            {
                "place": (
                    "infinity" if int(row[0]) == -1 else str(int(row[0]))
                ),
                "allowed_subspace_basis_columns_in_global_s_squareclasses": binary_basis_columns(
                    local_allowed_subspaces[index]
                ),
                "allowed_subspace_dimension_in_global_s_squareclasses": int(row[1]),
                "norm_subspace_intersection_dimension_for_this_place_alone": int(row[2]),
                "selmer_candidate_dimension_after_deleting_this_place": int(row[3]),
                "matrix_rank_after_deleting_this_place": int(row[4]),
                "rank_drop_from_full_matrix": (
                    full_local_condition_matrix_rank - int(row[4])
                ),
                "ambient_local_kummer_dimension": int(local_metadata[index][0]),
                "computed_local_kummer_image_dimension": int(local_metadata[index][1]),
                "localized_global_s_squareclass_image_dimension": int(local_metadata[index][2]),
                "elliptic_bad_place": (
                    int(row[0]) != -1
                    and str(int(row[0])) in elliptic_bad_local_data
                ),
                "auxiliary_descent_place": (
                    int(row[0]) != -1
                    and int(row[0]) in auxiliary_descent_primes
                ),
                "component_group_data": (
                    None
                    if int(row[0]) == -1
                    else elliptic_bad_local_data.get(str(int(row[0])))
                ),
            }
            for index, row in enumerate(local_audit)
        ],
    },
    "basis": basis,
    "bnf_checkpoint_sha256": meta["bnf_checkpoint_sha256"],
    "local_conditions_completed": True,
    "worker_seconds": elapsed,
}
Path(OUTPUT_PATH).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
stage("selmer_worker", "complete", output=OUTPUT_PATH)
'''


QUOTIENT_COVER_WORKER = r'''from __future__ import annotations
import json
from pathlib import Path
import sys
import time

from sage.all import PolynomialRing, QQ, ZZ, pari

payload = json.loads(Path(INPUT_PATH).read_text())
sys.path.insert(0, payload["cas_directory"])
from build_bnf_free_two_covers import cover_for, multiply_mod_cubic
meta = json.loads(Path(payload["bnf_metadata"]).read_text())
selmer = json.loads(Path(payload["selmer_result"]).read_text())
pari.allocatemem(int(payload["pari_stack_bytes"]))
curve_f = pari(meta["cubic"])
field_f = pari(meta["field_cubic"])
nf = pari.read(meta["bnf_checkpoint"])
if not bool(pari.bnfcertify(nf)):
    raise ArithmeticError("reloaded BNF failed certification")
cubic_coefficients = [QQ(value) for value in meta["cubic_coefficients_ascending"]]
ring = PolynomialRing(QQ, names=("u", "v", "w", "z"))

def square_root_q(value):
    value = QQ(value)
    if value < 0 or not ZZ(value.numerator()).is_square() or not ZZ(value.denominator()).is_square():
        return None
    return QQ(ZZ(value.numerator()).sqrt()) / QQ(ZZ(value.denominator()).sqrt())

def search_cover(alpha, bound):
    for u in range(-bound, bound + 1):
      for v in range(-bound, bound + 1):
       for w in range(-bound, bound + 1):
        if u == v == w == 0:
            continue
        constant, theta, theta2 = multiply_mod_cubic(
            alpha,
            multiply_mod_cubic([QQ(u), QQ(v), QQ(w)], [QQ(u), QQ(v), QQ(w)], cubic_coefficients),
            cubic_coefficients,
        )
        if theta2 != 0:
            continue
        z = square_root_q(-theta)
        if z is None or z == 0:
            continue
        x_value = constant / (z*z)
        y_value = square_root_q(
            cubic_coefficients[0] + cubic_coefficients[1]*x_value
            + cubic_coefficients[2]*x_value*x_value + x_value*x_value*x_value
        )
        if y_value is not None:
            return [str(QQ(u)), str(QQ(v)), str(QQ(w)), str(z)], [str(x_value), str(y_value)]
    return None, None

basis_alphas = []
basis_field_alphas = []
for row in selmer["basis"]:
    coeff = [QQ(value) for value in row["alpha_coefficients"]]
    basis_alphas.append(pari.Mod(coeff[0] + coeff[1]*pari("x") + coeff[2]*pari("x")**2, curve_f))
    field_coeff = [QQ(value) for value in row["field_alpha_coefficients"]]
    basis_field_alphas.append(pari.Mod(field_coeff[0] + field_coeff[1]*pari("x") + field_coeff[2]*pari("x")**2, field_f))
quotient_basis = [[int(bit) for bit in row] for row in payload["quotient_basis"]]
quotient_dimension = len(quotient_basis)
total_nonzero_classes = 2**quotient_dimension - 1
constructed_count = min(total_nonzero_classes, int(payload["quotient_class_limit"]))
classes = []
for mask in range(1, constructed_count + 1):
    quotient_bits = [(mask >> index) & 1 for index in range(quotient_dimension)]
    selmer_bits = [0] * len(basis_alphas)
    for bit, row in zip(quotient_bits, quotient_basis):
        if bit:
            selmer_bits = [left ^ right for left, right in zip(selmer_bits, row)]
    alpha = pari.Mod(1, curve_f)
    field_alpha = pari.Mod(1, field_f)
    for bit, basis_alpha, basis_field_alpha in zip(selmer_bits, basis_alphas, basis_field_alphas):
        if bit:
            alpha *= basis_alpha
            field_alpha *= basis_field_alpha
    lifted = pari.lift(alpha)
    coefficients = [QQ(str(lifted.polcoef(degree))) for degree in range(3)]
    norm = QQ(str(pari.nfeltnorm(nf, field_alpha)))
    norm_root = square_root_q(norm)
    if norm_root is None:
        raise ArithmeticError("a quotient representative has nonsquare norm")
    search_started = time.monotonic()
    witness, image = search_cover(coefficients, int(payload["search_bound"]))
    classes.append({
        "quotient_class_integer": mask,
        "quotient_bits": quotient_bits,
        "selmer_bits": selmer_bits,
        "alpha_coefficients": [str(value) for value in coefficients],
        "norm": str(norm),
        "norm_square_root": str(norm_root),
        "cover": cover_for(coefficients, cubic_coefficients, ring),
        "blind_search": {
            "generic_points_supplied": True,
            "exceptional_points_supplied": False,
            "coefficient_bound": int(payload["search_bound"]),
            "status": "point_found" if witness else "no_point_within_bound",
            "cover_point": witness,
            "transformed_elliptic_point": image,
            "seconds": time.monotonic() - search_started,
        },
    })

result = {
    "schema": "elliptic-curves.elkies-2026-relative-2selmer-quotient-covers.v1",
    "case_id": payload["case_id"],
    "quotient_dimension": quotient_dimension,
    "quotient_basis": quotient_basis,
    "nonzero_quotient_class_count": total_nonzero_classes,
    "constructed_class_count": constructed_count,
    "enumeration_complete": constructed_count == total_nonzero_classes,
    "quotient_class_limit": int(payload["quotient_class_limit"]),
    "exceptional_points_supplied": False,
    "classes": classes,
}
Path(OUTPUT_PATH).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(f"ELKIESR17CHECKREL2|case={payload['case_id']}|stage=quotient_covers|status=complete|constructed={constructed_count}|total={total_nonzero_classes}", flush=True)
'''


EMBED_WORKER = r'''from __future__ import annotations
import json
from pathlib import Path
import sys

from sage.all import QQ, pari, prime_range
from sage.env import SAGE_EXTCODE

payload = json.loads(Path(INPUT_PATH).read_text())
sys.path.insert(0, payload["cas_directory"])
from run_fermigier_rank20_auxiliary_fingerprints import prime_local_rows
meta = json.loads(Path(payload["bnf_metadata"]).read_text())
selmer = json.loads(Path(payload["selmer_result"]).read_text())
pari.allocatemem(int(payload["pari_stack_bytes"]))
bnf = pari.read(meta["bnf_checkpoint"])
if not bool(pari.bnfcertify(bnf)):
    raise ArithmeticError("reloaded BNF failed certification")
nf = bnf
f = pari(meta["field_cubic"])
theta = pari(meta["curve_theta_in_field"])
simon = Path(SAGE_EXTCODE) / "pari" / "simon"
for name in ("ellQ.gp", "ell.gp", "qfsolve.gp", "resultant3.gp"):
    pari.read(simon / name)

def rank(rows):
    pivots = {}
    for row in rows:
        value = sum((int(bit)&1)<<index for index,bit in enumerate(row))
        while value:
            pivot=value.bit_length()-1
            if pivot in pivots:value^=pivots[pivot]
            else:pivots[pivot]=value;break
    return len(pivots)

def solve(vectors,target):
    width=len(target); basis={}
    for index,vector in enumerate(vectors):
        value=sum((int(bit)&1)<<column for column,bit in enumerate(vector)); combination=1<<index
        while value:
            pivot=value.bit_length()-1
            if pivot in basis:value^=basis[pivot][0];combination^=basis[pivot][1]
            else:basis[pivot]=(value,combination);break
    value=sum((int(bit)&1)<<column for column,bit in enumerate(target)); combination=0
    while value:
        pivot=value.bit_length()-1
        if pivot not in basis:return None
        value^=basis[pivot][0];combination^=basis[pivot][1]
    return [(combination>>index)&1 for index in range(len(vectors))]

def transform_point(values):
    point = pari([QQ(values[0]), QQ(values[1])])
    for change in meta["point_change_sequence"]:
        point = pari.ellchangepoint(point, pari([QQ(value) for value in change]))
    return point

basis_alphas = []
for row in selmer["basis"]:
    coeff = [QQ(value) for value in row["field_alpha_coefficients"]]
    basis_alphas.append(pari.Mod(coeff[0] + coeff[1]*pari("x") + coeff[2]*pari("x")**2, f))
points = [transform_point(row) for row in payload["points"]]
point_alphas = [pari.Mod(point[0], f) - theta for point in points]
all_alphas = basis_alphas + point_alphas
signatures = [[] for _ in all_alphas]
prime_records = []
coordinates = [None] * len(points)
bad_discriminant = int(pari.poldisc(f))
nfissquare = pari("nfissquare")
for q in prime_range(3, int(payload["auxiliary_prime_bound"]) + 1):
    q = int(q)
    if bad_discriminant % q == 0:
        continue
    local, places = prime_local_rows(pari, nf, all_alphas, q)
    for index, row in enumerate(local):
        signatures[index].extend(int(bit) for bit in row)
    basis_signatures = signatures[:len(basis_alphas)]
    basis_rank = rank(basis_signatures)
    prime_records.append({"prime": q, "place_count": len(places), "basis_signature_rank": basis_rank})
    if basis_rank < len(basis_alphas):
        continue
    for point_index, alpha in enumerate(point_alphas):
        if coordinates[point_index] is not None:
            continue
        candidate = solve(basis_signatures, signatures[len(basis_alphas)+point_index])
        if candidate is None:
            continue
        ratio = alpha
        for bit, basis_alpha in zip(candidate, basis_alphas):
            if bit:
                ratio /= basis_alpha
        if bool(nfissquare(nf, ratio)):
            coordinates[point_index] = candidate
    if all(row is not None for row in coordinates):
        break
if any(row is None for row in coordinates):
    raise ArithmeticError("auxiliary signatures did not embed every supplied point")

result = {
    "schema": "elliptic-curves.elkies-2026-relative-2selmer-point-embedding.v1",
    "case_id": payload["case_id"],
    "role": payload["role"],
    "point_count": len(points),
    "point_selmer_rows": coordinates,
    "basis_signature_rank": rank(signatures[:len(basis_alphas)]),
    "auxiliary_primes": prime_records,
    "global_square_verification_completed": True,
    "public_exceptional_points_supplied": payload["role"] == "held-out-exceptional",
}
Path(OUTPUT_PATH).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(f"ELKIESR17CHECKREL2|case={payload['case_id']}|stage=point_embedding|status=complete|role={payload['role']}|points={len(points)}",flush=True)
'''


def validate_bnf_cache(
    meta_path: Path, checkpoint_path: Path, case: Any, bnf_flag: int
) -> bool:
    if not meta_path.exists() or not checkpoint_path.exists():
        return False
    try:
        meta = json.loads(meta_path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return (
        meta.get("schema")
        == "elliptic-curves.elkies-2026-relative-2selmer-bnf-checkpoint.v2"
        and meta.get("case_id") == case.case_id
        and meta.get("global_minimal_model") == [str(value) for value in case.model]
        and meta.get("bnf_flag") == bnf_flag
        and meta.get("bnfcertify_completed") is True
        and isinstance(meta.get("field_cubic"), str)
        and isinstance(meta.get("curve_theta_in_field"), str)
        and isinstance(meta.get("field_generator_in_curve_field"), str)
        and meta.get("bnf_checkpoint") == str(checkpoint_path)
        and meta.get("bnf_checkpoint_sha256") == file_sha256(checkpoint_path)
    )


def xor_rows(rows: list[list[int]]) -> list[int]:
    if not rows:
        return []
    result = [0] * len(rows[0])
    for row in rows:
        result = [left ^ int(right) for left, right in zip(result, row)]
    return result


def extend_standard(initial: list[list[int]], width: int) -> list[list[int]]:
    if f2_rank(initial) != len(initial):
        raise ArithmeticError("known point Selmer rows are not independent")
    rows = list(initial)
    extension = []
    current = len(initial)
    for index in range(width):
        unit = [1 if column == index else 0 for column in range(width)]
        if f2_rank([*rows, unit]) > current:
            rows.append(unit)
            extension.append(unit)
            current += 1
        if current == width:
            break
    if current != width:
        raise ArithmeticError("failed to extend the generic image to the Selmer basis")
    return extension


def combine_results(
    case: Any,
    selmer: dict[str, Any],
    generic: dict[str, Any],
    exceptional: dict[str, Any] | None,
    quotient_covers: dict[str, Any] | None = None,
) -> dict[str, Any]:
    dimension = int(selmer["two_selmer_dimension"])
    generic_rows = [[int(bit) for bit in row] for row in generic["point_selmer_rows"]]
    if len(generic_rows) != GENERIC_RANK or f2_rank(generic_rows) != GENERIC_RANK:
        raise ArithmeticError("generic MW17 image is not full rank")
    quotient_basis = extend_standard(generic_rows, dimension)
    aligned = [*generic_rows, *quotient_basis]
    exceptional_rows = [] if exceptional is None else [
        [int(bit) for bit in row] for row in exceptional["point_selmer_rows"]
    ]
    exceptional_quotient_rows = []
    for row in exceptional_rows:
        coordinates = f2_linear_combination_coefficients(aligned, row)
        if coordinates is None:
            raise ArithmeticError("an exceptional point missed the full Selmer basis")
        exceptional_quotient_rows.append(coordinates[GENERIC_RANK:])
    exceptional_rank = f2_rank(exceptional_quotient_rows)
    expected_exceptional = len(case.exceptional_points)
    if exceptional_rank != expected_exceptional:
        raise ArithmeticError("held-out exceptional quotient rank changed")
    residual = dimension - GENERIC_RANK
    if residual < expected_exceptional:
        raise ArithmeticError("Selmer dimension contradicts the known subgroup")
    rigid_character_quotient = None
    declared_rigid_rows = tuple(getattr(case, "rigid_quotient_rows", ()))
    if declared_rigid_rows:
        if any(len(row) != expected_exceptional for row in declared_rigid_rows):
            raise ArithmeticError("a rigid row has the wrong displayed-quotient width")
        rigid_residual_rows = [
            xor_rows(
                [
                    point_row
                    for bit, point_row in zip(coefficients, exceptional_quotient_rows)
                    if int(bit) & 1
                ]
            )
            for coefficients in declared_rigid_rows
        ]
        rigid_rank = f2_rank(rigid_residual_rows)
        if rigid_rank != len(declared_rigid_rows):
            raise ArithmeticError("the rigid-character image lost rank in the Selmer quotient")
        after_rigid_basis = extend_standard(rigid_residual_rows, residual)
        rigid_aligned_basis = [*rigid_residual_rows, *after_rigid_basis]
        displayed_rows_after_rigid = []
        for row in exceptional_quotient_rows:
            coordinates = f2_linear_combination_coefficients(rigid_aligned_basis, row)
            if coordinates is None:
                raise ArithmeticError("a displayed point missed the rigid-aligned basis")
            displayed_rows_after_rigid.append(coordinates[rigid_rank:])
        displayed_after_rigid_rank = f2_rank(displayed_rows_after_rigid)
        expected_nonrigid_rank = expected_exceptional - rigid_rank
        if displayed_after_rigid_rank != expected_nonrigid_rank:
            raise ArithmeticError("the displayed nonrigid quotient rank changed")
        rigid_character_quotient = {
            "rigid_direction_labels": list(
                getattr(case, "rigid_direction_labels", ())
            ),
            "rigid_rows_in_displayed_P18_through_P29_coordinates": [
                list(map(int, row)) for row in declared_rigid_rows
            ],
            "rigid_rows_in_residual_selmer_coordinates": rigid_residual_rows,
            "rigid_image_dimension": rigid_rank,
            "dimension_after_quotienting_rigid_plane": residual - rigid_rank,
            "basis_after_quotienting_rigid_plane_in_residual_selmer_coordinates": after_rigid_basis,
            "displayed_point_rows_after_quotienting_rigid_plane": displayed_rows_after_rigid,
            "displayed_nonrigid_image_dimension": displayed_after_rigid_rank,
            "canonical_coordinate_complement_point_labels": list(
                getattr(case, "rigid_complement_point_labels", ())
            ),
            "additional_dimension_beyond_all_twenty_nine_known_points": (
                residual - expected_exceptional
            ),
            "interpretation": (
                "When the residual Selmer dimension is twelve, this quotient is the "
                "ten-dimensional displayed nonrigid block. Any larger dimension is "
                "additional Selmer space beyond the twenty-nine known points."
            ),
        }
    quotient_basis_records = [
        {
            "quotient_basis_index": index + 1,
            "selmer_bits": row,
            "explicit_basis_cover_index": row.index(1) + 1,
        }
        for index, row in enumerate(quotient_basis)
    ]
    cover_rows = [] if quotient_covers is None else quotient_covers["classes"]
    if quotient_covers is not None and quotient_covers["quotient_basis"] != quotient_basis:
        raise ArithmeticError("quotient-cover worker used a different quotient basis")
    point_masks = [
        sum((int(bit) & 1) << index for index, bit in enumerate(row))
        for row in exceptional_quotient_rows
    ]
    realized_masks = {0}
    for point_mask in point_masks:
        realized_masks |= {value ^ point_mask for value in tuple(realized_masks)}
    labeled_covers = []
    for cover in cover_rows:
        class_integer = int(cover["quotient_class_integer"])
        labeled_covers.append(
            {
                "quotient_class_integer": class_integer,
                "quotient_bits": cover["quotient_bits"],
                "known_exceptional_point_indices": [
                    index + 1 for index, point_mask in enumerate(point_masks)
                    if point_mask == class_integer
                ],
                "known_exceptional_subgroup_realizes_class": class_integer in realized_masks,
                "blind_search_status": cover["blind_search"]["status"],
            }
        )
    return {
        "total_two_selmer_dimension": dimension,
        "generic_image_dimension": GENERIC_RANK,
        "quotient_dimension": residual,
        "quotient_basis": quotient_basis_records,
        "exceptional_selmer_rows": exceptional_rows,
        "exceptional_quotient_rows": exceptional_quotient_rows,
        "exceptional_quotient_rank": exceptional_rank,
        "rigid_character_quotient": rigid_character_quotient,
        "known_realized_quotient_classes_including_zero": 2**exceptional_rank,
        "classes_not_realized_by_known_exceptional_subgroup": 2**residual - 2**exceptional_rank,
        "unexplained_quotient_dimension": residual - exceptional_rank,
        "quotient_cover_enumeration_complete": None if quotient_covers is None else quotient_covers["enumeration_complete"],
        "quotient_cover_classification": labeled_covers,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE_PYTHON)
    parser.add_argument("--bnf-timeout", type=float, default=3600)
    parser.add_argument("--selmer-timeout", type=float, default=3600)
    parser.add_argument("--embedding-timeout", type=float, default=600)
    parser.add_argument("--cover-timeout", type=float, default=3600)
    parser.add_argument("--rss-limit-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--pari-stack-bytes", type=int, default=4_000_000_000)
    parser.add_argument("--bnf-tech", default="default")
    parser.add_argument("--bnf-flag", type=int, choices=(0, 1), default=1)
    parser.add_argument("--pari-debug", type=int, default=1)
    parser.add_argument("--simon-verbose", type=int, default=1)
    parser.add_argument("--random-seed", type=int, default=20260902)
    parser.add_argument("--raw-basis-search-bound", type=int, default=8)
    parser.add_argument("--quotient-search-bound", type=int, default=8)
    parser.add_argument("--quotient-class-limit", type=int, default=4095)
    parser.add_argument("--auxiliary-prime-bound", type=int, default=5000)
    parser.add_argument("--rebuild-bnf", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.bnf_tech == "default":
        bnf_tech = []
    else:
        try:
            parts = args.bnf_tech.split(",")
            if len(parts) != 3:
                raise ValueError
            bnf_tech = [float(parts[0]), float(parts[1]), int(parts[2])]
        except ValueError:
            parser.error("--bnf-tech must be 'default' or c1,c2,nrpid")
    if min(args.bnf_timeout, args.selmer_timeout, args.embedding_timeout, args.cover_timeout) <= 0:
        parser.error("worker timeouts must be positive")
    if args.quotient_class_limit <= 0:
        parser.error("--quotient-class-limit must be positive")
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != INPUT_SCHEMA:
        raise SystemExit("unexpected input manifest")
    sage_python = shutil.which(str(args.sage_python))
    if sage_python is None:
        raise SystemExit(f"Sage Python unavailable: {args.sage_python}")
    versions = query_open_source_versions(sage_python)
    authoritative = load_authoritative_cases()
    selected = selected_manifest_cases(manifest, set(args.case), args.controls_only)
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    runs = []
    for manifest_case in selected:
        case = authoritative[manifest_case["case_id"]]
        case_dir = args.cache_dir / case.case_id
        case_dir.mkdir(parents=True, exist_ok=True)
        bnf_meta = case_dir / "bnf.json"
        bnf_checkpoint = case_dir / "bnf.bin"
        if args.rebuild_bnf:
            bnf_meta.unlink(missing_ok=True)
            bnf_checkpoint.unlink(missing_ok=True)
        record: dict[str, Any] = {"case_id": case.case_id, "parameter": case.parameter, "role": case.role}
        if validate_bnf_cache(bnf_meta, bnf_checkpoint, case, args.bnf_flag):
            record["bnf"] = {"outcome": "reused_certified_checkpoint", "metadata": str(bnf_meta), "checkpoint": str(bnf_checkpoint)}
        else:
            bnf_meta.unlink(missing_ok=True)
            bnf_checkpoint.unlink(missing_ok=True)
            payload = {
                "case_id": case.case_id,
                "global_minimal_model": [str(value) for value in case.model],
                "factor_hint_primes": list(factor_hints(case.case_id)),
                "pari_stack_bytes": args.pari_stack_bytes,
                "pari_debug": args.pari_debug,
                "bnf_tech": bnf_tech,
                "bnf_flag": args.bnf_flag,
                "bnf_checkpoint": str(bnf_checkpoint),
            }
            record["bnf"] = supervise_source(
                sage_python, BNF_WORKER, payload, bnf_meta, case_dir / "bnf.log",
                timeout=args.bnf_timeout, rss_limit_bytes=args.rss_limit_bytes,
            )
        if not validate_bnf_cache(bnf_meta, bnf_checkpoint, case, args.bnf_flag):
            record.update({"status": "INCOMPLETE_BNF", "selmer": None, "generic_embedding": None, "quotient_covers": None, "exceptional_embedding": None, "classification": None})
            runs.append(record)
            continue

        selmer_path = case_dir / "selmer.json"
        if selmer_path.exists() and not args.overwrite:
            raise FileExistsError(selmer_path)
        selmer_path.unlink(missing_ok=True)
        selmer_payload = {
            "case_id": case.case_id,
            "cas_directory": str(CAS),
            "bnf_metadata": str(bnf_meta),
            "pari_stack_bytes": args.pari_stack_bytes,
            "pari_debug": args.pari_debug,
            "random_seed": args.random_seed,
            "simon_verbose": args.simon_verbose,
            "raw_basis_search_bound": args.raw_basis_search_bound,
        }
        selmer_source = SELMER_WORKER.replace("SIMON_FUNCTION", repr(SIMON_GP_FUNCTION))
        record["selmer"] = supervise_source(
            sage_python, selmer_source, selmer_payload, selmer_path, case_dir / "selmer.log",
            timeout=args.selmer_timeout, rss_limit_bytes=args.rss_limit_bytes,
        )
        if record["selmer"]["outcome"] != "completed":
            record.update({"status": "INCOMPLETE_SELMER", "generic_embedding": None, "quotient_covers": None, "exceptional_embedding": None, "classification": None})
            runs.append(record)
            continue

        generic_path = case_dir / "generic-embedding.json"
        generic_payload = {
            "case_id": case.case_id,
            "cas_directory": str(CAS),
            "role": "generic-mw17",
            "bnf_metadata": str(bnf_meta),
            "selmer_result": str(selmer_path),
            "points": [[str(x), str(y)] for x, y in case.generic_points],
            "pari_stack_bytes": args.pari_stack_bytes,
            "auxiliary_prime_bound": args.auxiliary_prime_bound,
        }
        record["generic_embedding"] = supervise_source(
            sage_python, EMBED_WORKER, generic_payload, generic_path, case_dir / "generic-embedding.log",
            timeout=args.embedding_timeout, rss_limit_bytes=args.rss_limit_bytes,
        )
        if record["generic_embedding"]["outcome"] != "completed":
            record.update({"status": "INCOMPLETE_GENERIC_EMBEDDING", "quotient_covers": None, "exceptional_embedding": None, "classification": None})
            runs.append(record)
            continue

        selmer_result = json.loads(selmer_path.read_text())
        generic_result = json.loads(generic_path.read_text())
        quotient_basis = extend_standard(
            [[int(bit) for bit in row] for row in generic_result["point_selmer_rows"]],
            int(selmer_result["two_selmer_dimension"]),
        )
        quotient_path = case_dir / "quotient-covers.json"
        quotient_payload = {
            "case_id": case.case_id,
            "cas_directory": str(CAS),
            "bnf_metadata": str(bnf_meta),
            "selmer_result": str(selmer_path),
            "quotient_basis": quotient_basis,
            "pari_stack_bytes": args.pari_stack_bytes,
            "search_bound": args.quotient_search_bound,
            "quotient_class_limit": args.quotient_class_limit,
        }
        record["quotient_covers"] = supervise_source(
            sage_python, QUOTIENT_COVER_WORKER, quotient_payload,
            quotient_path, case_dir / "quotient-covers.log",
            timeout=args.cover_timeout, rss_limit_bytes=args.rss_limit_bytes,
        )
        if record["quotient_covers"]["outcome"] != "completed":
            record.update({"status": "INCOMPLETE_QUOTIENT_COVERS", "exceptional_embedding": None, "classification": None})
            runs.append(record)
            continue

        exceptional_result = None
        if case.exceptional_points:
            exceptional_path = case_dir / "exceptional-embedding.json"
            exceptional_payload = {
                "case_id": case.case_id,
                "cas_directory": str(CAS),
                "role": "held-out-exceptional",
                "bnf_metadata": str(bnf_meta),
                "selmer_result": str(selmer_path),
                "points": [[str(x), str(y)] for x, y in case.exceptional_points],
                "pari_stack_bytes": args.pari_stack_bytes,
                "auxiliary_prime_bound": args.auxiliary_prime_bound,
            }
            record["exceptional_embedding"] = supervise_source(
                sage_python, EMBED_WORKER, exceptional_payload, exceptional_path, case_dir / "exceptional-embedding.log",
                timeout=args.embedding_timeout, rss_limit_bytes=args.rss_limit_bytes,
            )
            if record["exceptional_embedding"]["outcome"] != "completed":
                record.update({"status": "INCOMPLETE_EXCEPTIONAL_EMBEDDING", "classification": None})
                runs.append(record)
                continue
            exceptional_result = json.loads(exceptional_path.read_text())
        else:
            record["exceptional_embedding"] = None
        classification = combine_results(
            case,
            selmer_result,
            generic_result,
            exceptional_result,
            json.loads(quotient_path.read_text()),
        )
        record["classification"] = classification
        record["status"] = "COMPLETE_CERTIFIED_CHECKPOINTED_TWO_SELMER"
        runs.append(record)

    complete = sum(row["status"].startswith("COMPLETE") for row in runs)
    output = {
        "schema": OUTPUT_SCHEMA,
        "status": "COMPLETE_ALL_SELECTED_CHECKPOINTED_DESCENTS" if complete == len(runs) else "INCOMPLETE_ONE_OR_MORE_CHECKPOINTED_DESCENTS",
        "backend": {"license": "open_source", "sage_version": versions["sage"], "pari_version": versions["pari"], "method": "certified persisted BNF + Simon field-Selmer/local conditions + explicit intersections of quadrics"},
        "input_manifest": {"path": str(args.manifest), "sha256": file_sha256(args.manifest)},
        "parameters": {
            "selected_case_ids": [row["case_id"] for row in runs],
            "bnf_timeout_seconds": args.bnf_timeout,
            "selmer_timeout_seconds": args.selmer_timeout,
            "embedding_timeout_seconds": args.embedding_timeout,
            "cover_timeout_seconds": args.cover_timeout,
            "rss_limit_bytes": args.rss_limit_bytes,
            "pari_stack_bytes": args.pari_stack_bytes,
            "bnf_tech": bnf_tech,
            "bnf_flag": args.bnf_flag,
            "random_seed": args.random_seed,
            "raw_basis_search_bound": args.raw_basis_search_bound,
            "quotient_search_bound": args.quotient_search_bound,
            "quotient_class_limit": args.quotient_class_limit,
            "auxiliary_prime_bound": args.auxiliary_prime_bound,
        },
        "completed_case_count": complete,
        "runs": runs,
        "claim_boundary": [
            "No Selmer result is recorded without a full certified BNF, all Simon local conditions, and exact point squareclass embeddings.",
            "The Selmer worker receives no known rational points; exceptional points are loaded only after blind raw-basis and quotient-cover searches finish.",
            "A bounded cover-search miss is not evidence of a nontrivial Tate-Shafarevich class.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"{PROTOCOL}|stage=complete|selected={len(runs)}|completed={complete}|status={output['status']}|output={args.output}", flush=True)


if __name__ == "__main__":
    main()
