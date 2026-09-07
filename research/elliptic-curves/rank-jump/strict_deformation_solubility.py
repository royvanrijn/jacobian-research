#!/usr/bin/env python3
"""Strict zero-localization spaces and CT cross-obstructions in six retained deformations."""
import argparse
from itertools import combinations
from pathlib import Path
import retrospective as r
import local_collision as lc
import affine_selmer as af

PROTOCOL=Path(__file__).with_name("STRICT_DEFORMATION_SOLUBILITY_PROTOCOL.json")
OUTPUT=r.OUT/"rank_jump_strict_deformation_solubility_v1.json"


def strict_constraints(case):
    rows=[];places=[]
    for local in case["local"]:
        signatures=local["class_signature_rows"][:20]
        columns=r.transpose(signatures)
        start=len(rows);rows.extend(columns)
        places.append({"place":local["place"],"coordinate_count":len(columns),
                       "row_indices":list(range(start,len(rows))),"constraint_rows":columns})
    return rows,places


def cross_report(basis,u,local_source):
    old=next(x for x in local_source["rows"] if int(x["parameter_u"])==u)
    W=old["W_u_basis"]
    CT=next(x["matrix"] for x in local_source["ct"] if x["u"]==u)
    coordinates=[lc.coordinates(x,W) for x in basis]
    cross=[r.pack(lc.pairing(w,1<<j,CT) for j in range(len(W))) for w in coordinates]
    conditions=r.transpose([[(x>>j)&1 for j in range(len(W))] for x in cross]) if basis else []
    kernel=lc.orthogonal(conditions,len(basis))
    witness=None
    if any(cross):
        i=next(i for i,v in enumerate(cross) if v)
        j=next(j for j in range(len(W)) if cross[i]>>j&1)
        witness={"strict_basis_index":i,"strict_anchor_mask":basis[i],
                 "partner_inherited_basis_index":j,"partner_anchor_mask":W[j],"CT_value":1,
                 "strict_coordinates_in_inherited_basis":coordinates[i]}
    return {"u":u,"inherited_selmer_dimension":len(W),"strict_basis_in_inherited_coordinates":coordinates,
            "cross_pairing_rows":cross,"cross_pairing_rank":r.rank(cross),
            "necessary_soluble_dimension":len(kernel),
            "necessary_soluble_anchor_basis":[lc.lift(w,basis) for w in kernel],
            "first_obstructed_strict_class":witness}


def calculate():
    source=r.read(lc.INPUT);raw=r.read(af.INPUT)
    assert raw["bindings"]==af.bindings()
    us=r.read(PROTOCOL)["parameters"]
    cases={c["u"]:c for c in raw["cases"]}
    rows={};local=[]
    for u in us:
        rows[u],places=strict_constraints(cases[u])
        old=next(x for x in source["rows"] if int(x["parameter_u"])==u)
        assert old["all_local_kummer_images_complete"]
        expected={x["prime"] for x in old["finite_local_conditions"]}|{"infinity"}
        assert {x["place"] for x in places}==expected
        K=lc.orthogonal(rows[u],20)
        local.append({"u":u,"places":places,"strict_zero_localization_constraint_rank":r.rank(rows[u]),
                      "strict_anchor_basis":K,"strict_dimension":len(K),
                      "CT_cross_report":cross_report(K,u,source)})
    common=lc.orthogonal([x for u in us for x in rows[u]],20)
    pairs=[]
    for u,v in combinations(us,2):
        K=lc.orthogonal(rows[u]+rows[v],20)
        reports=[cross_report(K,t,source) for t in (u,v)]
        pairkernel=lc.intersection(*(x["necessary_soluble_anchor_basis"] for x in reports))
        pairs.append({"u":u,"v":v,"common_strict_basis":K,"common_strict_dimension":len(K),
                      "CT_cross_reports":reports,"common_necessary_soluble_basis":pairkernel})
    return {"single_deformations":local,"all_six_common_strict_basis":common,
            "all_six_common_strict_dimension":len(common),
            "all_six_CT_cross_reports":[cross_report(common,u,source) for u in us],
            "pairs":pairs,
            "summary":{"single_strict_dimensions":[x["strict_dimension"] for x in local],
                       "single_CT_cross_ranks":[x["CT_cross_report"]["cross_pairing_rank"] for x in local],
                       "common_all_six_dimension":len(common),
                       "pairs_with_nonzero_strict_space":sum(bool(x["common_strict_basis"]) for x in pairs)},
            "boundary":"Strict means zero localization at every retained bad place, 2 and infinity, not merely membership in both local point images. Nonzero CT cross-pairing certifies nonrationality; zero cross-pairing does not certify a point."}


def build(check=False):
    paths=(Path(__file__),PROTOCOL,af.INPUT,lc.INPUT)
    report={"schema":"rank-jump.strict-deformation-solubility.v1",
            "bindings":{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},**calculate()}
    if check:
        assert r.read(OUTPUT)==report;print("PASS strict localization and CT cross-pairing accounting")
    else:r.write_new(OUTPUT,report)
    print(report["summary"])
    for row in report["single_deformations"]:
        print(row["u"],row["strict_anchor_basis"],row["CT_cross_report"]["first_obstructed_strict_class"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=("build","check"))
    build(p.parse_args().mode=="check")
