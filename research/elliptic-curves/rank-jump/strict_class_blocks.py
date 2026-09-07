#!/usr/bin/env python3
"""Strict local kernels and unramified quadratic character blocks."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
import remaining_bad_primes as rem

OUTPUT = r.OUT/"rank_jump_strict_class_blocks_v1.json"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),rem.INPUT,rem.OUTPUT,r.INPUT,rem.bad.INPUT,rem.dy.INPUT)}


def local_rows(index):
    old = r.read(rem.bad.INPUT)["cases"][index]
    two = r.read(rem.dy.INPUT)["cases"][index]
    more = r.read(rem.INPUT)["cases"][index]["local"]
    assert more.get("status") != "UNKNOWN"
    rows = old["local"]+two["local"]+more["local"]
    return [p for p in rows if p["prime"] in (2,"infinity") or
            p["local_reduction"].get("conductor_exponent",0)>0]


def kernel(rows,n):
    columns = []
    metadata = []
    for row in rows:
        for j in range(len(row["point_signature_rows"][0])):
            columns.append(r.pack(s[j] for s in row["point_signature_rows"]))
            metadata.append([row["prime"],j])
    return lc.orthogonal(columns,n),columns,metadata


def calculate(index):
    source = r.read(rem.bad.INPUT)["cases"][index]
    profile,_,_ = r.characterize(rem.bad.cases()[index])
    assert profile["independent_input_indices"] == source["selected_input_indices"]
    n,m = source["witness_dimension"],source["generic_dimension"]
    assert profile["certified_independent_subgroup_rank_exact"] == n
    rows = local_rows(index)
    full,columns,metadata = kernel(rows,n)
    generic = lc.orthogonal([c&((1<<m)-1) for c in columns],m)
    extension = []
    span = list(generic)
    for v in full:
        if r.rank(span+[v]) > len(span):
            span.append(v)
            extension.append(v)
    assert r.rank(span) == len(full)
    assert all(not ((v&c).bit_count()%2) for v in span for c in columns)
    joint = rem.bad.characterize({**source,"local":rows})
    residual = joint["joint_quotient_support_dimension"]
    assert len(extension) == n-m-residual
    # Linear local functionals annihilating the marked generic image.
    # Their values on all known points give the exact joint obstruction matrix.
    point_rows = [r.pack((c>>i)&1 for c in columns) for i in range(n)]
    annihilator = lc.orthogonal(point_rows[:m],len(columns))
    checks = [r.pack((a&v).bit_count()%2 for v in point_rows[m:]) for a in annihilator]
    selected_checks = []
    images = []
    for a,v in zip(annihilator,checks):
        if r.rank(images+[v]) > len(images):
            images.append(v)
            selected_checks.append({"local_character_mask":a,"exceptional_value_mask":v})
    assert len(selected_checks) == residual
    complete = r.read(rem.INPUT)["cases"][index]["factor"]["factorization_complete"]
    return {"id":source["id"],"case_index":index,"generic_dimension":m,
            "witness_dimension":n,"all_bad_places_complete":complete,
            "tested_places":[p["prime"] for p in rows],
            "joint_generic_rank":m-len(generic),"joint_witness_rank":n-len(full),
            "joint_exceptional_rank":residual,"generic_strict_kernel_masks":generic,
            "witness_strict_kernel_masks":full,"relative_strict_lift_masks":extension,
            "local_character_coordinates":metadata,"joint_quotient_checks":selected_checks,
            "ordinary_class_group_two_rank_lower_bound":len(full) if complete else None,
            "generic_unramified_character_dimension":len(generic) if complete else None,
            "relative_unramified_character_dimension":len(extension) if complete else None,
            "unramified_compositum_degree":2**len(full) if complete else None,
            "scope":"Characters of the ordinary ideal class group, not an independently computed basis of ideal-half classes. Incomplete rows have only a tested-place kernel."}


def build(check=False):
    result = {"schema":"rank-jump.strict-class-blocks.v1","bindings":bindings(),
              "rows":[calculate(i) for i in range(6)],
              "boundary":"Independent known point classes define unramified quadratic extensions only on rows with complete bad-place coverage. Class-group ranks are lower bounds. This is retrospective incidence evidence, not a solubility criterion or prospective selector."}
    if check:
        assert r.read(OUTPUT) == result
        print("PASS strict-kernel and class-block algebra")
    else:
        r.write_new(OUTPUT,result)
    for row in result["rows"]:
        print(row["id"],"local quotient",row["joint_exceptional_rank"],
              "class rank at least",row["ordinary_class_group_two_rank_lower_bound"],
              "relative characters",row["relative_unramified_character_dimension"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",choices=("build","check"))
    args = parser.parse_args()
    build(args.mode == "check")
