#!/usr/bin/env python3
"""Finite Artin pairing on a point-blind prime-ideal dictionary."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
import strict_class_blocks as strict
import remaining_bad_primes as rem

PROTOCOL = Path(__file__).with_name("STRICT_ARTIN_PROTOCOL.json")
OUTPUT = r.OUT/"rank_jump_strict_artin_v1.json"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),PROTOCOL,strict.OUTPUT,rem.INPUT,rem.bad.INPUT,r.INPUT)}


def ideal_dictionary(index):
    # This routine consumes no generic or exceptional point coordinates.
    source = r.read(rem.bad.INPUT)["cases"][index]
    B,A,_,_ = map(int,source["integral_cubic_ascending"])
    disc = -4*A**3-27*B*B
    bound = r.read(PROTOCOL)["limits"]["rational_prime_bound"]
    ideals = []
    for p in r.primes(bound+1):
        if disc%p == 0:
            continue
        for a in range(p):
            if (a*a*a+A*a+B)%p == 0:
                ideals.append({"prime":p,"theta_residue":a,"ideal_norm":p,
                               "label":"(p, theta-a) in the maximal order"})
    return ideals


def calculate(index):
    frozen = ideal_dictionary(index)
    old = r.read(rem.bad.INPUT)["cases"][index]
    block = r.read(strict.OUTPUT)["rows"][index]
    assert block["all_bad_places_complete"]
    source = rem.bad.cases()[index]
    _,allpoints = r.short(source["model"],source["generic_points"]+source["points"])
    d = r.F(old["elliptic_scaling_d"])
    points = [[str(r.F(allpoints[i][0])*d*d),str(r.F(allpoints[i][1])*d**3)]
              for i in old["selected_input_indices"]]
    B,A,_,_ = old["integral_cubic_ascending"]
    model = ["0","0","0",A,B]
    masks = block["generic_strict_kernel_masks"]+block["relative_strict_lift_masks"]
    n = len(masks)
    assert r.rank(masks) == n == block["ordinary_class_group_two_rank_lower_bound"]
    point_columns = []
    for ideal in frozen:
        characters = [r.point_signature(model,P,[(ideal["prime"],[ideal["theta_residue"]])]) for P in points]
        point_columns.append(r.pack(characters))
    columns = [r.pack((mask&c).bit_count()%2 for mask in masks) for c in point_columns]
    selected = []
    independent = []
    for j,col in enumerate(columns):
        if r.rank(independent+[col]) > len(independent):
            selected.append(j)
            independent.append(col)
    dual = []
    if len(selected) == n:
        for i in range(n):
            coordinates = lc.coordinates(1<<i,independent)
            word = sum(1<<selected[j] for j in range(n) if coordinates>>j&1)
            assert lc.lift(word,columns) == 1<<i
            dual.append(word)
    return {"case_index":index,"id":old["id"],"integral_cubic_ascending":old["integral_cubic_ascending"],
            "tested_places":block["tested_places"],"ideal_dictionary":frozen,
            "generic_character_count":len(block["generic_strict_kernel_masks"]),
            "character_masks":masks,"point_evaluation_columns":point_columns,
            "artin_columns":columns,"artin_matrix_rows":[[(c>>i)&1 for c in columns] for i in range(n)],
            "artin_rank":len(selected),"selected_ideal_indices":selected,"dual_ideal_words":dual,
            "full_dual_basis":len(selected)==n,"ordinary_and_S_class_two_rank_lower_bound":len(selected),
            "relative_dual_ideal_count":len(block["relative_strict_lift_masks"]) if len(selected)==n else None,
            "boundary":"Exact Artin evaluations of retained unramified characters. Ideal pool is point-blind; characters and the certified rank lower bound are retrospective."}


def build(check=False):
    data = {"schema":"rank-jump.strict-artin.v1","bindings":bindings(),
            "rows":[calculate(i) for i in r.read(PROTOCOL)["cases"]]}
    if check:
        assert r.read(OUTPUT) == data
        print("PASS exact Artin matrix and dual ideal words")
    else:
        r.write_new(OUTPUT,data)
    for row in data["rows"]:
        print(row["id"],"ideals",len(row["ideal_dictionary"]),"rank",row["artin_rank"],
              "generic",row["generic_character_count"],"relative",row["relative_dual_ideal_count"],
              "selected",[(row["ideal_dictionary"][i]["prime"],row["ideal_dictionary"][i]["theta_residue"])
                           for i in row["selected_ideal_indices"]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",choices=("build","check"))
    args = parser.parse_args()
    build(args.mode=="check")
