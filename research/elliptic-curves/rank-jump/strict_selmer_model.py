#!/usr/bin/env python3
"""Exact numeric specialization of the strict Selmer / S-class-group sequence."""
import argparse
from pathlib import Path
import retrospective as r
import strict_artin as art
import strict_class_blocks as strict
import remaining_bad_primes as rem

OUTPUT = r.OUT/"rank_jump_strict_selmer_model_v1.json"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__),art.OUTPUT,strict.OUTPUT,rem.OUTPUT)}


def calculate():
    rows = []
    support = r.read(rem.OUTPUT)["rows"]
    blocks = r.read(strict.OUTPUT)["rows"]
    for arithmetic in r.read(art.OUTPUT)["rows"]:
        i = arithmetic["case_index"]
        block = blocks[i]
        local = support[i]["bad_places_and_two_real"]
        assert block["all_bad_places_complete"]
        n,m = block["witness_dimension"],block["generic_dimension"]
        k = len(block["witness_strict_kernel_masks"])
        a = local["joint_witness_image_dimension"]
        ell = local["full_product_point_image_dimension"]
        assert n == k+a
        assert ell-a == 1
        missing_characters = art.lc.orthogonal(arithmetic["artin_columns"],k)
        missing_point_masks = [art.lc.lift(v,arithmetic["character_masks"]) for v in missing_characters]
        assert len(missing_characters) == k-arithmetic["artin_rank"]
        rows.append({"id":block["id"],"case_index":i,"witness_rank":n,"generic_rank":m,
                     "known_quotient_rank":n-m,"known_strict_dimension":k,
                     "generic_strict_dimension":len(block["generic_strict_kernel_masks"]),
                     "witness_local_image_dimension":a,"local_product_dimension":ell,
                     "unknown_local_image_dimension_bound":[0,ell-a],
                     "unknown_S_class_excess_definition":"epsilon = dim_F2 Cl(O_K,S)/2 - known_strict_dimension >= 0",
                     "full_Selmer_dimension_formula":f"{n} + epsilon + b",
                     "full_relative_Selmer_dimension_formula":f"{n-m} + epsilon + b",
                     "b_definition":f"dim localization(Sel_2) - {a}, an integer in [0,{ell-a}]",
                     "ideal_pool_rank":arithmetic["artin_rank"],
                     "characters_unseparated_by_ideal_pool":missing_characters,
                     "unseparated_character_point_masks":missing_point_masks,
                     "claim_boundary":"epsilon and b are not computed. Formula concerns the full 2-Selmer group, not the full rational rank or the rational solubility of its unknown classes."})
    return {"schema":"rank-jump.strict-selmer-model.v1","bindings":bindings(),"rows":rows,
            "theorem":"For Q and an irreducible two-division cubic, with S containing 2, infinity and all bad places, Sel_2^S(E) is canonically Hom(Cl(O_K,S),F2). See STRICT_SELMER_AND_ARTIN_BLOCKS.md for proof and exact sequences."}


def build(check=False):
    data = calculate()
    if check:
        assert r.read(OUTPUT) == data
        print("PASS strict Selmer dimension accounting")
    else:
        r.write_new(OUTPUT,data)
    for row in data["rows"]:
        print(row["id"],"Selmer =",row["full_Selmer_dimension_formula"],
              "unseparated",row["characters_unseparated_by_ideal_pool"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode",choices=("build","check"))
    args = parser.parse_args()
    build(args.mode=="check")
