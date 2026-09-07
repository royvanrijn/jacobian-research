#!/usr/bin/env python3
"""Exact split class factors and unresolved left/right Artin kernels."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
import half_ideal_artin as art
import half_ideal_artin_completion as complete

OUTPUT=r.OUT/"rank_jump_half_ideal_class_blocks_v1.json"


def independent_indices(vectors):
    selected=[];basis=[]
    for i,v in enumerate(vectors):
        if r.rank(basis+[v])>len(basis):
            selected.append(i);basis.append(v)
    return selected


def calculate(row,raw):
    n=row["character_dimension"];g=row["generic_strict_dimension"]
    columns=row["columns"]
    ci=independent_indices(columns)
    # The generic columns are independent in all three retained matrices.
    assert list(range(g))==ci[:g]
    restricted_rows=[r.pack((columns[j]>>i)&1 for j in ci) for i in range(n)]
    ri=independent_indices(restricted_rows)
    d=len(ci);assert len(ri)==d
    subcols=[r.pack((columns[j]>>i)&1 for i in ri) for j in ci]
    assert r.rank(subcols)==d
    dual_words=[]
    for i in range(d):
        coefficients=lc.coordinates(1<<i,subcols)
        word=r.pack(bool(coefficients>>ci.index(j)&1) if j in ci else False for j in range(n))
        assert r.pack((lc.lift(word,columns)>>k)&1 for k in ri)==1<<i
        dual_words.append(word)
    left=lc.orthogonal(columns,n)
    right=lc.orthogonal([r.pack(x) for x in row["matrix_rows"]],n)
    masks=raw["character_point_masks"]
    return {"case_index":row["case_index"],"id":row["id"],
            "known_strict_dimension":n,"generic_strict_dimension":g,
            "split_elementary_S_class_factor_dimension":d,
            "generic_dimension_inside_selected_factor":g,
            "relative_dimension_inside_selected_factor":d-g,
            "selected_half_ideal_indices":ci,"selected_character_indices":ri,
            "dual_half_ideal_words":dual_words,
            "remaining_S_class_two_rank":"%d + epsilon"%(n-d),
            "left_character_kernel_masks":left,
            "left_character_kernel_point_masks":[lc.lift(x,masks) for x in left],
            "right_half_ideal_kernel_masks":right,
            "right_half_ideal_kernel_point_masks":[lc.lift(x,masks) for x in right],
            "complete_known_block_splits":d==n,
            "boundary":"C=H_selected direct sum ker(selected characters), H_selected elementary of the displayed dimension. epsilon=c_S-known_strict_dimension is uncomputed. A right kernel is only invisible to retained characters, not certified zero in C/2."}


def build(check=False):
    rows=r.read(complete.OUTPUT)["rows"]
    originals={x["case_index"]:x for x in r.read(art.INPUT)["cases"]}
    paths=(Path(__file__),complete.OUTPUT,art.INPUT)
    report={"schema":"rank-jump.half-ideal-class-blocks.v1",
            "bindings":{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
            "rows":[calculate(row,originals[row["case_index"]]) for row in rows]}
    if check:
        assert r.read(OUTPUT)==report;print("PASS explicit split factors and two separate Artin kernels")
    else:r.write_new(OUTPUT,report)
    for row in report["rows"]:
        print(row["id"],"factor",row["split_elementary_S_class_factor_dimension"],
              "relative",row["relative_dimension_inside_selected_factor"],
              "remaining",row["remaining_S_class_two_rank"])


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=("build","check"))
    build(p.parse_args().mode=="check")
