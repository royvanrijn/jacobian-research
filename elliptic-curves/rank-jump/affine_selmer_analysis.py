#!/usr/bin/env python3
"""Proof-level consequences of six affine local systems and four CT entries."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
import affine_selmer as af
import affine_ct as ct
from cubic_bridge import Cubic

OUTPUT=r.OUT/"rank_jump_affine_selmer_analysis_v1.json"


def extension_rank(matrix,radical_pairings):
    n=len(matrix)
    rank=r.rank(map(r.pack,matrix))
    radical=lc.orthogonal(map(r.pack,matrix),n)
    assert len(radical)==len(radical_pairings)
    return rank+2*int(any(radical_pairings))


def build(check=False):
    local=r.read(af.OUTPUT);pairs=r.read(ct.OUTPUT);old=r.read(lc.INPUT)
    assert local["bindings"]==af.bindings() and pairs["bindings"]==ct.bindings()
    A,B=map(r.F,old["anchor"]["short_model_ainvariants"][3:])
    K=Cubic(A,B);rows=[]
    for f in local["cases"]:
        u=f["u"];Bct=next(x["matrix"] for x in old["ct"] if x["u"]==u)
        oldbasis=next(x["W_u_basis"] for x in old["rows"] if int(x["parameter_u"])==u)
        radicals=[lc.lift(v,oldbasis) for v in lc.orthogonal(map(r.pack,Bct),len(Bct))]
        row={"u":u,"inherited_selmer_dimension":len(Bct),
             "inherited_CT_rank":r.rank(map(r.pack,Bct)),
             "affine_locally_admissible":f["affine_solution"]["consistent"]}
        if not f["affine_solution"]["consistent"]:
            word=f["affine_solution"]["inconsistent_row_combination"]
            assert lc.lift(word,f["all_constraint_rows"])==1<<20
            row.update({"affine_transport_excluded_for_all_inherited_representatives":True,
                        "contradiction_places":[p["place"] for p in f["places"]
                                               if any(word>>j&1 for j in p["global_row_indices"])],
                        "maximum_independent_anchor_classes_transported":0})
        else:
            evidence=[next(p for p in pairs["pairs"] if p["u"]==u and p["old_radical_mask"]==m)
                      for m in radicals]
            assert all("pair" in p for p in evidence)
            values=[p["pair"]["value"] for p in evidence]
            rank=extension_rank(Bct,values)
            n=len(Bct)+1
            row.update({"affine_anchor_particular":f["affine_solution"]["particular_anchor_mask"],
                        "old_radical_masks":radicals,"new_class_pairings_with_old_radical":values,
                        "enlarged_selmer_dimension":n,"enlarged_CT_rank":rank,
                        "enlarged_restricted_radical_dimension":n-rank,
                        "affine_coset_contains_restricted_CT_radical_class":not any(values),
                        "maximum_independent_anchor_classes_transported":len(radicals)+1 if not any(values) else 0,
                        "point_or_Sha_status":"UNKNOWN"})
            if u==-1:
                D=1+A-B;x=1+A;y=D
                assert y*y==x**3-2*A*x*x+(A-3*B+A*A)*x+B+A*B+B*B
                alpha=K.sub(K.theta,K.square(K.theta))
                gamma=K.add(K.one,K.theta)
                eta=K.scale(gamma,D)
                assert K.mul(K.sub(K.scalar(x),alpha),K.square(gamma))==eta
                row.update({"eta_rational_point":[str(x),str(y)],
                            "eta_Kummer_square_correction":list(map(str,gamma)),
                            "eta_rational_class_independent_of_inherited_space":True,
                            "maximum_independent_anchor_classes_transported":len(radicals),
                            "point_or_Sha_status":"eta is a certified rational class; other radical directions UNKNOWN"})
        rows.append(row)
    bindings={str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
              for p in (af.OUTPUT,ct.OUTPUT,lc.INPUT,Path(__file__).resolve())}
    out={"schema":"rank-jump.affine-selmer-analysis.v1","bindings":bindings,"cases":rows,
         "scope":"Ranks refer to the displayed Selmer subspaces and restricted CT forms. No full Selmer dimension or curve rank upper bound."}
    if check:
        assert r.read(OUTPUT)==out
        print("PASS affine Selmer and CT consequences")
    else:r.write_new(OUTPUT,out)
    for row in rows:
        print(row["u"],row["affine_locally_admissible"],
              row.get("enlarged_restricted_radical_dimension"),row["maximum_independent_anchor_classes_transported"])


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("build","check"))
    build(parser.parse_args().mode=="check")
