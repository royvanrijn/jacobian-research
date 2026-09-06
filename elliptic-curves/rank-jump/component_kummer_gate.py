#!/usr/bin/env python3
"""Exact additive-component Kummer obstruction, with a same-class control."""
import argparse
import itertools
from pathlib import Path
import retrospective as r
import triple_character as previous
from cubic_bridge import Cubic

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"COMPONENT_KUMMER_GATE_PROTOCOL.json"
INPUT=previous.INPUT
OUTPUT=r.OUT/"rank_jump_component_kummer_gate_v1.json"


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes())
            for p in (Path(__file__),PROTOCOL,INPUT,HERE/"cubic_bridge.py")}


def point_add(A,B,P,Q):
    x,y=P;z,w=Q
    assert y*y==x**3+A*x+B and w*w==z**3+A*z+B
    if x==z:
        assert y==w and y
        slope=(3*x*x+A)/(2*y)
    else:slope=(w-y)/(z-x)
    intercept=y-slope*x
    xx=slope*slope-x-z;yy=-y+slope*(x-xx)
    assert yy*yy==xx**3+A*xx+B
    return (xx,yy),(intercept,slope,r.F(0))


def build(check=False):
    rows=[]
    for row in r.read(INPUT)["cases"]:
        A,B=map(r.F,row["model"][3:]);points=[tuple(map(r.F,P)) for P in row["generic_points"]]
        K=Cubic(A,B);betas=[K.sub(K.scalar(P[0]),K.theta) for P in points]
        infinity=K.scale(K.theta,-B)
        assert K.norm(infinity)==B**4
        samples=[]
        for p in r.primes(r.read(PROTOCOL)["limits"]["odd_character_primes_max"]):
            roots=r.roots_at(str(A),str(B),p)
            if roots is None:continue
            try:
                values=[r.mod(P[0],p) for P in points];bb=r.mod(B,p)
            except ValueError:continue
            for root in roots:
                residues=[(x-root)%p for x in values]+[(-bb*root)%p]
                if not all(residues):continue
                samples.append({"p":p,"root":root,"bits":[int(pow(v,(p-1)//2,p)==p-1) for v in residues]})
        signatures=[sum(v["bits"][i]<<j for j,v in enumerate(samples)) for i in range(4)]
        assert r.rank(signatures[:3])==3
        pair_witnesses=[]
        for i,j in itertools.combinations(range(4),2):
            witness=next((v for v in samples if v["bits"][i]!=v["bits"][j]),None)
            pair_witnesses.append({"indices":[i,j],"unequal_squareclasses":witness is not None,"witness":witness})
        ranks=[]
        for mask in (3,5,6,7):
            indices=[i for i in range(3) if mask>>i&1]
            i,j=indices[:2]
            witness=next(v for v in pair_witnesses if v["indices"]==[i,j])
            assert witness["unequal_squareclasses"]
            ranks.append({"mask":mask,"arithmetic_generic_rank":0,"obstructing_indices":[i,j],
                          "witness":witness["witness"]})
        controls=[]
        for j in (1,2):
            P,Q=points[0],points[j]
            middle,line1=point_add(A,B,P,Q)
            target,line2=point_add(A,B,middle,Q)
            sigma=K.mul(line2,K.inverse(line1))
            beta_target=K.sub(K.scalar(target[0]),K.theta)
            assert beta_target==K.mul(betas[0],K.square(sigma))
            controls.append({"combination":[1,2 if j==1 else 0,2 if j==2 else 0],
                             "point":list(map(str,target)),"intermediate_point":list(map(str,middle)),
                             "first_chord":list(map(str,line1)),"second_chord":list(map(str,line2)),
                             "square_root_of_class_ratio":list(map(str,sigma))})
        derived=[points[0][0]]+[r.F(v["point"][0]) for v in controls]
        assert all(derived) and len(set(derived))==3
        rows.append({"id":row["id"],"samples":samples,"signature_vectors":signatures,
                     "generic_class_rank":3,"class_index_3":"[-B*theta], the infinity additive condition",
                     "pairwise_squareclass_witnesses":pair_witnesses,"mixed_characters":ranks,
                     "pair_base_arithmetic_generic_rank":2,"triple_base_arithmetic_generic_rank":3,
                     "production_curve_rank":"UNKNOWN",
                     "same_class_control":{"source_points":[0,1,2],"derived_points":controls,
                                           "index_in_original_three_point_lattice":4,
                                           "finite_additive_component_gate":"PASS",
                                           "mixed_triple_rational_section_exists":"UNKNOWN"}})
    permutations=list(itertools.permutations(range(3)))
    def perm(mask,p):return sum(((mask>>i)&1)<<p[i] for i in range(3))
    even=[0,3,5,6]
    invariants=[v for v in even if all(perm(v,p)==v for p in permutations)]
    assert invariants==[0]
    flags=[v for v in range(8) if all(perm(v,p)==v for p in permutations)]
    assert flags==[0,7]
    out={"schema":"rank-jump.component-kummer-gate.v1","bindings":bindings(),
         "component_module":{"even_subset_invariants":invariants,"I2_component_flag_invariants":flags,
                             "kappa_valuation_rows":[[0,1,1],[1,0,1],[1,1,0]]},
         "cases":rows,
         "boundary":"Exact arithmetic ranks for the new function fields only. Equal additive conditions are necessary, not sufficient, for a rational mixed section."}
    if check:assert out==r.read(OUTPUT);print("PASS component-Kummer gate replay")
    else:r.write_new(OUTPUT,out)
    for row in rows:print(row["id"],"class rank",row["generic_class_rank"],"pair rank 2, triple rank 3")


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=["build","check"]);args=p.parse_args()
    build(args.mode=="check")
