#!/usr/bin/env python3
"""Norm-square rank gate and a coefficient-only nonsplit structural control."""
import argparse
from itertools import permutations
from pathlib import Path
import retrospective as r
import local_collision as lc
import torsion_difference as td
from cubic_bridge import Cubic

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/"NORM_SQUARE_RANK_PROTOCOL.json"
OUTPUT=r.OUT/"rank_jump_norm_square_rank_v1.json"


def even_sign_module():
    states={td.encode_state(k,v):k|(v<<2) for k in range(4) for v in range(4)}
    matrices=[];actions=[]
    for perm in permutations(range(3)):
        g=(perm[0]+1,perm[1]+1)
        for signs in range(8):
            if signs.bit_count()%2:continue
            c=td.sign_matrix(signs);cg=td.compose(c,g)
            M=(g[0],g[1],cg[0]|(g[0]<<2),cg[1]|(g[1]<<2))
            for subset,state in states.items():
                assert states[td.permute_subset(subset,perm,signs)]==td.act(M,state)
            matrices.append(M)
            actions.append({"permutation":list(perm),"sign_mask":signs,"matrix_columns":list(M)})
    assert len(set(matrices))==24
    assert all(td.compose(a,b) in matrices for a in matrices for b in matrices)
    elementary=[tuple((1<<i) if j==k else 0 for j in range(4))
                for k in range(4) for i in range(4)]
    equations=[]
    for M in matrices:
        images=[td.add(td.compose(X,M),td.compose(M,X)) for X in elementary]
        for j in range(4):
            for i in range(4):
                equations.append(r.pack((X[j]>>i)&1 for X in images))
    commutant=lc.orthogonal(equations,16)
    pack=lambda M:sum(v<<(4*j) for j,v in enumerate(M))
    identity=(1,2,4,8);nilpotent=(0,0,1,2)
    assert commutant==lc.canonical([pack(identity),pack(nilpotent)])
    idempotents=[]
    for mask in range(1<<len(commutant)):
        word=lc.lift(mask,commutant)
        M=tuple((word>>(4*j))&15 for j in range(4))
        if td.compose(M,M)==M:idempotents.append(word)
    assert sorted(idempotents)==[0,pack(identity)]
    return {"actions":actions,"root_action_checks":384,"group_order":24,
            "commutant_dimension":2,"commutant_basis_packed":commutant,
            "idempotents_packed":sorted(idempotents),"module":"INDECOMPOSABLE_NONSPLIT",
            "abstract_group":"S4"}


def calculate():
    anchor=r.read(lc.INPUT)["anchor"]
    A,B=map(r.F,anchor["short_model_ainvariants"][3:])
    disc=-4*A**3-27*B**2
    assert A<0<B and disc>0
    K=Cubic(A,B);u=-A/B
    gamma=K.sub(K.one,K.scale(K.theta,u))
    assert K.norm(gamma)==1+A*u*u+B*u**3==1
    assert gamma==K.scale(K.mul(K.square(K.theta),K.theta),-1/B)
    base_a2=A;base_a6=B*B
    assert (-B)**2==(-A)**3+base_a2*(-A)**2+base_a6
    # Doubling (0,B) on y^2=x^3+A*x^2+B^2 has slope zero.
    assert -A/B==u
    return {"control":{"u":str(u),"A":str(A),"B":str(B),"D":1,
                       "cubic_discriminant":str(disc),"gamma_coordinates":list(map(str,gamma)),
                       "gamma_identity":"gamma=-theta^3/B",
                       "root_sign_counts":{"negative":1,"positive":2},
                       "gamma_sign_counts":{"negative":2,"positive":1},
                       "gamma_is_square_in_totally_real_splitting_field":False,
                       "relative_Kummer_degree_over_splitting_field":4,
                       "Jacobian_two_torsion_field_degree_over_Q":24,
                       "auxiliary_curve_ainvariants":["0",str(A),"0","0",str(B*B)],
                       "auxiliary_point":["0",str(B)],"auxiliary_double":[str(-A),str(-B)],
                       "specialized_Mordell_Weil_rank":"UNKNOWN"},
            "surface":{"twist":"D(u)=1+A*u^2+B*u^3",
                       "twisted_discriminant":"16*(-4*A^3-27*B^2)*D(u)^8",
                       "twisted_c4":"16*D(u)^2*(A^2*u^2-9*B*u-3*A)",
                       "finite_fibres":["I2*","I2*","I2*"],"infinity":"I0",
                       "Euler_number":24,"holomorphic_Euler_characteristic":2,
                       "geometric_Neron_Severi_rank":20,"trivial_lattice_rank":20,
                       "geometric_generic_Mordell_Weil_rank":0,
                       "norm_base_geometric_generic_rank":1,
                       "norm_base_arithmetic_generic_rank_for_retained_anchor":0,
                       "new_generic_rank_from_norm_base_change":0},
            "boundary":r.read(PROTOCOL)["boundary"]}


def build(check=False):
    paths=(Path(__file__),PROTOCOL,lc.INPUT,td.OUTPUT,HERE/"torsion_difference.py")
    data={"schema":"rank-jump.norm-square-rank.v1",
          "bindings":{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
          "even_sign_module":even_sign_module(),**calculate()}
    if check:assert r.read(OUTPUT)==data;print("PASS norm-square control and 384 root-action checks")
    else:r.write_new(OUTPUT,data)
    print("D-twist generic rank",data["surface"]["geometric_generic_Mordell_Weil_rank"])
    print("D=1 control torsion field degree",data["control"]["Jacobian_two_torsion_field_degree_over_Q"])


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("mode",choices=("build","check"))
    build(p.parse_args().mode=="check")
