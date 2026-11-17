#!/usr/bin/env python3
"""Independent rational verification of the89-variable cubic-linear certificate.

The determinant is proved through explicit cube identities, the checked
12-variable base determinant, companion block elimination and Sylvester's
identity. An89-variable polynomial determinant is not expanded.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import lcm
import hashlib,json
from verify import add,mul,variable,decode,evaluate
from verify_sparse12 import audit
HERE=Path(__file__).resolve().parent;DB=HERE.parent/"de-bondt-2610.24889"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mv(M,v):return [sum((Q(c)*q for c,q in zip(row,v)),Q(0)) for row in M]
def component_cube(v):return [q**3 for q in v]
def scale_polynomial(p,c):return {e:q*c for e,q in p.items() if q*c}
def main():
    path=HERE/"power89/witness.json";d=json.loads(path.read_text());n=d["dimension"]
    assert n==89 and d["source_sha256"]==sha(HERE/"build_power89.py")
    assert d["candidate_sha256"]==sha(HERE/"power-basis-anneal/candidate.json")
    assert d["base_source_sha256"]==sha(DB/"verify.py")
    candidate=json.loads((HERE/"power-basis-anneal/candidate.json").read_text())
    assert candidate["source_sha256"]==sha(HERE/"power_basis_anneal.py")
    L=d["forms"];rawB=[[Q(c) for c in row] for row in d["raw_B"]]
    assert L==candidate["forms"] and rawB==[[Q(c) for c in row] for row in candidate["B"]]
    assert len(L)==n and all(len(row)==18 and all(type(c) is int for c in row) for row in L)
    assert len(rawB)==18 and all(len(row)==n for row in rawB)
    # Independently reconstruct scale11=2 from the previously pinned scale11=1 input.
    source=HERE/"de-bondt12-witness.json";base=json.loads(source.read_text())
    assert base["source_sha256"]==d["base_source_sha256"]
    initial=[decode(row,12) for row in base["coefficients"]]
    diagonal=[Q(1)]*12;diagonal[10]=Q(1,2);K=[]
    for i,p in enumerate(initial):
        row={}
        for e,c in p.items():
            factor=Q(1)/diagonal[i]
            for a,power in zip(diagonal,e):factor*=a**power
            row[e]=c*factor
        K.append(row)
    audit(K)
    quadratic=[{e:c for e,c in p.items() if sum(e)==2} for p in K]
    cubic=[{e:c for e,c in p.items() if sum(e)==3} for p in K]
    assert all(not p for p in cubic[5:])
    H=[]
    for i in range(12):
        p={}
        for e,c in quadratic[i].items():
            ee=list(e+(0,)*6);ee[17]=1;p[tuple(ee)]=c
        if i<5:
            ee=[0]*18;ee[12+i]=1;ee[17]=2;p[tuple(ee)]=Q(-1)
        H.append(p)
    H += [{e+(0,)*6:c for e,c in cubic[i].items()} for i in range(5)]+[{}]
    assert len(H)==18 and all(sum(e)==3 for p in H for e in p)
    tau3={(0,)*17+(3,):Q(1)}
    cubes=[]
    for row in L:
        lin={}
        for j,c in enumerate(row):
            if c:lin=add(lin,variable(j,18),Q(c))
        cubes.append(mul(mul(lin,lin),lin))
    for i in range(18):
        reconstructed={}
        for j,c in enumerate(rawB[i]):
            if c:reconstructed=add(reconstructed,cubes[j],c)
        assert reconstructed==(tau3 if i==17 else H[i])
    scale=d["scale"];assert type(scale) is int and scale==24
    assert scale==lcm(*(c.denominator for row in rawB[:17] for c in row))
    B=d["B"];assert len(B)==18 and all(len(row)==n and all(type(c) is int for c in row) for row in B)
    expectedB=[[c*(scale if i<17 else 1) for c in rawB[i]] for i in range(18)]
    assert all(c.denominator==1 for row in expectedB for c in row)
    assert B==[[int(c) for c in row] for row in expectedB]
    padding=d["padding_index"];assert 0<=padding<n
    assert L[padding]==[0]*17+[1] and B[17]==[int(j==padding) for j in range(n)]
    L0=[list(row) for row in L];L0[padding]=[0]*18
    A=d["A"];assert len(A)==n and all(len(row)==n and all(type(c) is int for c in row) for row in A)
    assert A==[[sum(L0[i][k]*B[k][j] for k in range(18)) for j in range(n)] for i in range(n)]
    assert not any(A[padding]) and any(c for row in A for c in row)
    shear=[Q(B[i][padding])-Q(i==17) for i in range(18)]
    assert shear==[Q(c) for c in d["target_shear"]] and shear[-1]==0
    for i in range(18):
        reconstructed={}
        for j in range(n):
            if j!=padding and B[i][j]:reconstructed=add(reconstructed,cubes[j],Q(B[i][j]))
        expected=add(scale_polynomial(H[i],Q(scale)),tau3,-shear[i])
        assert reconstructed==expected
    # A finite right-inverse certificate proves surjectivity of B.
    piv=d["right_inverse_columns"];inverse=[[Q(c) for c in row] for row in d["right_inverse_minor"]]
    assert len(piv)==len(set(piv))==18 and all(0<=j<n for j in piv)
    assert len(inverse)==18 and all(len(row)==18 for row in inverse)
    assert [[sum((Q(B[i][piv[k]])*inverse[k][j] for k in range(18)),Q(0)) for j in range(18)] for i in range(18)]==[[Q(i==j) for j in range(18)] for i in range(18)]
    def right_lift(z):
        values=mv(inverse,z);out=[Q(0)]*n
        for j,c in zip(piv,values):out[j]=c
        assert mv(B,out)==z
        return out
    points=[[Q(v)/a for v,a in zip(p,diagonal)] for p in base["points"]]
    assert points==[[Q(v) for v in p] for p in d["base_points"]]
    base_images=[[evaluate(f,p) for f in K] for p in points]
    assert points[0]!=points[1] and base_images[0]==base_images[1]
    lifted=[p+[-Q(scale)*evaluate(f,p) for f in cubic[:5]]+[Q(1,scale)] for p in points]
    assert lifted==[[Q(v) for v in p] for p in d["lifted_points"]]
    def companion(z):return [a+Q(scale)*evaluate(h,z) for a,h in zip(z,H)]
    images=[companion(z) for z in lifted];assert images[0]==images[1]
    u=right_lift(lifted[0]);trial=right_lift(lifted[1])
    image=[a+b for a,b in zip(u,component_cube(mv(A,u)))]
    v=[a-b for a,b in zip(image,component_cube(mv(A,trial)))]
    assert u==[Q(a) for a in d["u"]] and v==[Q(a) for a in d["v"]]
    assert image==[Q(a) for a in d["image"]]
    assert u!=v and mv(B,v)==lifted[1] and mv(A,v)==mv(A,trial)
    assert [a+b for a,b in zip(v,component_cube(mv(A,v)))]==image
    sheared=[a-b*lifted[0][-1]**3 for a,b in zip(images[0],shear)]
    assert mv(B,image)==sheared
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),
       arithmetic_sha256=sha(HERE/"verify.py"),base_determinant_checker_sha256=sha(HERE/"verify_sparse12.py"),
       witness_sha256=sha(path),primary_source_sha256=sha(HERE/"build_power89.py"),
       candidate_sha256=sha(HERE/"power-basis-anneal/candidate.json"),base_input_sha256=sha(source),
       dimension=n,previous_constructed_upper_bound=90,degree=3,coefficient_ring="Z",scale=scale,
       all18_cube_identities=True,all18_projection_identities=True,integer_matrix_factorization=True,
       right_inverse_verified=True,base12_full_polynomial_determinant=True,
       determinant_via_companion_and_sylvester=True,expanded89_variable_determinant=False,
       source_and_target_shear_retained=True,full_rational_collision_verified=True,
       distinct_points_verified=True,identity_linear_part=True,
       scope="Explicit integer89-variable cubic-linear Keller collision, improving the90-variable construction. No minimum dimension, priority, formalization or external review claim.")
    (HERE/"power89/validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

