#!/usr/bin/env python3
"""Independent rational-polynomial replay and full determinant of the reused cubic witness."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from verify import add,mul,scale,variable,evaluate,derivative,identity,madd,mmul,determinant3,decode
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def invert_constant(A):
    n=len(A);M=[list(row)+[Q(i==j) for j in range(n)] for i,row in enumerate(A)];det=Q(1)
    for col in range(n):
        pivot=next(i for i in range(col,n) if M[i][col])
        if pivot!=col:M[col],M[pivot]=M[pivot],M[col];det=-det
        c=M[col][col];det*=c;M[col]=[v/c for v in M[col]]
        for i in range(n):
            if i!=col and M[i][col]:
                c=M[i][col];M[i]=[v-c*w for v,w in zip(M[i],M[col])]
    return [row[n:] for row in M],det
def initial():
    # Independent full coefficient construction; no discovery helper is imported.
    return [
      {(1,0,0):Q(1),(2,1,0):Q(-3,2),(3,0,1):Q(-1,2)},
      {(0,1,0):Q(1),(1,0,1):Q(3),(2,1,1):Q(6),(3,2,1):Q(3),(1,2,0):Q(12),(2,3,0):Q(9)},
      {(0,0,1):Q(1),(1,1,1):Q(3),(2,2,1):Q(3),(3,3,1):Q(1),(0,2,0):Q(4),(1,3,0):Q(7),(2,4,0):Q(3)}]
def replay(raw):
    F=initial();points=[[Q(0),Q(-2),Q(-16)],[Q(1),Q(-2),Q(8)],[Q(-1),Q(1),Q(5)]]
    baseJ=[[derivative(f,j) for j in range(3)] for f in F]
    assert determinant3(baseJ,3)=={(0,0,0):Q(1)}
    assert [[evaluate(f,p) for f in F] for p in points]==[[Q(0),Q(-2),Q(0)]]*3
    stats={"carrier":0,"product":0,"linear":0}
    for op in raw["witness"]["moves"]:
        kind=op["kind"];assert kind in stats;stats[kind]+=1
        if kind=="carrier":
            n=len(F);e=tuple(op["exponent"]);assert len(e)==n and min(e)>=0
            values=[evaluate({e:Q(1)},p) for p in points]
            F=[{a+(0,):c for a,c in f.items()} for f in F]
            F.append(add(variable(n,n+1),{e+(0,):Q(1)},Q(-1)))
            points=[p+[v] for p,v in zip(points,values)]
        elif kind=="product":
            i,j=op["factors"];assert 0<=i<len(F) and 0<=j<len(F)
            targets=[row[0] for row in op["targets"]]
            assert len(set(targets))==len(targets) and not ({i,j}&set(targets))
            product=mul(F[i],F[j])
            for target,c in op["targets"]:
                assert 0<=target<len(F) and Q(c)
                F[target]=add(F[target],product,-Q(c))
        else:
            i=op["target"];j=op["factor"];assert i!=j and 0<=i<len(F) and 0<=j<len(F)
            F[i]=add(F[i],F[j],-Q(op["coefficient"]))
    n=len(F)
    assert F==[decode(row,n) for row in raw["witness"]["coefficients"]]
    assert points[:2]==[[Q(c) for c in p] for p in raw["witness"]["points"]]
    linear=[[f.get(tuple(int(j==k) for j in range(n)),Q(0)) for k in range(n)] for f in F]
    inverse,det=invert_constant(linear);assert det==1
    normalized=[]
    for row in inverse:
        p={}
        for c,f in zip(row,F):
            if c:p=add(p,f,c)
        normalized.append(p)
    return normalized,points,linear,stats
def main():
    data=HERE/"cubic13-witness.json";d=json.loads(data.read_text())
    assert d["source_sha256"]==sha(HERE/"construct13.py")
    assert d["discovery_sha256"]==sha(HERE/"reuse-witness.json")
    raw=json.loads((HERE/"reuse-witness.json").read_text())
    assert raw["status"]=="complete" and raw["source_sha256"]==sha(HERE/"search_reuse.py")
    for name,digest in raw["helper_sha256"].items():assert sha(HERE/name)==digest
    F,points,linear,stats=replay(raw);n=len(F)
    assert n==d["dimension"]==13
    assert F==[decode(row,n) for row in d["coefficients"]]
    assert linear==[[Q(c) for c in row] for row in d["tangent_matrix"]]
    degrees=[max(sum(e) for e in f) for f in F];assert max(degrees)==d["degree"]==3
    counts=[sum(sum(e)>=2 for e in f) for f in F];assert counts==d["nonlinear_counts"] and sum(counts)==47
    assert points==[[Q(c) for c in p] for p in d["points"]]
    images=[[evaluate(f,p) for f in F] for p in points]
    assert len({tuple(p) for p in points})==3 and images[0]==images[1]==images[2]==[Q(c) for c in d["common_image"]]
    J=[[derivative(f,j) for j in range(n)] for f in F];zero=(0,)*n
    assert all(not f.get(zero,0) for f in F)
    assert [[p.get(zero,Q(0)) for p in row] for row in J]==[[Q(i==j) for j in range(n)] for i in range(n)]
    A=[r[:3] for r in J[:3]];B=[r[3:] for r in J[:3]]
    C=[r[:3] for r in J[3:]];D=[r[3:] for r in J[3:]];I=identity(n-3,n);N=madd(D,I,Q(-1))
    power=I;inv=I;nilpotence=None
    for exponent in range(1,n-2):
        power=mmul(power,N)
        if not any(p for row in power for p in row):nilpotence=exponent;break
        inv=madd(inv,power,Q((-1)**exponent))
    assert nilpotence==d["carrier_nilpotence_exponent"]==3
    assert mmul(D,inv)==I and mmul(inv,D)==I
    S=madd(A,mmul(B,mmul(inv,C)),Q(-1))
    assert determinant3(S,n)=={zero:Q(1)}
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/"verify.py"),
        primary_source_sha256=sha(HERE/"construct13.py"),witness_sha256=sha(data),
        discovery_sha256=sha(HERE/"reuse-witness.json"),dimension=n,degree=3,
        nonlinear_counts=counts,nonlinear_total=47,collision_points=3,determinant_one=True,
        identity_derivative=True,all_variables_symbolic=True,full_polynomial_schur_determinant=True,
        carrier_nilpotence_exponent=nilpotence,elementary_moves=stats,
        base_determinant_verified=True,entire_transformation_replayed=True,
        scope="Explicit cubic dimension13 upper bound and47 displayed nonlinear occurrences; no minimum dimension or support claim; no change to the three-variable interval [10,48].")
    (HERE/"validation13.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

