#!/usr/bin/env python3
"""Independent exact audit of the 17-variable cubic map. Standard library only."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
import hashlib,json
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clean(p):return {e:c for e,c in p.items() if c}
def add(a,b,scale=Q(1)):
    out=dict(a)
    for e,c in b.items():out[e]=out.get(e,Q(0))+scale*c
    return clean(out)
def mul(a,b):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(x+y for x,y in zip(e,f));out[g]=out.get(g,Q(0))+c*d
    return clean(out)
def scale(p,c):return clean({e:d*c for e,d in p.items()})
def derivative(p,j):
    out={}
    for e,c in p.items():
        if e[j]:
            f=list(e);f[j]-=1;out[tuple(f)]=c*e[j]
    return out
def variable(j,n):
    e=[0]*n;e[j]=1;return {tuple(e):Q(1)}
def identity(n,nvars):
    one={(0,)*nvars:Q(1)}
    return [[one if i==j else {} for j in range(n)] for i in range(n)]
def madd(A,B,coefficient=Q(1)):
    return [[add(a,b,coefficient) for a,b in zip(r,s)] for r,s in zip(A,B)]
def mmul(A,B):
    out=[[{} for _ in B[0]] for _ in A]
    for i,row in enumerate(A):
        for k,p in enumerate(row):
            if p:
                for j,q in enumerate(B[k]):
                    if q:out[i][j]=add(out[i][j],mul(p,q))
    return out
def mzero(A):return not any(p for row in A for p in row)
def determinant3(A,nvars):
    out={}
    for p in permutations(range(3)):
        sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        term={(0,)*nvars:Q(sign)}
        for i in range(3):term=mul(term,A[i][p[i]])
        out=add(out,term)
    return out
def evaluate(p,point):
    value=Q(0)
    for e,c in p.items():
        term=c
        for a,n in zip(point,e):term*=a**n
        value+=term
    return value
def compose(p,sub,nvars):
    one={(0,)*nvars:Q(1)}
    powers=[]
    for j,q in enumerate(sub):
        stack=[one]
        for _ in range(max((e[j] for e in p),default=0)):stack.append(mul(stack[-1],q))
        powers.append(stack)
    out={}
    for e,c in p.items():
        term={(0,)*nvars:c}
        for j,n in enumerate(e):
            if n:term=mul(term,powers[j][n])
        out=add(out,term)
    return out
def decode(row,n):
    out={tuple(t["exponent"]):Q(t["coefficient"]) for t in row}
    assert len(out)==len(row) and all(len(e)==n and min(e)>=0 for e in out)
    assert all(out.values())
    return out
def external_rows(d):
    names=["x","y","z"]+[f"w{i}" for i in range(1,17)]
    rows=[]
    for row in d["components"]:
        out={}
        for text,c in row.items():
            e=[0]*19
            for atom in text.split("*"):
                bits=atom.split("^");e[names.index(bits[0])]+=int(bits[1]) if len(bits)==2 else 1
            assert tuple(e) not in out
            out[tuple(e)]=Q(c)
        rows.append(out)
    return rows
def verify_construction(witness,external):
    # Retained 17 variables, followed by complementary r,s.
    X=[variable(i,19) for i in range(19)]
    sub=X[:13]+[X[17],add(scale(X[17],Q(6)),X[13],Q(-1)),
                 X[14],X[18],add(scale(X[18],Q(3)),X[15],Q(-1)),X[16]]
    old=[compose(p,sub,19) for p in external]
    retained=old[:3]+old[3:13]+[add(scale(old[13],Q(6)),old[14],Q(-1)),
             old[15],add(scale(old[16],Q(3)),old[17],Q(-1)),old[18]]
    assert len(retained)==17
    assert all(e[17:]==(0,0) for p in retained for e in p)
    expected=[{e[:17]:c for e,c in p.items()} for p in retained]
    expected=[scale(expected[2],Q(1,2)),expected[1],expected[0]]+expected[3:]
    assert expected==witness
    # The removed output rows are r+y*w3 and s+y*w7.
    assert old[13]==add(X[17],mul(X[1],X[5]))
    assert old[16]==add(X[18],mul(X[1],X[9]))
    return True
def midpoint_identity():
    # Exact monomial check in Q[a,b] for every univariate monomial degree<=2.
    # Polarization in arbitrary coordinates gives the same identity for mixed quadratics.
    a=variable(0,2);b=variable(1,2);half=scale(add(a,b),Q(1,2));delta=add(a,b,Q(-1))
    assert add(mul(a,a),mul(b,b),Q(-1))==scale(mul(half,delta),Q(2))
    # Explicit mixed quadratic in four independent endpoint coordinates.
    a0,a1,b0,b1=[variable(i,4) for i in range(4)]
    lhs=add(mul(a0,a1),mul(b0,b1),Q(-1))
    rhs=add(mul(scale(add(a1,b1),Q(1,2)),add(a0,b0,Q(-1))),
            mul(scale(add(a0,b0),Q(1,2)),add(a1,b1,Q(-1))))
    assert lhs==rhs
def main():
    data=HERE/"cubic17-witness.json";d=json.loads(data.read_text())
    assert d["source_sha256"]==sha(HERE/"construct17.py")
    source=HERE/"sources/degree3-map-external.json"
    provenance=json.loads((HERE/"sources/provenance.json").read_text())
    assert d["external_data_sha256"]==sha(source)==provenance["sha256"]
    n=d["dimension"];assert n==17
    F=[decode(row,n) for row in d["coefficients"]]
    assert len(F)==n and max(sum(e) for p in F for e in p)==3
    counts=[sum(sum(e)>=2 for e in p) for p in F];assert counts==d["nonlinear_counts"] and sum(counts)==58
    verify_construction(F,external_rows(json.loads(source.read_text())))
    J=[[derivative(p,j) for j in range(n)] for p in F]
    origin=(0,)*n;assert all(not p.get(origin,0) for p in F)
    assert [[p.get(origin,Q(0)) for p in row] for row in J]==[[Q(i==j) for j in range(n)] for i in range(n)]
    A=[r[:3] for r in J[:3]];B=[r[3:] for r in J[:3]]
    C=[r[:3] for r in J[3:]];D=[r[3:] for r in J[3:]]
    I=identity(n-3,n);L=madd(D,I,Q(-1));power=I;inv=I;nilpotence=None
    for exponent in range(1,n-2):
        power=mmul(power,L)
        if mzero(power):nilpotence=exponent;break
        inv=madd(inv,power,Q((-1)**exponent))
    assert nilpotence==d["carrier_nilpotence_exponent"]==3
    assert mmul(inv,D)==I and mmul(D,inv)==I
    S=madd(A,mmul(B,mmul(inv,C)),Q(-1))
    det=determinant3(S,n);assert det=={origin:Q(1)}
    points=[[Q(q) for q in p] for p in d["points"]]
    assert len({tuple(p) for p in points})==3
    images=[[evaluate(f,p) for f in F] for p in points]
    assert images[0]==images[1]==images[2]==[Q(q) for q in d["common_image"]]
    midpoint_identity()
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),witness_sha256=sha(data),
       primary_source_sha256=sha(HERE/"construct17.py"),external_data_sha256=sha(source),
       dimension=n,degree=3,nonlinear_counts=counts,nonlinear_total=58,
       determinant_one=True,identity_derivative=True,collision_points=3,
       retained_rows_independent_of_removed_directions=True,
       removed_rows_are_identity_plus_retained_polynomials=True,
       full_polynomial_schur_determinant=True,carrier_nilpotence_exponent=3,
       quadratic_midpoint_identity=True,
       scope="Minimum ordinary degree for characteristic-zero Keller collisions is3 across all finite dimensions. Dimension17 is an upper bound for cubic witnesses; no minimum dimension or term count is proved.")
    (HERE/"validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt))
if __name__=="__main__":main()

