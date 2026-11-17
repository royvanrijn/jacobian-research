"""Replay the reusable-factor construction and normalize the resulting cubic map."""
from pathlib import Path
import hashlib,json
import sympy as s
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    raw=json.loads((HERE/"reuse-witness.json").read_text());w=raw["witness"]
    assert raw["source_sha256"]==sha(HERE/"search_reuse.py")
    for name,digest in raw["helper_sha256"].items():assert sha(HERE/name)==digest
    x,y,z=s.symbols("x y z");X=[x,y,z];u=1+x*y
    F=s.Matrix([x-s.Rational(3,2)*x*x*y-s.Rational(1,2)*x**3*z,
                y+3*x*u*u*z+3*x*y*y*(4+3*x*y),
                u**3*z+y*y*u*(4+3*x*y)])
    points=[[s.Rational(0),s.Rational(-2),s.Rational(-16)],
            [s.Rational(1),s.Rational(-2),s.Rational(8)],
            [s.Rational(-1),s.Rational(1),s.Rational(5)]]
    assert s.expand(F.jacobian(X).det())==1
    for move in w["moves"]:
        kind=move["kind"]
        if kind=="carrier":
            e=move["exponent"];assert len(e)==len(X) and min(e)>=0
            a=s.prod(q**n for q,n in zip(X,e));new=s.Symbol("a"+str(len(X)-2))
            for point in points:point.append(a.subs(dict(zip(X,point))))
            F=F.col_join(s.Matrix([new-a]));X.append(new)
        elif kind=="product":
            i,j=move["factors"];p=F[i]*F[j]
            for target,c in move["targets"]:
                assert target not in (i,j);F[target]=s.expand(F[target]-s.Rational(c)*p)
        elif kind=="linear":
            i=move["target"];j=move["factor"];assert i!=j
            F[i]=s.expand(F[i]-s.Rational(move["coefficient"])*F[j])
        else:raise ValueError(kind)
    expected=s.Matrix([sum(s.Rational(t["coefficient"])*s.prod(q**n for q,n in zip(X,t["exponent"])) for t in row) for row in w["coefficients"]])
    assert (F-expected).applyfunc(s.expand)==s.zeros(len(X),1)
    J0=F.jacobian(X).subs(dict.fromkeys(X,0));assert J0.det()==1
    F=(J0.inv()*F).applyfunc(s.expand)
    n=len(X);J=F.jacobian(X);assert J.subs(dict.fromkeys(X,0))==s.eye(n)
    A=J[:3,:3];B=J[:3,3:];C=J[3:,:3];D=J[3:,3:];N=D-s.eye(n-3)
    power=s.eye(n-3);inv=power;nilpotence=None
    for exponent in range(1,n-2):
        power=(power*(-N)).applyfunc(s.expand)
        if power==s.zeros(n-3):nilpotence=exponent;break
        inv=(inv+power).applyfunc(s.expand)
    assert nilpotence and (D*inv).applyfunc(s.expand)==s.eye(n-3)
    Schur=(A-B*inv*C).applyfunc(s.expand)
    determinant=s.expand(Schur.det(method="berkowitz"));assert determinant==1
    images=[[s.expand(q.subs(dict(zip(X,point)))) for q in F] for point in points]
    assert images[0]==images[1]==images[2]
    coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in s.Poly(q,*X).terms()] for q in F]
    counts=[sum(sum(t["exponent"])>=2 for t in row) for row in coefficients]
    keys={}
    for i,q in enumerate(F):
        for e,c in s.Poly(q,*X).terms():
            if sum(e)!=3:continue
            for j in range(3,n):
                if e[j]:
                    ee=list(e);ee[j]-=1
                    keys.setdefault((i,tuple(ee)),[s.Integer(0)]*(n-3))[j-3]+=c*e[j]
    radical=s.Matrix(list(keys.values())).nullspace()
    receipt=dict(status="exact-primary-pass",source_sha256=sha(Path(__file__)),discovery_sha256=sha(HERE/"reuse-witness.json"),
       dimension=n,degree=max(s.Poly(q,*X).total_degree() for q in F),variables=list(map(str,X)),
       coefficients=coefficients,nonlinear_counts=counts,nonlinear_total=sum(counts),
       points=[[str(q) for q in p] for p in points],common_image=list(map(str,images[0])),
       determinant="1",identity_derivative=True,carrier_nilpotence_exponent=nilpotence,
       tangent_matrix=[[str(q) for q in row] for row in J0.tolist()],
       cubic_carrier_radical=[list(map(str,v)) for v in radical],
       schur_complement=[[str(q) for q in row] for row in Schur.tolist()],
       boundary="Existence witness only; no minimum dimension or support claim.")
    (HERE/"cubic13-witness.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k not in ("coefficients","tangent_matrix","points","schur_complement")}))
if __name__=="__main__":main()

