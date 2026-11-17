#!/usr/bin/env python3
"""Independent SymPy replay; no producer or retained algebra imports."""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "dixmier-residue-transport"
X,Q,Z = sp.symbols("X Q Z")
def poly(expr):
    return sp.Poly(expr,X,Q,Z,domain=sp.QQ)
def load(rows):
    return poly(sum(sp.Rational(c)*X**m[0]*Q**m[1]*Z**m[2] for m,c in rows))
def dz(p,r=1):
    for _ in range(r):
        p=p.diff(Z)
    return p
def delta(p,r=1):
    for _ in range(r):
        p=poly(3*X**2)*p.diff(X)+poly(2-6*X*Q)*p.diff(Q)
    return p
def bidiff(p,q,r=1):
    return sum(((-1)**j*sp.binomial(r,j)*
                delta(dz(p,r-j),j)*dz(delta(q,r-j),j)
                for j in range(r+1)),poly(0))
def linear(S,T,u,v):
    return bidiff(u,T)+bidiff(S,v)

def main():
    # Verify the universal calculus identity BEFORE assuming {S,T}=1.
    s=sp.Function("s")(X,Q,Z)
    t=sp.Function("t")(X,Q,Z)
    I=sp.Function("I")(X,Q,Z)
    a=sp.Function("a")(X,Q)
    b=sp.Function("b")(X,Q)
    def deriv(f):
        return a*sp.diff(f,X)+b*sp.diff(f,Q)
    def bracket(f,g):
        return sp.diff(f,Z)*deriv(g)-deriv(f)*sp.diff(g,Z)
    lhs=bracket(I*sp.diff(s,Z),t)+bracket(s,I*sp.diff(t,Z))
    rhs=sp.diff(I,Z)*bracket(s,t)+I*sp.diff(bracket(s,t),Z)
    assert sp.expand(lhs-rhs)==0

    receipt=json.loads((HERE/"receipt.json").read_text())
    for name,digest in receipt["hashes"].items():
        assert hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest()==digest,name
    data=json.loads((SOURCE/"inputs.json").read_text())
    witness=json.loads((HERE/"escape_witness.json").read_text())
    S,T,s2,t2=(load(data[k]) for k in ("S","T","s2","t2"))
    assert bidiff(S,T)==poly(1)
    for k,p in [("S",S),("T",T),("s",s2),("t",t2)]:
        assert load(witness[k])==p
    assert linear(S,T,s2,t2)+bidiff(S,T,3)/24==poly(0)
    defect=(bidiff(s2,t2)+bidiff(s2,T,3)/24+
            bidiff(S,t2,3)/24+bidiff(S,T,5)/1920)
    assert defect==load(witness["O5"])
    primitive=defect.integrate(Z)
    u,v=-S.diff(Z)*primitive,-T.diff(Z)*primitive
    assert u==load(witness["u"]) and v==load(witness["v"])
    assert linear(S,T,u,v)+defect==poly(0)
    lam={tuple(m):sp.Rational(c) for m,c in
         json.loads((SOURCE/"combined_cocycle.json").read_text())["functional"]}
    def period(p):
        return sum(c*lam.get(m,0) for m,c in p.terms())
    assert period(defect)==1
    assert period(linear(S,T,u,v))==-1
    result={"claim":"QFL1","status":"PASS",
            "generic_commuting_derivation_identity":True,
            "retained_lower_equation":True,
            "independently_rebuilt_order_five_defect":True,
            "independently_rebuilt_escape_corrections":True,
            "full_pair_relation_through_order_five":True,
            "residues":["1","-1"],
            "imports":"SymPy and standard library only; no producer or legacy algebra.",
            "coverage_boundary":"All-orders induction and the full-frame all-rank theorem are written proofs; not formally verified or externally reviewed.",
            "hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in [HERE/"check_independent.py",HERE/"receipt.json",HERE/"escape_witness.json"]}}
    (HERE/"independent.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__":
    main()
