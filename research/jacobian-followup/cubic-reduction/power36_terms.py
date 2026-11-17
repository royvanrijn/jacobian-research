#!/usr/bin/env python3
"""Change the companion tensor to KCUB3's36-term base; bounded termwise polarization."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json
from power_search import X,sp,homogeneous18,linear,decomposition,DB
from polarize import search,build
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def make_terms():
    base=HERE/"cubic12-integer-witness.json";d=json.loads(base.read_text())
    K=[sp.Add(*(sp.Rational(t["coefficient"])*sp.prod(x**e for x,e in zip(X,t["exponent"])) for t in row)) for row in d["coefficients"]]
    H,_=homogeneous18(K);H[-1]=X[-1]**3
    y=sp.symbols("u1:4");z=[3*X[1]-X[6],X[1]+X[9],3*X[1]-2*X[10]]
    grouped=list(H)
    replacements=[(11,X[0]*y[1]*X[17]),(12,-X[0]**2*y[0]),(13,12*X[0]*y[0]*y[1]),
                  (14,2*X[0]*X[2]*y[2]+4*X[0]*X[8]*y[0])]
    for i,t in replacements:
        grouped[i]=t+sp.expand(grouped[i]-t.subs(dict(zip(y,z)),simultaneous=True))
    assert all(sp.expand(a-b.subs(dict(zip(y,z)),simultaneous=True))==0 for a,b in zip(H,grouped))
    vs=[tuple(int(sp.expand(f).coeff(x)) for x in X) for f in list(X)+z];ts=[]
    for row,f in enumerate(grouped):
        for mon,co in sp.Poly(f,*(X+y)).terms():
            if not co:continue
            factors=[v for v,e in zip(vs,mon) for _ in range(e)]
            assert len(factors)==3;ts.append((row,Q(co),factors))
    return ts,H

def main():
    out=HERE/"power36-termwise";out.mkdir(exist_ok=True);ts,H=make_terms()
    print("TERMS",len(ts),"BASELINE",len(build(ts,[0]*len(ts))[0]),flush=True)
    n,choices,forms,B=search(ts,restarts=50,seed=261118)
    cubes=[dict(sp.Poly(linear(f)**3,*X).terms()) for f in forms]
    for i,h in enumerate(H):
        result={}
        for p,c in zip(cubes,B[i]):
            for e,a in p.items():result[e]=result.get(e,0)+c*a
        result={e:c for e,c in result.items() if c}
        assert result=={e:Q(c) for e,c in sp.Poly(h,*X).terms() if c}
    payload=dict(status="exact-cube-identity",source_sha256=sha(Path(__file__)),
      base_sha256=sha(HERE/"cubic12-integer-witness.json"),size=n,forms=forms,
      B=[[str(c) for c in row] for row in B],choices=choices,padding_index=forms.index((0,)*17+(1,)),
      seed=261118,restarts=50,scope="Bounded termwise search on the changed KCUB3 companion. No smaller Keller map without full determinant and collision construction.")
    (out/"candidate.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in payload.items() if k not in ["forms","B","choices"]}),flush=True)
if __name__=="__main__":main()
