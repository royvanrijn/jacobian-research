"""Try quadratic carrier-target corrections after a single graph slice."""
from pathlib import Path
import hashlib,json,time
import sympy as sp
HERE=Path(__file__).resolve().parent
def main():
    d=json.loads((HERE/"cubic17-witness.json").read_text())
    X=sp.symbols(d["variables"]);n=len(X)
    R=[sp.Poly(sum(sp.Rational(t["coefficient"])*sp.prod(v**e for v,e in zip(X,t["exponent"])) for t in row),*X) for row in d["coefficients"]]
    log=[];t0=time.monotonic()
    for j in range(3,n):
        f=sp.expand(R[j].as_expr()-X[j])
        if sp.degree(f,X[j])>0 or sp.Poly(f,*X).total_degree()>2:continue
        kept=[i for i in range(n) if i!=j];Y=[X[i] for i in kept]
        S=[sp.Poly(R[i].as_expr().subs(X[j],-f),*Y) for i in kept]
        if max(p.total_degree() for p in S[3:])>3:continue
        pairs=[(a,b) for a in range(3,n-1) for b in range(a,n-1)]
        basis=[S[a]*S[b] for a,b in pairs]
        high=sorted({e for p in basis+S[:3] for e,c in p.terms() if sum(e)>3})
        M=sp.Matrix([[p.coeff_monomial(e) for p in basis] for e in high])
        rhs=sp.Matrix([[S[i].coeff_monomial(e) for i in range(3)] for e in high])
        rr,piv=M.row_join(rhs).rref()
        ok=all(q<len(basis) for q in piv)
        row=dict(removed=d["variables"][j],raw_degree=max(p.total_degree() for p in S),high_equations=len(high),quadratic_products=len(basis),correctable=ok,seconds=round(time.monotonic()-t0,3))
        log.append(row);print(json.dumps(row),flush=True)
        if ok:
            solutions=sp.zeros(len(basis),3)
            for a,b in enumerate(piv):
                for i in range(3):solutions[b,i]=rr[a,len(basis)+i]
            for i in range(3):
                S[i]-=sum((basis[a]*solutions[a,i] for a in range(len(basis))),sp.Poly(0,*Y))
            assert max(p.total_degree() for p in S)<=3
            result=dict(source_dimension=17,dimension=16,removed=d["variables"][j],variables=list(map(str,Y)),
                coefficients=[[dict(exponent=list(e),coefficient=str(c)) for e,c in p.terms()] for p in S],
                quadratic_corrections=[[dict(pair=list(pairs[a]),coefficient=str(solutions[a,i])) for a in range(len(basis)) if solutions[a,i]] for i in range(3)])
            (HERE/"cubic16-candidate.json").write_text(json.dumps(result,indent=2)+"\n")
            break
    (HERE/"quadratic-slice-repair.json").write_text(json.dumps(dict(status="complete",source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=log),indent=2)+"\n")
if __name__=="__main__":main()

