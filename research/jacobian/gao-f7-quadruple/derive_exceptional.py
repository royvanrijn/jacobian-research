"""Construct the full D=0 source projection; no acceptance or exclusion decision."""
from pathlib import Path
import json,time
import sympy as S
from flint import fmpz_mpoly_ctx,fmpz_poly
HERE=Path(__file__).resolve().parent;START=time.monotonic()
s,a,b,z=S.symbols('s a b z');D=3*s**4-2*s
source=json.loads((HERE/'quadruple-input.json').read_text())
f=[S.rem(S.Poly(S.sympify(row['polynomial']),s),S.Poly(D,s)).as_expr() for row in source['equations']]
f=[S.Poly(S.expand(q),s,a,b).clear_denoms()[1].primitive()[1].as_expr() for q in f]
first,second=[S.Poly(q,b) for q in f];assert first.degree()==second.degree()==1
result=first.nth(1)*second.nth(0)-second.nth(1)*first.nth(0)
H=S.rem(S.Poly(S.expand(result),s),S.Poly(D,s)).as_expr();H=S.Poly(S.expand(H),a,s).clear_denoms()[1].primitive()[1].as_expr()
ctx=fmpz_mpoly_ctx.get(['s','z']);SN,ZN=ctx.gens();carrier=3*SN**4-2*SN
shifted=S.Poly(S.expand(H.subs(a,z-s)),s,z);native=ctx.from_dict({ex:int(c) for ex,c in shifted.terms()});R=native.resultant(carrier,'s');terms=R.to_dict();assert all(ex[0]==0 for ex in terms)
coeff=[0]*(1+max(ex[1] for ex in terms))
for (_,j),c in terms.items():coeff[j]=int(c)
P=fmpz_poly(coeff);P=P//P.content();sf=P//P.gcd(P.derivative())
data={'variables':['s','a','b'],'carrier':str(D),'equations':[str(q) for q in f],'linear_resultant':str(H),'projection_coordinate':'z=a+s','projection_coefficients':list(map(int,P.coeffs())),'squarefree_projection':list(map(int,sf.coeffs())),'scope':'Necessary finite projection on D=0 only; no count or exclusion registered.'}
(HERE/'exceptional-source.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'status':'EXCEPTIONAL_FINITE_SOURCE_DERIVED','projection_degree':P.degree(),'seconds':time.monotonic()-START,'exclusion':False}))
