"""Derive the residual septic at a putative fivefold root; no exclusion decision."""
from pathlib import Path
import hashlib,json,time
import sympy as S
from flint import fmpz_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic()
ss,aa,bb=S.symbols('s a b');G2=aa*aa+aa+ss*bb+ss**3;G3=2*ss*G2+aa+ss**3+ss*ss*bb;G4=3*ss*ss*G2+bb+ss**4
X1=S.expand(-S.det(S.Matrix([G2,G3,G4]).jacobian([ss,aa,bb])))
X=[X1,S.expand(G2+ss*X1),S.expand(G3+ss*ss*X1),S.expand(G4+ss**3*X1)]
ctx=fmpz_mpoly_ctx.get(['t','s','a','b']);t,s,a,b=ctx.gens()
def native(expression):return ctx.from_dict({(0,)+ex:int(c) for ex,c in S.Poly(expression,ss,aa,bb).terms()})
Y1,Y2,Y3,Y4=map(native,X)
bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4;at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt
R=at*at+at+t*bt+t**3+Y1*t-Y2
Q,remainder=divmod(R,(t-s)**5)
assert (t-s)**5*Q+remainder==R
raw=json.loads((HERE.parent/'gao-f7-fibres/sixfold-input.json').read_text());f=[S.sympify(q) for q in raw['polynomials']]
# Divide in the shifted variable: the lower Taylor coefficients are exactly
# the three necessary derivative polynomials, up to nonzero constants.
x=S.symbols('x')
remexpr=sum(int(c)*x**ex[0]*ss**ex[1]*aa**ex[2]*bb**ex[3] for ex,c in remainder.to_dict().items())
shift=S.Poly(S.expand(remexpr.subs(x,x+ss)),x)
assert shift.nth(0)==shift.nth(1)==0
for j in range(2,5):
 ratio=S.cancel(shift.nth(j)/f[j-2]);assert ratio.is_Rational and ratio!=0
rows=[]
for j in range(8):
 terms={ex[1:]:int(c) for ex,c in Q.to_dict().items() if ex[0]==j}
 expr=sum(c*ss**ex[0]*aa**ex[1]*bb**ex[2] for ex,c in terms.items())
 rows.append({'coefficient':j,'polynomial':str(expr),'terms':len(terms)})
data={'variables':['s','a','b'],'factorization_on_fivefold_locus':'R(t;X)=(t-s)^5*Q7(t); f2=f3=f4=0','source':'../gao-f7-fibres/sixfold-input.json','equations':[{'derivative':j+2,'polynomial':str(f[j])} for j in range(3)],'quotient_coefficients':rows,'scope':'Necessary source equations and exact monic quotient only. Neither residual squarefreeness nor exclusion of a three-point type is yet proved.'}
(HERE/'quintuple-input.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({'status':'RESIDUAL_SEPTIC_DERIVED','quotient_terms':sum(r['terms'] for r in rows),'source_terms':[len(S.Poly(q,ss,aa,bb).terms()) for q in f[:3]],'seconds':time.monotonic()-START,'source_sha256':hashlib.sha256((HERE.parent/'gao-f7-fibres/sixfold-input.json').read_bytes()).hexdigest(),'exclusion':False}))
