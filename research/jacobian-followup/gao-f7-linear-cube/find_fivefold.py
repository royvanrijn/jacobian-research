"""Bounded finite-field search for a nonsingular fivefold-root Hensel seed.

An empty search is inconclusive. A seed must be separately replayed.
"""
from pathlib import Path
import json,math,time
import sympy as S
HERE=Path(__file__).resolve().parent;START=time.monotonic()
s,a,b,t=S.symbols('s a b t')
f=[S.sympify(q) for q in json.loads((HERE.parent/'gao-f7-fibres/sixfold-input.json').read_text())['polynomials']]
def table(expression):
 q=S.Poly(expression,b,a,s);bd=q.degree(b);ad=q.degree(a);sd=q.degree(s)
 result=[[[0]*(sd+1) for _ in range(ad+1)] for _ in range(bd+1)]
 for (i,j,k),c in q.terms():result[i][j][k]=int(c)
 return result
tables=[table(q) for q in f[:3]]
def evaluate(coefficients,x,p):
 out=0
 for c in reversed(coefficients):out=(out*x+c)%p
 return out
def evaluate_b_coefficients(tab,S0,A0,p):
 return [evaluate([evaluate(q,S0,p) for q in rows],A0,p) for rows in tab]
G2=a*a+a+s*b+s**3;G3=2*s*G2+a+s**3+s*s*b;G4=3*s*s*G2+b+s**4
X1=S.expand(-S.det(S.Matrix([G2,G3,G4]).jacobian([s,a,b])));X=[X1,G2+s*X1,G3+s*s*X1,G4+s**3*X1]
J=S.Matrix(f[:3]).jacobian([s,a,b]);tested=0;candidates=0
for p in [17,19,23,31]:
 for S0 in range(p):
  for A0 in range(p):
   coefficients=[evaluate_b_coefficients(tab,S0,A0,p) for tab in tables]
   for B0 in range(p):
    tested+=1
    if any(evaluate(q,B0,p) for q in coefficients):continue
    candidates+=1;point={s:S0,a:A0,b:B0};minor=int(J.subs(point).det())%p
    if minor==0 or int(f[3].subs(point))%p==0:continue
    target=[int(q.subs(point))%p for q in X]
    Y1,Y2,Y3,Y4=target;bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4;at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt
    R=S.Poly(at*at+at+t*bt+t**3+Y1*t-Y2,t,modulus=p);factor=S.Poly((t-S0)**5,t,modulus=p)
    Q,rem=R.div(factor)
    if not rem.is_zero or int(Q.eval(S0))%p==0 or Q.gcd(Q.diff()).degree()!=0:continue
    result={'status':'HENSEL_SEED_CANDIDATE','prime':p,'variables':['s','a','b'],'point':[S0,A0,B0],'jacobian_mod':minor,'target_mod':target,'quotient_coefficients_ascending_mod':[int(Q.nth(j))%p for j in range(8)],'tested_points':tested,'vanishing_candidates':candidates,'seconds':time.monotonic()-START,'boundary':'Only a seed. Separate exact source, Jacobian, quotient and squarefreeness replay plus the written Hensel/algebraicity argument establish a characteristic-zero occurrence.'}
    (HERE/'fivefold-seed.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));raise SystemExit(0)
print(json.dumps({'status':'NO_SEED_WITHIN_BOUNDS','tested':tested,'candidates':candidates,'seconds':time.monotonic()-START,'exclusion':False}))
