"""Check whether every retained root is an actual fivefold point over Q.

Positive-identity construction only; no new registered claim from this script.
"""
from pathlib import Path
import json,time,resource
import sympy as S
from flint import fmpq_poly,nmod_poly
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
HERE=Path(__file__).resolve().parent;START=time.monotonic();cert=json.loads((HERE/'certificate.json').read_text());source=json.loads((HERE/'quintuple-input.json').read_text());old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text())
s,a,b=S.symbols('s a b');h=fmpq_poly(cert['generic_retained']);co=cert['cofactor'];num=fmpq_poly(co['a_cofactor']);den=fmpq_poly(co['constant_cofactor']);g,u,v=den.xgcd(h);assert g.degree()==0
A=num*u/g[0]%h
print(json.dumps({'stage':'exact a recovered','seconds':time.monotonic()-START,'coefficient_bits':max(max(int(c.p).bit_length(),int(c.q).bit_length()) for c in A.coeffs())}),flush=True)
def evaluate(expression,values):
 variables=[s,a] if len(values)==2 else [s,a,b];q=S.Poly(expression,*variables,domain=S.ZZ);terms={ex:int(c) for ex,c in q.terms()}
 def recurse(data,level):
  if not data:return fmpq_poly([])
  if level==0:
   coeff=[0]*(max(ex[0] for ex in data)+1)
   for ex,c in data.items():coeff[ex[0]]=c
   return fmpq_poly(coeff)%h
  answer=fmpq_poly([])
  for j in range(max(ex[level] for ex in data),-1,-1):
   sub={ex[:level]:c for ex,c in data.items() if ex[level]==j};answer=(answer*values[level]+recurse(sub,level-1))%h
  return answer
 return recurse(terms,len(values)-1)
S0=fmpq_poly([0,1]);L=S.Poly(S.sympify(old['linear_equations'][0]),b);P=evaluate(L.nth(1),[S0,A]);Q=evaluate(L.nth(0),[S0,A]);g,u,v=P.xgcd(h);assert g.degree()==0;B=-Q*u/g[0]%h
print(json.dumps({'stage':'exact b recovered','seconds':time.monotonic()-START,'coefficient_bits':max(max(int(c.p).bit_length(),int(c.q).bit_length()) for c in B.coeffs())}),flush=True)
passed=[]
for row in source['equations']:
 value=evaluate(S.sympify(row['polynomial']),[S0,A,B]);assert not value,row['derivative'];passed.append(row['derivative']);print(json.dumps({'stage':'exact source identity','derivative':row['derivative'],'seconds':time.monotonic()-START}),flush=True)
p=cert['prime'];hp=nmod_poly(cert['generic_retained'],p);assert hp.degree()==117;g=hp.gcd(hp.derivative())
print(json.dumps({'status':'POSITIVE_IDENTITIES_CHECKED_OVER_Q','derivatives':passed,'degree':h.degree(),'squarefree_mod_prime_gcd_degree':g.degree(),'constant_nonzero':h[0]!=0,'seconds':time.monotonic()-START,'registered':False}),flush=True)
(HERE/'positive-rur-probe-receipt.json').write_text(json.dumps({'status':'POSITIVE_IDENTITIES_CHECKED_OVER_Q','derivatives':passed,'degree':h.degree(),'squarefree_mod_prime_gcd_degree':g.degree(),'constant_nonzero':h[0]!=0,'seconds':time.monotonic()-START,'boundary':'Separate proof checker and a characteristic-zero squarefreeness certificate are still required before any exact118 count is registered.'},indent=2)+'\n')
