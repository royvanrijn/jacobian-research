#!/usr/bin/env python3
"""Exact elementary quadratic-conjugacy screen; all parameter roots retained."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb
import hashlib,json,time
import sympy as s
from verify import decode,add,mul
from optimize12 import serialize
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def family(F,i,j,k):
 n=12;rows=[]
 for f in F:
  out={}
  for e,c in f.items():
   for a in range(e[i]+1):
    ee=list(e)+[a];ee[i]-=a;ee[j]+=a;ee[k]+=a;ee=tuple(ee);out[ee]=out.get(ee,Q(0))+c*comb(e[i],a)
  rows.append({e:c for e,c in out.items() if c})
 product=mul(rows[j],rows[k]);product={e[:-1]+(e[-1]+1,):c for e,c in product.items()};rows[i]=add(rows[i],product,Q(-1))
 grouped=[]
 for row in rows:
  by={}
  for e,c in row.items():by.setdefault(e[:-1],{})[e[-1]]=c
  grouped.append(by)
 return grouped

def main():
 source=HERE/'cubic12-integer-witness.json';d=json.loads(source.read_text());F=[decode(p,12) for p in d['coefficients']];out=HERE/'quadratic-conjugacy12';out.mkdir(exist_ok=True);started=time.monotonic();t=s.Symbol('t');cases=[];best=36;cache={}
 def poly(p):return s.Poly(sum(s.Rational(c)*t**a for a,c in p.items()),t,domain=s.QQ)
 for i in range(12):
  for j,k in combinations_with_replacement([a for a in range(12) if a!=i],2):
   fam=family(F,i,j,k);high=[(row,e,p) for row,f in enumerate(fam) for e,p in f.items() if sum(e)>3];gcd=None;selected=[]
   for row,e,p in high:
    pp=poly(p);gcd=pp.monic() if gcd is None else s.gcd(gcd,pp).monic();selected.append(dict(row=row,exponent=list(e),coefficient=[[a,str(c)] for a,c in sorted(p.items())]))
    if len(gcd.terms())==1:break
   row=dict(i=i,j=j,k=k,high_coefficient_count=len(high),selected_high_coefficients=selected)
   if gcd is not None:
    assert gcd.eval(0)==0;row['high_gcd']=[[e[0],str(c)] for e,c in gcd.terms()];roots=set(gcd.ground_roots())-{s.Rational(0)};row['nonzero_rational_cubic_parameters']=list(map(str,sorted(roots)))
   else:
    row['high_gcd']=[];roots={s.Rational(-2),s.Rational(-1),s.Rational(-1,2),s.Rational(1,2),s.Rational(1),s.Rational(2)}
    for f in fam:
     for e,p in f.items():
      if sum(e)<2:continue
      pp=poly(p)
      if pp.degree()>0:roots.update(pp.ground_roots())
    roots.discard(s.Rational(0));row['rational_support_cancellation_parameters']=list(map(str,sorted(roots)))
   local=36
   for root in roots:
    rr=Q(root);G=[{e:v for e,p in f.items() if (v:=sum((c*rr**a for a,c in p.items()),Q(0)))} for f in fam];assert max(sum(e) for f in G for e in f)<=3;score=sum(sum(e)>=2 for f in G for e in f);local=min(local,score)
    if score<best:
     best=score;payload=dict(status='candidate',source_sha256=sha(Path(__file__)),input_sha256=sha(source),move=dict(i=i,j=j,k=k,t=str(root)),nonlinear_total=score,coefficients=serialize(G),scope='Requires independent conjugacy/determinant/collision replay before acceptance.');(out/'candidate.json').write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n');print('BEST',best,payload['move'],flush=True)
   row['line_minimum']=local;cases.append(row)
  receipt=dict(status='running',source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases,completed_axes=i+1,lines=len(cases),best=best,seconds=round(time.monotonic()-started,3),scope='Single elementary quadratic conjugacies H^-1 F H with H_i=x_i+t*x_j*x_k, i distinct from j,k, j<=k. Exact high-coefficient GCD handles all complex cubic parameters; support optimization on all-cubic lines covers rational cancellation parameters. No multistage/multimonomial/global optimum claim.');(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('AXIS',i+1,'LINES',len(cases),'BEST',best,flush=True)
 receipt['status']='complete';(out/'checkpoint.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('DONE',len(cases),'allcubic',sum(not r['high_gcd'] for r in cases),'nonzero-parameters',sum(bool(r.get('nonzero_rational_cubic_parameters')) for r in cases),'BEST',best,flush=True)
if __name__=='__main__':main()
