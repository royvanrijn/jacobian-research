#!/usr/bin/env python3
"""Produce complete translated jets and two short selected-minor units; no solver."""
from pathlib import Path
from math import comb
import sympy as s
import json,hashlib,time
B=Path(__file__).resolve().parent;O=B/'raw13-affine-jet'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,vs):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(p),*vs).terms() if c]
def main():
 started=time.monotonic();O.mkdir(exist_ok=True);src=B/'raw13-specialization/input.json';atoms=json.loads(src.read_text())['atoms'];a,b,c=s.symbols('a b c');vs=(a,b,c);models={};matrices={}
 for degree in range(2,8):
  exps=[(i,j,degree-i-j) for i in range(degree,-1,-1) for j in range(degree-i,-1,-1)]
  M=s.Matrix([[s.expand(sum(s.Rational(z['coefficient'])*s.prod(comb(x,y)*t**(x-y) for x,y,t in zip(z['exponent'],e,vs)) for z in atoms if z['row']==row and all(x>=y for x,y in zip(z['exponent'],e)))) for e in exps] for row in range(3)])
  models[str(degree)]={'monomials':[list(e) for e in exps],'matrix':[[enc(z,vs) for z in M.row(i)] for i in range(3)]};matrices[degree]=M
 # Explicit two-variable units for complete rank3 coverage in a!=0,normalized to a=1.
 certificates={}
 M=matrices[2].subs(a,1);cols2=[(0,2,4),(1,2,4),(0,1,2),(0,2,3)];ms=[s.expand(M[:,cols].det()) for cols in cols2];P=s.cancel((ms[2]-s.Rational(3,4)*ms[3]-s.Rational(27,2))/(2*b+1));assert s.denom(P)==1
 mult=[4*P/1458,P/1458,s.Rational(2,27),-s.Rational(1,18)];assert s.expand(sum(x*y for x,y in zip(ms,mult)))==1
 certificates['degree2']={'columns':[list(k) for k in cols2],'minors':[enc(z,(b,c)) for z in ms],'multipliers':[enc(z,(b,c)) for z in mult]}
 M=matrices[3].subs(a,1);cols3=[(1,4,7),(1,6,7),(0,1,2)];ms=[s.expand(M[:,cols].det()) for cols in cols3];z,h=s.symbols('z h');inv=s.Rational(16,25)+s.Rational(256,125)*z;f3=(z-1)**2*((5-3*z)/2+h)**2/4;res=s.expand(inv*f3-1);r0=res.subs(h,0);qz=s.cancel(r0/z**2);qh=s.cancel((res-r0)/h);back={z:2*b+1,h:3*b+c-1};mult=[qz.subs(back)/648,-qh.subs(back)/288,inv.subs(back)/144];assert s.expand(sum(x*y for x,y in zip(ms,mult)))==1
 certificates['degree3']={'columns':[list(k) for k in cols3],'minors':[enc(z,(b,c)) for z in ms],'multipliers':[enc(z,(b,c)) for z in mult]}
 payload={'source_sha256':sha(Path(__file__)),'input_sha256':sha(src),'centre_variables':['a','b','c'],'degrees':models,'selected_minor_units':certificates,'scope':'All affine source centres,arbitrary invertible linear source/output frames. Rank proofs are written and finite algebra independently replayed.'};ip=O/'input.json';ip.write_text(json.dumps(payload,indent=2)+'\n')
 counts={k:sum(len(z) for z in v['multipliers']) for k,v in certificates.items()};cp={'status':'complete','source_sha256':sha(Path(__file__)),'input_sha256':sha(ip),'source_input_sha256':sha(src),'discovery_sha256':sha(O/'discovery.json'),'degrees':list(range(2,8)),'selected_minor_unit_counts':counts,'unit_multiplier_terms':sum(counts.values()),'new_solver_queries':0,'preliminary_primary_queries':3,'primary_answers_not_used_as_proof':True,'seconds':time.monotonic()-started};(O/'checkpoint.json').write_text(json.dumps(cp,indent=2)+'\n');print('COMPLETE six full translated jet matrices,two selected-minor units',counts,'no solver in producer',flush=True)
if __name__=='__main__':main()
