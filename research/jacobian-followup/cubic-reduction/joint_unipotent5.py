#!/usr/bin/env python3
"""Complete five-parameter joint linear conjugacy coefficient model."""
from pathlib import Path
import json,hashlib,time
from collections import Counter
import sympy as sp
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def encode(p,variables):return [{'exponent':list(e),'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*variables).terms()) if c]
def main():
 start=time.monotonic();folder=HERE/'joint-unipotent5';folder.mkdir(exist_ok=True);path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());X=sp.symbols('x1:13');P=sp.symbols('u v w z t');u,v,w,z,t=P
 F=sp.Matrix([sum(sp.Rational(a['coefficient'])*sp.prod(x**k for x,k in zip(X,a['exponent'])) for a in row) for row in raw['coefficients']]);S=sp.eye(12)
 for i,j,p in [(1,6,u),(1,9,v),(2,8,w),(2,7,z+w*t),(8,7,t)]:S[i,j]=p
 A=S-sp.eye(12);assert A**3==sp.zeros(12);Sinv=sp.eye(12)-A+A**2;assert S*Sinv==Sinv*S==sp.eye(12)
 sub=dict(zip(X,list(Sinv*sp.Matrix(X))));R=[sp.expand(p) for p in F.subs(sub,simultaneous=True)];G=[sp.expand(p) for p in S*sp.Matrix(R)];slots=[];factorlist=[];inventory=Counter()
 for i,p in enumerate(G):
  for e,c in sorted(sp.Poly(p,*X).terms()):
   if sum(e)<2:continue
   coefficient,factors=sp.factor_list(c,*P);atomids=[]
   for f,k in factors:
    f=sp.Poly(f,*P).monic().as_expr();fc=sp.Poly(c,*P)
    if f not in factorlist:factorlist.append(f)
    atomids.append({'factor':factorlist.index(f),'power':int(k)})
   # Recover the scalar after normalizing each factor to monic.
   product=sp.prod(factorlist[a['factor']]**a['power'] for a in atomids);scalar=sp.cancel(c/product);assert not scalar.free_symbols
   inventory[tuple((a['factor'],a['power']) for a in atomids)]+=1
   slots.append({'output':i+1,'source_exponent':list(e),'parameter_polynomial':encode(c,P),'factor_scalar':str(scalar),'factors':atomids})
 source_columns={}
 for i,p in enumerate(R):
  for e,c in sp.Poly(p,*X).terms():
   if sum(e)>=2:source_columns.setdefault(e,{})[i]=c
 source_records=[];source_inventory=Counter()
 for e,cs in sorted(source_columns.items()):
  gcd=sp.Poly(next(iter(cs.values())),*P)
  for c in cs.values():gcd=gcd.gcd(sp.Poly(c,*P))
  gcd=gcd.monic().as_expr();quotients={i:sp.cancel(c/gcd) for i,c in cs.items()};mons=sorted(set().union(*(set(sp.Poly(q,*P).monoms()) for q in quotients.values()))|{(0,)*5});qs=list(quotients.values());B=sp.Matrix([[sp.Poly(q,*P).coeff_monomial(m) for q in qs] for m in mons]);b=sp.Matrix([sp.Integer(not any(m)) for m in mons]);sol=sp.linsolve((B,b));unit=None
  if sol is not sp.EmptySet and sol:
   vector=next(iter(sol));free=set().union(*(a.free_symbols for a in vector));vector=[a.subs({v:0 for v in free}) for a in vector];assert sp.expand(sum(a*q for a,q in zip(vector,qs)))==1;unit=[{'output':i+1,'coefficient':str(a)} for i,a in zip(quotients,vector) if a]
  _,fs=sp.factor_list(gcd,*P);ids=[]
  for f,k in fs:
   f=sp.Poly(f,*P).monic().as_expr()
   if f not in factorlist:factorlist.append(f)
   ids.append({'factor':factorlist.index(f),'power':int(k)})
  scalar=sp.cancel(gcd/sp.prod(factorlist[a['factor']]**a['power'] for a in ids));source_inventory[tuple((a['factor'],a['power']) for a in ids)]+=1
  source_records.append({'source_exponent':list(e),'gcd_polynomial':encode(gcd,P),'factor_scalar':str(scalar),'factors':ids,'quotients':{str(i+1):encode(q,P) for i,q in quotients.items()},'quotient_unit':unit})
 print('SOURCE COLUMNS',len(source_records),'unresolved',sum(a['quotient_unit'] is None for a in source_records),flush=True)
 print('SOURCE INVENTORY',[(count,ids) for ids,count in sorted(source_inventory.items())],flush=True)
 values=[0,0,0,0,3];special=dict(zip(P,values));Sp=S.subs(special);Ip=Sinv.subs(special);B=sp.eye(12)
 for i,j,c in [(1,6,-3),(1,9,-12),(2,8,-4),(2,7,-20),(8,7,3),(10,9,3)]:B[i,j]=c
 Rpoint=sp.Matrix([sp.expand(p.subs(special)) for p in R]);examples=[]
 for label,out in [('identity-linear-t3',sp.Matrix(G).subs(special)),('free-output-t3',B*Rpoint)]:
  out=[sp.expand(p) for p in out];nonlinear=[sum(sum(e)>=2 for e,c in sp.Poly(p,*X).terms()) for p in out];linear=[sum(sum(e)==1 for e,c in sp.Poly(p,*X).terms()) for p in out]
  Aout=Sp if label=='identity-linear-t3' else B;points=[[sp.Rational(v) for v in point] for point in raw['points']];pp=[list(Sp*sp.Matrix(point)) for point in points];images=[[p.subs(dict(zip(X,point))) for p in out] for point in pp];assert images[0]==images[1]
  examples.append({'label':label,'parameters':list(map(str,values)),'target_matrix':[[str(v) for v in Aout.row(i)] for i in range(12)],'source_matrix':[[str(v) for v in Ip.row(i)] for i in range(12)],'coefficients':[encode(p,X) for p in out],'points':[[str(v) for v in point] for point in pp],'common_image':[str(v) for v in images[0]],'nonlinear_counts':nonlinear,'linear_counts':linear,'nonlinear_total':sum(nonlinear),'linear_total':sum(linear),'total':sum(nonlinear)+sum(linear),'nonlinear_source_columns':len(set().union(*(set(e for e,c in sp.Poly(p,*X).terms() if sum(e)>=2) for p in out)))})
 assert [a['nonlinear_total'] for a in examples]==[39,34]
 result={'status':'complete-symbolic-family','source_sha256':sha(Path(__file__)),'input_sha256':sha(path),'parameter_names':list(map(str,P)),'matrix_parameterization':'S27=u,S2,10=v,S3,9=w,S3,8=z+w*t,S9,8=t; other diagonal1','source_variables':raw['variables'],'source_matrix':[[encode(Sinv[i,j],P) for j in range(12)] for i in range(12)],'target_matrix':[[encode(S[i,j],P) for j in range(12)] for i in range(12)],'coefficients':[encode(p,X+P) for p in G],'nonlinear_slots':slots,'right_source_coefficients':[encode(p,X+P) for p in R],'source_columns':source_records,'examples':examples,'source_factor_inventory':[{'count':count,'factors':[{'factor':j,'power':k} for j,k in ids]} for ids,count in sorted(source_inventory.items())],'factors':[{'id':j,'expression':str(f),'polynomial':encode(f,P)} for j,f in enumerate(factorlist)],'factor_inventory':[{'count':count,'factors':[{'factor':j,'power':int(k)} for j,k in ids]} for ids,count in sorted(inventory.items())],'generic_nonlinear_total':len(slots),'generic_degree':max(sum(e) for p in G for e,c in sp.Poly(p,*X).terms()),'seconds':round(time.monotonic()-start,3),'scope':'All scalar parameters, fixed twelve-variable input. No minimum or improved witness accepted without independent replay.'};(folder/'model.json').write_text(json.dumps(result,indent=2)+'\n')
 print('MODEL',len(slots),'slots',len(factorlist),'factors',result['seconds'],'seconds',flush=True)
 for j,f in enumerate(factorlist):print('FACTOR',j,f,flush=True)
 print('INVENTORY',[(count,ids) for ids,count in sorted(inventory.items())],flush=True)
if __name__=='__main__':main()
