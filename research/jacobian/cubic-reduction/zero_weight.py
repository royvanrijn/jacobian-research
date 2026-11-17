#!/usr/bin/env python3
"""Complete degree-seven one-zero/opposite-weight models and principal faces."""
from pathlib import Path
from itertools import product
from fractions import Fraction
import sympy as s
import json,hashlib,time
HERE=Path(__file__).resolve().parent;OUT=HERE/'zero-weight';OUT.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':list(e),'coefficient':str(c)} for e,c in s.Poly(s.expand(p),*V).terms() if c]
def main():
 start=time.monotonic();x,y,z,t=s.symbols('x y z t');cases=[(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,3),(2,5),(3,4)];records=[];cp=OUT/'checkpoint.json'
 def checkpoint(status):cp.write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'degree_cap':7,'pairs':[list(p) for p in cases],'models':records,'seconds':time.monotonic()-start},indent=2)+'\n')
 checkpoint('running')
 for p,q in cases:
  weights=[0,p,-q];basis=[(i,e) for i in range(3) for e in product(range(8),repeat=3) if 2<=sum(e)<=7 and p*e[1]-q*e[2]==weights[i]];names=list(s.symbols('c0:'+str(len(basis))));A=x;B=s.Integer(1);C=s.Integer(1);module=[];zero=[]
  for k,((i,e),a) in enumerate(zip(basis,names)):
   ey=e[1]-int(i==1);ez=e[2]-int(i==2);assert ey%q==0 and ez%p==0 and ey//q==ez//p;j=ey//q;term=a*x**e[0]*t**j
   if i==0:A+=term
   elif i==1:B+=term
   else:C+=term
   module.append({'output':i,'x_degree':e[0],'t_degree':j,'parameter':k})
   if j==0:zero.append(k)
  assert len(zero)==18;J=s.expand(s.diff(A,x)*(B*C+t*(q*C*s.diff(B,t)+p*B*s.diff(C,t)))-s.diff(A,t)*t*(q*C*s.diff(B,x)+p*B*s.diff(C,x)));F=[s.expand(A.subs(t,y**q*z**p)),s.expand(y*B.subs(t,y**q*z**p)),s.expand(z*C.subs(t,y**q*z**p))];kill={names[j]:0 for j in zero};forms=[s.expand(a.subs(kill)) for a in [A,B,C]];slopes=sorted({Fraction(z['x_degree']-int(z['output']==0),z['t_degree']) for z in module if z['t_degree'] and z['x_degree']-int(z['output']==0)>0});faces=[]
  for slope in slopes:
   r,ss=slope.numerator,slope.denominator;kept=[m['parameter'] for m in module if m['t_degree'] and Fraction(m['x_degree']-int(m['output']==0),m['t_degree'])<=slope];top=[m['parameter'] for m in module if m['t_degree'] and Fraction(m['x_degree']-int(m['output']==0),m['t_degree'])==slope];faces.append({'slope':str(slope),'auxiliary_weights':[q*ss,-r,0],'ray_exponent':[r,q*ss,p*ss],'retained_parameters':kept,'top_parameters':top})
  d={'status':'complete-primary-model','source_sha256':sha(Path(__file__)),'pair':[p,q],'weights':weights,'basis':[{'output':i,'exponent':list(e)} for i,e in basis],'parameter_names':list(map(str,names)),'module_terms':module,'forms':[enc(a,[x,t]+names) for a in [A,B,C]],'source_outputs':[enc(a,[x,y,z]+names) for a in F],'quotient_jacobian':enc(J,[x,t]+names),'native_generators':[enc(a,names) for a in s.Poly(J-1,x,t).coeffs()],'pure_x_parameters':zero,'boundary_reduced_forms':[enc(a,[x,t]+names) for a in forms],'positive_principal_faces':faces};ip=OUT/('p%dq%d.json'%(p,q));ip.write_text(json.dumps(d,indent=2)+'\n');records.append({'pair':[p,q],'input':ip.name,'input_sha256':sha(ip),'native_atoms':len(basis),'boundary_removed':len(zero),'remaining_atoms':len(basis)-len(zero),'native_generators':len(d['native_generators']),'principal_faces':len(faces)});checkpoint('running');print('MODEL',p,q,len(basis),'->',len(basis)-18,'faces',len(faces),flush=True)
 checkpoint('complete');print('COMPLETE',len(records),'models',sum(z['native_atoms'] for z in records),'atoms',sum(z['principal_faces'] for z in records),'faces',flush=True)
if __name__=='__main__':main()
