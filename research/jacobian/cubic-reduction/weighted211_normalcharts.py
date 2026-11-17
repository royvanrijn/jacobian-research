#!/usr/bin/env python3
"""Checkpointed exact discovery for all binary power-form normal charts."""
from pathlib import Path
import sympy as sp
import json,hashlib,time,signal
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,vs):return [{'exponent':[int(k) for k in e],'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*vs).terms()) if c]
class TimeCap(Exception):pass
def alarm(*args):raise TimeCap()
def main():
 a,b,c=sp.symbols('a b c');q,s,t,r,rho=sp.symbols('q s t r rho');q0,q1,q2,h=sp.symbols('q0 q1 q2 h');start=time.monotonic();records=[];folder=HERE/'weighted211-normalcharts';folder.mkdir(exist_ok=True);signal.signal(signal.SIGALRM,alarm)
 cases=[('bb',b,b),('bc',b,c),('b-bc',b,b+c),('cc',c,c),('cb',c,b)]
 def checkpoint(status):
  (folder/'checkpoint.json').write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'cases':records,'seconds':round(time.monotonic()-start,3),'scope':'Exact primary discovery only. Full binary normalization and independent unit replay required before exclusion.'},indent=2)+'\n')
 for family in ['quartic-first','mixed']:
  for label,L,M in cases:
   if time.monotonic()-start>45:checkpoint('total-time-cap');return
   if family=='quartic-first':
    P=[q,s,t,r,rho];Q=q*L**2;S=s*L*M**2;U=t*L**2*M**3;T=sp.Integer(0);sat=rho*q*s*t-1
   else:
    P=[q0,q1,q2,s,t,h,r,rho];Q=q0*b*b+q1*b*c+q2*c*c;S=s*L**2*M;T=t*L*M**3;U=h*L**3*M**2;sat=rho*s*t*h-1
   R=sp.integrate(sp.expand(-sp.diff(S,c)-2*Q),b)+r*c**3;F=[a+a*a*Q+a**3*T,b+a*R+a*a*U,c+a*S];J=sp.Matrix(F).jacobian([a,b,c]);det=sp.expand(J.det()-1);terms=sorted(sp.Poly(det,a,b,c).terms());gens=[value for e,value in terms if value];keys=[list(e) for e,value in terms if value];gens.append(sat);keys.append('saturation');casepath=folder/(family+'-'+label+'.json')
   model={'family':family,'chart':label,'parameter_names':list(map(str,P)),'source_sha256':sha(Path(__file__)),'binary_lines':[str(L),str(M)],'Q':enc(Q,[b,c]+P),'S':enc(S,[b,c]+P),'T':enc(T,[b,c]+P),'U':enc(U,[b,c]+P),'R':enc(R,[b,c]+P),'outputs':[enc(p,[a,b,c]+P) for p in F],'generator_monomials':keys,'generators':[enc(g,P) for g in gens]};casepath.write_text(json.dumps(model,indent=2)+'\n');rec={'family':family,'chart':label,'parameters':len(P),'generators':len(gens),'model':casepath.name,'model_sha256':sha(casepath),'status':'running-groebner'};records.append(rec);checkpoint('running')
   signal.alarm(10)
   try:
    tick=time.monotonic();G=sp.groebner(gens,*P,order='grevlex');rec['seconds']=round(time.monotonic()-tick,3);rec['status']='primary-unit' if list(G)==[1] else 'primary-proper-ideal';rec['basis']=[enc(p,P) for p in G];print(family,label,rec['status'],rec['seconds'],[str(p) for p in G][:10],flush=True)
   except TimeCap:rec['status']='groebner-timeout';print(family,label,'TIMEOUT',flush=True)
   finally:signal.alarm(0)
   checkpoint('running')
 checkpoint('complete')
if __name__=='__main__':main()
