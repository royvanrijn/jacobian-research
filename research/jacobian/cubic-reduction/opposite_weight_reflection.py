#!/usr/bin/env python3
"""Complete opposite-magnitude primitive weighted models throughdegree7."""
from pathlib import Path
from math import gcd
import sympy as sp
import json,hashlib,time,signal
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':[int(k) for k in e],'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*V).terms()) if c]
class TimeCap(Exception):pass
def alarm(*args):raise TimeCap()
def main():
 a,b,c,u,v,rho=sp.symbols('a b c u v rho');folder=HERE/'opposite-weight-reflection';folder.mkdir(exist_ok=True);pairs=[(m,s-m) for s in range(3,9) for m in range(2,s) if gcd(m,s-m)==1];models=[];cases=[];start=time.monotonic();signal.signal(signal.SIGALRM,alarm)
 def checkpoint(status):(folder/'checkpoint.json').write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'models':models,'cases':cases,'seconds':round(time.monotonic()-start,3)},indent=2)+'\n')
 def build(m,q,support,ksupport):
  P=[];forms=[]
  for prefix in ['f','g','h']:
   vs=sp.symbols(prefix+'0:'+str(len(support)));P.extend(vs);forms.append(sp.Integer(1)+sum(p*u**i*v**j for p,(i,j) in zip(vs,support)))
  ks=sp.symbols('k0:'+str(len(ksupport)));P.extend(ks);k=sum(p*v**j for p,j in zip(ks,ksupport));f,g,h=forms;sub={u:a*c,v:a**q*b**m};F=[sp.expand(a*f.subs(sub,simultaneous=True)),sp.expand(b*g.subs(sub,simultaneous=True)),sp.expand(c*h.subs(sub,simultaneous=True)+a**(q-1)*b**m*k.subs(v,a**q*b**m))];return P,forms,k,F
 for m,q in pairs:
  s=m+q;support=[(i,j) for j in range(3) for i in range(4) if 0<2*i+s*j<=6];ksupport=[j for j in range(3) if s-1+s*j<=7];P,forms,k,F=build(m,q,support,ksupport);f,g,h=forms;J=sp.Matrix(F).jacobian([a,b,c]);boundary=sp.expand(J.subs(b,0).det());fu,gu,hu=[p.subs(v,0) for p in forms];bp=sp.expand(gu*sp.diff(u*fu*hu,u));assert boundary==sp.expand(bp.subs(u,a*c))
  rsupport=[e for e in support if e[1]];RP,rforms,rk,RF=build(m,q,rsupport,ksupport);rf,rg,rh=rforms;R=u*rh+v*rk;br=lambda x,y:sp.diff(x,u)*sp.diff(y,v)-sp.diff(x,v)*sp.diff(y,u);dictionary=sp.expand(rg*(rf*sp.diff(R,u)+R*sp.diff(rf,u)+q*v*br(R,rf))+m*v*(rf*br(R,rg)+R*br(rf,rg)));actual=sp.expand(sp.Matrix(RF).jacobian([a,b,c]).det());assert actual==sp.expand(dictionary.subs({u:a*c,v:a**q*b**m},simultaneous=True));delta=sp.Poly(dictionary-1,u,v);terms=sorted((e,p) for e,p in delta.terms() if p);name='m%dq%d'%(m,q);mp=folder/(name+'.json');payload={'status':'complete-primary-model','source_sha256':sha(Path(__file__)),'weights':[m,-q,-m],'support':[list(e) for e in support],'k_support':ksupport,'parameter_names':list(map(str,P)),'outputs':[enc(p,[a,b,c]+P) for p in F],'boundary_product':enc(bp,[u]+P),'boundary_determinant':enc(boundary,[a,c]+P),'reduced_support':[list(e) for e in rsupport],'reduced_parameter_names':list(map(str,RP)),'reduced_outputs':[enc(p,[a,b,c]+RP) for p in RF],'invariant_dictionary':enc(dictionary,[u,v]+RP),'invariant_generators':[{'invariant_exponent':list(e),'polynomial':enc(p,RP)} for e,p in terms]};mp.write_text(json.dumps(payload,indent=2)+'\n');models.append({'name':name,'input':mp.name,'input_sha256':sha(mp),'full_atoms':len(P),'reduced_atoms':len(RP),'invariant_equations':len(terms)});checkpoint('running');print('MODEL',name,len(P),'->',len(RP),'equations',len(terms),flush=True)
  if s not in [3,4]:continue
  for j in range(len(rsupport)):
   cname=name+'-f'+str(j);cp=folder/(cname+'.json');eqs=[p for e,p in terms]+[rho*RP[j]-1];data={'source_sha256':sha(Path(__file__)),'model_sha256':sha(mp),'model':mp.name,'parameter_names':list(map(str,RP+[rho])),'nonzero_first_index':j,'generators':[enc(p,RP+[rho]) for p in eqs]};cp.write_text(json.dumps(data,indent=2)+'\n');rec={'name':cname,'input':cp.name,'input_sha256':sha(cp),'status':'running'};cases.append(rec);checkpoint('running');signal.alarm(10)
   try:
    tick=time.monotonic();G=sp.groebner(eqs,*(RP+[rho]),order='grevlex');rec.update(status='primary-unit' if list(G)==[1] else 'primary-proper-ideal',seconds=round(time.monotonic()-tick,3),basis=[enc(p,RP+[rho]) for p in G]);print('FIRST',cname,rec['status'],rec['seconds'],flush=True)
   except TimeCap:rec['status']='timeout';print('TIMEOUT',cname,flush=True)
   finally:signal.alarm(0)
   checkpoint('running')
 checkpoint('complete')
if __name__=='__main__':main()
