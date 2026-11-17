#!/usr/bin/env python3
"""Full small-ratio binary normal charts, exact primary discovery only."""
from pathlib import Path
import sympy as sp
import json,hashlib,time,signal
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,vs):return [{'exponent':[int(k) for k in e],'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*vs).terms()) if c]
class TimeCap(Exception):pass
def alarm(*args):raise TimeCap()
def main():
 a,b,c=sp.symbols('a b c');q,s,t,r,rho=sp.symbols('q s t r rho');q0,q1,q2,h=sp.symbols('q0 q1 q2 h');start=time.monotonic();records=[];folder=HERE/'equal-negative-smallweights';folder.mkdir(exist_ok=True);signal.signal(signal.SIGALRM,alarm)
 five=[('bb',b,b),('bc',b,c),('b-bc',b,b+c),('cc',c,c),('cb',c,b)]
 def checkpoint(status):
  (folder/'checkpoint.json').write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'cases':records,'seconds':round(time.monotonic()-start,3),'scope':'Primary17 normal charts for primitive ratios2/1 and1/2. Independent coefficient/unit replay and complete normalization proof required.'},indent=2)+'\n')
 for pw,qw in [(2,1),(1,2)]:
  for family in ['quartic-first','mixed']:
   cases=[('b',b,b),('c',c,c)] if (pw,qw,family)==(1,2,'mixed') else five
   for label,L,M in cases:
    if time.monotonic()-start>60:checkpoint('total-time-cap');return
    if family=='quartic-first':
     P=[q,s,t,r,rho];T=sp.Integer(0);sat=rho*q*s*t-1
     if pw==2:Q=q*L**2;S=s*L*M**2;U=t*L**2*M**3
     else:Q=q*L;S=s*M**2;U=t*L**2*M
    else:
     if pw==2:P=[q0,q1,q2,s,t,h,r,rho];Q=q0*b*b+q1*b*c+q2*c*c;S=s*L**2*M;T=t*L*M**3;U=h*L**3*M**2
     else:P=[q0,q1,s,t,h,r,rho];Q=q0*b+q1*c;S=s*L**2;T=t*L**2;U=h*L**3
     sat=rho*s*t*h-1
    R=sp.integrate(sp.expand(-sp.diff(S,c)-(qw+1)*Q),b)+r*c**(pw+1);F=[a+a**(qw+1)*Q+a**(2*qw+1)*T,b+a**qw*R+a**(2*qw)*U,c+a**qw*S];det=sp.expand(sp.Matrix(F).jacobian([a,b,c]).det()-1);terms=sorted(sp.Poly(det,a,b,c).terms());gens=[value for e,value in terms if value];keys=[list(e) for e,value in terms if value];gens.append(sat);keys.append('saturation');name='p%dq%d-%s-%s'%(pw,qw,family,label);casepath=folder/(name+'.json')
    model={'name':name,'weights':[pw,-qw,-qw],'family':family,'chart':label,'parameter_names':list(map(str,P)),'source_sha256':sha(Path(__file__)),'binary_lines':[str(L),str(M)],'Q':enc(Q,[b,c]+P),'S':enc(S,[b,c]+P),'T':enc(T,[b,c]+P),'U':enc(U,[b,c]+P),'R':enc(R,[b,c]+P),'outputs':[enc(p,[a,b,c]+P) for p in F],'generator_monomials':keys,'generators':[enc(g,P) for g in gens]};casepath.write_text(json.dumps(model,indent=2)+'\n');rec={'name':name,'parameters':len(P),'generators':len(gens),'model':casepath.name,'model_sha256':sha(casepath),'status':'running-groebner'};records.append(rec);checkpoint('running')
    signal.alarm(10)
    try:
     tick=time.monotonic();G=sp.groebner(gens,*P,order='grevlex');rec['seconds']=round(time.monotonic()-tick,3);rec['status']='primary-unit' if list(G)==[1] else 'primary-proper-ideal';rec['basis']=[enc(p,P) for p in G];print(name,rec['status'],rec['seconds'],[str(p) for p in G][:8],flush=True)
    except TimeCap:rec['status']='groebner-timeout';print(name,'TIMEOUT',flush=True)
    finally:signal.alarm(0)
    checkpoint('running')
 checkpoint('complete')
if __name__=='__main__':main()
