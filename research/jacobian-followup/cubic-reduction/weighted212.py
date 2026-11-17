#!/usr/bin/env python3
"""Complete weight(2,-1,-2) degree7 model and boundary-reduced discovery."""
from pathlib import Path
import json,hashlib,time,signal
import sympy as sp
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def enc(p,V):return [{'exponent':[int(k) for k in e],'coefficient':str(c)} for e,c in sorted(sp.Poly(p,*V).terms()) if c]
class TimeCap(Exception):pass
def alarm(*args):raise TimeCap()
def main():
 a,b,c,u,v=sp.symbols('a b c u v');support=[(1,0),(0,1),(2,0),(1,1),(3,0),(0,2)];params=[];forms=[]
 for prefix in ['f','g','h']:
  P=sp.symbols(prefix+'0:6');params.extend(P);forms.append(1+sum(p*u**i*v**j for p,(i,j) in zip(P,support)))
 k0,k1=sp.symbols('k0 k1');params.extend([k0,k1]);f,g,h=forms;sub={u:a*c,v:a*b*b};F=[sp.expand(a*f.subs(sub,simultaneous=True)),sp.expand(b*g.subs(sub,simultaneous=True)),sp.expand(c*h.subs(sub,simultaneous=True)+b*b*(k0+k1*a*b*b))];J=sp.Matrix(F).jacobian([a,b,c]);boundary=sp.expand(J.subs(b,0).det());fu,gu,hu=[p.subs(v,0) for p in forms];expected=sp.expand((gu*sp.diff(u*fu*hu,u)).subs(u,a*c));assert boundary==expected
 folder=HERE/'weighted212';folder.mkdir(exist_ok=True);full={'status':'complete-primary-weighted-model','source_sha256':sha(Path(__file__)),'parameter_names':list(map(str,params)),'weights':[2,-1,-2],'invariant_support':[list(e) for e in support],'outputs':[enc(p,[a,b,c]+params) for p in F],'boundary_determinant':enc(boundary,[a,c]+params),'boundary_product':enc(gu*sp.diff(u*fu*hu,u),[u]+params),'nonlinear_atoms':20};(folder/'full-model.json').write_text(json.dumps(full,indent=2)+'\n')
 A,B,C,D,E,G,H,I,K,L,M,rho=sp.symbols('A B C D E G H I K L M rho');P=[A,B,C,D,E,G,H,I,K,L,M];fr=1+A*v+B*u*v+C*v*v;gr=1+D*v+E*u*v+G*v*v;hr=1+H*v+I*u*v+K*v*v;FR=[sp.expand(a*fr.subs(sub,simultaneous=True)),sp.expand(b*gr.subs(sub,simultaneous=True)),sp.expand(c*hr.subs(sub,simultaneous=True)+b*b*(L+M*a*b*b))];delta=sp.expand(sp.Matrix(FR).jacobian([a,b,c]).det()-1);invariant={}
 for e,coef in sp.Poly(delta,a,b,c).terms():
  assert e[1]%2==0 and e[0]==e[2]+e[1]//2
  key=(e[2],e[1]//2);assert key not in invariant;invariant[key]=coef
 gens=[v for e,v in sorted(invariant.items())];model={'status':'complete-primary-boundary-reduced-model','source_sha256':sha(Path(__file__)),'parameter_names':list(map(str,P)),'outputs':[enc(p,[a,b,c]+P) for p in FR],'invariant_generators':[{'invariant_exponent':list(e),'polynomial':enc(p,P)} for e,p in sorted(invariant.items())],'nonlinear_atoms':11,'scope':'Boundary forces all invariant-u-only coefficients to zero over fields. Independent boundary and complete-source replay required; no counterexample or exclusion accepted from this primary model.'};(folder/'reduced-model.json').write_text(json.dumps(model,indent=2)+'\n');print('MODEL20 ->11; invariant equations',len(gens),flush=True)
 signal.signal(signal.SIGALRM,alarm);start=time.monotonic();records=[]
 def checkpoint(status):(folder/'checkpoint.json').write_text(json.dumps({'status':status,'source_sha256':sha(Path(__file__)),'model_sha256':sha(folder/'reduced-model.json'),'cases':records,'seconds':round(time.monotonic()-start,3)},indent=2)+'\n')
 for label,saturated in [('A',A),('B',B),('C',C)]:
  cp=folder/('first-'+label+'.json');eqs=gens+[rho*saturated-1];data={'source_sha256':sha(Path(__file__)),'parameter_names':list(map(str,P+[rho])),'model_sha256':sha(folder/'reduced-model.json'),'generators':[enc(p,P+[rho]) for p in eqs],'nonzero_first_coefficient':label};cp.write_text(json.dumps(data,indent=2)+'\n');record={'case':label,'input':cp.name,'input_sha256':sha(cp),'status':'running'};records.append(record);checkpoint('running');signal.alarm(10)
  try:
   tick=time.monotonic();basis=sp.groebner(eqs,*(P+[rho]),order='grevlex');record.update(status='primary-unit' if list(basis)==[1] else 'primary-proper-ideal',seconds=round(time.monotonic()-tick,3),basis=[enc(p,P+[rho]) for p in basis]);print('FIRST',label,record['status'],record['seconds'],[str(p) for p in basis][:10],flush=True)
  except TimeCap:record['status']='timeout';print('TIMEOUT',label,flush=True)
  finally:signal.alarm(0)
  checkpoint('running')
 checkpoint('complete')
if __name__=='__main__':main()
