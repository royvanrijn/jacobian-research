#!/usr/bin/env python3
"""Independent sparse-weight model audit and complete(1,-2,-3) exclusion."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from verify import add,mul,variable,derivative,determinant3
from verify_power_single_compress import decode,scale
from weighted_verify_utils import substitution,expression,equal_generators,rank,coefficient_buckets,power
HERE=Path(__file__).resolve().parent
WEIGHTS=[(1,-1,-3),(2,-1,-3),(1,-2,-3)]
CASES={'both':(-1,-1),'mu-zero':(0,-1),'nu-zero':(-1,0),'both-zero':(0,0)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def main():
 started=time.monotonic();root=HERE/'weighted-sparse7';models={};gates=[];certs=[]
 for weights in WEIGHTS:
  tag='_'.join(map(str,weights));ip=root/tag/'input.json';d=load(ip);free=d['variables'];m=len(free);names=list(d['substitution']);n=len(names);assert d['source_sha256']==sha(HERE/'weighted_sparse7_search.py') and d['weights']==list(weights) and d['max_degree']==7 and d['raw_nonlinear_bound']==12;expected=set()
  for row,target in enumerate(weights):
   for a in range(8):
    for b in range(8):
     for c in range(8):
      if 2<=a+b+c<=7 and sum(w*e for w,e in zip(weights,[a,b,c]))==target:expected.add((row,(a,b,c)))
  assert {(a['row'],tuple(a['exponent'])) for a in d['atoms']}==expected;N=3+n;F=[variable(i,N) for i in range(3)]
  for atom in d['atoms']:
   e=list(atom['exponent'])+[0]*n;e[3+names.index(atom['coefficient'])]=1;F[atom['row']]=add(F[atom['row']],{tuple(e):Q(1)})
  det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),const(1,N),Q(-1));native=list(coefficient_buckets(det,3).values());linear=[p for p in native if p and max(sum(e) for e in p)<=1];vals=[expression(d['substitution'][name],free) for name in names]
  for i,name in enumerate(names):
   if name in free:assert vals[i]==variable(free.index(name),m)
  assert all(not substitution(p,vals,m) for p in linear);assert rank([[p.get(tuple(int(i==j) for i in range(n)),Q(0)) for j in range(n)] for p in linear],n)==n-m==d['linear_rank'];derived=[substitution(p,vals,m) for p in native];prepared=[decode(p,m) for p in d['equations']];equal_generators(derived,prepared)
  first=[a for a in d['atoms'] if a['row']==0];assert all(a['exponent'][0]>=2 and sum(a['exponent'][1:])>0 for a in first);assert len({tuple(a['exponent'][1:]) for a in first})==len(first)
  # Exact a=0 tail restriction: identity plus only a possible b^(q/p) term.
  atzero=[]
  for row in [1,2]:
   p0={e:c for e,c in F[row].items() if e[0]==0};expected0=variable(row,N)
   for atom in d['atoms']:
    if atom['row']==row and atom['exponent'][0]==0:
     assert row==2 and atom['exponent'][2]==0 and atom['exponent'][1]*weights[1]==weights[2];e=list(atom['exponent'])+[0]*n;e[3+names.index(atom['coefficient'])]=1;expected0=add(expected0,{tuple(e):Q(1)})
   assert p0==expected0;atzero.append(True)
  gates.append(dict(weights=list(weights),input_sha256=sha(ip),nonlinear_slots=n,first_slots=len(first),negative_plane_triangular=True,first_nonlinearity_forces_collision='Written proof using the a=1 zero and triangular a=0 fibre. No collision variables in model.'))
  models[weights]=dict(d=d,free=free,prepared=prepared,F=F)
 model=models[1,-2,-3];d=model['d'];V=model['free'];assert len(V)==13 and d['linear_rank']==0;names=V;nn=2+len(V);qn=['u','w']+V;f=expression('1+c1*u+c2*w+c3*u**2',qn);g=expression('u+c4*w+c5*u**2+c6*u*w+c7*w**2+c8*u**3',qn);h=expression('w+c9*u**2+c10*u*w+c11*w**2+c12*u**3+c13*u**2*w',qn)
 fu,fw=derivative(f,0),derivative(f,1);gu,gw=derivative(g,0),derivative(g,1);hu,hw=derivative(h,0),derivative(h,1);J=add(add(scale(mul(g,add(mul(fu,hw),mul(fw,hu),Q(-1))),Q(2)),scale(mul(h,add(mul(gu,fw),mul(gw,fu),Q(-1))),Q(3))),mul(f,add(mul(gu,hw),mul(gw,hu),Q(-1))));buckets=coefficient_buckets(add(J,const(1,nn),Q(-1)),2);equal_generators(list(buckets.values()),model['prepared'])
 degree_pairs=[(m,n) for m in range(5) for n in range(5) if m+n>0 and (m+n==1 or 2*n==3*m)];assert degree_pairs==[(0,1),(1,0),(2,3)];tt=variable(0,1);GG=expression('t**2+1',['t']);HH=expression('t**3/3+t/2',['t']);assert add(scale(mul(GG,derivative(HH,0)),Q(2)),scale(mul(derivative(GG,0),HH),Q(-3)))==const(1,1)
 slices=root/'1_-2_-3/slices';cache={}
 def sliced(case,mode):
  if (case,mode) in cache:return cache[case,mode]
  folder=slices/case/mode;inp=load(folder/'input.json');P=inp.get('parameters',[])+inp['variables'];width=len(P);mu,nu=CASES[case];vals=[const(mu,width),const(nu,width),variable(P.index('c3'),width) if mode in ['field','boundary'] else {}]+[variable(P.index(name),width) for name in V[3:]];base=[substitution(p,vals,width) for p in model['prepared']];base=[p for p in base if p]
  if mode=='boundary':
   nt=['t']+P;nnt=len(nt);coords=[variable(0,nnt),expression('1-t+c3*t**2' if case=='both' else '1+c3*t**2',nt)]+[const(mu,nnt),const(nu,nnt),variable(nt.index('c3'),nnt)]+[variable(nt.index(name),nnt) for name in V[3:]];G,H=[substitution(p,coords,nnt) for p in [g,h]];assert max(e[0] for e in G)<=4 and max(e[0] for e in H)<=4;boundary=[]
   for poly,limit in [(G,2),(H,3)]:boundary += [p for (k,),p in coefficient_buckets(poly,1).items() if k>limit]
   equal_generators(boundary,[decode(p,width) for p in inp['boundary_equations']]);base+=boundary
   # At f=0 the quotient Jacobian is 2G Hprime-3Gprime H.
   wr=add(scale(mul(G,derivative(H,0)),Q(2)),scale(mul(derivative(G,0),H),Q(-3)));assert wr==substitution(J,coords,nnt)
  equal_generators(base,[decode(p,width) for p in inp['equations']]);cache[case,mode]=(P,base);return P,base
 for case,target in [('both','c3**7*(c3+810/2021)**2'),('mu-zero','c3**6')]:
  P,derived=sliced(case,'boundary');folder=slices/case/'boundary';cp=folder/'certificate.json';cert=load(cp);width=len(P);assert cert['source_sha256']==sha(HERE/'weighted_123_lift.py') and cert['input_sha256']==sha(folder/'input.json') and cert['singular_source_sha256']==sha(folder/'lift.sing') and cert['variables']==P;gen=[decode(p,width) for p in cert['generators']];equal_generators(derived,gen);t=decode(cert['target'],width);assert t==expression(target,P);mult=[decode(p,width) for p in cert['multipliers']];assert len(gen)==len(mult);res={}
  for a,b in zip(gen,mult):res=add(res,mul(a,b))
  assert res==t;certs.append(dict(case=case,mode='boundary',input_sha256=sha(folder/'input.json'),certificate_sha256=sha(cp),target=target,multiplier_monomials=sum(map(len,mult)),max_multiplier_degree=max(sum(e) for p in mult for e in p)))
 for case,mode in [('both','zero'),('mu-zero','zero'),('nu-zero','zero'),('both','exception')]:
  folder=slices/case/mode;inp=load(folder/'input.json');P=inp['variables'];width=len(P)
  if mode=='zero':_,gen=sliced(case,mode)
  else:
   old,derived=cache[case,'boundary'];vals=[const(-Q(810,2021),width) if name=='c3' else variable(P.index(name),width) for name in old];gen=[substitution(p,vals,width) for p in derived];assert inp['specialization']=={'c3':'-810/2021'} and inp['input_sha256']==sha(slices/'both/boundary/input.json');equal_generators(gen,[decode(p,width) for p in inp['equations']]);gen=[decode(p,width) for p in inp['equations']]
  cp=folder/'certificate.json';cert=load(cp);assert cert['source_sha256']==sha(HERE/'weighted_123_unit.py') and cert['input_sha256']==sha(folder/'input.json') and cert['singular_source_sha256']==sha(folder/'lift.sing');mult=[decode(p,width) for p in cert['multipliers']];assert len(gen)==len(mult);res={}
  for a,b in zip(gen,mult):res=add(res,mul(a,b))
  assert res==const(1,width);certs.append(dict(case=case,mode=mode,input_sha256=sha(folder/'input.json'),certificate_sha256=sha(cp),target='1',multiplier_monomials=sum(map(len,mult)),max_multiplier_degree=max(sum(e) for p in mult for e in p)))
 # Algebra of the nu=0 univariate branch. A root of phi forces G,H linear
 # in w because degree_w<=2 and 2degH=3degG cannot cancel leading terms.
 XX=['w','phi','phip','A','B','C','D','Ap','Bp','Cp','Dp'];nx=len(XX);GG=expression('A*w+B',XX);HH=expression('C*w+D',XX);jj=add(add(scale(mul(expression('phip',XX),mul(GG,derivative(HH,0))),Q(2)),scale(mul(expression('phip',XX),mul(HH,derivative(GG,0))),Q(-3))),mul(expression('phi',XX),add(mul(expression('Ap*w+Bp',XX),derivative(HH,0)),mul(derivative(GG,0),expression('Cp*w+Dp',XX)),Q(-1))));assert {tuple([0]+list(e[1:])):c for e,c in jj.items() if e[0]==1}==expression('-phip*A*C+phi*(Ap*C-A*Cp)',XX)
 X=['t','mu','P','A','B','C','D','Ap','Bp','Cp','Dp'];xn=len(X);phi=expression('1+mu*t+P*t**2',X);phip=derivative(phi,0);jc1=add(scale(mul(phip,expression('A*C',X)),Q(-1)),mul(phi,expression('Ap*C-A*Cp',X)));# phi*C=k*A is the integrated coefficient equation.
 # After A=phi*C/k, the constant equation factors C times a polynomial.
 Y=['phi','phip','A','B','C','D','Bp','Dp','k'];yn=len(Y);a=expression('phi*C',Y);j0=add(mul(expression('phip',Y),add(scale(expression('B*C*k',Y),Q(2)),scale(mul(a,expression('D',Y)),Q(-3)))),mul(expression('phi',Y),add(expression('Bp*C*k',Y),mul(a,expression('Dp',Y)),Q(-1))));assert j0==expression('C*(k*(2*phip*B+phi*Bp)-3*phi*phip*D-phi**2*Dp)',Y)
 Z=['t','mu','P'];zz=len(Z);Phi=expression('1+mu*t+P*t**2',Z);integral=expression('t+mu*t**2/2+P*t**3/3',Z);assert derivative(integral,0)==Phi;assert max(e[0] for e in integral)==3 and max(e[0] for e in power(Phi,2,zz))==4
 # f=1: g,h have ordinary degree<=3 and the pinned plane theorem fromKCUB15
 # applies. Their inverse (U,V) has zero constants and V no first-target
 # linear term. Weighted degrees(2,3) then clear a^-2 and a^-3 globally.
 nq=['u','w']+V[3:];nqn=len(nq);sub=[variable(0,nqn),variable(1,nqn)]+[{}]*3+[variable(nq.index(name),nqn) for name in V[3:]];gf,hf=[substitution(p,sub,nqn) for p in [g,h]];assert max(e[0]+e[1] for p in [gf,hf] for e in p)<=3;assert {e:c for e,c in gf.items() if e[0]+e[1]==1}==expression('u+c4*w',nq) and {e:c for e,c in hf.items() if e[0]+e[1]==1}==variable(1,nqn)
 impath=HERE/'weighted-zero7/sources/plane-quartic-import.json';imp=load(impath);assert imp['doi']=='10.18910/8846' and imp['sha256']==sha(HERE/'weighted-zero7/sources/nakai-baba-1977.pdf') and not imp['reproved_here'];prior=load(HERE/'weighted-zero7/validation.json');assert prior['status']=='pass' and prior['source_sha256']==sha(HERE/'verify_weighted_zero7.py')
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),gates=gates,certificates=certs,certificate_count=6,multiplier_monomials=sum(c['multiplier_monomials'] for c in certs),max_multiplier_degree=max(c['max_multiplier_degree'] for c in certs),all_three_universal_models_reconstructed=True,automatic_collision_gate_restrictions_checked=True,quotient_jacobian_coefficients=[2,3],all_six_identities_multiplied=True,necessary_boundary_equations_reconstructed=True,balanced_degree23_control_checked=True,exceptional_parameter_values_exactly_resolved=True,univariate_integral_divisibility_algebra_checked=True,plane_degree_bound=3,plane_theorem_import_sha256=sha(impath),plane_theorem_reproved=False,complete_weights=[1,-2,-3],ordinary_degree_bound=7,all_maps_automorphisms=True,other_two_weight_searches_closed=False,seconds=round(time.monotonic()-started,3),scope='Every identity-linear weight(1,-2,-3) Keller map of ordinary degree<=7 is a polynomial automorphism. Nonlinear first components are excluded by four normalized coefficient cases, six exact identities and a written univariate antiderivative-divisibility argument; the f=1 branch uses the attributed plane degree<=3 consequence of the pinned Nakai-Baba theorem and a polynomial weighted inverse. The general automatic-collision lemma is written proof. The(1,-1,-3) and(2,-1,-3) bounded searches remain UNKNOWN/timeouts, not exclusions; global minima remain open.')
 (root/'1_-2_-3/validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
