#!/usr/bin/env python3
"""Independent exact audit of the remaining degree-seven negative-weight sector.

The proof uses written curve/divisibility lemmas and the explicitly pinned
Nakai-Baba plane low-degree theorem. Neither is proof-assistant formalized.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from verify import add,mul,variable,derivative,determinant3
from verify_power_single_compress import decode,scale
from weighted_verify_utils import substitution,expression,equal_generators,coefficient_buckets,power
from verify_weighted_unmarked7 import laurent,clear_monomial_denominator,linear_audit
HERE=Path(__file__).resolve().parent
CASES={'mu-zero':(0,-1),'nu-zero':(-3,0),'both-zero':(0,0)}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def jac(f,g,h):
 fu,fw=derivative(f,0),derivative(f,1);gu,gw=derivative(g,0),derivative(g,1);hu,hw=derivative(h,0),derivative(h,1)
 return add(add(mul(g,add(mul(fu,hw),mul(fw,hu),Q(-1))),scale(mul(h,add(mul(gu,fw),mul(gw,fu),Q(-1))),Q(2))),mul(f,add(mul(gu,hw),mul(gw,hu),Q(-1))))
def strip(poly,k):return {(0,)+e[1:]:c for e,c in poly.items() if e[0]==k}
def main():
 started=time.monotonic();root=HERE/'weighted-zero7';models={};cases=[];branch_cache={};boundary_cache={}
 for case,(mu,nu) in CASES.items():
  source=root/case/'input.json';d=load(source);names=d['variables'];n=len(names);assert n==21 and d['source_sha256']==sha(HERE/'weighted_zero7.py') and d['input_sha256']==sha(HERE/'weighted-unmarked7/input.json') and d['fixed_first_coefficients']==[mu,nu];atoms=d['atoms'];expected=set()
  for row,target in enumerate([1,-1,-2]):
   for a in range(8):
    for b in range(8):
     for c in range(8):
      if 2<=a+b+c<=7 and a-b-2*c==target:expected.add((row,(a,b,c)))
  assert len(atoms)==23 and {(a['row'],tuple(a['exponent'])) for a in atoms}==expected;N=3+n;F=[variable(i,N) for i in range(3)]
  for atom in atoms:
   e=list(atom['exponent'])+[0]*n
   if atom['coefficient'] in names:e[3+names.index(atom['coefficient'])]=1;c=Q(1)
   else:
    c=Q(atom['coefficient']);assert atom['row']==0 and ((atom['exponent']==[2,1,0] and c==mu) or (atom['exponent']==[3,0,1] and c==nu))
   F[atom['row']]=add(F[atom['row']],{tuple(e):c})
  det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),const(1,N),Q(-1));buckets=coefficient_buckets(det,3);prepared=[decode(p,n) for p in d['equations']];equal_generators(list(buckets.values()),prepared);assert all(a==b+2*c for a,b,c in buckets)
  qnames=['u','w']+names;nn=2+n;u,w=variable(0,nn),variable(1,nn);qgh=[add(add(const(1,nn),scale(u,Q(mu))),scale(w,Q(nu))),u,w]
  for atom in atoms:
   if 'quotient_exponent' not in atom:continue
   a,b=atom['quotient_exponent'];qgh[atom['row']]=add(qgh[atom['row']],mul(variable(qnames.index(atom['coefficient']),nn),mul(power(u,a,nn),power(w,b,nn))))
  for poly,text in zip(qgh,d['quotients']):assert poly==expression(text,qnames)
  qj=add(jac(*qgh),const(1,nn),Q(-1));assert qj=={(b,c)+e:co for (a,b,c),p in buckets.items() for e,co in p.items()};models[case]=dict(d=d,names=names,qgh=qgh,qj=qj,prepared=prepared,F=F)
  print('RECONSTRUCTED',case,flush=True)
 Z=['u','w','R','P','Q','T','mu','nu'];zz=len(Z);ff=expression('1+mu*u+nu*w+R*u**3+P*u**2+Q*u*w+T*w**2',Z);yy=expression('2*T*w+Q*u+nu',Z);cc=expression('(Q*u+nu)**2-4*T*(R*u**3+P*u**2+mu*u+1)',Z);assert add(power(yy,2,zz),cc,Q(-1))==scale(mul(variable(Z.index('T'),zz),ff),Q(4)) and derivative(ff,1)==yy
 tails=[name for name in names if name[0]!='f']
 def remap(newnames,mode):
  z={name:variable(i,len(newnames)) for i,name in enumerate(newnames)};zero={};R,P,Qp,T=[z.get(name,zero) for name in ['f1','f2','f3','f4']]
  if mode=='graph':vals=[R,P,Qp,zero]
  elif mode=='conic':vals=[zero,P,Qp,T]
  elif mode=='cubicgraph':vals=[R,P,zero,zero]
  elif mode=='quadraticgraph':vals=[zero,P,zero,zero]
  elif mode=='parabola':vals=[zero,expression('sigma**2*f4',newnames),expression('2*sigma*f4',newnames),T]
  elif mode=='linear':vals=[zero]*4
  else:raise AssertionError(mode)
  mapping=dict(zip(['f1','f2','f3','f4'],vals));return [mapping[name] if name in mapping else z[name] for name in names]
 def homcomp(poly,values,den,width):
  degree=max(e[0]+e[1] for e in poly);lift={e+(degree-e[0]-e[1],):c for e,c in poly.items()};return substitution(lift,values+[den],width),degree
 def branch(case,mode):
  if (case,mode) in branch_cache:return branch_cache[case,mode]
  model=models[case];f,g,h=model['qgh'];ip=root/case/mode/'input.json';inp=load(ip);assert inp['input_sha256']==sha(root/case/'input.json');P=inp['parameters']+inp['variables'];width=len(P);ren=remap(P,mode);base=[substitution(p,ren,width) for p in model['prepared']];base=[p for p in base if p];boundary=[]
  if mode not in ['graph','conic']:
   nt=['t']+P;ntn=len(nt);den=const(1,ntn)
   if mode=='cubicgraph':U=expression('t',nt);B=expression('1+f2*t**2+f1*t**3',nt)
   elif mode=='quadraticgraph':U=expression('t',nt);B=expression('1+f2*t**2',nt)
   elif mode=='linear':U=expression('t' if case=='mu-zero' else '1/3',nt);B=expression('1' if case=='mu-zero' else 't',nt)
   elif case=='mu-zero':U=expression('-1+t-f4*t**2',nt);B=expression('sigma*t-sigma*(-1+t-f4*t**2)',nt);den=expression('sigma',nt)
   else:U=expression('(1+f4*t**2)/3',nt);B=expression('t-sigma*(1+f4*t**2)/3',nt)
   coords=[U,B]+remap(nt,mode)
   for poly,limit in [(g,1),(h,2)]:
    val,degree=homcomp(poly,coords,den,ntn);boundary += [p for (k,),p in coefficient_buckets(val,1).items() if k>limit]
   vf=['t','v']+P;nv=len(vf);dd=const(1,nv)
   if mode=='parabola' and case=='mu-zero':Uf=expression('v-1+t-f4*t**2',vf);Bf=expression('sigma*t-sigma*(v-1+t-f4*t**2)',vf);dd=expression('sigma',vf);expected_jac=scale(dd,Q(-1))
   elif mode=='parabola':Uf=expression('(1+f4*t**2-v)/3',vf);Bf=expression('t-sigma*(1+f4*t**2-v)/3',vf);expected_jac=const(Q(1,3),nv)
   elif mode=='linear' and case=='nu-zero':Uf=expression('(1-v)/3',vf);Bf=expression('t',vf);expected_jac=const(Q(1,3),nv)
   else:Uf=expression('t',vf);Bf=expression({'cubicgraph':'1+f2*t**2+f1*t**3-v','quadraticgraph':'1+f2*t**2-v','linear':'1-v'}[mode],vf);expected_jac=const(-1,nv)
   ff,degree=homcomp(f,[Uf,Bf]+remap(vf,mode),dd,nv);assert ff==mul(power(dd,degree,nv),variable(1,nv));assert add(mul(derivative(Uf,0),derivative(Bf,1)),mul(derivative(Uf,1),derivative(Bf,0)),Q(-1))==expected_jac;equal_generators(boundary,[decode(p,width) for p in inp['boundary_equations']]);boundary_cache[case,mode]=boundary
  equal_generators(base+boundary,[decode(p,width) for p in inp['equations']]);branch_cache[case,mode]=(P,base+boundary);return P,base+boundary
 def record(cp,ip,mults,case,mode,target):
  cases.append(dict(case=case,mode=mode,target=target,input_sha256=sha(ip),certificate_sha256=sha(cp),multiplier_monomials=sum(map(len,mults)),max_multiplier_degree=max(sum(e) for p in mults for e in p)));print('REPLAYED',case,mode,flush=True)
 def field(case,mode,targettext):
  P,derived=branch(case,mode);width=len(P);folder=root/case/mode;cp=folder/'certificate.json';cert=load(cp);assert cert['source_sha256']==sha(HERE/'weighted_zero7_lift.py') and cert['input_sha256']==sha(folder/'input.json') and cert['singular_source_sha256']==sha(folder/'lift.sing') and cert['variables']==P;gen=[decode(p,width) for p in cert['generators']];equal_generators(derived,gen);target=decode(cert['target'],width);assert target==expression(targettext,P);mults=[decode(p,width) for p in cert['multipliers']];assert len(gen)==len(mults);result={}
  for a,b in zip(gen,mults):result=add(result,mul(a,b))
  assert result==target;record(cp,folder/'input.json',mults,case,mode,targettext)
 for c,m,target in [('mu-zero','graph','f3**16'),('mu-zero','conic','f4**2*(f2*f4-f3**2/4)**7'),('mu-zero','cubicgraph','f2*f1**3'),('mu-zero','quadraticgraph','f2**4'),('mu-zero','parabola','sigma**5*f4**9'),('nu-zero','graph','f3**5'),('nu-zero','conic','(f2*f4-f3**2/4)**7'),('both-zero','graph','f3**5'),('both-zero','conic','(f2*f4-f3**2/4)**3')]:field(c,m,target)
 def unit(case,mode,ip,gen):
  inp=load(ip);width=len(inp['variables']);cp=ip.parent/'certificate.json';cert=load(cp);assert cert['source_sha256']==sha(HERE/'weighted_zero7_unit.py') and cert['input_sha256']==sha(ip) and cert['singular_source_sha256']==sha(ip.parent/'lift.sing');mults=[decode(p,width) for p in cert['multipliers']];assert len(gen)==len(mults);result={}
  for a,b in zip(gen,mults):result=add(result,mul(a,b))
  assert result==const(1,width);record(cp,ip,mults,case,mode,'1')
 for case,mode in [('mu-zero','purecubic'),('nu-zero','parabola')]:
  base='cubicgraph' if mode=='purecubic' else mode;P,orig=branch(case,base);linear=boundary_cache[case,base];ip=root/case/'reduced'/mode/'input.json';inp=load(ip);assert inp['source_sha256']==sha(HERE/'weighted_zero7_reduce.py') and inp['input_sha256']==sha(root/case/base/'input.json');N=inp['variables'];width=len(N)
  if mode=='purecubic':
   old=P;P=[name for name in old if name!='f2'];vals=[{} if name=='f2' else variable(P.index(name),len(P)) for name in old];orig=[substitution(p,vals,len(P)) for p in orig];linear=[substitution(p,vals,len(P)) for p in linear]
  allowed={'f1'} if mode=='purecubic' else {'f4'};assert linear_audit(linear,P,N,inp['linear_solution'],allowed)==inp['linear_rank'];vals=[laurent(inp['linear_solution'][name],N) if name in inp['linear_solution'] else variable(N.index(name),width) for name in P];expected=[clear_monomial_denominator(substitution(p,vals,width),width) for p in orig];sat='f1' if mode=='purecubic' else 'f4**2';assert expression(inp['saturation'],N)==expression(sat,N);expected.append(expression('1-branch_inv*('+sat+')',N));gen=[decode(p,width) for p in inp['equations']];equal_generators(expected,gen);unit(case,mode,ip,gen)
 P,gen=branch('mu-zero','linear');unit('mu-zero','linear',root/'mu-zero/linear/input.json',gen)
 # Highest-w branch for f=1-3u. Boundary gives G2=K*v,H3=0.
 model=models['nu-zero'];f,g,h=model['qgh'];P,orig=branch('nu-zero','linear');width=len(P);ip=root/'nu-zero/high/generic/input.json';inp=load(ip);assert inp['source_sha256']==sha(HERE/'weighted_zero7_high.py') and inp['input_sha256']==sha(root/'nu-zero/linear/input.json');nt=['x','v']+P;ntn=len(nt);lookup={(a['row'],tuple(a['quotient_exponent'])):a['coefficient'] for a in model['d']['atoms'] if 'quotient_exponent' in a};kn=lookup[1,(0,2)];gn=lookup[1,(1,2)];hn=lookup[2,(0,3)];vals=remap(nt,'linear')
 for name,value in [(gn,scale(variable(nt.index(kn),ntn),Q(-3))),(hn,{})]:vals[names.index(name)]=value
 coords=[expression('(1-v)/3',nt),expression('x',nt)]+vals;G,H=[substitution(p,coords,ntn) for p in [g,h]];K=variable(nt.index(kn),ntn);vv=variable(1,ntn);assert strip(G,2)==mul(K,vv) and max(e[0] for e in H)<=2;G1=strip(G,1);H1,H2=strip(H,1),strip(H,2);assert max(e[1] for e in G1)<=2 and max(e[1] for e in H1)<=3 and max(e[1] for e in H2)<=1
 EE=add(mul(add(G,mul(vv,derivative(G,1))),derivative(H,0)),mul(derivative(G,0),add(scale(H,Q(2)),mul(vv,derivative(H,1)))),Q(-1));assert substitution(add(model['qj'],const(1,2+len(names))),coords,ntn)==scale(EE,Q(-3));assert strip(EE,3)==scale(mul(mul(K,power(vv,2,ntn)),derivative(H2,1)),Q(-2));rhs=mul(vv,add(add(scale(mul(derivative(G1,1),H2),Q(2)),mul(G1,derivative(H2,1)),Q(-1)),scale(mul(K,add(H1,mul(vv,derivative(H1,1)))),Q(-2))));assert strip(EE,2)==rhs
 atone=[const(1,ntn) if i==1 else variable(i,ntn) for i in range(ntn)];assert substitution(H1,atone,ntn)==const(1,ntn)
 # Reconstruct, rather than trust, all five high-branch linear constraints.
 cnames=['u','w']+P;cn=len(cnames);ren=remap(cnames,'linear');nh=substitution(h,[variable(0,cn),variable(1,cn)]+ren,cn);bucket=coefficient_buckets({(e[1],e[0])+e[2:]:co for e,co in nh.items()},1);linear=list(boundary_cache['nu-zero','linear'])
 for wpower,limit in [(2,0),(1,1)]:
  pu=bucket.get((wpower,),{})
  linear += [co for (u,),co in coefficient_buckets(pu,1).items() if u>limit]
 equal_generators(linear,[decode(p,len(P)) for p in inp['linear_equations']]);N=inp['variables'];nw=len(N);claim={name:text for name,text in inp['substitution'].items() if name not in N};linear_audit(linear,P,N,claim,set());sv=[expression(inp['substitution'][name],N) for name in P];expected=[substitution(p,sv,nw) for p in orig];extra='1-branch_inv*('+inp['nonzero']+')';expected.append(expression(extra,N));gen=[decode(p,nw) for p in inp['equations']];equal_generators(expected,gen);unit('nu-zero','high',ip,gen)
 # Audit the global-degree reductions in both written univariate lemmas.
 model=models['both-zero'];f,g,h=model['qgh'];lemma_receipts=[]
 for mode in ['parallel','univariate-u']:
  ip=root/'lemmas'/mode/'input.json';inp=load(ip);P=inp['parameters']+inp['variables'];nt=['t','l']+P;ntn=len(nt);U=expression('t' if mode=='parallel' else 'l',nt);B=expression('l-sigma*t' if mode=='parallel' else 't',nt);coords=[U,B]+remap(nt,'linear');GH=[substitution(poly,coords,ntn) for poly in [g,h]];linear=[]
  for poly,limit in zip(GH,[1,2]):
   for (t,l),co in coefficient_buckets(poly,2).items():
    if t>limit:assert l<2;linear.append(co)
  equal_generators(linear,[decode(p,len(P)) for p in inp['linear_equations']]);N=inp['parameters']+inp['free_tail_variables'];assert linear_audit(linear,P,N,inp['linear_solution'],set())==inp['linear_rank'];ns=['t','l']+N;nns=len(ns);coefs=[expression(inp['linear_solution'][name],ns) if name in inp['linear_solution'] else variable(ns.index(name),nns) for name in P];expand=[variable(0,nns),variable(1,nns)]+coefs;G,H=[substitution(poly,expand,nns) for poly in GH];assert max(e[0] for e in G)<=1 and max(e[0] for e in H)<=2;A,B=strip(G,1),strip(G,0);C,D,E=strip(H,2),strip(H,1),strip(H,0);assert max(e[1] for e in A)<=2 and max(e[1] for e in C)<=1;atz=[{} if i==1 else variable(i,nns) for i in range(nns)];assert not substitution(B,atz,nns) and not substitution(E,atz,nns)
  if mode=='parallel':assert substitution(D,atz,nns)==expression('-sigma',ns) and substitution(A,atz,nns)==expression('1-sigma*'+lookup[1,(0,1)],ns)
  else:assert substitution(D,atz,nns)==const(1,nns)
  lemma_receipts.append(dict(mode=mode,input_sha256=sha(ip),linear_rank=inp['linear_rank'],bounds_checked=True))
 # Universal coefficient and factorization identities behind the written proof.
 JN=['t','A','B','C','D','E','Ap','Bp','Cp','Dp','Ep','phi','phip'];jn=len(JN);G=expression('A*t+B',JN);H=expression('C*t**2+D*t+E',JN);Gu=expression('Ap*t+Bp',JN);Hu=expression('Cp*t**2+Dp*t+Ep',JN);J=add(mul(expression('phip',JN),add(mul(G,derivative(H,0)),scale(mul(H,derivative(G,0)),Q(-2)))),mul(expression('phi',JN),add(mul(Gu,derivative(H,0)),mul(derivative(G,0),Hu),Q(-1))));assert strip(J,2)==expression('-phi*(A*Cp-2*Ap*C)',JN);assert strip(J,1)==expression('-phip*A*D+phi*(Ap*D-A*Dp+2*Bp*C)+2*phip*B*C',JN)
 X=['phi','phip','k','B','Bp','D','E','Ep'];xn=len(X);a=laurent('phi*D/k',X);j=add(mul(expression('phip',X),add(expression('B*D',X),scale(mul(a,expression('E',X)),Q(-2)))),mul(expression('phi',X),add(expression('Bp*D',X),mul(a,expression('Ep',X)),Q(-1))));assert j==laurent('D*(phip*B+phi*Bp-(2*phi*phip*E+phi**2*Ep)/k)',X)
 # Constant f1 branch: degree<=4 plane Keller map; source/target rescalings
 # give its polynomial inverse globally because inverse Q has no linear p term.
 model=models['both-zero'];f,g,h=model['qgh'];N=['u','w']+tails;nn=len(N);coords=[variable(0,nn),variable(1,nn)]+remap(N,'linear');fg=[substitution(p,coords,nn) for p in [f,g,h]];assert fg[0]==const(1,nn) and max(e[0]+e[1] for p in fg[1:] for e in p)<=4;assert jac(*fg)==add(mul(derivative(fg[1],0),derivative(fg[2],1)),mul(derivative(fg[1],1),derivative(fg[2],0)),Q(-1));linear_g={e:c for e,c in fg[1].items() if e[0]+e[1]==1};linear_h={e:c for e,c in fg[2].items() if e[0]+e[1]==1};assert linear_g==expression('u+'+lookup[1,(0,1)]+'*w',N) and linear_h==variable(1,nn)
 importpath=root/'sources/plane-quartic-import.json';imp=load(importpath);assert imp['doi']=='10.18910/8846' and imp['sha256']==sha(root/'sources/nakai-baba-1977.pdf') and imp['primary_statement_visually_checked'] and not imp['reproved_here'];prior=load(HERE/'weighted-unmarked7/validation.json');assert prior['status']=='pass' and prior['source_sha256']==sha(HERE/'verify_weighted_unmarked7.py')
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),linear_audit_sha256=sha(HERE/'verify_weighted_unmarked7.py'),models={c:sha(root/c/'input.json') for c in CASES},cases=cases,lemma_cases=lemma_receipts,cleared_parameter_identities=9,unit_identities=4,multiplier_monomials=sum(c['multiplier_monomials'] for c in cases),max_multiplier_degree=max(c['max_multiplier_degree'] for c in cases),universal_three_variable_determinants_reconstructed=True,general_mu_nu_completed_square_checked=True,collision_equations_assumed=0,all_boundary_constraints_reconstructed=True,linear_parameterizations_certified_with_allowed_monomial_pivots=True,all_certificates_multiplied=True,highest_w_identities_and_bounds_checked=True,written_univariate_lemma_source_bounds_and_algebra_checked=True,constant_first_quotient_reduces_to_plane_degree_at_most4=True,plane_degree4_theorem_import_sha256=sha(importpath),plane_theorem_reproved=False,prior_nonzero_sector_validation_sha256=sha(HERE/'weighted-unmarked7/validation.json'),weighted_counterexample_minimum_degree=7,weighted_counterexample_minimum_terms_through_degree7=13,seconds=round(time.monotonic()-started,3),scope='Together withKCUB14, the full identity-linear weight(1,-1,-2) degree<=7 class consists of polynomial automorphisms withF1=a or diagonal conjugates of the established13-term map. Zero first-coefficient branches with nonconstant first quotient have no Keller completions. Written smooth-cubic/univariate degree-divisibility arguments and an explicit imported Nakai-Baba plane-degree4 theorem complement thirteen replayed exact identities. No global simplicity minimum is asserted.')
 (root/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
