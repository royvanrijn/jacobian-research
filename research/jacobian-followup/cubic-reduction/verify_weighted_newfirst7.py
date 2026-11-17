#!/usr/bin/env python3
"""Independent Fraction replay for the full marked 23-atom degree-seven class.

The geometric pole-order lemma is a written proof, not a formal proof-system
artifact. This audit checks its algebraic chart identities and all finite
coefficient certificates used for the complementary branches.
"""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from verify import add,mul,variable,derivative,determinant3
from verify_power_single_compress import decode,scale
from weighted_verify_utils import substitution,expression,equal_generators,rank,coefficient_buckets,power
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def main():
 started=time.monotonic();source=HERE/'weighted-newfirst7/input.json';d=load(source);free=d['variables'];m=len(free);assert m==19 and d['degree']==7 and d['source_sha256']==sha(HERE/'weighted_newfirst.py');atoms=d['atoms'];assert len(atoms)==23
 expected=set()
 for row,target in enumerate([1,-1,-2]):
  for a in range(8):
   for b in range(8):
    for c in range(8):
     if 2<=a+b+c<=7 and a-b-2*c==target:expected.add((row,(a,b,c)))
 assert {(a['row'],tuple(a['exponent'])) for a in atoms}==expected
 names=list(d['substitution']);n=len(names);assert n==21;index={name:i for i,name in enumerate(names)};N=3+n;F=[variable(i,N) for i in range(3)]
 for atom in atoms:
  e=list(atom['exponent'])+[0]*n
  if atom['coefficient'] in index:e[3+index[atom['coefficient']]]=1;c=Q(1)
  else:c=Q(atom['coefficient']);assert atom['row']==0 and c in [Q(-3),Q(-1)]
  F[atom['row']]=add(F[atom['row']],{tuple(e):c})
 det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),const(1,N),Q(-1));buckets=coefficient_buckets(det,3);assert all(e[0]==e[1]+2*e[2] for e in buckets)
 def point(f):
  out={}
  for e,c in f.items():out=add(out,{e[3:]:c*(Q(-1)**e[0])*(Q(4)**e[2])})
  return out
 original=list(buckets.values())+[point(F[0]),add(point(F[1]),const(1,n),Q(-1))];linear=[p for p in original if p and max(sum(e) for e in p)<=1];assert all(sum(e)==1 for p in linear for e in p)
 values=[expression(d['substitution'][name],free) for name in names]
 for i,name in enumerate(names):
  if name in free:assert values[i]==variable(free.index(name),m)
 assert all(not substitution(p,values,m) for p in linear);assert rank([[p.get(tuple(int(k==j) for k in range(n)),Q(0)) for j in range(n)] for p in linear],n)==n-m==2
 native=[substitution(p,values,m) for p in original];prepared=[decode(p,m) for p in d['equations']];equal_generators(native,prepared)
 # Reconstruct f,g,h in u=ab,w=a^2c directly from every source atom.
 qnames=['u','w']+free;nn=2+m;u,w=[variable(i,nn) for i in range(2)];qf=[add(add(const(1,nn),scale(u,Q(-3))),scale(w,Q(-1))),u,w]
 for atom in atoms:
  if 'quotient_exponent' not in atom:continue
  b,c=atom['quotient_exponent'];coef=expression(d['substitution'][atom['coefficient']],qnames);qf[atom['row']]=add(qf[atom['row']],mul(coef,mul(power(u,b,nn),power(w,c,nn))))
 f,g,h=qf;fu,fw=[derivative(f,i) for i in range(2)];gu,gw=[derivative(g,i) for i in range(2)];hu,hw=[derivative(h,i) for i in range(2)]
 qj=add(add(mul(g,add(mul(fu,hw),mul(fw,hu),Q(-1))),scale(mul(h,add(mul(gu,fw),mul(gw,fu),Q(-1))),Q(2))),mul(f,add(mul(gu,hw),mul(gw,hu),Q(-1))));qj=add(qj,const(1,nn),Q(-1));nativeq={}
 for (a,b,c),poly in buckets.items():
  for e,co in substitution(poly,values,m).items():nativeq=add(nativeq,{(b,c)+e:co})
 assert qj==nativeq
 assert f==expression('1-3*u-w+(f2-4*f3+16*f4)*u**3+f2*u**2+f3*u*w+f4*w**2',qnames)
 # Completed-square chart for the nonzero cubic/quadratic branch.
 gnames=['u','w','P','Q','T','R'];gn=len(gnames);gf=expression('1-3*u-w+R*u**3+P*u**2+Q*u*w+T*w**2',gnames);gy=expression('2*T*w+Q*u-1',gnames);gc=expression('(Q*u-1)**2-4*T*(R*u**3+P*u**2-3*u+1)',gnames)
 assert add(power(gy,2,gn),gc,Q(-1))==scale(mul(variable(4,gn),gf),Q(4));assert derivative(gf,1)==gy
 # Degree arithmetic in the pole proof: k>=1 and 3k-1<=2k force k=1;
 # at infinity a zero order 3m+1 exceeds the total pole degree 2m.
 assert expression('(3*k-1)-2*k', ['k'])==expression('k-1',['k']);assert expression('(3*m+1)-2*m',['m'])==expression('m+1',['m'])
 out=HERE/'weighted-newfirst7/parameter-field';tails=[name for name in free if name not in ['f2','f3','f4']];cases=[];boundary_native={}
 def remap(target_names,mode):
  z={name:variable(i,len(target_names)) for i,name in enumerate(target_names)};zero={}
  if mode=='conic':mapping={'f2':expression('4*f3-16*f4',target_names),'f3':z['f3'],'f4':z['f4']}
  elif mode=='graph':mapping={'f2':z['f2'],'f3':z['f3'],'f4':zero}
  elif mode=='cubicgraph':mapping={'f2':z['f2'],'f3':zero,'f4':zero}
  elif mode=='parabolic':mapping={'f2':scale(z['f4'],Q(16)),'f3':scale(z['f4'],Q(8)),'f4':z['f4']}
  else:raise AssertionError(mode)
  return [mapping[name] if name in mapping else z[name] for name in free]
 for mode in ['conic','graph','cubicgraph-boundary','parabolic-boundary']:
  folder=out/mode;ip=folder/('input-v0.json' if mode=='parabolic-boundary' else 'input.json');inp=load(ip);certpath=folder/'certificate.json';cert=load(certpath);P=cert['variables'];width=len(P);base=mode.split('-')[0];rename=remap(P,base);derived=[substitution(p,rename,width) for p in prepared];derived=[p for p in derived if p]
  if mode.endswith('-boundary'):
   assert inp['input_sha256']==sha(source) and inp['source_sha256']==sha(out/'boundary-source-v0.py')
   nt=['t']+P;ntn=len(nt);vnames=['t','v']+P;vn=len(vnames)
   if base=='cubicgraph':U=expression('t',nt);B=expression('1-3*t+f2*t**2+f2*t**3',nt);Uf=expression('t',vnames);Bf=expression('1-3*t+f2*t**2+f2*t**3-v',vnames)
   else:U=expression('-1+t-f4*t**2',nt);B=expression('t-4*(-1+t-f4*t**2)',nt);Uf=expression('v-1+t-f4*t**2',vnames);Bf=expression('t-4*(v-1+t-f4*t**2)',vnames)
   assert add(mul(derivative(Uf,0),derivative(Bf,1)),mul(derivative(Uf,1),derivative(Bf,0)),Q(-1))==const(-1,vn)
   assert substitution(f,[Uf,Bf]+remap(vnames,base),vn)==variable(1,vn)
   coords=[U,B]+remap(nt,base);G,H=[substitution(p,coords,ntn) for p in [g,h]];assert not substitution(f,coords,ntn)
   wr=add(add(mul(G,derivative(H,0)),scale(mul(derivative(G,0),H),Q(-2))),const(1,ntn),Q(-1));assert wr==substitution(qj,coords,ntn)
   boundary=[]
   for poly,limit in [(G,1),(H,2)]:boundary += [p for (k,),p in coefficient_buckets(poly,1).items() if k>limit]
   equal_generators(boundary,[decode(p,width) for p in inp['boundary_equations']]);derived+=boundary;boundary_native[base]=(P,derived)
  else:assert inp['source_sha256']==sha(out/'source-v0.py') and inp['input_sha256']==sha(source)
  equal_generators(derived,[decode(p,width) for p in inp['equations']]);assert cert['input_sha256']==sha(ip) and cert['source_sha256']==sha(HERE/'weighted_newfirst7_lift.py') and cert['singular_source_sha256']==sha(folder/'lift.sing');assert cert['variables']==inp['parameters']+inp['variables']
  generators=[decode(p,width) for p in cert['generators']];assert generators==[decode(p,width) for p in inp['equations']];multipliers=[decode(p,width) for p in cert['multipliers']];assert len(generators)==len(multipliers);target=decode(cert['target'],width)
  expected_target={'conic':'f4**2*(f3-8*f4)**14','graph':'f3**16','cubicgraph-boundary':'f2**3','parabolic-boundary':'f4**8*(f4-7/4)'}[mode];assert target==expression(expected_target,P);result={}
  for a,b in zip(generators,multipliers):result=add(result,mul(a,b))
  assert result==target;cases.append(dict(mode=mode,input_sha256=sha(ip),certificate_sha256=sha(certpath),target=expected_target,multiplier_monomials=sum(map(len,multipliers)),max_multiplier_degree=max(sum(e) for p in multipliers for e in p)))
 # The sole nonzero parabolic exception is checked separately over Q.
 folder=out/'parabolic-exception';ip=folder/'input.json';inp=load(ip);cp=folder/'certificate.json';cert=load(cp);assert inp['specialization']=={'f4':'7/4'} and inp['source_sha256']==sha(HERE/'weighted_newfirst7_exception.py') and inp['input_sha256']==sha(out/'parabolic-boundary/input.json');P,derived=boundary_native['parabolic'];W=inp['variables'];width=len(W);rename=[const(Q(7,4),width)]+[variable(W.index(name),width) for name in P[1:]];expected=[substitution(p,rename,width) for p in derived];generators=[decode(p,width) for p in inp['generators']];equal_generators(expected,generators)
 assert cert['source_sha256']==sha(HERE/'weighted_newfirst7_exception.py') and cert['input_sha256']==sha(ip) and cert['singular_source_sha256']==sha(folder/'lift.sing');multipliers=[decode(p,width) for p in cert['multipliers']];assert len(multipliers)==len(generators);result={}
 for a,b in zip(generators,multipliers):result=add(result,mul(a,b))
 assert result==const(1,width);cases.append(dict(mode='parabolic-exception',input_sha256=sha(ip),certificate_sha256=sha(cp),target='1',multiplier_monomials=sum(map(len,multipliers)),max_multiplier_degree=max(sum(e) for p in multipliers for e in p)))
 # Exhaustive branch logic: RT!=0 impossible by the pole lemma; T0 forces Q0,
 # then P0; R0,T!=0 forces Q8T,P16T and the excluded nonzero exception7/4.
 # At P=Q=T=0 every extra first-output coefficient vanishes.
 zeros=[{} if name in ['f2','f3','f4'] else variable(i,m) for i,name in enumerate(free)]
 for name in ['f1','f2','f3','f4']:assert not substitution(expression(d['substitution'][name],free),zeros,m)
 tailpath=HERE/'weighted-tail7/validation.json';tail=load(tailpath);assert tail['status']=='pass' and tail['source_sha256']==sha(HERE/'verify_weighted_tail7.py') and tail['unique_nonlinear_total']==13;wide=load(HERE/'weighted-tail7/input.json');assert {(a['row'],tuple(a['exponent'])) for a in wide['tail_atoms']}=={(a['row'],tuple(a['exponent'])) for a in atoms if a['row']>0}
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),input_sha256=sha(source),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),tail_classification_sha256=sha(tailpath),nonlinear_atoms=23,coefficient_parameters=21,free_parameters=19,cases=cases,cleared_parameter_identities=4,exceptional_unit_identities=1,multiplier_monomials=sum(c['multiplier_monomials'] for c in cases),max_multiplier_degree=max(c['max_multiplier_degree'] for c in cases),universal_three_variable_determinant_reconstructed=True,quotient_jacobian_checked=True,complete_linear_parameterization_checked=True,completed_square_identity_checked=True,boundary_charts_and_wronskians_checked=True,boundary_constraints_reconstructed=True,all_used_certificates_multiplied=True,all_exceptional_branches_covered=True,new_first_coefficients_forced_zero=True,zero_slice_reduced_to_KCUB11=True,unique_nonlinear_total=13,geometric_pole_argument='Written proof in WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md; algebraic charts checked here. Not formalized.',seconds=round(time.monotonic()-started,3),scope='The full23-atom degree<=7 identity-linear weight(1,-1,-2) class with first a^2b/a^3c coefficients-3,-1 and collision(0,1,c0),(-1,1,4) contains only the established13-term representative. The smooth cubic boundary pole argument and five replayed identities cover every parameter branch over characteristic zero. No global raw13-term minimum or unrestricted degree7 classification is asserted.')
 (HERE/'weighted-newfirst7/validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
