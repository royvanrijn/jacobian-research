#!/usr/bin/env python3
"""Fixed generic-basis incidence gate; no exceptional specialization input."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import resource
import time
import sympy as S
from verify_q80_branch_trace_specialization import quotient, rank
from audit_q80_branch_trace_specialization import prime

ROOT=Path(__file__).resolve().parents[2]
SOURCES={
 'q80':'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json',
 'r17':'artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/input.json',
 'curve302':'artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json',
 'x1092':'artifacts/generated-results/elliptic-curves/x1092_class1_realization_compact_parent_v1.json'}
t=S.Symbol('t')
def pol(a):return S.Poly.from_list(list(map(S.Rational,reversed(a))),t,domain=S.QQ)
def coeff(f):return list(map(str,reversed(f.all_coeffs())))
def modpoly(f,p):return S.Poly.from_list([int(c.p)%p*pow(int(c.q),-1,p)%p for c in f.all_coeffs()],t,modulus=p)
def ev(cs,x,p):
 v=0
 for c in reversed(cs):v=(v*x+c)%p
 return v
def reduce(cs,p):return [c.numerator*pow(c.denominator,-1,p)%p for c in cs]
def infinity(n,d,weight):
 k=n.degree()-d.degree()-weight
 return 'O' if k>0 else str(n.LC()/d.LC()) if k==0 else '0'


def load(tag):
 raw=json.loads((ROOT/SOURCES[tag]).read_text())
 if tag=='q80':
  points=[tuple((pol(a[k]['numerator_coefficients_low_to_high']),pol(a[k]['denominator_coefficients_low_to_high'])) for k in ['X','Y']) for a in raw['sections']['records']]
  a1=a3=pol(['0'])
 elif tag=='r17':
  points=[((pol(a['x']),pol(['1'])),(pol(a['y']),pol(['1']))) for a in raw['basis']]
  a1=a3=pol(['0'])
 else:
  # Only the generic equation/basis fields are used, never the specialized coordinates or parameter.
  points=[tuple((pol(a[k]['numerator']),pol(a[k]['denominator'])) for k in [0,1]) for a in raw['basis_weierstrass_coordinates']]
  av=raw['a_invariants'];assert av[0]['denominator']==av[2]['denominator']==['1']
  a1=pol(av[0]['numerator']);a3=pol(av[2]['numerator'])
 assert len(points)==17 and a1.degree()<=0 and a3.degree()<=0
 return raw,points,a1,a3


def incidence(tag):
 raw,points,a1,a3=load(tag);loci=[];ordinates=[]
 for i,((xn,xd),(yn,yd)) in enumerate(points):
  assert S.gcd(xn,xd).degree()==S.gcd(yn,yd).degree()==0
  if xd.degree()>0:loci.append((f'pole-{i}',1<<i,xd.sqf_part().monic()))
  zn=2*yn*xd+(a1*xn+a3*xd)*yd;zd=yd*xd
  ordinates.append(zn.exquo(S.gcd(zn,zd)).sqf_part().monic())
 for i,j in combinations(range(17),2):
  xn,xd=points[i][0];un,ud=points[j][0]
  n=xn*ud-un*xd;d=xd*ud;assert not n.is_zero
  loci.append((f'pair-{i}-{j}',(1<<i)|(1<<j),n.exquo(S.gcd(n,d)).sqf_part().monic()))
 # All polynomials are monic and p-integral, making the modular coprimality gate sound.
 assert all(f.LC()==1 for _,_,f in loci) and all(f.LC()==1 for f in ordinates)
 for p in range(131,998):
  if not prime(p):continue
  try:reduced=[modpoly(f,p) for _,_,f in loci];oy=[modpoly(f,p) for f in ordinates]
  except ValueError:continue
  break
 else:raise AssertionError('No common integral reduction prime in fixed pool')
 joint=[];mixed=[]
 for i,j in combinations(range(len(loci)),2):
  if S.gcd(reduced[i],reduced[j]).degree()==0:continue
  g=S.gcd(loci[i][2],loci[j][2]);assert g.degree()==0
  joint.append([loci[i][0],loci[j][0]])
 for j,y in enumerate(ordinates):
  for i,(_,_,f) in enumerate(loci):
   if S.gcd(reduced[i],oy[j]).degree()==0:continue
   g=S.gcd(f,y);assert g.degree()==0
   mixed.append([loci[i][0],j])
 xs=[infinity(*a[0],4) for a in points];ys=[infinity(*a[1],6) for a in points]
 # Constant a1,a3 terms disappear in the minimal infinity ordinate chart.
 assert 'O' not in xs and 'O' not in ys and '0' not in ys and len(set(xs))==17
 return {'source':SOURCES[tag],'source_sha256':sha256((ROOT/SOURCES[tag]).read_bytes()).hexdigest(),
         'loci':[{'label':name,'parity_mask':mask,'polynomial':coeff(f)} for name,mask,f in loci],
         'torsion_ordinate_polynomials':[coeff(f) for f in ordinates],
         'prime':p,'joint_modular_survivors_rejected_exactly':joint,'mixed_modular_survivors_rejected_exactly':mixed,
         'infinity_x':xs,'infinity_y':ys,'exact_joint_incidence_count':0,'exact_mixed_incidence_count':0},(raw,points,loci)


def poles(packet):
 raw,points,loci=packet;model=raw['weierstrass_model'];candidates=[]
 for name,mask,g in loci:
  if not name.startswith('pole-'):continue
  i=int(name.split('-')[1]);xn,xd=points[i][0];yn,yd=points[i][1]
  assert g.degree() in [1,2] and g.is_irreducible
  assert xd.monic()==g**2 and yd.monic()==g**3
  assert S.gcd(g,xn).degree()==S.gcd(g,yn).degree()==0
  candidates.append({'label':name,'expected_parity':mask,'q':list(map(F,coeff(g))), 'rows':[],'places':[]})
 assert len(candidates)==12
 polys=[list(map(F,model[k])) for k in ['A_coefficients_low_to_high','B_coefficients_low_to_high']]
 for P in points:
  for n,d in P:polys.extend([list(map(F,coeff(n))),list(map(F,coeff(d)))])
 for p in range(5,998):
  if not prime(p):continue
  try:ff=[reduce(f,p) for f in polys]
  except ValueError:continue
  cache={}
  for c in candidates:
   if rank(c['rows'])==16:continue
   try:q=reduce(c['q'],p)
   except ValueError:continue
   if len(q)==3 and (q[1]**2-4*q[0]*q[2])%p==0:continue
   for z in range(p):
    if ev(q,z,p):continue
    if z not in cache:
     a,b=ev(ff[0],z,p),ev(ff[1],z,p)
     if (4*a**3+27*b*b)%p==0:cache[z]=None;continue
     pts=[];ok=True
     for i in range(17):
      xn,xd,yn,yd=[ev(f,z,p) for f in ff[2+4*i:6+4*i]]
      if xd==0:
       if xn==0:ok=False;break
       pts.append(None)
      elif yd==0:ok=False;break
      else:
       x,y=xn*pow(xd,-1,p)%p,yn*pow(yd,-1,p)%p
       assert (y*y-x**3-a*x-b)%p==0;pts.append((x,y))
     if not ok:cache[z]=None;continue
     labels,dim,order,doubles=quotient(a,b,p)
     rr=[sum(((labels[P]>>j)&1)<<i for i,P in enumerate(pts)) for j in range(dim)]
     cc=[]
     for e in range(p):
      if (e**3+a*e+b)%p:continue
      vals=[1 if P is None else (P[0]-e)%p if P[0]!=e else (3*e*e+a)%p for P in pts]
      assert all(vals)
      cc.append(sum((pow(v,(p-1)//2,p)==p-1)<<i for i,v in enumerate(vals)))
     assert rank(cc)==rank(rr)==rank(cc+rr)
     cache[z]={'p':p,'t':z,'A':a,'B':b,'rows':rr,'character_rows':cc,'points':pts,'group_order':order,'double_subgroup_order':doubles}
    place=cache[z]
    if place is None:continue
    assert all((row&c['expected_parity']).bit_count()%2==0 for row in place['rows'])
    c['places'].append(place);c['rows']+=place['rows']
    if rank(c['rows'])==16:break
  if all(rank(c['rows'])==16 for c in candidates):break
 for c in candidates:
  assert rank(c['rows'])==16
  c['q']=list(map(str,c['q']));c['matrix_rank']=16
  c['kernel']=[0,c['expected_parity']]
  c['no_branch_two_torsion_witness']=next(r for r in c['places'] if r['group_order']%2==1)
 return candidates


if __name__=='__main__':
 resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));start=time.monotonic()
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--record',type=Path,required=True);args=p.parse_args()
 result={'schema':'basis-branch-incidence-v1','parents':{}}
 for tag in SOURCES:
  result['parents'][tag],packet=incidence(tag)
  if tag=='q80':result['q80_pole_kernels']=poles(packet)
 result.update(script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),software={'python':'3','sympy':S.__version__},limits={'cpu_seconds':20,'memory_bytes':1024**3},elapsed_seconds=time.monotonic()-start,
               helper_hashes={n:sha256(Path(__file__).with_name(n).read_bytes()).hexdigest() for n in ['verify_q80_branch_trace_specialization.py','audit_q80_branch_trace_specialization.py']},
               boundary='Only the four retained generic bases, their individual pole relations, pairwise signed equalities and their own nonzero two-torsion values. The twelve Q80 pole fields have exact one-dimensional specialization kernels and no nonzero residue-field two-torsion. No other MW combination, branch divisor or positive construction is excluded.')
 args.record.parent.mkdir(parents=True,exist_ok=True)
 with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps({'parent_loci':{k:len(v['loci']) for k,v in result['parents'].items()},'pole_kernel_ranks':[c['matrix_rank'] for c in result['q80_pole_kernels']],'elapsed_seconds':result['elapsed_seconds']}))
