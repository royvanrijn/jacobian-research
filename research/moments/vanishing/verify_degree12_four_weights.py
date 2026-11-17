#!/usr/bin/env python3
"""Exact degree12 sector ideals and full-module quartic/prefix controls;not a full model solver."""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,comb
from itertools import product
import sys,json,hashlib

ROOT=Path(__file__).resolve().parent
import verify_haar_shear_witnesses as S
import verify_octahedral_minimal_type as M

def deriv12(poly,up):
 out={}
 for a,c in poly.items():
  for i,r in enumerate(a):
   j=i-1 if up else i+1
   if r and 0<=j<=12:
    b=list(a);b[i]-=1;b[j]+=1;b=tuple(b)
    out[b]=out.get(b,Q(0))+c*r*(i if up else 12-i)
 return {a:c for a,c in out.items() if c}

def full_moment(m,ends,middle,end_product):
 seed={}
 for r in range(m//2+1):
  a=[0]*13;a[ends]=a[12-ends]=r;a[6]=m-2*r
  seed[tuple(a)]=Q(factorial(m)*middle**(m-2*r)*end_product**r,factorial(r)**2*factorial(m-2*r))
 out={};low=seed
 for k in range(6*m+1):
  high=low
  for _ in range(k):high=deriv12(high,True)
  out=M.add(out,M.scale(high,Q((-1)**k,factorial(k)*factorial(k+1))))
  low=deriv12(low,False)
 assert not low and not deriv12(out,True) and not deriv12(out,False)
 return out

def moment_sector(m,basis):
 out={}
 powers=[[S.power(v,r) for r in range(m+1)] for v in basis]
 for a in product(range(m+1),repeat=4):
  if sum(a)!=m or sum(w*e for w,e in zip([3,1,-1,-3],a)):continue
  poly={(0,0,0,0):Q(1)}
  for i,e in enumerate(a):poly=S.mul(poly,powers[i][e])
  c=Q(factorial(m),factorial(a[0])*factorial(a[1])*factorial(a[2])*factorial(a[3]))*S.haar(poly)
  if c:out[a]=c
 return out

def sector_columns(m):
 u={(1,0,0,1):Q(1)};v={(0,1,1,0):Q(1)}
 z={(1,0,3,0):Q(1),(0,3,0,1):Q(1)}
 columns=[];names=[]
 for ez in range(m//4+1):
  for eu in range(m//2+1):
   ev=(m-4*ez-2*eu)//2
   if ev<0:continue
   column=S.mul(S.mul(S.power(u,eu),S.power(v,ev)),S.power(z,ez))
   columns.append(column);names.append((eu,ev,ez))
 return names,columns


def dim_invariants(degree):
 counts={(0,0):1}
 for j in range(13):
  out={}
  for (d,w),c in counts.items():
   for r in range(degree-d+1):out[(d+r,w+j*r)]=out.get((d+r,w+j*r),0)+c
  counts=out
 return counts.get((degree,6*degree),0)-counts.get((degree,6*degree-1),0)

def groups(name):
 if name=='2O':
  field=S.Cyclotomic4((1,0,0,0));z,o,zz=field.z,field.one,field.zero
  i=field.power(z,2);half=lambda a:field.scale(a,Q(1,2))
  gens=[(z,zz,zz,field.scale(field.power(z,3),-1)),(zz,o,field.scale(o,-1),zz),
        (half(field.add(o,i)),half(field.add(o,i)),half(field.add(field.scale(o,-1),i)),half(field.add(o,field.scale(i,-1))))]
  expected=48
 else:
  field=S.Cyclotomic4((1,1,1,1));z,o,zz=field.z,field.one,field.zero
  z2,z3,z4=(field.power(z,k) for k in (2,3,4))
  root5=field.add(o,field.scale(field.add(z,z4),2))
  a=field.scale(field.mul(field.add(z4,field.scale(z,-1)),root5),Q(1,5))
  b=field.scale(field.mul(field.add(z2,field.scale(z3,-1)),root5),Q(1,5))
  gens=[(field.scale(z3,-1),zz,zz,field.scale(z2,-1)),(zz,o,field.scale(o,-1),zz),(a,b,b,field.scale(a,-1))]
  expected=120
 values,orders=S.matrix_group(field,gens,expected)
 assert orders==({1:1,2:9,3:8,4:6} if name=='2O' else {1:1,2:15,3:20,5:24})
 return field,gens,values

def trim(p):
 while p and not p[-1]:p.pop()
 return p

def rem(p,f):
 p=trim(p.copy())
 while len(p)>=len(f):
  at=len(p)-len(f);c=p[-1]/f[-1]
  for j,x in enumerate(f):p[at+j]-=c*x
  trim(p)
 return p

def poly_gcd(a,b):
 while b:a,b=b,rem(a,b)
 return [x/a[-1] for x in a]

def algebraic_prefix_control(kappa,lam,delta,expressions):
 f=[kappa**3,0,0,-lam,0,0,Q(1)]
 assert poly_gcd(f,[0,0,-3*lam,0,0,Q(6)])==[1]
 def add(a,b):
  out=[Q(0)]*max(len(a),len(b))
  for i,x in enumerate(a):out[i]+=x
  for i,x in enumerate(b):out[i]+=x
  return rem(out,f)
 def mul(a,b):
  if not a or not b:return []
  out=[Q(0)]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]+=x*y
  return rem(out,f)
 def power(a,n):
  out=[Q(1)]
  for _ in range(n):out=mul(out,a)
  return out
 b=[Q(0),Q(1)];c=[Q(0),Q(0),lam/kappa**2,Q(0),Q(0),-1/kappa**2]
 v=mul(b,c);z=add(power(b,3),power(c,3))
 assert v==[kappa] and z==[lam]
 images={}
 for m,terms in expressions.items():
  value=[]
  for eu,ev,ez,coefficient in terms:
   part=mul(power(v,ev),power(z,ez))
   value=add(value,[coefficient*x for x in part])
  images[m]=value
 assert images[2]==[] and images[4]==[] and images[6]==[delta]
 return {'sextic_coefficients_ascending':[str(x) for x in f],
         'squarefree':True,'b_inverse_and_c_relations':True,
         'first_five_pure_moments_zero':True,'sixth_moment':str(delta),
         'all_six_algebraic_roots_give_valid_witnesses':True}

def main():
 frozen=json.loads((ROOT/'degree12_four_weight_relations.json').read_text())['cases']
 vars13=[{tuple(int(i==j) for j in range(13)):Q(1)} for i in range(13)]
 f={i:M.scale(vars13[i],Q(comb(12,i))) for i in range(13)}
 i2=M.trans(f,12,f,12,12)[0]
 g4=M.trans(f,12,f,12,10);g8=M.trans(f,12,f,12,8)
 A=M.trans(g4,4,g4,4,4)[0];B=M.trans(g8,8,g8,8,8)[0]
 i2square=M.mul(i2,i2);columns=[i2square,A,B]
 assert [dim_invariants(d) for d in range(5)]==[1,0,1,1,3]
 rows={}
 for name,coeffs,seed in [('2O',{2:Q(1),6:Q(-2),10:Q(1)},(2,-2,1)),('2I',{1:Q(1),6:Q(11),11:Q(-1)},(1,11,-1))]:
  expected=frozen[name];field,gens,values=groups(name)
  S.invariant_form(field,gens,coeffs,12)
  assert S.invariant_dimensions(field,values,12)[12]==1
  p={(12-i,i,0,0):c for i,c in coeffs.items()};basis=[p]
  for j in range(1,13):basis.append(S.scale(S.lower(basis[-1]),Q(1,j)))
  assert basis[-1]==S.compact_star(p) and not S.lower(basis[-1])
  sector=[basis[j] for j in [0,4,8,12]]
  for j,b in zip([0,4,8,12],sector):
   assert all(sum(a)==12 and a[0]+a[1]-a[2]-a[3]==12-2*j for a in b)
   assert (12-2*j)%8==4
  expressions={};coefficient_count=0
  for m in [2,4,6]:
   hm=moment_sector(m,sector);names,cc=sector_columns(m)
   solved=M.solve_columns(cc,hm)
   exact={(eu,ev,ez):c for (eu,ev,ez),c in zip(names,solved)}
   target={(t['u'],t['v'],t['z']):Q(t['coefficient']) for t in expected['moments'][str(m)]}
   assert exact==target
   expressions[m]=[(*key,c) for key,c in exact.items()]
   coefficient_count+=len(hm)
  kappa=Q(expected['kappa']);lam=Q(expected['lambda']);delta=Q(expected['delta']);gamma=Q(expected['gamma'])
  assert kappa==Q(-1,495) and gamma!=0 and delta!=0
  for m,terms in expressions.items():
   normal=sum(c*kappa**ev*lam**ez for eu,ev,ez,c in terms)
   assert normal==(delta if m==6 else 0)
  full2=full_moment(2,*seed);full4=full_moment(4,*seed)
  restriction={tuple(a[i] for i in [0,4,8,12]):c for a,c in full4.items() if all(a[i]==0 for i in range(13) if i not in [0,4,8,12])}
  assert restriction==moment_sector(4,sector)
  coordinates=M.solve_columns(columns,full4)
  assert [str(c) for c in coordinates]==expected['full_quartic_coordinates'] and coordinates[2]!=0
  h2scale=M.solve_columns([i2],full2)[0];assert h2scale!=0
  assert len(full4)==expected['full_quartic_terms']
  mixed=S.haar(S.mul(basis[12],S.power(basis[4],3)))
  assert mixed==gamma/4
  for m in range(1,11):assert not any(12*j-10*(m-j)==0 for j in range(m+1))
  eleventh=comb(11,5)*S.haar(S.mul(S.power(basis[0],5),S.power(basis[11],6)))
  assert str(eleventh)==expected['eleventh_moment_p0_plus_p11'] and eleventh!=0
  algebraic=algebraic_prefix_control(kappa,lam,delta,expressions)
  rows[name]={'sector_parameter_moment_coefficients_checked':coefficient_count,
    'kappa':str(kappa),'lambda':str(lam),'delta':str(delta),'gamma':str(gamma),
    'full_quartic_moment_terms':len(full4),'full_quartic_coordinates':[str(c) for c in coordinates],
    'full_quartic_invariant_dimension':3,'moment_ideal_invariant_quartic_dimension':2,
    'explicit_missing_quartic':'A=((f,f)_10,(f,f)_10)_4,in the primitive convention',
    'sharp_sector_set_and_ordinary_prefix':6,'sector_pure_locus':'a=b=0 or c=d=0',
    'sector_ordinary_ideal':'(v+u/495,z-lambda*u^2,u^3)',
    'sharp_uniform_marker_coefficient':'1/4','marker_p12_times_p4_cubed':str(mixed),
    'algebraic_sharp_prefix_control':algebraic,
    'full_V12_prefix_lower_bound':11,'first_ten_moments_p0_plus_p11_zero_by_weight':True,
    'eleventh_moment_p0_plus_p11':str(eleventh)}
  print(name,'sector prefix6,delta',delta,'full prefix lower11',eleventh,flush=True)
 receipt={'completed':True,'cases':rows,'symmetric_power_invariant_dimensions_degrees_0_to_4':[1,0,1,1,3],
  'all_order_sector_proof':'Weyl-fixed torus invariants C[u,v,z] and graded normal form in the written proof',
  'full_ordinary_moment_ideal_equals_full_invariant_nullcone_ideal':False,
  'full_V12_pure_locus_classified':False,'any_global_exceptional_MZ_decision':False,
  'independent_replay':False,'external_review':False,'formal_verification':False,
  'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'relations_sha256':hashlib.sha256((ROOT/'degree12_four_weight_relations.json').read_bytes()).hexdigest(),
  'imported_arithmetic_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['verify_haar_shear_witnesses.py','verify_octahedral_minimal_type.py']}}
 (ROOT/'degree12_four_weight_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
