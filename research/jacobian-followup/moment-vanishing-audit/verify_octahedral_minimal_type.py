#!/usr/bin/env python3
"""Exact universal octavic moment certificates; classical ring/HSOP inputs are written."""
from pathlib import Path
import hashlib
import sys
from collections import defaultdict
from fractions import Fraction as Q
from math import factorial, comb, gcd, lcm
import json

N = 8
ZERO = (0,)*9
ONE = {ZERO: Q(1)}
VARS = [{tuple(int(i == j) for j in range(9)): Q(1)} for i in range(9)]

def add(*ps):
 out = defaultdict(Q)
 for p in ps:
  for m,c in p.items():out[m] += c
 return {m:c for m,c in out.items() if c}

def scale(p,c):return {m:v*c for m,v in p.items() if v*c}

def mul(p,q):
 out=defaultdict(Q)
 for a,c in p.items():
  for b,d in q.items():out[tuple(x+y for x,y in zip(a,b))] += c*d
 return {m:c for m,c in out.items() if c}

def power(p,n):
 out=ONE
 for _ in range(n):out=mul(out,p)
 return out

def deriv(p,up):
 out=defaultdict(Q)
 for a,c in p.items():
  for i,r in enumerate(a):
   j=i-1 if up else i+1
   if r and 0<=j<=N:
    b=list(a);b[i]-=1;b[j]+=1
    out[tuple(b)]+=c*r*(i if up else N-i)
 return {m:c for m,c in out.items() if c}

def moment(m):
 seed={}
 for j in range(m//2+1):
  a=[0]*9;a[0]=a[8]=j;a[4]=m-2*j
  seed[tuple(a)]=Q(factorial(m)*14**(m-2*j),factorial(j)**2*factorial(m-2*j))
 out={};lower=seed
 for k in range(4*m+1):
  top=lower
  for _ in range(k):top=deriv(top,True)
  out=add(out,scale(top,Q((-1)**k,factorial(k)*factorial(k+1))))
  lower=deriv(lower,False)
 assert not lower and not deriv(out,True) and not deriv(out,False)
 return out

def primitive_covariant(f):
 vals=[c for a in f.values() for c in a.values()]
 den=lcm(*(c.denominator for c in vals));num=gcd(*(int(c*den) for c in vals))
 return {i:scale(a,Q(den,num)) for i,a in f.items() if a}

def fall(a,k):return factorial(a)//factorial(a-k) if a>=k else 0

def trans(a,da,b,db,r):
 out={}
 for i,pa in a.items():
  for j,pb in b.items():
   index=i+j-r
   coefficient=sum((-1)**k*comb(r,k)*fall(da-i,r-k)*fall(i,k)*fall(db-j,k)*fall(j,r-k) for k in range(r+1))
   if coefficient:
    out[index]=add(out.get(index,{}),scale(mul(pa,pb),coefficient))
 return primitive_covariant(out)

def invariants():
 f={i:scale(VARS[i],Q(comb(8,i))) for i in range(9)}
 g=trans(f,8,f,8,4);k=trans(f,8,f,8,6)
 h=trans(k,4,k,4,2);m=trans(f,8,k,4,4);nn=trans(f,8,h,4,4);pp=trans(g,8,k,4,4);qq=trans(g,8,h,4,4)
 pairs={2:(f,8,f,8,8),3:(f,8,g,8,8),4:(k,4,k,4,4),5:(m,4,k,4,4),6:(k,4,h,4,4),7:(m,4,h,4,4),8:(pp,4,h,4,4),9:(nn,4,h,4,4),10:(qq,4,h,4,4)}
 return {d:trans(*args)[0] for d,args in pairs.items()}

def partitions(d,minimum=2):
 if not d:yield []
 for x in range(minimum,min(d,10)+1):
  for rest in partitions(d-x,x):yield [x]+rest

def solve_columns(columns,target):
 rows=[];pivots={}
 keys=sorted(set(target).union(*(set(c) for c in columns)))
 for key in keys:
  row=[c.get(key,Q(0)) for c in columns]+[target.get(key,Q(0))]
  for col,pivot in sorted(pivots.items()):
   if row[col]:row=[a-row[col]*b for a,b in zip(row,pivot)]
  first=next((i for i,x in enumerate(row[:-1]) if x),None)
  if first is None:assert not row[-1]
  else:
   row=[x/row[first] for x in row];pivots[first]=row
 assert len(pivots)==len(columns)
 result=[Q(0)]*len(columns)
 for col,row in sorted(pivots.items(),reverse=True):result[col]=row[-1]-sum(row[i]*result[i] for i in range(col+1,len(columns)))
 assert add(*(scale(p,c) for p,c in zip(columns,result)))==target
 return result


ROOT=Path(__file__).resolve().parent

def evaluate(poly,values):
 return sum((c*product_value(a,values) for a,c in poly.items()),Q(0))

def product_value(exponents,values):
 result=Q(1)
 for e,v in zip(exponents,values):result*=Q(v)**e
 return result

def coefficient_hash(poly):
 data=[[list(a),str(c)] for a,c in sorted(poly.items())]
 return hashlib.sha256(json.dumps(data,separators=(',',':')).encode()).hexdigest()

def invariant_dimensions(max_degree):
 counts={(0,0):1}
 for j in range(9):
  out=defaultdict(int)
  for (d,w),c in counts.items():
   for r in range(max_degree-d+1):out[(d+r,w+j*r)]+=c
  counts=out
 return [counts.get((d,4*d),0)-counts.get((d,4*d-1),0) for d in range(max_degree+1)]

def direct_controls(hmoments):
 import verify_haar_shear_witnesses as S
 field=S.Cyclotomic4((1,0,0,0));z,o,zz=field.z,field.one,field.zero
 i=field.power(z,2);half=lambda a:field.scale(a,Q(1,2))
 generators=[
  (z,zz,zz,field.scale(field.power(z,3),-1)),
  (zz,o,field.scale(o,-1),zz),
  (half(field.add(o,i)),half(field.add(o,i)),
   half(field.add(field.scale(o,-1),i)),half(field.add(o,field.scale(i,-1))))]
 group,orders=S.matrix_group(field,generators,48)
 assert orders=={1:1,2:9,3:8,4:6}
 dimensions=S.invariant_dimensions(field,group,12)
 assert dimensions==[1]+[0]*7+[1]+[0]*3+[1]
 S.invariant_form(field,generators,{0:Q(1),4:Q(14),8:Q(1)},8)
 p={(8,0,0,0):Q(1),(4,4,0,0):Q(14),(0,8,0,0):Q(1)}
 basis=[p]
 for j in range(1,9):basis.append(S.scale(S.lower(basis[-1]),Q(1,j)))
 assert not S.lower(basis[-1])
 for j,b in enumerate(basis):
  assert not S.omega(b)
  assert S.raise_(b)==(S.scale(basis[j-1],9-j) if j else {})
  assert all(sum(mon)==8 and mon[0]+mon[1]-mon[2]-mon[3]==8-2*j for mon in b)
 identity_values=[sum(c for (a,b,c0,d),c in poly.items() if b==0 and c0==0) for poly in basis]
 assert identity_values==[1,0,0,0,14,0,0,0,1]
 point=[0]*9;point[4]=1
 pure=S.power(basis[4],0);control_values=[]
 for m in range(1,11):
  pure=S.mul(pure,basis[4]);exact=S.haar(pure)
  assert exact==evaluate(hmoments[m],point),(m,exact)
  control_values.append(str(exact))
 point=[0]*9;point[0]=point[7]=1
 prefix_values=[evaluate(hmoments[m],point) for m in range(1,8)]
 assert prefix_values[:6]==[0]*6
 direct7=35*S.haar(S.mul(S.power(basis[0],3),S.power(basis[7],4)))
 assert direct7==prefix_values[6]==Q(-8811708416,154040315)
 q=S.compact_star(p);assert q==basis[8]
 mixed=S.haar(S.mul(q,S.power(basis[3],4)))
 assert mixed==Q(-67436544,230945)
 point=[0]*9;point[3]=1
 assert all(evaluate(hmoments[m],point)==0 for m in range(1,11))
 return {
  'binary_octahedral_order':48,'fixed_dimensions_degrees_0_to_12':dimensions,
  'basis_ladder_harmonicity_weight_and_identity_evaluations':True,
  'direct_matrix_haar_checks_of_middle_basis_moments':control_values,
  'sharp_set_prefix_witness':'P_0+P_7,where P_j=F^j p/j!',
  'sharp_set_prefix_first_seven_moments':[str(v) for v in prefix_values],
  'sharp_set_prefix_seventh_direct_integral':str(direct7),
  'sharp_uniform_marker':'Q=P_8, P=P_3, m=4, degree(Q)=8',
  'sharp_uniform_marker_moment':str(mixed),
  'imported_matrix_arithmetic_sha256':hashlib.sha256((ROOT/'verify_haar_shear_witnesses.py').read_bytes()).hexdigest()}

def main():
 frozen=json.loads((ROOT/'octahedral_moment_relations.json').read_text())
 inv=invariants();hmoments={1:moment(1)};assert hmoments[1]=={}
 dims=invariant_dimensions(10)
 assert dims==[1,0,1,1,2,2,4,4,7,8,12]
 rows={}
 for d in range(2,11):
  h=moment(d);hmoments[d]=h
  parts=list(partitions(d));columns=[]
  for part in parts:
   column=ONE
   for degree in part:column=mul(column,inv[degree])
   columns.append(column)
  solved=solve_columns(columns,h)
  expected={tuple(t['invariant_degrees']):Q(t['coefficient']) for t in frozen['relations'][str(d)]}
  assert len(parts)==dims[d] and {tuple(part):c for part,c in zip(parts,solved)}==expected
  assert expected[(d,)]!=0
  reconstructed=add(*(scale(column,expected[tuple(part)]) for part,column in zip(parts,columns)))
  assert reconstructed==h and not deriv(inv[d],True) and not deriv(inv[d],False)
  rows[d]={'moment_monomials':len(h),'moment_coefficients_sha256':coefficient_hash(h),
           'primitive_invariant_monomials':len(inv[d]),
           'primitive_invariant_sha256':coefficient_hash(inv[d]),
           'full_invariant_space_dimension':dims[d],
           'all_product_columns_independent':True,
           'coefficient_of_new_generator':str(expected[(d,)]),
           'universal_polynomial_identity':True}
  print('degree',d,'terms',len(h),'basis',dims[d],'new coefficient',expected[(d,)],flush=True)
 # Omitting k+1 from the projector denominator leaves the V8 zero vector alive.
 seed=scale(VARS[4],Q(14));bad={};lower=seed
 for k in range(5):
  raised=lower
  for _ in range(k):raised=deriv(raised,True)
  bad=add(bad,scale(raised,Q((-1)**k,factorial(k)**2)))
  lower=deriv(lower,False)
 assert bad==seed and deriv(bad,True) and moment(1)=={}
 receipt={
  'completed':True,'universal_moment_relations':rows,
  'symmetric_power_invariant_dimensions_0_to_10':dims,
  'all_relations_verified_coefficientwise':True,
  'projector_denominator_mutation_detected':True,
  'direct_controls':direct_controls(hmoments),
  'least_pure_only_set_prefix':7,
  'least_full_ordinary_ideal_prefix':10,
  'sharp_uniform_mixed_zero_offset_coefficient':'1/2',
  'octahedral_MZ_counterexample_matrix_degree_lower_bound':12,
  'classical_HSOP_and_ring_generation_imported':True,
  'classical_theorems_formally_verified_here':False,
  'global_octahedral_MZ_decided':False,'icosahedral_lowest_module_decided':False,
  'independent_replay':False,'external_review':False,'formal_verification':False,
  'python_version':sys.version.split()[0],
  'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'frozen_relations_sha256':hashlib.sha256((ROOT/'octahedral_moment_relations.json').read_bytes()).hexdigest()}
 (ROOT/'octahedral_minimal_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps({k:v for k,v in receipt.items() if k!='universal_moment_relations'},indent=2))

if __name__=='__main__':main()
