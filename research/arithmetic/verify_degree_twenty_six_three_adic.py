#!/usr/bin/env python3
"""Exact tangent and first-Newton controls for degree26 odd Q3 pencils.

The accompanying note proves the infinite p-adic construction.
This checker verifies the rank-eight minor and an actual unrestricted
paired-critical lift modulo nine.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,platform
import sympy as sp

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--write-receipt",type=Path)
 args=parser.parse_args()
 assert sp.__version__=="1.14.0"
 dependency=Path(__file__).with_name("verify_degree_twenty_six_reciprocal_local.py")
 spec=importlib.util.spec_from_file_location("local_certificate",dependency)
 m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 N=[1,1,1,1,0,0,1,1,0,0,1,0,1,1]
 T=[1,0,0,1,1,1,2,1,1,1,0,0,1]
 assert T==T[::-1]
 assert sum(N)%3==0 and sum(j*c for j,c in enumerate(N))%3==0
 expanded=[0]*14
 for j,c in enumerate(T):
  expanded[j]-=c;expanded[j+1]+=c
 assert [(2*j-1)*c%3 for j,c in enumerate(N)]==[x%3 for x in expanded]
 assert m.gcd(N,T)==[1]
 invu=[(-T[i])%3 for i in range(1,13)]
 assert sum(T)%3!=0
 derivative=[i*T[i]%3 for i in range(1,13)]
 assert m.gcd(T,m.trim(derivative))==[1]
 assert m.substitute(invu,invu,T)==[0,1]
 tau_columns=[];power=[1]
 for j in range(12):
  tau_columns.append(power+[0]*(12-len(power)))
  power=m.mul(power,invu,T)
 anti_matrix=[[(int(i==j)-tau_columns[j][i])%3 for j in range(12)]for i in range(12)]
 assert m.rank(anti_matrix)==6 and m.rank(anti_matrix[:6])==6
 W=m.mul([2],m.mul(m.rem(N,T),invu,T),T)
 columns=[]
 for j in range(13):
  value=m.mul(W,[0]*j+[1],T)
  anti=m.add(value,m.substitute(value,invu,T),-1)
  columns.append([1,j%3]+anti+[0]*(12-len(anti)))
 J=[[columns[j][i]for j in range(13)]for i in range(8)]
 variables=[0,1,2,3,4,5,6,8]
 minor=sp.Matrix([[J[i][j]for j in variables]for i in range(8)])
 assert int(minor.det())==16 and m.rank(J)==8

 # Arithmetic in A/9A, with literal reciprocal T and tau(t)=t^-1.
 mod=9
 def trim(a):
  a=[int(x)%mod for x in a]
  while len(a)>1 and a[-1]==0:a.pop()
  return a
 def add(a,b,sign=1):
  c=[0]*max(len(a),len(b))
  for i,x in enumerate(a):c[i]+=x
  for i,x in enumerate(b):c[i]+=sign*x
  return trim(c)
 def rem(a):
  a=trim(a[:])
  while len(a)>12:
   j=len(a)-13;q=a[-1]
   for k,x in enumerate(T):a[j+k]=(a[j+k]-q*x)%mod
   a=trim(a)
  return a
 def mul(a,b):
  c=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):c[i+j]+=x*y
  return rem(c)
 def evaluate(a,value):
  ans=[0]
  for x in reversed(a):ans=add(mul(ans,value),[x])
  return ans
 u=sp.symbols("u")
 def inverse(value):
  v=sp.Poly(sum(int(x)*u**j for j,x in enumerate(value)),u,modulus=3)
  tp=sp.Poly(sum(int(x)*u**j for j,x in enumerate(T)),u,modulus=3)
  inv=sp.invert(v,tp)
  first=[int(inv.nth(j))%3 for j in range(12)]
  result=mul(first,add([2],mul(value,first),-1))
  assert mul(value,result)==[1]
  return result
 invt=[(-T[i])%9 for i in range(1,13)]
 assert mul([0,1],invt)==[1]
 def target(N):
  R=[(2*j-1)*c%9 for j,c in enumerate(N)]
  derivative=[j*R[j]%9 for j in range(1,len(R))]
  U=add([0,1],mul(evaluate(R,[0,1]),inverse(evaluate(derivative,[0,1]))),-1)
  assert evaluate(R,U)==[0]
  value=mul(mul(evaluate(N,U),evaluate(N,U)),inverse(U))
  anti=add(value,evaluate(value,invt),-1)
  anti+= [0]*(12-len(anti))
  G=[sum(N)%9,sum(j*c for j,c in enumerate(N))%9]+anti[:6]
  return G,U,anti
 G,U,anti=target(N)
 assert all(x%3==0 for x in G)
 rhs=[(-x//3)%3 for x in G]
 equations=[[int(minor[i,j])%3 for j in range(8)]+[rhs[i]]for i in range(8)]
 rows,pivots,free=m.affine(equations,8)
 assert not free
 correction=[0]*8
 for row,j in zip(rows,pivots):correction[j]=row[-1]
 lifted=N[:]
 for j,c in zip(variables,correction):lifted[j]=(lifted[j]+3*c)%9
 G,U,anti=target(lifted)
 assert not any(G)and not any(anti)
 S=[0]*13;partial=0
 for j in reversed(range(1,14)):
  partial=(partial+(2*j-1)*lifted[j])%9;S[j-1]=partial
 assert all(any((S[j]-pow(lam,6-j,9)*S[12-j])%9 for j in range(6))
            for lam in(1,2,4,5,7,8))
 report={
  "schema_version":1,"kind":"exact_three_adic_hensel_controls",
  "python":platform.python_version(),"sympy":sp.__version__,
  "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  "local_checker_sha256":hashlib.sha256(dependency.read_bytes()).hexdigest(),
  "anti_involution_rank":6,"first_six_anti_coordinates_have_rank":6,
  "constraint_rank":8,"solved_coefficient_indices":variables,
  "free_coefficient_indices":[7,9,10,11,12],
  "minor":[[int(x) for x in row] for row in minor.tolist()],"minor_determinant":16,"minor_determinant_mod3":1,
  "mod9_numerator_coefficients":lifted,"mod9_critical_root_in_A":U,
  "mod9_paired_critical_values_invariant":True,
  "mod9_residual_critical_polynomial_is_not_reciprocal_for_any_unit_lambda":True,
  "scope":"Finite rank/minor and modulo9 controls. The written Newton proof gives the five-parameter Z3 family, full regular A26 over Q3(s), and its exclusion from the good reciprocal stratum. No rational coefficient descent or full-A26 arithmetic specialization over Q3 is asserted.",
  "independent_replay":False,"formal_verification":False,
  "external_review":False,"result":"PASS"
 }
 output=json.dumps(report,indent=2,sort_keys=True)+"\n"
 if args.write_receipt:args.write_receipt.write_text(output)
 print(output,end="")
if __name__=="__main__":main()

