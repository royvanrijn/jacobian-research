#!/usr/bin/env python3
"""Exact finite witnesses for S24 j-monodromy of four retained MW17 parents."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import resource
import time
import sympy as S
from audit_basis_branch_incidence import load, pol, coeff, reduce, SOURCES
from audit_q80_branch_trace_specialization import prime

ROOT=Path(__file__).resolve().parents[2]
t,z=S.symbols('t z')

def invariants(tag):
 d,_,_,_=load(tag)
 if tag=='q80':
  m=d['weierstrass_model'];return -48*pol(m['A_coefficients_low_to_high']),-864*pol(m['B_coefficients_low_to_high'])
 if tag=='r17':return -48*pol(d['A']),-864*pol(d['B'])
 av=d['a_invariants'];assert all(a['denominator']==['1'] for a in av)
 a1,a2,a3,a4,a6=[pol(a['numerator']) for a in av]
 b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6
 return b2*b2-24*b4,-b2**3+36*b2*b4-216*b6


def mult(A,B,p):
 n=len(A);return [[sum(A[i][k]*B[k][j] for k in range(n))%p for j in range(n)] for i in range(n)]


def witness(tag):
 c4,c6=invariants(tag);skipped=[]
 for p in range(131,998):
  if not prime(p):continue
  try:c=reduce(list(map(F,coeff(c4))),p);e=reduce(list(map(F,coeff(c6))),p)
  except ValueError:skipped.append([p,'nonintegral invariant coefficients']);continue
  u=S.Poly.from_list(c[::-1],t,modulus=p);v=S.Poly.from_list(e[::-1],t,modulus=p)
  D=u**3-v**2;H=3*u.diff()*v-2*u*v.diff()
  if [u.degree(),v.degree(),D.degree(),H.degree()]!=[8,12,24,18]:skipped.append([p,'degree drop']);continue
  if any(S.gcd(f,f.diff()).degree()>0 for f in [u,v,D,H]) or S.gcd(H,u*v*D).degree()>0 or S.gcd(u,v).degree()>0:
   skipped.append([p,'repeated or overlapping divisors']);continue
  N=1728*u**3
  assert (N.diff()*D-N*D.diff()+1728*u*u*v*H).is_zero
  H=H.monic();den_inverse=S.invert(D,H);J=(N*den_inverse).rem(H)
  assert (D*den_inverse).rem(H)==S.Poly(1,t,modulus=p)
  cols=[];x=S.Poly(t,t,modulus=p);w=J
  for i in range(18):cols.append([int(w.nth(k))%p for k in range(18)]);w=(w*x).rem(H)
  M=[[cols[j][i] for j in range(18)] for i in range(18)]
  R=S.Poly(S.Matrix(M).charpoly(z).as_expr(),z,modulus=p)
  if S.gcd(R,R.diff()).degree()>0:skipped.append([p,'colliding critical values']);continue
  assert R.degree()==18 and R.eval(0)%p and R.eval(1728)%p
  # Independent characteristic polynomial from modular matrix traces and Newton identities.
  power=[[int(i==j) for j in range(18)] for i in range(18)];traces=[];newton=[1]
  for k in range(1,19):
   power=mult(power,M,p);traces.append(sum(power[i][i] for i in range(18))%p)
   newton.append(-pow(k,-1,p)*sum(newton[k-i]*traces[i-1] for i in range(1,k+1))%p)
  assert newton==[int(c)%p for c in R.all_coeffs()]
  remainder=S.Poly(0,t,modulus=p)
  for c in newton:remainder=(remainder*J+int(c)).rem(H)
  assert remainder.is_zero
  return {'source':SOURCES[tag],'source_sha256':sha256((ROOT/SOURCES[tag]).read_bytes()).hexdigest(),
          'prime':p,'skipped_primes':skipped,'c4':c,'c6':e,'D':coeff(D),'critical_polynomial':coeff(H),
          'inverse_denominator_mod_critical':coeff(den_inverse),'j_mod_critical':coeff(J),
          'multiplication_matrix':M,'critical_value_polynomial':coeff(R),'matrix_power_traces':traces,
          'newton_coefficients_high_to_low':newton,'degrees':[8,12,24,18,18],
          'ramification_indices':{'over0':[3]*8,'over1728':[2]*12,'other_distinct_values':[2]*18},
          'total_ramification':8*2+12+18,'infinity_unramified':True}
 raise AssertionError('No witness in fixed prime pool')


if __name__=='__main__':
 resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3));start=time.monotonic()
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--record',type=Path,required=True);args=parser.parse_args()
 data={'schema':'mw17-j-monodromy-v1','parents':{tag:witness(tag) for tag in SOURCES}}
 assert all(p['total_ramification']==46 for p in data['parents'].values())
 data.update(script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
             helper_hashes={name:sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ['audit_basis_branch_incidence.py','audit_q80_branch_trace_specialization.py','verify_q80_branch_trace_specialization.py']},
             elapsed_seconds=time.monotonic()-start,limits={'cpu_seconds':20,'memory_bytes':1024**3},software={'sympy':S.__version__},
             boundary='Finite polynomial, matrix and ramification premises. S24 monodromy, the degree23 different-parameter obstruction and the symmetry/isogeny consequences use the written proof. No two-gain construction or general twist rank exclusion.')
 args.record.parent.mkdir(parents=True,exist_ok=True)
 with args.record.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps({'primes':{k:v['prime'] for k,v in data['parents'].items()},'elapsed_seconds':data['elapsed_seconds']}))
