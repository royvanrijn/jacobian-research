"""Exact bounded controls only; the all-index argument is in the TeX."""
from pathlib import Path
from math import factorial
import hashlib,json
HERE=Path(__file__).resolve().parent
tests=[]
def check(name,ok):
 if not ok:raise AssertionError(name)
 tests.append(name)
def vp(x):
 if not x:return 10**9
 x=abs(x);v=0
 while x%7==0:x//=7;v+=1
 return v
def s7(n):
 s=0
 while n:s,n=s+n%7,n//7
 return s
def raise_poly(a,c):
 out=[0]*(len(a)+1)
 for k,x in enumerate(a):out[k]+=(c-k)*x;out[k+1]+=x
 return out
for n in [3,8,20,48,145,342]:
 e=[(s7(n)-a+4)//6 for a in range(6)]
 f=[1]
 for j in range(n):f=raise_poly(f,n+6*j+5)
 for r in range(50):
  denval=vp(factorial(n))+vp(factorial(r))
  for a in range(6):
   for j in range(a,len(f),6):
    val=vp(f[j])+(j-a)//6-denval
    check("actual_binomial_Q_content",val>=e[a])
  f=raise_poly(f,n+5-r)
# Factorial valuation proof of integrality of every formal gauge coefficient.
for j in range(1,2001):
 ell,r=divmod(j,6)
 for epsilon in [0,1]:
  value=ell+epsilon-vp(factorial(j))
  check("formal_gauge_coefficient_integral",value>=0)
  check("formal_gauge_digit_identity",value==epsilon+(s7(j)-r)//6)
# Explicit surviving-edge coefficients prevent restricted-series convergence.
for depth in range(1,6):
 j=7**depth
 check("unit_gauge_tail_coefficient",(j-1)//6-vp(factorial(j))==0)
# All finite residue cut patterns: exactly one 7-weighted link, cycle product7v.
for digit_sum in range(1,201):
 e=[(digit_sum-a+4)//6 for a in range(6)]
 weights=[1+e[5]-e[0]]+[e[a-1]-e[a]for a in range(1,6)]
 check("one_exact_cycle_valuation",sum(weights)==1 and set(weights)<={0,1})
# Corruption control: removing the required seven from a unit tail would fail.
check("reject_wrong_factorial_normalization",(7-1)//6-vp(factorial(7))-1<0)
pins={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
 for p in [Path(__file__),HERE/"GAP_SIX_INTEGER_VALUED_RAISING_CALCULUS.tex",
 HERE/"GAP_SIX_DIGIT_NORMALIZATION_AND_FINITE_GAUGE.tex"]}
receipt={"status":"PASS","controls":len(tests),
"scope":"Only bounded coefficient/valuation controls, not a finite proof of the all-index theorem or any noncollision claim.",
"source_sha256":pins,"independent_whole_proof_replay":False,
"formal_verification":False,"external_review":False}
(HERE/"gap_six_integer_valued_raising_calculus_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"status":"PASS","controls":len(tests)}))
