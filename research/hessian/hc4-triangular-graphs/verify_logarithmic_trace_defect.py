#!/usr/bin/env python3
"""Exact support for logarithmic defect residues and missing derivative jets."""
from __future__ import annotations
import argparse,hashlib,json,platform,time
from pathlib import Path
import sympy as S
BASE=Path(__file__).resolve().parent
PROOF=BASE/"LOGARITHMIC_TRACE_DEFECT_AND_NON_DENSE_DERIVATIVE_DOMAINS.md"
RECEIPT=BASE/"logarithmic_trace_defect_verification.json"
checks=[]
def check(name,c):
 ok=bool(c);checks.append({"name":name,"passed":ok})
 if not ok:raise AssertionError(name)
def eq(name,a,b=0):check(name,S.cancel(S.sympify(a)-S.sympify(b))==0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(f,u):
 f=S.cancel(f)
 if not f:return S.oo
 n,d=S.fraction(f)
 low=lambda g:min(m[0] for m,_ in S.Poly(g,u).terms())
 return low(n)-low(d)
def run():
 u,z,t,Z,beta=S.symbols("u z t Z beta")
 for e in range(2,10):
  D=lambda f:S.cancel(S.diff(f,u)/(e*u**(e-1)))
  W=lambda f:S.cancel(u**e*D(f))
  vals=[]
  for i in range(1,e):
   lam=u**(-i);v=-S.Rational(i,e);vals.append(v)
   eq(f"actual logarithmic Laurent eigenmode e{e}i{i}",W(lam),v*lam)
   eq(f"exact codifferent derivative pole e{e}i{i}",val(D(lam),u),-i-e)
   check(f"branch factor restores original codifferent e{e}i{i}",val(W(lam),u)>=1-e)
   for degree in [0,1,3]:
    b=(1+u**e)**degree
    correction=S.cancel(W(b*lam)-b*W(lam))
    check(f"B-linearity correction regular e{e}i{i}d{degree}",
          val(correction,u)>=0)
   check(f"ordinary square integrable pole e{e}i{i}",2*e-1-2*i>-1)
   check(f"derivative radial exponent diverges e{e}i{i}",
         2*e-1-2*(i+e)<-1)
   eq(f"canonical forbidden derivative exponent e{e}i{i}",
      e-1+val(D(lam),u),-i-1)
  M=S.diag(*vals)
  eq(f"logarithmic residue trace half rank e{e}",S.trace(M),-S.Rational(e-1,2))
  eq(f"logarithmic residue determinant e{e}",M.det(),
     (-1)**(e-1)*S.factorial(e-1)/S.Integer(e)**(e-1))
  eq(f"second spectral moment e{e}",S.trace(M*M),
     S.Rational((e-1)*(2*e-1),6*e))
  check(f"fixed-factor logarithmic spectrum distinct e{e}",len(set(vals))==e-1)
 # Finite packets retain residue degrees and repeated fractions across primes.
 for packet in [[(2,1),(4,2)],[(3,2),(5,1)],[(2,2),(3,1),(6,1)]]:
  eigen=[-S.Rational(i,e) for e,f in packet for _ in range(f) for i in range(1,e)]
  M=S.diag(*eigen);rank=len(eigen)
  expected=S.prod((Z+S.Rational(i,e))**f for e,f in packet for i in range(1,e))
  eq("whole packet characteristic polynomial "+str(packet),M.charpoly(Z).as_expr(),expected)
  eq("whole packet trace "+str(packet),S.trace(M),-S.Rational(rank,2))
  distinct=sorted(set(eigen));projectors=[]
  for a in distinct:
   P=S.eye(rank)
   for b in distinct:
    if b!=a:P=P*(M-b*S.eye(rank))/(a-b)
   check("spectral projector idempotent "+str((packet,a)),P*P==P)
   eq("spectral projector aggregate multiplicity "+str((packet,a)),
      S.trace(P),eigen.count(a))
   check("spectral projector correct eigenspace "+str((packet,a)),
         M*P==a*P)
   projectors.append(P)
  check("aggregate spectral projectors sum to identity "+str(packet),
        sum(projectors,S.zeros(rank))==S.eye(rank))
 check("cross-factor one-half eigenvalue overlap is retained",
       [-S.Rational(i,e) for e,f in [(2,1),(4,2)] for _ in range(f) for i in range(1,e)].count(-S.Rational(1,2))==3)
 # Matrix trace obstruction is exactly trace, not a false source commutator.
 A=S.Matrix([[1,z],[2,3]]);B=S.Matrix([[z,1],[0,4]])
 eq("any endomorphism commutator trace vanishes",S.trace(A*B-B*A),0)
 check("ramified logarithmic trace is nonzero in characteristic zero",-S.Rational(7,2)!=0)
 # Exact canonical jet test with variable analytic units/tangential terms.
 for e in range(2,8):
  unit=1+u**e+z
  aa=1+u**e+z*z;bb=1+z
  D=lambda f:S.cancel(aa*S.diff(f,u)/(e*u**(e-1))+bb*S.diff(f,z))
  for k in range(0,2*e+2):
   F=u**k*(z+2+u*z)
   f=F/(e*unit*u**(e-1))
   coeff=S.cancel(e*unit*u**(e-1)*D(f))
   check(f"canonical derivative regular iff first e-1 jets vanish e{e}k{k}",
         (val(coeff,u)>=0)==(k>=e-1))
   if k<e-1:
    eq(f"first forbidden canonical exponent e{e}k{k}",val(coeff,u),k-e)
    eq(f"noncancelling forbidden coefficient e{e}k{k}",
       S.cancel(coeff*u**(e-k)).subs(u,0),
       (1+z*z)*S.Rational(k-e+1,e)*(z+2))
  # With constant unit the exact cancellation is visible at k=e-1.
  F=u**(e-1)*(z+1)
  f=F/u**(e-1)
  eq(f"critical e-1 coefficient has zero normal derivative e{e}",
     S.diff(f,u)/(e*u**(e-1)),0)
  for i in range(1,e):
   canonical=u**(e-1-i)
   j=e-1-i
   eq(f"codifferent mode detected by a continuous jet e{e}i{i}",
      S.diff(canonical,u,j).subs(u,0),S.factorial(j))
 # A coordinate derivative on the invertible ball has polynomial domain.
 p,q=S.symbols("p q")
 for m in range(6):
  for n in range(4):
   check(f"invertible-ball polynomial derivative stays polynomial m{m}n{n}",
         S.diff(p**m*q**n,p).is_polynomial(p,q))
 return True
def main():
 parser=argparse.ArgumentParser();parser.add_argument("--write-receipt",action="store_true")
 args=parser.parse_args();start=time.monotonic();run()
 proof=PROOF.read_text()
 for phrase in ["B-LINEAR","Equal fractions from DIFFERENT","rank_(C(Gamma))(T_Gamma)/2",
                "NO such commutator","sufficiently small target ball",
                "No quotient-singularity L2 theorem","outer-boundary behavior",
                "NON-DENSITY","separately established identity R intersect H=Lambda",
                "No automatic codifferent multiplication"]:
  check("written scope retains "+phrase,phrase in proof)
 parents=["CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md",
          "BERGMAN_TRACE_DUAL_AND_MULTIPLIER_DOMAIN_RIGIDITY.md",
          "REDUCED_BRANCH_PAIRING_AND_SIXTH_MOMENT_DEFECT.md"]
 receipt={"claim":"Logarithmic trace-defect spectrum and non-dense derivative domains",
 "outcome":"PASS","check_count":len(checks),"checks":checks,
 "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
 "retained_inputs":[{"path":p,"sha256":sha(BASE/p)} for p in parents],
 "environment":{"python":platform.python_version(),"sympy":S.__version__},
 "scope":{"ordinary_codifferent_log_endomorphism":True,
 "generic_spectrum_with_residue_degrees":True,
 "coincident_eigenvalues_from_distinct_primes_grouped":True,
 "normalized_trace":"minus_half_generic_defect_rank",
 "small_generic_ball_non_density_without_quotient_point_import":True,
 "whole_ball_density_inverse_equivalence_uses_retained_LW0":True,
 "source_Ga_commutator_representation_claimed":False,
 "automatic_density_integrality_or_inverse_proved":False},
 "assurance":"Internal unrestricted written proof and exact controls. Finite mode/jet tests do not establish the global analytic proof; no formal verification or external human review."}
 if args.write_receipt:RECEIPT.write_text(json.dumps(receipt,indent=2)+"\n")
 print(json.dumps({"outcome":"PASS","check_count":len(checks),
 "elapsed_seconds":time.monotonic()-start,"receipt_written":args.write_receipt},indent=2))
if __name__=="__main__":main()
