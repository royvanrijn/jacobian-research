"""Exact replay for the boundary-degree invariant and fourth-power gap.

Universal conclusions rely on the retained induction and mixed-rank
proof, not the finite controls here. Arithmetic uses exact integers
and Fractions, including actual ramified quadratic normalization.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial,prod,inf
import hashlib,json,platform

HERE=Path(__file__).resolve().parent
PROOF=HERE/"LOWER_TWO_FOURTH_POWER_GAP_SFC.tex"
IMPORTS={
 "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex":
 "9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07",
 "LOWER_TWO_PRIME_GAP_SFC.tex":
 "8453b5b0ccd220573ef106fa26790e1f691885970b8c2116f036ed893ab2c200",
}
for name,digest in IMPORTS.items():
 assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,(name,"changed exact import")
counts={}
def vp(x,p):
 x=F(x)
 if not x:return inf
 a,b=abs(x.numerator),x.denominator;v=0
 while a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v
def bval(p,s):return (p**s-1)//(p-1)
def trim(A):
 A=A.copy()
 while len(A)>1 and A[-1]==0:A.pop()
 return A
def add(A,B,sgn=1):
 C=[F(0)]*max(len(A),len(B))
 for i,a in enumerate(A):C[i]+=a
 for i,b in enumerate(B):C[i]+=sgn*b
 return trim(C)
def scale(A,c):return trim([F(c)*a for a in A])
def mul(A,B):
 C=[F(0)]*(len(A)+len(B)-1)
 for i,a in enumerate(A):
  for j,b in enumerate(B):C[i+j]+=a*b
 return trim(C)
def shift(A,k):return [F(0)]*k+A
def linproduct(values):
 A=[F(1)]
 for c in values:A=mul(A,[F(c),F(1)])
 return A
def apply_A(P,c):
 Q=[F(0)]*(len(P)+1)
 for k,a in enumerate(P):
  Q[k]+=(c-k)*a;Q[k+1]+=a
 return trim(Q)
def eval_operator(P,c):
 Q=[F(0)]
 for a in reversed(P):
  Q=apply_A(Q,c);Q[0]+=a
 return trim(Q)
def cyclic(P,q,u):
 C=[F(0)]*q
 for k,a in enumerate(P):C[k%q]+=a*u**(k//q)
 return C

# Complete small boundary-saturation Hasse-derivative controls.
hasse=0;basis=0;transfers=0
for p0,s0 in ((2,2),(2,3),(3,2),(5,2)):
 q=p0**s0;qp=p0**(s0-1)
 assert (bval(p0,s0)-s0)-(qp-1)==bval(p0,s0-1)-(s0-1)
 transfers+=1
 C=[0]*q;C[q-1]=1;C[qp-1]=-1
 assert next(j for j,c in enumerate(C) if c%p0)==qp-1
 for e in range(1,p0):
  for r in range(qp):
   h=sum(comb(j,r)*C[j]*e**(j-r) for j in range(r,q))
   assert h%p0==0
   hasse+=1
 # The unit-basis interface after seed reduction, for every scalar residue.
 for c in range(p0):
  V=[1]
  for r in range(qp):
   assert len(V)==r+1 and V[-1]%p0==1
   V=[int(x)%p0 for x in apply_A(V,c)]
   V=trim(V)
  basis+=1
counts["boundary_saturated_Hasse_derivatives"]=hasse
counts["complete_small_scalar_unit_basis_cases"]=basis
counts["exact_stage_constant_transfer_controls"]=transfers

# Actual arithmetic in Q(sqrt(p)); integer and half-integer valuations
# cannot cancel one another, so the valuation formula is exact.
def qa(a=0,b=0):return F(a),F(b)
def qm(x,y,p):return (x[0]*y[0]+p*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qi(x,p):
 den=x[0]*x[0]-p*x[1]*x[1]
 assert den
 return x[0]/den,-x[1]/den
def qs(x,c):return x[0]*c,x[1]*c
def qv(x,p):return min(vp(x[0],p),vp(x[1],p)+F(1,2))
fractional=[]
for p0,s0 in ((2,3),(3,2),(5,2),(2,4)):
 q=p0**s0;qp=p0**(s0-1);Bs=bval(p0,s0)-s0
 C={0:qa(p0**(Bs+1)),1:qa(0,1),2:qa(1),q+2:qa(0,1)}
 r0=min(r for r,c in C.items() if qv(c,p0)==0)
 assert r0<=qp-1
 M,rstar=min((qv(c,p0)+r,r) for r,c in C.items())
 assert M==F(3,2)
 lam=qs(C[rstar],p0**rstar)
 assert qv(lam,p0)==M
 Q={r:qm(qs(c,p0**r),qi(lam,p0),p0) for r,c in C.items()}
 assert all(qv(c,p0)>=0 for c in Q.values())
 assert min(qv(c,p0) for c in Q.values())==0
 assert max(r for r,c in Q.items() if qv(c,p0)==0)<=qp-1
 assert all(qv(c,p0)>0 for r,c in Q.items() if r>=qp)
 assert qv(Q[0],p0)>bval(p0,s0-1)-(s0-1)
 # A coefficient of higher ordinary degree is present but reduces to zero.
 assert max(C)>q-1 and max(r for r,c in C.items() if qv(c,p0)==0)<=q-1
 # Dividing by an element of strictly larger valuation is invalid.
 bad={r:qs(c,F(1,p0)) for r,c in Q.items()}
 assert min(qv(c,p0) for c in bad.values())<0
 fractional.append({"p":p0,"s":s0,"ordinary_degree":max(C),
                    "actual_minimum":str(M),"first_unit_index":r0,
                    "next_constant_valuation":str(qv(Q[0],p0))})
counts["actual_fractional_coefficient_normalizations"]=len(fractional)
counts["over_normalization_rejected"]=len(fractional)

# All-stage sharpness family, controlled on exact full operators.
sharp=[]
for p0,s0 in ((2,1),(2,2),(2,3),(3,1),(3,2),(5,1),(5,2)):
 q=p0**s0;u=(-1)**q*factorial(q)
 chi=linproduct([-j for j in range(1,q+1)])
 chi[0]-=u
 assert chi[0]==0
 P=chi[1:]
 assert len(P)==q and P[-1]==1
 H=sum((F(1,j) for j in range(1,q+1)),F(0))
 assert P[0]==(-1)**(q-1)*factorial(q)*H
 assert vp(P[0],p0)==bval(p0,s0)-s0
 assert vp(u,p0)==bval(p0,s0)
 operator=apply_A(eval_operator(P,q),q)
 expected=[F(-u)]+[F(0)]*(q-1)+[F(1)]
 assert operator==expected
 assert not any(cyclic(operator,q,u))
 sharp.append({"p":p0,"s":s0,"q":q,"constant_valuation":vp(P[0],p0),
               "radius_valuation":vp(u,p0),"actual_cyclic_collision":True})
 if q==4:
  assert P==list(map(F,[-50,35,-10,1])) and u==24
  corrupted=P.copy();corrupted[0]+=1
  assert any(cyclic(apply_A(eval_operator(corrupted,q),q),q,u))
counts["whole_operator_sharpness_controls"]=len(sharp)
counts["corrupted_sharpness_control_rejected"]=1

def Bfamily(n,q,gamma):
 B=[F(1)];out={0:B}
 for r in range(n+q-1):
  C=[F(0)]*(len(B)+1)
  for k,a in enumerate(B):
   C[k]+=(gamma+q*r+q-1-k)*a
   C[k+1]+=a
  B=trim(C)
  if r+1>=n:out[r+1]=B
 return out
relations=[];root_coefficients=0
specs=[
 (2,2,1,[F(1),F(2),F(3)]),
 (2,2,2,[F(1,3),F(2,3),F(5,3)]),
 (2,2,5,[F(0),F(1),F(1)]),
 (2,3,4,[F(1)]+[F(0)]*5+[F(1)]),
 (3,2,3,[F(1)]+[F(0)]*6+[F(1)]),
]
for p0,s0,n,a in specs:
 q=p0**s0;b=bval(p0,s0);Bs=b-s0;gamma=2*n+1
 assert len(a)==q-1 and all(vp(x,p0)>=0 for x in a) and min(vp(x,p0) for x in a)==0
 kap=F(q,(q+2)*((q+2)*n+q+1));u=p0**b
 assert vp(n,p0)>b-2*s0 and vp(kap,p0)>-s0
 P=[F(0)]
 for i,ai in enumerate(a,1):
  Ri=linproduct([q*j for j in range(n,n+i)])
  P=add(P,scale(Ri,ai))
 P[0]+=a[-1]*kap*u
 assert all(vp(c,p0)>=0 for c in P)
 assert min(vp(c,p0) for c in P)==0
 assert max(j for j,c in enumerate(P) if vp(c,p0)==0)<=q-1
 assert vp(P[0],p0)>Bs
 assert vp(P[0],p0)>=min(s0+vp(n,p0),vp(kap,p0)+b)
 fam=Bfamily(n,q,gamma)
 left=[F(0)]
 for i,ai in enumerate(a,1):left=add(left,scale(fam[n+i],ai))
 left=add(left,scale(shift(fam[n],q),a[-1]*kap))
 right=eval_operator(P,gamma+q-1)
 for ell in range(n):
  right=add(apply_A(right,gamma+q-1),scale(right,q*ell))
 difference=add(left,right,-1)
 expected=scale(add(shift(fam[n],q),scale(fam[n],u),-1),a[-1]*kap)
 assert difference==expected
 assert cyclic(left,q,u)==cyclic(right,q,u)
 for k in range(n+1):
  coefficient=comb(n,k)*prod(range(gamma,gamma+q*k))
  assert vp(coefficient,p0)>=b*k
  root_coefficients+=1
 if q==4:assert vp(kap,2)==1
 relations.append({"p":p0,"s":s0,"n":n,
                   "relation_constant_valuation":vp(P[0],p0),
                   "boundary_threshold":Bs,"actual_whole_operator_identity":True})
counts["literal_relation_polynomials_and_cyclic_transfer_controls"]=len(relations)
counts["literal_Gamma_root_bound_coefficients"]=root_coefficients
# A failed index condition is not promoted to a counterexample.
assert not (vp(1,2)>bval(2,3)-6)
counts["outside_prime_power_index_condition_control"]=1

receipt={
 "kind":"internal exact replay; all-index claims rely on the written induction",
 "passed":True,
 "scope":{
  "boundary_invariant":"all p prime,s,n>=1 in the exact closed AR class, degbarP<=p^s-1 and vP0>b_s-s",
  "ordinary_degree_unrestricted":True,
  "strict_bound_sharp_at_every_stage":True,
  "ordinary_support_two_six_all_indices":True,
  "prime_power_index_condition":"q=p^s>=3, v_p(n)>b_s-2s and v_p(kappa_q)>-s",
  "all_ordinary_degree_six_binomials":False,
  "unrestricted_binomial_SFC_or_general_FC_SFC":False,
  "boundary_sharpness_is_SFC_counterexample":False,
 },
 "controls":counts,
 "fractional_normalization_controls":fractional,
 "sharpness_controls":sharp,
 "relation_controls":relations,
 "exact_imports":{
  "stage_change":"ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex lem:artin-hasse-closure, exact class and displayed block identities",
  "mixed_rank":"LOWER_TWO_PRIME_GAP_SFC.tex thm:mixed, necessary direction at gap4 and general prime powers",
 },
 "python":platform.python_version(),
 "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in [PROOF,Path(__file__),*(HERE/name for name in IMPORTS)]},
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "command":"python research/factorial/image-rigidity/lower_two_fourth_power_gap_sfc.py",
}
(HERE/"lower_two_fourth_power_gap_sfc_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
