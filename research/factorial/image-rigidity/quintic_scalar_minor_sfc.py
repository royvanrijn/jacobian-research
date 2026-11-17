"""Exact replay for the scalar-minor ordinary quintic proof.

The all-index theorem is established by the written integrality,
reduction, injectivity and real-field descent arguments. This verifier
checks complete reduced-operator cases and full literal minor controls.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import platform

HERE=Path(__file__).resolve().parent
PROOF=HERE/"QUINTIC_SCALAR_MINOR_SFC.tex"
IMPORTS={
 "QUINTIC_THREE_ADIC_PARITY_REDUCTION.tex":
  "97e54c56e407a6ab7866497bbe673eed73c017bbdd0971229f4f8f6912dfd1b5",
 "lower_exponent_initial_and_boundary.tex":
  "e66e02d0354de5cd92d33eef54efea15c387b9f6a0c1b62ed379295e13f07d45",
 "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex":
  "9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07",
 "ALL_UNIVARIATE_TWO_TERM_SFC.tex":
  "112d995972906353784af4f9ceb2b03eb7e71ce3e7bd52adc5a44c6869aeadc0",
 "UNIFORM_PADE_PARABOLA_SFC.tex":
  "ddb8155dcce6308db9d7fdcbc971898b9725e9ac23c83362d6c45a63f4177cc5",
}
for name,digest in IMPORTS.items():
 assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,(name,"changed exact import")
counts={}
def matmul(A,B):
 return [[sum(A[i][k]*B[k][j] for k in range(3))%3 for j in range(3)] for i in range(3)]
def matvec(A,v):
 return [sum(A[i][j]*v[j] for j in range(3))%3 for i in range(3)]
def det3(A):
 return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
        -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
        +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))%3
def minor(C,B,ij):
 i,j=ij
 return (C[i]*B[j]-C[j]*B[i])%3
for gamma_bar in range(3):
 A=[[0]*3 for _ in range(3)]
 for j in range(3):
  A[j][j]=(gamma_bar-1-j)%3
  if j<2:A[j+1][j]=1
 assert matmul(matmul(A,A),A)==A
 powers=[[1,0,0]]
 for j in range(4):powers.append(matvec(A,powers[-1]))
 K=[[powers[j][i] for j in range(3)] for i in range(3)]
 assert det3(K)==1
 for parity in (1,2):
  B,C,D=powers[parity:parity+3]
  assert D==B
  rows=next(ij for ij in combinations(range(3),2) if minor(C,B,ij))
  # The two coefficients of det[C,D+X B] agree and are nonzero:
  # this is a whole polynomial statement, not testing only F3 points.
  assert minor(C,D,rows)==minor(C,B,rows)!=0
  if gamma_bar==0:
   assert minor(C,B,(0,1))==(1 if parity==1 else 2)
counts["complete_reduced_gamma_and_positive_power_parity_cases"]=6
counts["complete_reduced_affine_minor_coefficients"]=12
counts["fixed_row_01_parity_sign_cases"]=2

def trim(A):
 A=A.copy()
 while len(A)>1 and A[-1]==0:A.pop()
 return A
def add(A,B,sign=1):
 C=[0]*max(len(A),len(B))
 for i,a in enumerate(A):C[i]+=a
 for i,b in enumerate(B):C[i]+=sign*b
 return trim(C)
def mul(A,B):
 C=[0]*(len(A)+len(B)-1)
 for i,a in enumerate(A):
  for j,b in enumerate(B):C[i+j]+=a*b
 return trim(C)
def scale_variable(A,q):
 return [a*q**j for j,a in enumerate(A)]
def Bfamily(n,gamma):
 B=[1];out={0:B}
 for r in range(n+2):
  C=[0]*(len(B)+1)
  for k,a in enumerate(B):
   C[k]+=(gamma+3*r+2-k)*a
   C[k+1]+=a
  B=trim(C)
  if r+1 in (n,n+1,n+2):out[r+1]=B
 return out
def cyclic(B,row):
 return trim(B[row::3] or [0])
def affine_check(F,epsilon):
 assert F[0]%3==F[1]%3==epsilon
 assert all(a%3==0 for a in F[2:])
controls=[];all_coefficients=0;bad_seed=None
for n in (1,10,19,28,37,46,55,64):
 assert n%9==1
 gamma=2*n+1;q=5*(5*n+4)//3
 assert q%3==0
 family=Bfamily(n,gamma)
 b0,b1=[cyclic(family[n],row) for row in (0,1)]
 c0,c1=[cyclic(family[n+1],row) for row in (0,1)]
 d0,d1=[cyclic(family[n+2],row) for row in (0,1)]
 cd=add(mul(c0,d1),mul(c1,d0),-1)
 cb=add(mul(c0,b1),mul(c1,b0),-1)
 F=add(scale_variable(cd,q),[0]+scale_variable(cb,q))
 # A separately assembled version substitutes before multiplying columns.
 C0,C1,D0,D1,B0,B1=[scale_variable(A,q) for A in (c0,c1,d0,d1,b0,b1)]
 F_direct=add(mul(C0,add(D1,[0]+B1)),mul(C1,add(D0,[0]+B0)),-1)
 assert F==F_direct
 epsilon=1 if n%2 else 2
 affine_check(F,epsilon)
 all_coefficients+=len(F)
 controls.append({"n":n,"gamma":gamma,"q":q,"minor_degree":len(F)-1,
                  "constant_linear_reduction":epsilon,
                  "every_higher_coefficient_divisible3":True,
                  "substitute_before_and_after_multiplication_agree":True})
 if n==10:bad_seed=(F.copy(),epsilon)
counts["full_literal_minors"]=len(controls)
counts["literal_minor_coefficients_checked"]=all_coefficients
F,epsilon=bad_seed
corruptions=[]
bad=F.copy();bad[1]+=3-epsilon;corruptions.append(bad)
bad=F.copy();bad[2]+=1;corruptions.append(bad)
for bad in corruptions:
 try:affine_check(bad,epsilon)
 except AssertionError:pass
 else:raise AssertionError("corrupted critical affine-reduction input accepted")
counts["corrupted_linear_or_higher_coefficient_controls"]=2
receipt={
 "kind":"internal exact replay; universal proof is the retained TeX source",
 "passed":True,
 "scope":{
  "ordinary_support_two_five_all_indices":True,
  "ordinary_univariate_binomials_degree_at_most_five":True,
  "all_common_root_conjugates_argument":"real-field descent in characteristic zero",
  "all_index_basis":"integer cyclic columns and affine reduction of one fixed row-pair minor",
  "degree_six_frontier":"equivalence to supports{2,6},{3,6},{4,6} at all terminalN>=6",
  "three_degree_six_families_proved":False,
  "all_ordinary_univariate_binomials":False,
  "all_positive_real_Gamma_shapes":False,
  "unrestricted_FC_or_SFC":False,
 },
 "controls":counts,
 "full_minor_controls":controls,
 "exact_imports":{
  "mixed_rank":"QUINTIC_THREE_ADIC_PARITY_REDUCTION.tex lem:mixed-input; original lower_exponent_initial_and_boundary.tex thm:gap-three-offset-two-rank",
  "ordinary_degree_five_classification":"QUINTIC_THREE_ADIC_PARITY_REDUCTION.tex cor:degree-five-ray",
  "ordinary_terminal_at_most_five":"lower_exponent_initial_and_boundary.tex thm:ordinary-fifth",
  "constant_plus_power":"ALL_UNIVARIATE_TWO_TERM_SFC.tex thm:all",
  "lower_one_prime_gap_five":"ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex cor:prime-power-lower-one",
  "adjacent_exponents":"UNIFORM_PADE_PARABOLA_SFC.tex thm:main",
 },
 "python":platform.python_version(),
 "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in [PROOF,Path(__file__),*(HERE/name for name in IMPORTS)]},
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "command":"python research/factorial/image-rigidity/quintic_scalar_minor_sfc.py",
}
(HERE/"quintic_scalar_minor_sfc_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps(receipt,indent=2))
