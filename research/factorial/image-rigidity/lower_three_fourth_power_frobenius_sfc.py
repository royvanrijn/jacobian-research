"""Exact internal replay for the ordinary {3,7} Frobenius minor.

The written proof supplies the all-index quantifiers and restricted-series
integrality. Actual AH reconstruction below checks exact constant/linear
outer coefficients; it does not approximate the proof by an index prefix.
"""
from pathlib import Path
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
PROOF=HERE/"LOWER_THREE_FOURTH_POWER_FROBENIUS_SFC.tex"
IMPORTS={
 "LOWER_THREE_FOURTH_POWER_RESIDUE_REDUCTION.tex":
 "8d0806d61de2a8af064d1d499ed3dec6b3e3e95dbe8f4cd4a2796009911d6972",
 "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex":
 "9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07",
}
for name,h in IMPORTS.items():
 assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,(name,"changed exact input")
counts={}
def count(k,n=1):counts[k]=counts.get(k,0)+n
t,Y,n0,a0,w0=S.symbols('t Y n a w')
c0=3*n0+4;r0=4*n0+1
R1=t+4*n0;R2=(t+4*n0)*(t+4*n0+4)
lam=(4*n0+5)*(1-c0)
T3=(t+4*n0)*(t-1)*(t+n0+1)
assert S.expand(R2*(t-c0)-lam*R1-T3)==0
count("generic_fixed_column_basis_identities")
assert S.expand(R2*(t-c0)-lam*R1-(t+4*n0)*(t-1)*(t+n0+2))!=0
count("corrupted_column_basis_rejections")
alpha,beta,alpha3,beta3=S.symbols('alpha beta alpha3 beta3')
T=[R1,a0*R2-w0,T3]
pred0=[r0,a0*r0*(r0+4)-w0+4*a0*alpha,4*(5*n0+3)*alpha+8*alpha3]
pred1=[1,2*a0*(r0+2+beta),r0*(n0+2)+2*(5*n0+3)*beta+4*beta3]
for p,b0,b1 in zip(T,pred0,pred1):
 q=S.Poly(S.expand(p.subs(t,1+2*Y)),Y)
 got0=q.nth(0)+q.nth(2)*alpha+q.nth(3)*alpha3
 got1=(q.nth(1)+q.nth(2)*beta+q.nth(3)*beta3)/2
 assert S.expand(got0-b0)==0
 assert S.expand(got1-b1)==0
 count("generic_normalized_Taylor_coordinates",2)
AA,TT,VV,OO=S.symbols('A T V O')
av,tv,nv,ov=2*AA,1+2*TT,1+2*VV,1+2*OO
mat=S.Matrix([[1,av-w0,0],[1,0,tv],[2*nv,-w0*ov,nv]])
delta=S.Poly(S.expand(mat.det()-(av*nv+w0*(tv*ov-nv))),AA,TT,VV,OO,w0)
assert all(int(x)%4==0 for x in delta.coeffs())
count("generic_even_determinant_mod4_coefficients",len(delta.coeffs()))
odd=S.Matrix([[1,a0-w0,0],[1,0,1],[0,w0,0]])
assert all(int(x)%2==0 for x in S.Poly(S.expand(odd.det()-w0),w0,a0).coeffs())
count("generic_odd_minor_reduction")
d0=(n0+2)/2
assert S.expand(T3.subs(t,1+2*Y)-4*Y*(Y+d0)*(r0+2*Y))==0
count("generic_even_column_divisibility_identity")
# The even residue class gives all constants used in the linear coefficient.
h=S.symbols('h',integer=True);nn=64*h+54
for expression,modulus in [
 ((3*nn+4)/2-3,4),((3*nn+2)/2-2,4),
 ((nn+2)/2,4),(4*nn+1-1,4),(nn+1-3,4)]:
 assert all(S.denom(x)==1 and int(x)%modulus==0 for x in S.Poly(S.expand(expression),h).coeffs())
 count("generic_even_class_congruences")
z,w=S.symbols('z w')
series=S.Poly(S.series(S.exp(-z-z*z/2-z**4/4),z,0,8).removeO(),z)
coeff=[series.nth(i) for i in range(8)]
def op(c):
 A=S.zeros(4)
 for j in range(4):
  A[j,j]=c-j
  if j+1<4:A[j+1,j]=1
  if j+2<4:A[j+2,j]=-2
 return A
def seed(j):
 return S.Matrix([coeff[j+2*r]*(-2)**r for r in range(4)])
def prod_action(A,step,n,v):
 for r in range(n):
  v=(A+step*r*S.eye(4))*v
 return v
def r2(A,n):
 I=S.eye(A.rows)
 return (A+4*n*I)*(A+(4*n+4)*I)
def pol(A,n,a):
 I=S.eye(A.rows)
 return [A+4*n*I,a*r2(A,n),(A+4*n*I)*(A-I)*(A+(n+1)*I)]
def cyc01(f,external=None):
 # exact constant and linear coefficients in outer w, with inner x²=2w.
 lo=S.Matrix(f[:2]);hi=2*S.Matrix(f[2:4])
 if external is not None:hi+=S.Matrix(external[:2])
 return lo,hi
def val2(q):
 q=S.Rational(q)
 if not q:return None
 a,b=int(q.p),int(q.q)
 return (abs(a)&-abs(a)).bit_length()-(b&-b).bit_length()
def mod2(q):
 q=S.Rational(q);assert val2(q) is None or val2(q)>=0
 return int(q.p)%2*pow(int(q.q),-1,2)%2
def det01(M0,M1):
 d0=M0.det();d1=S.Integer(0)
 for j in range(3):
  M=M0.copy();M[:,j]=M1[:,j];d1+=M.det()
 return S.cancel(d0),S.cancel(d1)
results=[]
for n in [13,45,54,77,109,118,141,173,182,205,237,246]:
 c=3*n+4
 jr=c%2;jn=1-jr
 cr=(c-jr)//2;cn=(c-jn-1)//2
 LR,LN=op(cr),op(cn)
 a=S.Rational(7*(7*n+5)*(7*n+6),32)
 assert val2(a)>=0
 hR,hN=seed(jr),seed(jn)
 AN=S.eye(4)+2*LN
 hat=prod_action(AN,4,n,hN)
 v0,v1=cyc01(hat),cyc01(LN*hat)
 V0=S.Matrix.hstack(v0[0],v1[0]);V1=S.Matrix.hstack(v0[1],v1[1])
 I0=V0.inv();I1=-I0*V1*I0
 Fs=pol(AN,n,a)
 N0,N1=S.zeros(2,3),S.zeros(2,3)
 for j,f in enumerate(Fs):
  f0=f*hat
  external=-hat if j==1 else None
  y0,y1=cyc01(f0,external)
  N0[:,j]=I0*y0
  N1[:,j]=I0*y1+I1*y0
 N0[1,:]/=2;N1[1,:]/=2
 omega=prod_action(LR,2,n,hR)
 Fs=pol(2*LR,n,a)
 R0,R1=S.zeros(2,3),S.zeros(2,3)
 for j,f in enumerate(Fs):
  y0,y1=cyc01(f*omega,-omega if j==1 else None)
  R0[:,j]=y0;R1[:,j]=y1
 row=0 if cr%2 else 1
 M0=S.Matrix.vstack(N0,R0[row,:]);M1=S.Matrix.vstack(N1,R1[row,:])
 if n%2==0:
  M0[:,2]/=2;M1[:,2]/=2
 for v in list(M0)+list(M1):
  assert val2(v) is None or val2(v)>=0
  count("actual_normalized_block_coefficient_integrality")
 if n%2:
  P0=S.Matrix([[1,mod2(a),0],[1,0,1],[0,0,0]])
  P1=S.Matrix([[0,1,0],[0,0,0],[0,1,0]])
  for actual,wanted in zip(list(M0)+list(M1),list(P0)+list(P1)):
   assert mod2(actual)==int(wanted)
   count("actual_odd_block_low_coefficient_reductions")
 else:
  def coordinates(f):
   y0,y1=cyc01(f)
   return I0*y0,I0*y1+I1*y0
  v20,v21=coordinates(LN*LN*hat)
  v30,v31=coordinates(LN*LN*LN*hat)
  rr=4*n+1;dd=(n+2)//2
  tau0=2*v30[1]+(rr+2*dd)*v20[1]+rr*dd
  tau1=2*v31[1]+(rr+2*dd)*v21[1]
  om0,om1=cyc01(omega);nu0,nu1=cyc01(LR*omega)
  P0=S.Matrix([[1,a,0],[1,0,tau0],[2*nu0[0],0,nu0[0]]])
  P1=S.Matrix([[0,-1,0],[0,0,tau1],[2*nu1[0],-om0[0],nu1[0]]])
  for actual,wanted in zip(list(M0)+list(M1),list(P0)+list(P1)):
   vv=val2(actual-wanted);assert vv is None or vv>=2
   count("actual_even_block_low_coefficients_mod4")
  assert v20[1]==2*cn-1
  assert nu0[0]==cr*om0[0]
  assert val2(nu1[0]) is None or val2(nu1[0])>=1
  assert mod2((tau0-1)/2)==0
  assert mod2((tau0*om0[0]-nu0[0])/2)==1
  count("actual_even_constant_and_derivative_identities",5)
 g0,g1=det01(M0,M1)
 if n%2:
  assert mod2(g0)==0 and mod2(g1)==1
 else:
  assert val2(g0) is None or val2(g0)>=1
  assert val2(g1)==1
  assert mod2(g0/2)==mod2(a/2)
 count("actual_normalized_minor_linear_coefficient_cases")
 results.append({"n":n,"cR":cr,"cN":cn,"v_a":val2(a),
 "normalized_matrix_low_coefficients_integral":True,
 "v_G0":val2(g0),"v_G1":val2(g1),"claimed_derivative_unit":True})

# Independent finite Artin--Hasse recurrence and coefficient integrality.
for j,q in enumerate(coeff):
 assert val2(q) is None or val2(q)>=0
 count("actual_AH_inverse_coefficient_integrality")
 if j:
  expr=j*q+sum(coeff[j-2**r] for r in range(3) if 2**r<=j)
  assert expr==0
  count("actual_AH_inverse_Euler_recurrence")
receipt={
 "passed":True,
 "kind":"internal exact replay; universal proof is the retained TeX",
 "scope":{
  "ordinary_support_three_seven_all_indices":True,
  "arbitrary_positive_Gamma_shape_extension":False,
  "all_ordinary_degree_seven_binomials":False,
  "unrestricted_binomial_or_FC_SFC":False,
  "actual_AH_controls":"exact first two outer coefficients; inner degrees0..3 suffice because Euler operations never lower inner degree",
  "unbounded_integrality_and_injectivity":"written restricted-series and conjugacy proof",
 },
 "controls":counts,"total_checks":sum(counts.values()),"actual_AH_inverse_first_eight":[str(x) for x in coeff],
 "actual_normalized_block_profiles":results,
 "source_sha256":{
  PROOF.name:hashlib.sha256(PROOF.read_bytes()).hexdigest(),
  Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),**IMPORTS,
 },
 "exact_imports":{
  "remaining_classes":"LOWER_THREE_FOURTH_POWER_RESIDUE_REDUCTION.tex cor:sixty-one-sixty-four",
  "rank_identity":"same source lem:mixed-matrix",
  "strict_radius":"same source proof of cor:fifteen-sixteenths",
  "stage_and_seed_closure":"ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex lem:artin-hasse-closure",
 },
 "python":platform.python_version(),"sympy":S.__version__,
 "independent_replay":False,"formal_verification":False,"external_review":False,
 "priority_claim":False,"novelty_claim":False,
 "command":".venv/bin/python research/factorial/image-rigidity/lower_three_fourth_power_frobenius_sfc.py",
}
out=HERE/"lower_three_fourth_power_frobenius_sfc_verification.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"passed":True,"controls":counts,"receipt":str(out)},indent=2))
