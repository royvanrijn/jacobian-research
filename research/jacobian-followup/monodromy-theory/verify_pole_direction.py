"""Exact controls for PROJECTIVE_POLE_RIGIDITY.md; no bounded scalar search."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as S
from pole_direction import wave_tensor,leading_direction,conjugate_affine_pair
from rational_wave import rational_chart
from polynomial_wave import polynomial_chart
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
w,z,t,k,alpha,u=S.symbols("w z t kappa alpha u")
n,m,ell=S.symbols("n m ell",integer=True)
indicial=alpha**2*(m+1)*(m+2)+2*(n+1)*alpha*(m+1)+n*(n+1)
inverse=n*(n+1)*u*u+2*(n+1)*(m+1)*u+(m+1)*(m+2)
check("universal inverse-multiplier indicial equation",S.expand(indicial.subs(alpha,1/u)*u*u-inverse)==0)
check("universal negative discriminant formula",S.expand(S.discriminant(inverse,u)-4*(n+1)*(m+1)*(m+1-n))==0)
check("universal normalized trace",S.cancel(-S.Poly(inverse,u).nth(1)/(2*S.Poly(inverse,u).nth(2))+(m+1)/n)==0)
check("universal projective-degree inequality gap",S.expand(n-((n-2*ell)+ell+1)-(ell-1))==0)
# Full wave-to-ODE identity for nonpolynomial projective ratios and varied lambda.
for nv,R,L in [(4,t*t/(t*t+1),t*t+2*t+3),(6,(t**3+t+1)/(t**3+1),t**4+1),
                (3,k*t,t+1)]:
 rho=R.subs(t,z/w);lam=w**nv*L.subs(t,z/w)
 wave=S.diff(lam,w,2)+2*S.diff(rho*lam,w,z)+S.diff(rho*rho*lam,z,2)
 delta=R-t
 ode=S.diff(delta*delta*L,t,2)+2*(nv+1)*S.diff(delta*L,t)+nv*(nv+1)*L
 check("exact projective ODE identity n="+str(nv),S.cancel(wave-w**(nv-2)*ode.subs(t,z/w))==0)
# Algebraic residue identity: trace multiplication by Q/F' in k[t]/F.
residue_controls=[]
for P,Q in [(t*t,t*t+1),(t**3+t+1,t**3+1),(2*t+1,t+3)]:
 F=P-t*Q;dv=S.degree(F,t)
 check("fixed-point polynomial simple degree="+str(dv),S.gcd(F,S.diff(F,t))==1)
 uu=S.rem(Q*S.invert(S.diff(F,t),F,t),F,t)
 matrix=S.Matrix(dv,dv,lambda i,j:S.Poly(S.rem(uu*t**j,F,t),t).nth(i))
 check("complete algebraic fixed-point residue sum degree="+str(dv),S.trace(matrix)==-1)
 residue_controls.append(dict(R=str(P/Q),degree=int(dv-1),residue_sum=-1))
# Actual homogeneous solutions with poles in rotated directions.
T=S.Matrix([[1,1],[1,2]])
homogeneous=[]
for rv,sv in [(2,0),(4,0),(3,1),(2,3)]:
 minimum=rv*(rv-1)+2*k*(rv-1)*(sv+1)+k*k*(sv+1)*(sv+2)
 def red(e):
  num,den=S.fraction(S.cancel(e))
  return S.rem(S.rem(num,minimum,k)*S.invert(S.rem(den,minimum,k),minimum,k),minimum,k).expand()
 scalar=w**(rv-2)*z**sv
 vector=conjugate_affine_pair(w,k*z,T,S.zeros(2,1),w,z)
 old=T*S.Matrix([w,z])
 As=S.expand(T.det()*scalar.subs({w:old[0],z:old[1]},simultaneous=True))
 check(f"homogeneous rotated tensor wave {rv},{sv}",red(wave_tensor(As,vector[0],vector[1],w,z))==0)
 R=red((vector[1]/vector[0]).subs({w:1,z:t}))
 D=S.cancel(R-t)
 check(f"homogeneous ratio not fixed at infinity {rv},{sv}",S.limit(R,t,S.oo)!=S.oo)
 # The fixed directions are the two old coordinate axes.
 sums=[]
 for point in [-S.Rational(1,2),-S.Integer(1)]:
  ui=red(1/S.diff(D,t).subs(t,point));sums.append(ui)
  multiplicity=sv if point==-S.Rational(1,2) else rv-2
  # rational trace is the average under k -> conjugate root.
  conjugate=-2*(rv-1)/S.Integer(sv+2)-k
  avg=red((ui+ui.subs(k,conjugate))/2)
  check(f"homogeneous local trace {rv},{sv},{point}",S.cancel(avg+S.Rational(multiplicity+1,rv+sv))==0)
 check(f"homogeneous global residue sum {rv},{sv}",red(sum(sums)+1)==0)
 homogeneous.append(dict(r=rv,s=sv,tensor_degree=rv+sv,projective_pair_degree=1,scalar_degree=rv+sv-2))
# Necessary leading tests reject all scalar degrees for these entire pairs.
rejected=[]
for b,c in [(w*w+z*z+1,w*z+w+1),(w**3+1,z**3+w+1),
             (w*w+z*z+1,2*(w*w+z*z)+w),
             (w*(w+z)+1,k*z*(w+z)+w),(w+z,z)]:
 result=leading_direction(b,c,w,z)
 check("all-scalar leading rejection "+str(b),result["rejected"] is True)
 rejected.append(dict(b=str(b),c=str(c),reason=result["reason"]))
for b,c in [(1,w*w),(w,w*w+k*z),(w*w,k*w*z+1)]:
 result=leading_direction(S.sympify(b),c,w,z)
 check("known chart passes necessary leading test "+str(b),result["rejected"] is False)
# Affine direction reduction by a left eigenvector; no wave existence assumed.
for M,b0 in [(S.Matrix([[1,1],[0,2]]),S.Matrix([3,4])),
              (S.Matrix([[1,1],[0,1]]),S.Matrix([1,2])),
              (S.Matrix([[0,1],[0,0]]),S.Matrix([0,1]))]:
 eigen=M.T.eigenvects()[0][2][0]
 row=S.Matrix([list(eigen)])
 second=S.Matrix([[1,0]]) if eigen[1]!=0 else S.Matrix([[0,1]])
 B=row.col_join(second)
 old=B.inv()*S.Matrix([w,z])
 v=(B*(M*old+b0)).applyfunc(S.expand)
 check("affine eigen-coordinate denominator "+str(M.tolist()),S.diff(v[0],z)==0)
 if v[0]!=0:
  numerator,denominator=S.fraction(S.cancel(v[1]/v[0]))
  check("affine ratio in k(w)[z] "+str(M.tolist()),not denominator.has(z))
 else:check("affine zero direction is constant",v[0]==0)
# Exhaustive coefficient identities for the two quadratic top types.
aa,bb,cc,dd,ee,ff=S.symbols("a b0 c d0 e f")
qwave=wave_tensor(1,w*w+aa*z+bb,k*w*z+cc*z+dd,w,z)
expected=2*(k*k+3*k+6)*w*w+4*cc*(k+1)*w+4*aa*(k+1)*z+2*bb*(k+2)+2*cc*cc
check("full quadratic eigen-direction wave coefficients",S.expand(qwave-expected)==0)
b1=S.symbols("b1")
qwave=wave_tensor(1,aa*w+b1*z+bb,w*w+cc*w+ee*z+ff,w,z)
check("full quadratic constant-direction coefficients",S.expand(qwave-4*b1*w-2*(aa*aa+b1*cc+aa*ee+ee*ee))==0)
# The three surviving actions, using previously proved complete classifiers.
f5=polynomial_chart(w*w,[(0,1)],w=w,z=z)
check("quadratic constant scalar A5 case",f5["N"]==5 and f5["geometric_group"]=="A5")
f6=rational_chart(4,0,1/w,w,z,k)
check("quadratic constant scalar A6 case",f6["N"]==6 and f6["geometric_group"]=="A6")
check("quadratic A6 primitive pair",S.gcd(w*w,k*w*z+1)==1)
minimum=k*k+3*k+6
check("quadratic A6 exact tensor wave",S.rem(wave_tensor(1,w*w,k*w*z+1,w,z),minimum,k)==0)
f8=rational_chart(2,0,w*w,w,z,k)
check("quadratic constant scalar A8 case",f8["N"]==8 and f8["geometric_group"]=="A8")
check("quadratic A8 exact tensor wave",S.rem(wave_tensor(1,w,w*w+k*z,w,z),k*k+k+1,k)==0)
# Radial direction itself admits no positive homogeneous scalar.
for degree in [0,1,2,3]:
 scalar=sum(w**j*z**(degree-j) for j in range(degree+1))
 check("radial tensor exclusion scalar degree="+str(degree),S.expand(wave_tensor(scalar,w,z,w,z)-(degree+2)*(degree+3)*scalar)==0)
model=dict(homogeneous_normalizations="constant direction or two linear eigen-directions; primitive degree>=2 impossible for all scalars",
 homogeneous_controls=homogeneous,residue_controls=residue_controls,
 rejected_entire_pairs=rejected,
 leading_filter_boundary="necessary only; passing is not a polynomial wave or Keller realization",
 affine_primitive_pair="coordinate reduction to completed chart; full groups natural alternating/symmetric",
 quadratic_primitive_pair_constant_scalar=dict(full_groups=["A5","A6","A8"],complete=True),
 remaining_general_two_variable_pole_class_complete=False,abstract_PSL2_11_actions_decided=False,full_action_classification_complete=False)
(ROOT/"pole_direction_model.json").write_text(json.dumps(model,indent=2)+"\n")
files=["PROJECTIVE_POLE_RIGIDITY.md","pole_direction.py","verify_pole_direction.py","pole_direction_model.json",
 "rational_wave.py","polynomial_wave.py"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
 assurance="Written universal projective-residue trace proof, homogeneous classification, leading obstruction and complete quadratic constant-scalar reduction. Exact controls are not formal verification or independent whole-proof review.",
 goal_complete=False)
(ROOT/"pole_direction_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"])))
