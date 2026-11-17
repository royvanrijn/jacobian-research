"""One-use exact controls for Fuchsian witnesses and scoped Gaussian bridges."""
from pathlib import Path
import hashlib,json,platform
import sympy as s
import gaussian_gfunction_regular_core as G
import gaussian_saddle_arithmetic as A

p=Path(__file__).resolve().parent;family='gaussian_gfunction_regular_core'
outputs=[p/(family+tail+'.json') for tail in ['_certificate','_runtime','_verification']]
assert all(not q.exists() for q in outputs),'Sealed outputs exist; do not rerun.'
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
dump=lambda x:json.dumps(x,indent=2,ensure_ascii=False)+'\n'
t=G.t;e=s.Symbol('e');z,w,x=s.symbols('z w x')
elliptic=[-s.Rational(1,4),1-2*t,t*(1-t)]
fixtures=[('algebraic_quadratic',[-s.Rational(1,2),1-t]),('Witt_unit',[-4*t,1-4*t**2]),('Witt_basis_lclm',[0,-(1+8*t**2),t*(1-4*t**2)]),('abstract_elliptic',elliptic),('abstract_logarithm',[0,-1,1-t])]
models=[]
for name,op in fixtures:
 cert=G.fuchsian_certificate(op)
 for c in cert['finite_witnesses']:
  q=s.sympify(c['factor'],locals={'t':t});v=s.sympify(c['quotient'],locals={'t':t})
  assert s.expand(op[c['derivative_order']]-q**c['required_power']*v)==0
 models.append(dict(name=name,model=cert))

parameters=[]
for name,op in [('moving_elliptic',[-s.Rational(1,4)*e,1-2*e*t,t*(1-e*t)]),('colliding_poles',[1,t*(t-e)]),('rational_parameter',[1/e,1-t])]:
 cert=G.fuchsian_certificate(op,(e,))
 guard=s.sympify(cert['specialization_guard'],locals={'e':e})
 assert guard.subs(e,0)==0
 special=[]
 for val in [1,2,-1]:
  if guard.subs(e,val)!=0:
   c=G.fuchsian_certificate([s.sympify(v).subs(e,val) for v in op]);special.append(dict(value=val,order=c['order']))
 parameters.append(dict(name=name,model=cert,good_specializations=special))

indicial=[]
for point,want in [(0,s.Symbol('rho')**2),(1,s.Symbol('rho')**2),(s.oo,-(s.Symbol('rho')-s.Rational(1,2))**2)]:
 v=G.indicial_at(elliptic,point);assert s.expand(v-want)==0
 indicial.append(dict(point=str(point),polynomial=str(v)))

# The elliptic G-function is an ABSTRACT control, not a Gaussian realization.
coeff=[s.binomial(2*k,k)**2/s.Integer(16)**k for k in range(10)]
f=sum(c*t**k for k,c in enumerate(coeff))
assert all((c*16**k).is_Integer for k,c in enumerate(coeff))
assert G.derivative_recurrence(elliptic,f,7)==0
for k in range(9):assert (k+1)**2*coeff[k+1]==(k+s.Rational(1,2))**2*coeff[k]

def wick_witt(poly):
 ans=0
 for (aa,bb),c in s.Poly(s.expand(poly),z,w).terms():
  if aa==bb:ans+=c*s.factorial(aa)
 return s.expand(ans)

P=w*(1+z**2);Delta=1-4*t**2;f0=1/s.sqrt(Delta);g0=(1-s.sqrt(Delta))/(2*t)
original=[]
for Q in [s.Integer(1),z,z**2,w]:
 v=sum(wick_witt(Q*P**k)*t**k/s.factorial(k) for k in range(9))
 if Q==1:expected=f0
 elif Q==z:expected=(f0-1)/(2*t)
 elif Q==z**2:expected=f0*g0**2
 else:
  # Insertion W is the derivative in an external Z translation.
  # Implicit equation g=t(1+(g+u)^2), J=1-2t(g+u).
  # d_u(1/J) at u0 = 2t/J^3.
  expected=2*t/Delta**s.Rational(3,2)
 assert s.series(expected-v,t,0,9).removeO().expand()==0
 original.append(dict(phase=str(P),marker=str(Q),order=8,period=str(s.expand(v)),exact_period=str(expected)))
assert s.cancel(g0-t*(1+g0**2))==0

# Whole-class boundaries use existing all-index proofs, not finite convergence.
linear_f=A.period(x,1,(x,),8)
linear_H=A.jet(s.exp(-t**2/2)*linear_f,8)
assert linear_H==1
quartic=[s.factorial2(4*k-1)/s.factorial(k) if k else s.Integer(1) for k in range(7)]
assert all((quartic[k]*2**k).is_Integer for k in range(7))
ratios=[]
for k in range(6):
 actual=s.cancel(quartic[k+1]/quartic[k]);expected=s.Rational((4*k+1)*(4*k+3),k+1)
 assert actual==expected;ratios.append(str(actual))
# Kloosterman k2 original marker h-1 has a_(m+1)/a_m=(2m+1)^2.
kloo=[s.factorial2(2*k-1)**2 if k else s.Integer(1) for k in range(6)]
assert all(kloo[k+1]==(2*k+1)**2*kloo[k] for k in range(5))

rejects=[]
for name,fn in [('exponential_infinity',lambda:G.fuchsian_certificate([-t,1])),('finite_irregular_pole',lambda:G.fuchsian_certificate([-1,(1-t)**2])),('order_zero',lambda:G.fuchsian_certificate([1])),('zero_leading_operator',lambda:G.fuchsian_certificate([1,0])),('negative_infinity_balance',lambda:G.fuchsian_certificate([1,0,t]))]:
 try:fn()
 except ValueError:rejects.append(name)
 else:raise AssertionError(name)
packet=dict(date='2026-11-17',universal_claims_use_written_proofs=True,classical_ACK_independently_reproved=False,convergence_decided_from_jets=False,abstract_elliptic_Gaussian_realization_claimed=False,cofinite_prime_nilpotence_claimed=False,vanishing_pcurvature_claimed=False,field_safety_implies_purity_claimed=False,Fuchsian_models=models,parameter_models=parameters,elliptic_indicial=indicial,abstract_elliptic_coefficients=list(map(str,coeff)),original_Witt_controls=original,linear_normalized_control=str(linear_H),quartic_divergence_exact_ratio='(4*k+1)*(4*k+3)/(k+1)',quartic_ratio_controls=ratios,pure_Kloosterman_divergence_exact_ratio='(2*m+1)**2',pure_Kloosterman_coefficients=list(map(str,kloo)),rejections=rejects)
imports=['gaussian_saddle_arithmetic.py','CRITICAL_NORMALIZATION_AND_RANK_ONE_ARITHMETIC.md','gaussian_saddle_arithmetic_verification.json','GAUSSIAN_KLOOSTERMAN_QUOTIENTS_AND_FAITHFUL_SYMMETRIC_IMAGES.md','ALGEBRAIC_FIXED_POINT_LIFTS_AND_UNIVERSAL_GAUSSIAN_IMAGES.md']
runtime=dict(python=platform.python_version(),sympy=s.__version__,sealed_imports_sha256={nm:sha(p/nm) for nm in imports})
outputs[0].write_text(dump(packet));outputs[1].write_text(dump(runtime))
receipt=dict(completed=True,kind='exact finite Fuchsian witnesses and boundary controls;universal Gaussian/ACK deductions are written proofs',canonical_proof_sha256=sha(p/'GAUSSIAN_G_FUNCTIONS_AND_GLOBAL_REGULAR_CORE.md'),producer_sha256=sha(p/'gaussian_gfunction_regular_core.py'),checker_sha256=sha(Path(__file__)),certificate_sha256=sha(outputs[0]),runtime_sha256=sha(outputs[1]),exact_controls=dict(Fuchsian_models=len(models),parameter_models=len(parameters),elliptic_indicial=3,original_Witt_markers=len(original),linear_normalization=1,divergent_boundary_families=2,rejections=len(rejects)),independent_replay=False,formal_verification=False,external_review=False)
outputs[2].write_text(dump(receipt));print(dump(receipt))
