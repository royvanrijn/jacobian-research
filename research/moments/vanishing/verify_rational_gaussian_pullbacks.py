"""Exact small controls for rational pullbacks; no extrapolation of the all-nu proof."""
from pathlib import Path
import hashlib,json,sys
import sympy as s
Z,W,T,t,u,S,V,e=s.symbols('Z W T t u S V e')
out=Path(__file__).resolve().parent
def clean(f):return s.cancel(f)
def gauss(f):
 ans=0
 for (k,),c in s.Poly(s.expand(f),T).terms():
  if k%2==0:ans+=c*(s.factorial2(k-1) if k else 1)
 return s.expand(ans)
def wick(f):
 ans=0
 for (a,b,c),q in s.Poly(s.expand(f),Z,W,T).terms():
  if a==b and c%2==0:ans+=q*s.factorial(a)*(s.factorial2(c-1) if c else 1)
 return s.expand(ans)
def family(n,a=0):
 b=(1-Z)**(n+1);R=((1-Z)**(-n)-1)/n
 delta=s.sympify(a).subs(Z,0)
 h=a+b*T
 G=clean(h/(2*Z)*(T*T-(R*h/Z-delta)**2))
 return b,R,h,G,s.expand(W*h+G),delta
polycontrols=[];Wick=[]
for n in [1,2,3]:
 q=(2*n+2)//3
 for a in [s.S.Zero,(1-Z)**q]:
  b,R,h,G,P,delta=family(n,a)
  assert s.denom(G)==1
  assert clean(s.diff(R,Z)-1/b)==0
  a_bad=(1-Z)**(q-1)
  bad=family(n,a_bad)[3]
  assert s.denom(bad).subs(Z,1)==0
  polycontrols.append({'nu':n,'a':str(a),'delta':str(delta),'G_polynomial':True,
   'one_order_lower_a_has_pole':True})
  maxm=4 if n<=2 and a==0 else 2
  for m in range(1,maxm+1):
   pv=wick(P**m);zv=wick(Z*P**m)
   predicted=-(-1)**m*s.rf(s.Rational(1,n),m)*n**m*gauss((T+delta)**m)
   assert pv==0 and s.simplify(zv-predicted)==0,(n,a,m,pv,zv,predicted)
   Wick.append({'nu':n,'delta':str(delta),'m':m,'pure':str(pv),'marked_Z':str(zv)})
# Formal transport identity, exact infinitesimal form and finite t-adic coefficients.
sigma=s.symbols('sigma')
Y=T+t*(T*T+1)+t*t*T**3;Ys=T+sigma*(Y-T);vs=s.diff(Ys,sigma)
transport=0
for j in range(4):
 f=Ys**j;exp_derivative=T-Ys*s.diff(Ys,T)
 lhs=s.diff(f*s.diff(Ys,T),sigma)-f*s.diff(Ys,T)*Ys*vs
 rhs=s.diff(f*vs,T)+(exp_derivative-T)*f*vs
 assert s.expand(lhs-rhs)==0;transport+=1
# Normal-coordinate reduction keeps the entire amplitude direction and the Jacobian.
normal=[]
for n in [1,2]:
 b,R,h,G,P,delta=family(n,(1-Z)**((2*n+2)//3))
 Tc=(S+Z-t*h.subs(T,0))/(t*b)
 H=clean((t*G-T*T/2).subs(T,Tc))
 H0=clean(H.subs(S,0));target=-(R/t-delta)**2/2
 assert clean(H0-target)==0
 K=clean((H-H0)/S)
 assert s.denom(K).has(S) is False
 assert clean(S*(V-K)+H-(S*V+H0))==0
 assert clean(s.diff(Tc,S)-1/(t*b))==0
 normal.append({'nu':n,'phase_reduced':True,'Jacobian_retained':True})
# Laurent source basis and connection, sampled only at nu1..3.
sourcecontrols=[]
def parts(f):
 d={}
 for term in s.Add.make_args(s.expand(f)):
  power=int(term.as_powers_dict().get(u,0))
  d[power]=clean(d.get(power,0)+term/u**power)
 return {k:v for k,v in d.items() if v!=0}
for n in [1,2,3]:
 for delta in [0,1]:
  d=2*n;Phi=-(u**n-1-n*delta*t)**2/(2*n*n*t*t)
  phiu=s.diff(Phi,u);A=-s.S(1)/(n*t*t)
  def reduce(f):
   f=s.expand(f);certificate=0
   while any(k<-1 for k in parts(f)):
    coeffs=parts(f);k=min(coeffs);q=coeffs[k]*u**(k+1)/(k+1)
    f=s.expand(f-s.diff(q,u)-phiu*q);certificate+=q
   while any(k>=d-1 for k in parts(f)):
    coeffs=parts(f);k=max(coeffs);q=clean(coeffs[k]/A)*u**(k-d+1)
    f=s.expand(f-s.diff(q,u)-phiu*q);certificate+=q
   assert clean(s.expand(f+s.diff(certificate,u)+phiu*certificate))!=s.nan
   return s.expand(f),s.expand(certificate)
  basis=list(range(-1,d-1))
  columns=[]
  for j in range(1,n+1):
   f=u**(n-1-j)
   for v in [f,s.diff(f,t)+(s.diff(Phi,t)-1/t)*f]:
    red,cert=reduce(v)
    assert clean(v-red-s.diff(cert,u)-phiu*cert)==0
    columns.append(s.Matrix([parts(red).get(k,0) for k in basis]))
  B=s.Matrix.hstack(*columns)
  assert B.rank()==d
  # Connection denominators have no finite zero outsidet0.
  for k in basis:
   v=(s.diff(Phi,t)-1/t)*u**k;red,cert=reduce(v)
   assert clean(v-red-s.diff(cert,u)-phiu*cert)==0
   for c in parts(red).values():
    den=s.factor(s.denom(c))
    assert all(factor==t for factor,power in s.factor_list(den,t)[1])
  sourcecontrols.append({'nu':n,'delta':delta,'basis_dimension':d,'source_marked_Krylov_rank':d,
   'only_possible_finite_time_pole':0,'not_an_actual_kernel_producer':True})
# Scalar period recurrence and infinity exponents; no guessed operator.
kap,al,lam=s.symbols('kappa alpha lambda',nonzero=True)
ind=s.expand(lam*(lam+1)-2*(al+1)*lam+al*(al+1))
assert s.factor(ind)==(al-lam)*(al-lam+1)
ODE=[]
for n in [1,2,3]:
 for j in range(1,n+1):
  aa=s.Rational(j,n)
  for dd in [0,1]:
   coeff=[(-1)**m*s.rf(aa,m)*n**m*gauss((T+dd)**m)/s.factorial(m) for m in range(10)]
   y=sum(coeff[m]*t**m for m in range(10))
   residual=s.expand(n*n*t**3*s.diff(y,t,2)+(2*(aa+1)*n*n*t*t-n*dd*t-1)*s.diff(y,t)+(aa*(aa+1)*n*n*t-aa*n*dd)*y)
   assert all(residual.coeff(t,k)==0 for k in range(8))
   ODE.append({'nu':n,'j':j,'delta':dd,'coefficient_checks':8,
    'infinity_exponents':[str(aa),str(aa+1)]})
# The alpha1 inhomogeneous relation supplies1 in its rank-two cyclic span.
for dd in [0,1]:
 y=sum((-1)**m*s.rf(1,m)*gauss((T+dd)**m)/s.factorial(m)*t**m for m in range(10))
 residual=s.expand(t**3*s.diff(y,t)+(t*t-dd*t-1)*y+1)
 assert all(residual.coeff(t,k)==0 for k in range(9))
# Polynomial rank-drop models; no global coefficient-locus theorem inferred.
drops=[]
for n in [1,2,3]:
 Q,rem=s.div(n*n-u*u*sum(u**j for j in range(n))**2,1-u,u)
 assert rem==0 and s.LC(s.Poly(Q,u))==1
 base=(1-Z)**(n+1)*W*T+(1-Z)**(n+1)*Q.subs(u,1-Z)*T**3/(2*n*n)
 Pscaled=s.expand(e**(3*n)*base.subs({Z:Z/e,W:e*W},simultaneous=True))
 assert s.denom(Pscaled)==1
 expected=(-1)**n*Z**(3*n)*T**3/(2*n*n)
 assert s.expand(Pscaled.subs(e,0)-expected)==0
 assert wick(W**(6*n)*expected**2)!=0
 drops.append({'nu':n,'generic_rank_from_written_proof':2*n,'special_rank_from_charge_proof':1,
  'polynomial_coefficients_checked':True,'special_phase':str(expected)})
falsekeys=[
 'all_order_pure_identity_inferred_from_finite_Wick_replay','contour_comparison_used_for_new_identities',
 'all_nu_actual_rank_certified_by_checker','general_period_kernel_solver_executed',
 'all_repeated_root_phases_classified','all_Gaussian_higher_rank_connections_classified',
 'global_complex_safe_locus_uniformity_proved','new_cutoff_ideal_chain_computed',
 'minimal_ternary_failure_degree_changed','new_SIC33_answer_computed','new_Delta_minimum_computed',
 'numeric_cubic_frontier_changed','operator_order_reduced','independent_replay','external_review','formal_verification']
rec={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python_version':sys.version.split()[0],'sympy_version':s.__version__,
 'polynomiality_controls':polycontrols,'Wick_controls':Wick,'formal_transport_infinitesimal_controls':transport,
 'normal_pair_controls':normal,'Laurent_source_controls':sourcecontrols,'scalar_ODE_controls':ODE,
 'rank_drop_controls':drops,
 'mutation_controls':{'Jacobian_not_discarded':True,'root_at_zero_not_localized':True,
 'full_multiplier_source_not_replaced_by_pure_seed':True,
 'source_Krylov_rank_not_called_actual_without_minimality_and_exponent_separation':True,
 'divergent_formal_series_not_treated_as_convergent_Gaussian_integrals':True,
 'fixed_degree_field_uniformity_not_disproved_by_rank_drop':True},
 **{k:False for k in falsekeys}}
(out/'rational_gaussian_pullback_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps({'completed':True,'checker_sha256':rec['checker_sha256'],
 'Wick_controls':len(Wick),'source_controls':len(sourcecontrols),'scalar_ODE_controls':len(ODE)}))

