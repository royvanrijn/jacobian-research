"""Finite exact controls; no numerical contour or universal analytic certification."""
from pathlib import Path
import hashlib,json,sys
import sympy as s
t,x,Z,W,T,a,c,p,q,k,e=s.symbols('t x Z W T a c p q k e')
out=Path(__file__).resolve().parent
def expectation(poly,var):
 ans=0
 for (j,),coef in s.Poly(s.expand(poly),var).terms():
  if j%2==0:ans+=coef*(s.factorial2(j-1) if j else 1)
 return s.expand(ans)
def pure_arc(F,N):
 coeff=s.Poly(F,t);g=[s.S.One]
 for n in range(1,N+1):
  g.append(s.expand(sum(j*coeff.coeff_monomial(t**j)*g[n-j] for j in range(1,n+1))/n))
 return [expectation(h,x) for h in g]
good=[]
for B in [t,t+t*t,2*t*t-t**3]:
 F=B*x-B*B/2
 values=pure_arc(F,6);assert values==[1]+[0]*6
 good.append({'B':str(B),'checked_positive_indices':6})
bad=[]
for name,F,first in [('cubic',t*x**3,2),('quartic',t*x**4,1),('delayed_cubic',t**3*x**3,6),
 ('rational_time_boundary',t*x**3/(1+t)**2,2)]:
 if name=='rational_time_boundary':F=s.series(F,t,0,7).removeO()
 vals=pure_arc(F,6)
 actual=next(i for i,v in enumerate(vals[1:],1) if v!=0)
 assert actual==first
 bad.append({'name':name,'first_nonzero_in_control':actual,'value':str(vals[actual])})
# Finite source expansion for coefficientwise Gaussian IBP and the actual tangent.
F=t*x**3+t*t*x
coeff=s.Poly(F,t);g=[s.S.One]
for n in range(1,5):
 g.append(s.expand(sum(j*coeff.coeff_monomial(t**j)*g[n-j] for j in range(1,n+1))/n))
ibp=0
for marker in [1,x,x*x,x**3]:
 for n in range(4):
  # [t^n](q_x+Phi_x*q)exp(F), Phi_x=-x+F_x.
  h=(s.diff(marker,x)-x*marker)*g[n]
  h+=sum(s.diff(coeff.coeff_monomial(t**j),x)*marker*g[n-j] for j in range(1,n+1))
  assert expectation(h,x)==0;ibp+=1
# Normalized quadratic residue after the all-degree phase rigidity reduction.
B=p*t+k*c*t*t;C=q*t+k*a*t*t
defect=s.Poly(s.expand(C+B*B/2),t)
assert defect.coeff_monomial(t**4)==k*k*c*c/2
assert defect.coeff_monomial(t**2)==k*a+p*p/2
assert defect.coeff_monomial(t)==q
# Do not replace a radical field argument by equality of ordinary coefficient ideals.
assert s.expand((p*p+2*k*a)**2-4*k*a*(p*p+2*k*a)+4*k*k*a*a)==p**4
# Full shear identity for arbitrary a(Z), represented by three nontrivial choices.
shears=[]
for az in [2+3*Z,1+Z*Z,2+Z+2*Z**3]:
 cc=s.S(3);a0=az.subs(Z,0);u=s.cancel((az-a0)/(cc*Z))
 G=-(u*T+Z*u*u/2)*(az+cc*T)
 P=W*(az+cc*T)+G
 rem=s.expand((u*T-Z*u*u/2)*(a0+cc*T)+G.subs(T,T-Z*u))
 assert rem==0
 transported=s.expand(((W-u*T-Z*u*u/2)*(a0+cc*(T+Z*u)))-P)
 assert transported==0
 shears.append({'a':str(az),'u':str(u),'remainder_zero':True,'factorization_checked':True})
# Wick pure checks of the original sheared members, finite only.
def wick(poly):
 ans=0
 for (zz,ww,tt),coef in s.Poly(s.expand(poly),Z,W,T).terms():
  if zz==ww and tt%2==0:
   ans+=coef*s.factorial(zz)*(s.factorial2(tt-1) if tt else 1)
 return s.expand(ans)
wickchecks=0
for az in [1+Z,2+Z*Z]:
 cc=1;a0=az.subs(Z,0);u=s.cancel((az-a0)/(cc*Z))
 P=(W-u*T-Z*u*u/2)*(az+cc*T)
 for m in range(1,5):
  assert wick(P**m)==0;wickchecks+=1
assert wick(Z*W*(a+c*T))==a
assert wick(Z**2*(W*c*T)**2)==2*c*c
# Ordinary pure ideal on G=p*T: all orders follow from the written exponential.
vals=pure_arc(p*t*x,8)
assert vals[2]==p*p/2
for n in range(1,9):
 if vals[n]!=0:assert s.rem(vals[n],p*p,p)==0
assert all(s.rem(v.subs(p,e),e*e,e)==0 for v in vals[1:])
# Root disk control in exact complex rational arithmetic.
rootsets=[[s.S(-2),s.S(3)],[-1,s.I,2],[-2,-s.I,s.I,3]]
coverchecks=0;inside=0
for roots in rootsets:
 d=len(roots);aa=s.S(3)
 pol=aa*s.prod(x-r for r in roots)
 for xx in [s.Rational(i,2)+s.I*s.Rational(j,2) for i in range(-4,5) for j in range(-2,3)]:
  mod2=s.expand(pol.subs(x,xx)*s.conjugate(pol.subs(x,xx)))
  for R in [1,4]:
   if mod2<=R*R:
    dist2=min(s.expand((xx-r)*s.conjugate(xx-r)) for r in roots)
    assert dist2**d<=s.Rational(R*R,9);inside+=1
   coverchecks+=1
assert inside>0
# Exact monomial inverse-ray scaling; this does not verify arbitrary contour geometry.
scalings=0;h,y=s.symbols('h y',positive=True)
for d in range(2,8):
 assert s.expand(-h**d*(y/h)**d+y**d)==0
 assert s.diff(y/h,y)==1/h;scalings+=1
falsekeys=[
 'contour_length_estimate_certified_by_checker','action_separation_certified_by_checker',
 'rapid_decay_comparison_certified_by_checker','polynomial_arc_universal_proof_certified_by_checker',
 'constant_b_universal_classification_certified_by_checker','general_coefficient_ideal_solver_executed',
 'new_least_pure_prefix_computed','global_complex_safe_locus_uniformity_proved',
 'several_spectator_constant_b_classified','repeated_root_nonconstant_b_classified',
 'all_rational_time_Gaussian_arcs_classified','ordinary_pure_ideal_replaced_by_radical',
 'new_SIC33_answer_computed','new_Delta_minimum_computed','numeric_cubic_frontier_changed',
 'operator_order_reduced','independent_replay','external_review','formal_verification']
receipt={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'python_version':sys.version.split()[0],'sympy_version':s.__version__,
 'translated_Gaussian_arc_controls':good,'nonconstant_arc_controls':bad,'Gaussian_IBP_controls':ibp,
 'normalized_defect_coefficients':{str(j):str(defect.coeff_monomial(t**j)) for j in range(1,5)},
 'shear_factorization_controls':shears,'original_shear_Wick_controls':wickchecks,
 'normalized_sharpness_values':['a','2*c^2'],'nilpotent_line_ordinary_pure_ideal':'(p^2)',
 'root_disk_tests':coverchecks,'inside_root_disk_tests':inside,'inverse_ray_monomial_scalings':scalings,
 'mutation_controls':{'rank_one_not_assumed_in_new_phase_argument':True,
 'ordinary_measure_dx_and_normalization_retained':True,
 'formal_Wick_convergence_not_assumed':True,'moving_tangent_not_frozen':True,
 'field_square_steps_not_applied_to_nilpotents':True,
 'nonlinear_shear_cutoff_degree_loss_retained':True,
 'nonlinear_shear_not_used_for_ordinary_Laplacian_transfer':True},
 **{key:False for key in falsekeys}}
(out/'constant_b_all_degree_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','Gaussian_IBP_controls','root_disk_tests','inside_root_disk_tests','original_shear_Wick_controls']}))

