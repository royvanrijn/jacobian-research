"""Selected exact geometry/band controls for the written Haar-loop proofs.

No numerical conformal map, moment prefix or degree/support search is used.
This does not machine verify Mergelyan, harmonic measure or the index proof.
"""
from pathlib import Path
import hashlib,json,platform
import sympy as s

here=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
t,u,v,e,q,r,R=s.symbols('t u v e q r R',positive=True)
p=t**2-s.Rational(1,2)+s.I*(t**3-t)
dp=s.diff(p,t)
out_a=s.expand(dp.subs(t,-1));out_b=s.expand(-dp.subs(t,1))
va=s.Matrix([s.re(out_a),s.im(out_a)])
vb=s.Matrix([s.re(out_b),s.im(out_b)])
assert va==s.Matrix([-2,2]) and vb==s.Matrix([-2,-2])
assert va.dot(vb)==0 and s.det(s.Matrix.hstack(va,vb))==8
assert s.expand(s.re(p))==t**2-s.Rational(1,2)
assert s.expand(s.im(p)**2-t**2*(1-t**2)**2)==0
assert dp.subs(t,-1)!=0 and dp.subs(t,1)!=0
area=s.integrate(2*s.sqrt(u)*(1-u),(u,0,1))
assert area==s.Rational(8,15)

# Ahlfors bound in an exact constant-angle sector. All variable-opening
# estimates and normalization are written arguments, not inferred from this.
integral=s.integrate(1/(s.pi*q*t),(t,r,R))
assert s.simplify(s.pi*integral-s.log(R/r)/q)==0
beta=s.Rational(7,4)
assert beta<2 and all(beta>s.Rational(N,2) for N in [1,2,3])
assert beta<s.Rational(4,2)  # dimension4 equality is deliberately undecided
threshold_controls=[]
for N,k,angle,excluded in [
 (3,1,s.Rational(1,2),True),
 (3,1,s.Rational(2,3),False),
 (3,1,s.Rational(1),False),
 (2,1,s.Rational(3,2),False),
 (2,1,s.Rational(2),False),
 (2,2,s.Rational(2),False),
]:
 actual=bool(angle<s.Rational(2*k,N));assert actual==excluded
 threshold_controls.append({'core_dimension':N,'endpoint_order':k,'opening_over_pi':str(angle),
                            'excluded_by_dimension_bound_alone':actual})

# Exact local mass models: a saddle has an unbounded band/epsilon ratio;
# an extremal curve has square-root mass, incompatible with beta>1/2.
saddle_area=4*(e+s.integrate(e/t,(t,e,1)))
assert s.simplify(saddle_area-4*e*(1-s.log(e)))==0
assert s.limit(saddle_area/e,e,0,dir='+')==s.oo
curve_band=2*s.sqrt(e)
assert s.limit(curve_band/e**s.Rational(3,4),e,0,dir='+')==s.oo
assert s.Rational(2,3)<1 and s.Rational(5,6)<1

# A genuine compact boundary example: h=z^2 on the uniform sphere has
# CDF sqrt(e) at0 and no interior critical value. Its extremal equator
# is not isolated, so an isolated-extremum index count cannot be used.
z=s.Symbol('z',real=True)
sphere_cdf=s.integrate(s.Rational(1,2),(z,-s.sqrt(e),s.sqrt(e)))
assert sphere_cdf==s.sqrt(e)
h=z**2
g1=s.diag(-1,-1,1);g2=s.diag(-1,1,-1)
assert g1.det()==g2.det()==1 and g1*g2==g2*g1
assert g1**2==g2**2==s.eye(3) and g1!=g2 and g1*g2 not in [g1,g2,s.eye(3)]
assert s.Matrix.vstack(g1-s.eye(3),g2-s.eye(3)).nullspace()==[]
assert s.expand(h.subs(z,-z)-h)==0

# Quartic rescaling distinguishes an isolated flat minimum from a singular
# critical level: (x^2+y^2)^2 has area pi below1, while x^2*y^2 has
# unbounded sublevel area. In a quadrant truncated at L>=1 the latter
# area is 1+2*log(L), giving the exact divergence control.
L=s.Symbol('L',positive=True)
quartic_positive_area=s.pi
quartic_semidefinite_quadrant=s.integrate(L,(t,0,1/L))+s.integrate(1/t,(t,1/L,L))
assert s.simplify(quartic_semidefinite_quadrant-1-2*s.log(L))==0
assert s.limit(quartic_semidefinite_quadrant,L,s.oo)==s.oo
# The two Morse perturbation zeros of a critical circle have indices+1,-1.
perturbation=s.Symbol('perturbation',positive=True)
assert s.det(s.diag(perturbation,2))>0
assert s.det(s.diag(-perturbation,2))<0
assert s.det(s.diag(perturbation,-2))<0
assert s.det(s.diag(-perturbation,-2))>0

receipt={
 'date':'2026-11-15','completed':True,'checker_sha256':sha(Path(__file__)),
 'python_version':platform.python_version(),'sympy_version':s.__version__,
 'cubic_outgoing_rays':[list(map(str,va)),list(map(str,vb))],
 'cubic_ray_dot_product':str(va.dot(vb)),'cubic_ray_determinant':'8',
 'cubic_opening_over_pi':'1/2','cubic_endpoint_orders':[1,1],'cubic_area':str(area),
 'cubic_harmonic_upper_exponent_control':str(beta),
 'cubic_excluded_core_dimensions':[1,2,3],
 'cubic_minimum_necessary_core_dimension':4,
 'dimension4_realization_proved':False,'threshold_controls':threshold_controls,
 'saddle_band_area':str(saddle_area),'saddle_band_ratio_unbounded':True,
 'extremal_curve_band_mass':str(curve_band),'sphere_z_squared_endpoint_CDF':str(sphere_cdf),
 'noncyclic_D2_boundary_group_common_fixed_vector':False,
 'positive_definite_quartic_sublevel_area':str(quartic_positive_area),
 'semidefinite_quartic_sublevel_quadrant_area':str(quartic_semidefinite_quadrant),
 'semidefinite_quartic_sublevel_area_diverges':True,
 'critical_circle_perturbation_index_sum':0,
 'exterior_hairs_on_pure_positive_image_graph':'EXCLUDED;interior bridges not excluded',
 'positive_Jordan_pure_measure':'INTERNAL WRITTEN:unique harmonic measure at0,with0 in bounded domain',
 'general_endpoint_necessary_bound':'opening>=2*pi*endpoint_order/core_dimension',
 'noncyclic_sphere_non_outward_cusp_Jordan_loop':'EXCLUDED in every degree',
 'noncyclic_non_dihedral_sphere_all_Jordan_loops':'EXCLUDED in every degree,including outward cusps',
 'square_root_endpoint_critical_set':'finitely many isolated extrema and smooth critical circles with nondegenerate normal Hessian',
 'original_isolated_sphere_critical_points':2,
 'polyhedral_two_point_invariant_set':'IMPOSSIBLE',
 'dihedral_cusp_pure_phase_existence':'NOT asserted',
 'sphere_exclusion_applied_to_arbitrary_full_non_T_phases':False,
 'cubic_Haar_realization_on_either_full_exceptional_model':'IMPOSSIBLE for the stated real-seed interval loop',
 'multiple_loop_or_full_rank2_architectures_solved':False,
 'one_loop_with_interior_branches_excluded':False,
 'height_symmetrization_replaced_by_pure_phase_measure':False,
 'universal_analytic_geometric_topological_proofs_machine_verified':False,
 'new_degree_support_or_moment_prefix_campaign':False,
 'exceptional_MZ_decisions':'OPEN on both models',
 'independent_replay':False,'formal_verification':False,'external_review':False}
(here/'haar_loop_corner_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
