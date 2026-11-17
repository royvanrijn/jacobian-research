"""Selected exact controls for two written all-degree thin-phase notes.

No graph enumeration, moment prefix, support search or numerical Riemann map.
Universal topology, harmonic sweeping and degree assertions remain written.
"""
from pathlib import Path
import hashlib,json,platform
import sympy as s
here=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
t,v,u,beta=s.symbols('t v u beta',real=True)
x,y,z,w=s.symbols('x y z w')
G=s.groebner([x*w-y*z-1],x,y,z,w,order='grlex')
normal=lambda f:G.reduce(s.expand(f))[1]
degree=lambda f:int(s.Poly(normal(f),x,y,z,w).total_degree())
assert normal(x*w-y*z)==1 and degree(x*w-y*z)==0
f=x*w+y*z+x;g=(x+w)**2+y*z
assert degree(f*g)==degree(f)+degree(g)==4
h=(x+w)**2+(y+z)**2
assert degree(h)==2 and degree(h**3-2*h+7)==6
assert degree(h**2)==4
degree_gates=[]
for name,minh in [('2O',8),('2I',12)]:
 for branch,d in [('full',3),('original_T_invariant',4)]:
  D=d*minh
  assert (d-1)*minh==D-minh
  degree_gates.append({'model':name,'branch':branch,'minimum_normalized_degree':d,
                       'necessary_matrix_degree':D,'endpoint_seed_degree':minh,
                       'bad_marker_degree_at_endpoint_at_most':D-minh})

# Canonical cubic topology. Plane-affine coordinates are used only for topology.
cx=t**2;cy=t**3-beta*t
assert s.cancel((cx-cx.subs(t,v))/(t-v))==t+v
assert s.cancel((cy-cy.subs(t,v))/(t-v))==t**2+t*v+v**2-beta
assert s.expand((t**2+t*v+v**2-beta).subs(v,-t))==t**2-beta
a=s.Symbol('a',positive=True)
assert s.expand((1+a)**2-1)==a**2+2*a
assert s.gcd(2*t,3*t**2-1)==1

# The irreducible quartic loop has an inward continuation, unlike the cubic.
px=t**4-s.Rational(3,2)*t**2+s.Rational(1,4);py=t-t**3
assert (px+s.I*py).subs(t,-1)==(px+s.I*py).subs(t,1)==-s.Rational(1,4)
assert px.subs(t,0)==s.Rational(1,4) and py.subs(t,0)==0
assert s.gcd(s.diff(px,t),s.diff(py,t))==1
assert s.factor((u**2-s.Rational(3,2)*u)-(v**2-s.Rational(3,2)*v))==(u-v)*(2*u+2*v-3)/2
diff=s.cancel((u*(1-u)**2-v*(1-v)**2)/(u-v))
assert s.expand(diff-(u**2+u*v+v**2-2*u-2*v+1))==0
# On u+v=3/2, equality of squared imaginary coordinates requires uv=1/4.
assert s.expand(diff.subs(v,s.Rational(3,2)-u))==u**2-s.Rational(3,2)*u+s.Rational(1,4)
assert (3+s.sqrt(5))/4>1
raya=[s.diff(px,t).subs(t,-1),s.diff(py,t).subs(t,-1)]
rayb=[-s.diff(px,t).subs(t,1),-s.diff(py,t).subs(t,1)]
assert raya==[-1,-2] and rayb==[-1,2]
inner=[-raya[0],-raya[1]];assert inner==[1,2]
area=s.integrate((px*s.diff(py,t)-py*s.diff(px,t))/2,(t,-1,1))
assert area==s.Rational(12,35)
assert s.gcd(s.degree(px,t),s.degree(py,t))==1

# Herglotz lobe separator with shared boundary points+1,-1.
q=s.Symbol('q')
T=(1+q)/(1-q)+(-1+q)/(-1-q)
F=s.cancel(1/(1+T))
assert s.cancel(F-(1-q**2)/(3+q**2))==0
assert F.subs(q,1)==F.subs(q,-1)==0 and F.subs(q,0)==s.Rational(1,3)

# Exact positive harmonic kernel identity underlying the disk cap bound.
aa,bb,rr=s.symbols('aa bb rr',real=True)
ww=aa+s.I*bb;lam=1-rr
J=s.re((1+lam*ww)/(1-lam*ww)).expand(complex=True)
kernel=(1-lam**2*(aa**2+bb**2))/((1-lam*aa)**2+(lam*bb)**2)
assert s.simplify(J-kernel)==0
assert s.expand(1-(1-rr)**2-rr-rr*(1-rr))==0
# Squaring the triangle estimate gives denominator<=4r^2 on a radius-r cap.
assert s.expand((rr+rr)**2)==4*rr**2

receipt={
 'date':'2026-11-15','completed':True,'checker_sha256':sha(Path(__file__)),
 'python_version':platform.python_version(),'sympy_version':s.__version__,
 'filtered_degree_controls':{'determinant_class_degree':0,'product_degree':degree(f*g),
                            'real_seed_degree':degree(h),'cubic_composition_degree':degree(h**3-2*h+7)},
 'original_bad_marker_bounds':{'2O':'D-8','2I':'D-12'},'thin_phase_degree_gates':degree_gates,
 'cubic_loop_only_double_parameters':'plus/minus sqrt(beta),beta>0',
 'cubic_tail_topology':'all tails outside the node interval are exterior hairs',
 'plane_affine_changes_claimed_to_preserve_analytic_moments_or_angles':False,
 'quartic_loop_polynomial':str(px+s.I*py),'quartic_join_outgoing_rays':[list(map(str,raya)),list(map(str,rayb))],
 'quartic_inward_continuation_ray':list(map(str,inner)),
 'quartic_loop_area':str(area),'quartic_normalization_coordinate_degree_gcd':1,
 'quartic_interior_measure_swept_density_bound':'<=2*delta<1/2 for a sufficiently small delta<1/4',
 'quartic_abstract_interval_purity':'INTERNAL WRITTEN:exact harmonic sweeping balance in every positive order',
 'quartic_abstract_fixed_bad_marker_list':['t','t^2','t^3'],
 'quartic_which_bad_marker_identified':False,'quartic_actual_Haar_realization_proved':False,
 'pure_graph_hull':'INTERNAL WRITTEN:one closed Jordan domain containing0 in its interior',
 'inner_graph_mass_automatically_zero':False,
 'exact_sweeping_equation':'omega_0=nu_boundary+integral omega_z dnu_interior',
 'disk_boundary_cap_bound':'mass(B(xi,r))<=4r',
 'regular_outer_arc_seed_constraint_N_le_3':'no critical polynomial parameter can map there',
 'universal_topological_analytic_and_filtered_proofs_machine_verified':False,
 'new_degree_support_graph_or_moment_prefix_campaign':False,
 'original_matrix_bounds_apply_to_full_thin_phases':True,
 'torus_invariant32_48_gate_applied_to_full_non_T_phases':False,
 'either_exceptional_MZ_decision':'OPEN',
 'independent_replay':False,'formal_verification':False,'external_review':False}
(here/'graph_hull_degree_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
