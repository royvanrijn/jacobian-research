"""Exact algebraic controls; proper/Milnor/period arguments remain written proofs."""
from pathlib import Path
import hashlib,json
import sympy as s

a,b,c,z,w,u,t,r = s.symbols('a b c z w u t r')
q=a*a+b*b+c*c
cartan={a:(z*z-w*w)/2,b:(z*z+w*w)/(2*s.I),c:z*w}
assert s.expand(q.subs(cartan)) == 0
fO=z**8+14*z**4*w**4+w**8
gO=z*z*w*w*(z**4-w**4)**2
fI=z*w*(z**10+11*z**5*w**5-w**10)
gI=z**20-228*z**15*w**5+494*z**10*w**10+228*z**5*w**15+w**20
delta=lambda p:sum(s.diff(p,x,2) for x in (a,b,c))
hO=a**4+b**4+c**4-s.Rational(3,5)*q*q
R=2*a**5*c-20*a**3*b*b*c+10*a*b**4*c+11*c**6
hI=R-q*delta(R)/22+q*q*delta(delta(R))/(22*36)-q**3*delta(delta(delta(R)))/(22*36*42)
assert s.expand(delta(hO)) == 0 and s.expand(delta(hI)) == 0
assert s.expand(hO.subs(cartan)-fO/8) == 0
assert s.expand(hI.subs(cartan)-fI) == 0
assert all(v.is_real for v in s.Poly(s.expand(hI),a,b,c).coeffs())

squarefree=[];coprime=[]
for name,f,g,n in [('2O',fO,gO,8),('2I',fI,gI,12)]:
    F=s.Poly(f.subs(w,1),z);G=s.Poly(g.subs(w,1),z)
    assert s.gcd(F,F.diff()).degree() == 0
    infinity_order=n-F.degree()
    assert infinity_order in (0,1)
    assert s.gcd(F,G).degree() == 0
    assert not (f.subs({z:1,w:0}) == 0 and g.subs({z:1,w:0}) == 0)
    squarefree.append({'group':name,'binary_degree':n,'affine_derivative_gcd_degree':0,'infinity_zero_order':infinity_order})
    coprime.append({'group':name,'affine_pair_gcd_degree':0,'no_common_infinity_zero':True})

# Local graph Jacobian: v=t*u^D turns du wedge dv into u^D du wedge dt.
D,N=s.symbols('D N',integer=True,positive=True)
v=t*u**D
assert s.det(s.Matrix([[s.diff(u,u),s.diff(u,t)],[s.diff(v,u),s.diff(v,t)]])) == u**D
assert s.simplify(u**(-N)*u**D-u**(D-N)) == 0
threshold=[{'N':n,'D':d,'volume_order':d-n,'holomorphic':d>=n} for n,d in [(2,4),(2,6),(3,3),(3,24),(3,60),(3,2)]]
assert all(x['holomorphic'] for x in threshold[:-1]) and not threshold[-1]['holomorphic']

# An explicit projective chart of Q2 verifies its boundary volume pole order2.
AA,BB=s.symbols('AA BB')
# On AA^2+BB^2+1=u^2, AA_u=u/AA, AA_BB=-BB/AA.
chart_numerator=u*u/AA-BB*BB/AA-AA
assert s.expand(AA*chart_numerator-1).subs(AA**2,u*u-BB*BB-1) == 0

# The full-rank-one boundary is smooth and the leading invariant pencil has equal degrees.
x,y,zz,ww=s.symbols('x y zz ww')
det=x*ww-y*zz
assert [s.diff(det,j) for j in (x,y,zz,ww)] == [ww,-zz,-y,x]
assert 3*8 == 2*12 == 24 and 5*12 == 3*20 == 60
ratio_derivatives=[]
for name,f,g,A,B in [('2O',fO,gO,3,2),('2I',fI,gI,5,3)]:
    numerator=s.factor(A*s.diff(f,z)*g-B*f*s.diff(g,z))
    assert numerator != 0
    ratio_derivatives.append({'group':name,'right_ratio_nonconstant_exact':True})

# Original matrix leading forms survive: q=(xw-yz)^2 and the Cartan restrictions are nonzero.
hopf={a:x*zz-y*ww,b:s.I*(x*zz+y*ww),c:x*ww+y*zz}
assert s.expand(q.subs(hopf)-det**2) == 0
matrix_seed=[]
for name,h,d in [('2O',hO,4),('2I',hI,6)]:
    lifted=s.expand(h.subs(hopf))
    assert s.Poly(lifted,x,y,zz,ww).total_degree() == 2*d
    assert lifted.subs({x:1,y:2,zz:3,ww:6}) != 0
    matrix_seed.append({'group':name,'sphere_degree':d,'matrix_degree':2*d,'rank_one_leading_nonzero':True})

resolution_controls=[]
for multiplicities,zeros in [([2,3],[0,0]),([5,7,11],[0,1,2]),([1,1,1],[0,0,0])]:
    exponents=[s.Rational(k+1,m)-1 for k,m in zip(zeros,multiplicities)]
    assert all(alpha > -1 for alpha in exponents)
    resolution_controls.append({'multiplicities':multiplicities,'Jacobian_zero_orders':zeros,'powers':[str(x) for x in exponents],'strictly_above_minus_one':True})

receipt={
 'completed':True,'sympy_version':s.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'harmonic_seed_lifts_and_Cartan_restrictions_exact':True,'projective_squarefree_controls':squarefree,
 'binary_invariant_pair_projective_coprimality':coprime,'graph_Jacobian_and_volume_exponent_exact':True,
 'quadric_volume_threshold_controls':threshold,'sphere_chart_boundary_pole_order':2,
 'rank_one_projective_quadric_smoothness_gradient_exact':True,'invariant_right_ratio_controls':ratio_derivatives,
 'minimal_seed_matrix_degree_controls':matrix_seed,'resolution_exponent_controls':resolution_controls,
 'full_exceptional_gate':'Every nonzero pure original phase has singular leading zero divisor on the projective rank-one quadric',
 'nonempty_smooth_leading_invariant_degree_strata':'2O:24k;2I:60k;all k>=1;written incidence proof,not a coefficient search',
 'minimal_real_T_invariant_seed_families':'Only0 pure in C[h_O] and C[h_I] in every phase degree',
 'singular_leading_implies_purity_or_failure':False,'arbitrary_affine_maps_assumed_to_have_proper_holomorphic_models':False,
 'all_full_minimal_matrix_degree_phases_closed':False,'full_nonT_or_image_interior_problem_removed':False,
 'properness_Milnor_topology_resolution_and_period_estimates_machine_verified':False,
 'new_degree_support_density_or_prefix_campaign':False,'exceptional_MZ_decisions':'OPEN on both models',
 'independent_replay':False,'formal_verification':False,'external_review':False,
}
Path(__file__).with_name('proper_quadric_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
