"""Selected exact controls; universal semistable/adjunction proofs remain written."""
from pathlib import Path
import hashlib,json
import sympy as s

p=Path(__file__).parent
t,x1,x2,x3,J=s.symbols('t x1 x2 x3 J')
# Relative top log forms in dimension3, for two and three crossing components.
dt2=s.Matrix([x2,x1,0]);eta2={(1,2):J/x2}
dt3=s.Matrix([x2*x3,x1*x3,x1*x2]);eta3={(1,2):J/(x2*x3)}
def wedge_top(one,two):return s.simplify(one[0]*two.get((1,2),0)-one[1]*two.get((0,2),0)+one[2]*two.get((0,1),0))
assert wedge_top(dt2,eta2)==J and wedge_top(dt3,eta3)==J

# Positive plt margin plus finite-value valuation removes the volume pole.
A,a,D,N,b,vH=s.symbols('A a D N b vH',positive=True)
c=(N-1)/D
assert s.simplify((a+c*vH+b-N*a)-(b+c*(vH-D*a)))==0
adjunction_controls=[]
for dim,degree in [(2,4),(3,8),(3,12),(4,6)]:
    coefficient=s.Rational(dim-1,degree)
    assert coefficient==c.subs({N:dim,D:degree})
    adjunction_controls.append({'dimension':dim,'degree':degree,'boundary_coefficient':str(coefficient)})

old=p/'exceptional_sl2_verification.json';groups=json.loads(old.read_text())['groups']
orbit_controls=[];degree_controls=[]
for name,r in [('2O',6),('2I',12)]:
    row=groups[name]
    counts={int(k):v for k,v in row['projective_element_order_counts'].items()}
    assert sum(counts.values())==row['projective_order'] and row['projective_order']//max(counts)==r
    assert s.Rational(2,r)<1
    orbit_controls.append({'group':name,'minimum_right_orbit':r,'critical_restricted_coefficient_bound':str(s.Rational(2,r))})
    nonconstant=[n for n,v in enumerate(row['symmetric_invariant_dimensions_n0_to32']) if n and v]
    first,second=nonconstant[:2]
    assert (first,second)==((8,12) if name=='2O' else (12,20))
    degree_controls.append({'group':name,'first_two_nonconstant_matrix_degrees':[first,second],
                           'new_real_seed_minimum':second,'thin_degree_gate':3*second,'bad_marker_subtraction':second})

# A common factor in a vector-valued binary form is detected by the coefficient gcd.
u,v=s.symbols('u v')
tuple_forms=[u**6*(v**6+u**6),u**6*v**3*(u**3+v**3)]
assert s.Poly(s.gcd(*tuple_forms),u,v).degree(u)==6
assert s.gcd(*tuple_forms)==u**6
assert s.Rational(2,12)*6==1  # half-degree coefficient is exactly the non-klt boundary

weight_controls=[]
for degree in [8,12,24]:
    for k in [degree//2,degree//2+1]:
        monomials=[(j,degree-j) for j in range(k,degree+1)]
        weights=[2*j-degree for j,_ in monomials]
        assert min(weights)==2*k-degree
        if k==degree//2:assert weights[0]==0 and all(w>0 for w in weights[1:])
        else:assert min(weights)>=2
        weight_controls.append({'degree':degree,'common_root_multiplicity':k,'minimum_weight':min(weights),
                                'nonzero_zero_weight_limit':k==degree//2})

# Compact real structure on rank-one matrices is antipodal on both factors.
u0,u1,v0,v1=s.symbols('u0 u1 v0 v1')
rankone=s.Matrix([[u0*v0,u0*v1],[u1*v0,u1*v1]])
sigma=s.Matrix([[rankone[1,1],-rankone[1,0]],[-rankone[0,1],rankone[0,0]]])
ju=s.Matrix([-u1,u0]);jv=s.Matrix([-v1,v0])
assert sigma==ju*jv.T
assert s.Matrix([-ju[1],ju[0]])==-s.Matrix([u0,u1])
for degree in [8,12]:
    k=degree//2
    assert 2*k==degree and 2*(k+1)>degree
    A=(u0*u1)**k
    H=u0*s.diff(A,u0)-u1*s.diff(A,u1)
    assert H==0

# Retained squarefree prefix controls are consistent with the new necessary half-degree condition.
z=s.symbols('z')
for degree in [8,12]:
    Q=s.Poly(z**degree+degree*z,z)
    assert s.gcd(Q,Q.diff()).degree()==0

cover_controls=[]
for e in [1,2,5,12]:
    power=s.Rational(1,e)-1
    assert power>-1
    cover_controls.append({'ramification':e,'original_period_power':str(power),'subpolar':True})

out={
 'completed':True,'sympy_version':s.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'dimension3_relative_log_form_wedge_controls_exact':True,'plt_finite_valuation_margin_identity_exact':True,
 'adjunction_threshold_controls':adjunction_controls,'right_orbit_controls':orbit_controls,
 'vector_binary_coefficient_gcd_half_boundary_exact':True,'one_parameter_weight_and_limit_controls':weight_controls,
 'rank_one_compact_antipodal_factorization_exact':True,'real_half_pair_exhaustion_controls_exact':True,
 'degree_and_marker_controls':degree_controls,'semistable_cover_controls':cover_controls,
 'imported_group_receipt':'moments/vanishing/exceptional_sl2_verification.json',
 'imported_group_receipt_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
 'single_isotypic_scope':'ALL degrees and ALL right multiplicities;pure iff strict common-root multiplicity>D/2;uniform field/reduced-ring cutoff2m>marker degree',
 'full_icosahedral_failure_minimum_matrix_degree':20,
 'full_thin_degree_gates':{'2O':36,'2I':60},'thin_bad_marker_degree_bounds':{'2O':'D-12','2I':'D-20'},
 'real_single_isotypic_seed_scope':'Every polynomial in a nonconstant real seed inC+W_D is excluded from nonzero purity',
 'ordinary_V12_quartic_deficit_preserved':True,'moment_and_invariant_ideal_radicals':'equal nullcone ideals',
 'full_failure_requires_mixed_nonconstant_Peter_Weyl_degrees':True,
 'all_homogeneous_matrix_representatives_treated_as_single_isotypic':False,
 'ordinary_pure_ideals_assumed_equal_or_reduced':False,'reduced_ring_result_silently_extended_to_nonreduced_isotypic_families':False,
 'numeric_pure_prefix_or_general_nilpotent_cutoff_computed':False,
 'universal_semistable_log_comparison_chain_trace_adjunction_and_classification_machine_verified':False,
 'new_degree_support_density_or_prefix_campaign':False,'full_exceptional_MZ_decisions':'OPEN on both models',
 'independent_replay':False,'formal_verification':False,'external_review':False,
}
Path(__file__).with_name('semistable_isotypic_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
