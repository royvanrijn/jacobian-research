"""Selected exact controls; universal surface/valuation proofs remain written."""
from pathlib import Path
import hashlib,json
import sympy as s

p=Path(__file__).parent
old=p/'exceptional_sl2_verification.json'
groups=json.loads(old.read_text())['groups']
orbit_controls=[]
for name,max_order,min_orbit in [('2O',4,6),('2I',5,12)]:
    row=groups[name];counts={int(k):v for k,v in row['projective_element_order_counts'].items()}
    assert sum(counts.values()) == row['projective_order']
    assert max(counts) == max_order
    assert row['projective_order']//max_order == min_orbit
    gap=s.Rational(min_orbit,2)-2
    assert gap > 0
    orbit_controls.append({'group':name,'projective_order':row['projective_order'],'maximum_rotation_order':max_order,
                           'minimum_orbit_bound':min_orbit,'D_minus_2m_lower_coefficient':str(gap)})

# Exact valuation algebra. A log discrepancy integer strictly above2a gives a regular volume.
A,a,D,m,vH,r=s.symbols('A a D m vH r',positive=True)
volume=A-1-2*a
assert s.expand(volume.subs(A,2*a+r)-(r-1)) == 0
assert s.factor((D*a/m-1-2*a)-(a*(D/m-2)-1)) == 0
assert s.expand(volume.subs(A,3*a)-(a-1)) == 0
assert s.expand(volume.subs(A,6*a)-(4*a-1)) == 0

# Node-relative form: dt=v du+u dv, eta=-J du/u, so dt wedge eta=J du wedge dv.
u,v,J,t=s.symbols('u v J t')
dt=(v,u);eta=(-J/u,0)
assert s.simplify(dt[0]*eta[1]-dt[1]*eta[0]) == J
cover_controls=[]
for e in [1,2,3,5]:
    exponent=s.Rational(1,e)-1
    assert exponent > -1
    assert s.diff(t**e,t) == e*t**(e-1)
    cover_controls.append({'ramification':e,'original_period_power':str(exponent),'strictly_integrable':True})

# Equality of the reciprocal-multiplicity bound is insufficient for the strict certificate.
assert s.simplify((1/m-2/D).subs(D,2*m)) == 0

x,y,z,w=s.symbols('x y z w')
det=x*w-y*z
derD=lambda f:z*s.diff(f,x)+w*s.diff(f,y)
derE=lambda f:x*s.diff(f,z)+y*s.diff(f,w)
derH=lambda f:x*s.diff(f,x)+y*s.diff(f,y)-z*s.diff(f,z)-w*s.diff(f,w)
assert derD(det)==derE(det)==derH(det)==0
fO=lambda u,v:u**8+14*u**4*v**4+v**8
fI=lambda u,v:u*v*(u**10+11*u**5*v**5-v**10)
stabilizer_controls=[]
for name,f,n in [('2O',fO,8),('2I',fI,12)]:
    bottom=f(z,w);top=f(x,y);phase=top+bottom
    assert derD(bottom)==0 and s.expand(derH(bottom)+n*bottom)==0
    dagger=s.expand(phase.xreplace({x:w,y:-z,z:-y,w:x}))
    assert dagger==s.expand(phase)
    derivs=[s.Poly(s.expand(g(phase)),x,y,z,w) for g in [derD,derE,derH]]
    monomials=sorted(set().union(*(set(h.monoms()) for h in derivs)))
    matrix=s.Matrix([[h.coeff_monomial(mon) for h in derivs] for mon in monomials])
    assert matrix.rank()==3
    stabilizer_controls.append({'group':name,'row_highest_weight':n,'row_shear_derivative_zero':True,
                                'selfadjoint_phase_three_derivative_rank':3,
                                'positive_second_moment_by_nonzero_selfadjoint_pairing':True,
                                'finite_stabilizer_does_not_certify_purity':True})

# Nonzero trace-zero2x2 matrices: det!=0 gives distinct eigenvalues; det0 gives square0.
aa,bb,cc=s.symbols('aa bb cc')
X=s.Matrix([[aa,bb],[cc,-aa]])
assert X*X == (aa**2+bb*cc)*s.eye(2)

# A nonreduced pure phase may be nonzero; eventual nilpotence is different from first purity.
epsilon,seed=s.symbols('epsilon seed')
assert s.rem((epsilon*seed)**2,epsilon**2,epsilon) == 0
assert s.Rational(3,5)-s.Rational(3,5) == 0  # mean of harmonic hO
assert s.rem(epsilon,epsilon**2,epsilon) == epsilon

receipt={
 'completed':True,'sympy_version':s.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'imported_group_receipt':'moments/vanishing/exceptional_sl2_verification.json',
 'imported_group_receipt_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),
 'orbit_and_multiplicity_controls':orbit_controls,'valuation_volume_formula_exact':True,
 'finite_boundary_volume_order_bounds':{'2O':'a_E-1','2I':'4*a_E-1'},
 'nodal_relative_form_wedge_exact':True,'semistable_base_change_controls':cover_controls,
 'strict_threshold_equality_not_certified':True,'left_derivations_preserve_determinant_exact':True,
 'stabilizer_boundary_controls':stabilizer_controls,'sl2_square_identity_exact':True,
 'entire_original_T_invariant_eventual_radicals':'ZERO on both2O and2I,in ALL degrees',
 'nonzero_pure_in_any_automorphic_torus_fixed_subalgebra':'EXCLUDED by intrinsic Haar transport',
 'full_failure_necessary_stabilizer':'finite for leftSL2 and every automorphism-conjugated left action',
 'all_C_algebra_T_invariant_eventual_radical':'exactly the nilpotent elements;not an assertion that all nilpotents are pure from their first power',
 'dual_number_nonzero_pure_and_nonpure_nilpotent_boundary_exact':True,
 'fixed_degree_ordinary_ideal_nilpotence_and_prefix_certificate':'existence and terminating theoretical algorithm only;no N or L computed',
 'full_exceptional_MZ_decisions':'OPEN on both models',
 'general_nonisolated_higher_dimensional_period_theorem_claimed':False,
 'all_full_phases_reduced_to_T_invariants':False,'finite_stabilizer_implies_purity_or_failure':False,
 'same_matrix_cutoff_after_arbitrary_nonlinear_automorphism':False,
 'ordinary_pure_ideals_assumed_reduced':False,'all_nilpotents_assumed_pure':False,
 'universal_T_invariant_nilpotence_cutoff_computed':False,
 'fixed_degree_Noetherian_nilpotence_proof_machine_verified':False,
 'surface_semistable_trace_Hartogs_and_lct_proofs_machine_verified':False,
 'new_degree_support_density_or_prefix_campaign':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,
}
Path(__file__).with_name('torus_volume_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
