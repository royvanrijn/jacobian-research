"""Selected exact controls; universal residue and cycle proofs are written."""
from pathlib import Path
import hashlib,json
import sympy as s

u,v,z=s.symbols('u v z');a,b,c,d,e=s.symbols('a b c d e')
form=a/(u**2*v)+b/(u*v**2)+c/(u*v)+d/u**3+e/(u*v**3)
res=lambda h:s.residue(s.residue(h,u,0),v,0)
assert res(form)==c
beta_v=a/(u*v)+b/(u**2*v**3)
beta_u=c/(u*v)+d/(u**3*v**2)
assert res(s.diff(beta_v,u)-s.diff(beta_u,v))==0

# Graph cycles pair with graph cohomology independently of representatives.
inc=s.Matrix([[-1]*4,[1]*4])
R=s.Matrix([[1,0,0],[0,1,0],[0,0,1],[-1,-1,-1]])
assert inc*R==s.zeros(2,3) and R.rank()==3
alpha=s.Matrix([a,b]);direction=s.Matrix([1,2,3,4])
assert ((direction+inc.T*alpha).T*R-direction.T*R)==s.zeros(1,3)

# A fully cyclic algebraic unit and a noncyclic unit.
def wronskian(functions,derivative):
    rows=[list(functions)]
    for _ in range(len(functions)-1):rows.append([derivative(x) for x in rows[-1]])
    return s.simplify(s.Matrix(rows).det())
assert wronskian([1,z,z**2],lambda h:s.diff(h,z))==2
assert wronskian([1,z,2*z],lambda h:s.diff(h,z))==0
# The derivation on C(v), z=v^2, must be d/dz=(2v)^-1 d/dv.
derivative=lambda h:s.diff(h,v)/(2*v)
algebraic_wronskian=wronskian([v,v**3,v**5],derivative)
assert algebraic_wronskian!=0

# Genuine SL2 composite: two torus fibres, both ORIGINAL markers retained.
unit=s.Matrix([[-1/(2*v),1/(2*v)]])
marked=s.Matrix([[-s.Rational(1,2),-s.Rational(1,2)]])
assert wronskian(list(unit),derivative)==0
marker_matrix=unit.col_join(marked)
assert s.simplify(marker_matrix.det())==1/(2*v)
swap=s.Matrix([[0,1],[1,0]]);projection=(s.eye(2)+swap)/2
actual=s.Matrix([-1,0]);projected=projection*actual
assert projected==s.Matrix([-s.Rational(1,2),-s.Rational(1,2)])
assert unit*projected==s.zeros(1,1)
assert marked*projected==s.Matrix([[s.Rational(1,2)]])
assert s.simplify(1/(v*v-u*u)-(1/(v-u)+1/(v+u))/(2*v))==0
assert s.simplify(s.diff(-s.log(z-u*u)/2,u)-u/(z-u*u))==0
analytic_part=s.series(s.log((1+v)/(1-v))/(2*v),v,0,7).removeO().expand()
assert analytic_part==1+v**2/3+v**4/5+v**6/7
assert s.simplify(s.diff((s.log(z)-s.log(z-1))/2,z)+1/(2*z*(z-1)))==0

# All-degree safe-family unit vectors are constant multiples of 1/z.
coef=s.symbols('c0:3');family_unit=[x/z for x in coef]
assert wronskian(family_unit,lambda h:s.diff(h,z))==0
q=s.symbols('q',positive=True,integer=True)
for order,exponent in [(24,12),(60,30)]:
    assert s.cancel(2*(q*exponent)/order-2)==q-2

out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'higher_pole_double_residue_and_exact_form_controls_exact':True,
 'graph_pairing_coboundary_independence_exact':True,
 'unit_Wronskian_cyclicity_and_dependence_controls_exact':True,
 'algebraic_cover_Wronskian':str(algebraic_wronskian),
 'actual_SL2_square_composite_marker_matrix_determinant':'1/(2*v)',
 'actual_SL2_projected_cycle':[str(x) for x in projected],
 'actual_SL2_unit_integer_log_coefficient':'0',
 'actual_SL2_fixed_marker_integer_log_coefficient':'1/2',
 'actual_SL2_nonpure_unit_moments':'1/(2*m+1)',
 'actual_SL2_fixed_marker_moments':'1/(2*m+2)',
 'square_composite_half_power_analytic_part':str(analytic_part),
 'safe_family_unit_Wronskian_zero_exact':True,
 'safe_family_unit_invisible_dimension_lower_bound_in_k_over_e':'q - 2',
 'universal_finite_field_residue_and_marker_basis_proofs_machine_verified':False,
 'actual_exceptional_graph_coupling_computed_for_every_phase':False,
 'unit_cyclicity_universal_on_entire_exceptional_pure_locus':False,
 'primitive_unit_cyclicity_automatically_preserved_by_recomposition':False,
 'pure_primitive_factor_inferred':False,
 'actual_SL2_control_claimed_as_pure_or_exceptional_failure':False,
 'algebraic_log_coefficients_identified_with_algebraic_full_resolvents':False,
 'all_coefficient_ring_topological_promotion_claimed':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('haar_residue_log_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
