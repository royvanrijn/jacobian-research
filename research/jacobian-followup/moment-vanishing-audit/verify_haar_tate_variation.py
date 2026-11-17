"""Exact transcription controls, not a verifier for universal MHS arguments."""
from pathlib import Path
import hashlib,json
import sympy as s

# A rank-one nilpotent cannot be skew for a nondegenerate symmetric form.
N=s.Matrix([[0,1],[0,0]])
a,b,c=s.symbols('a b c');S=s.Matrix([[a,b],[b,c]])
solutions=s.solve(list(N.T*S+S*N),[a,b],dict=True)
assert solutions==[{a:0,b:0}]
assert S.det().subs(solutions[0])==0
# The alternating-pairing countercontrol retains the even-weight hypothesis.
A=s.Matrix([[0,1],[-1,0]])
assert A.det()==1 and N.T*A+A*N==s.zeros(2)

# Cohomological weights 4,5,6: the allowed Tate extension drops two weights.
# Both have N=0 on gr_W. The adjacent-weight candidate fails M=W's axiom.
weights=[4,5,6]
allowed=s.zeros(3);allowed[0,2]=1
forbidden=s.zeros(3);forbidden[0,1]=1
assert allowed**2==forbidden**2==s.zeros(3)
def lowers_two(M):
    return all(M[i,j]==0 or weights[i]<=weights[j]-2
               for i in range(3) for j in range(3))
assert lowers_two(allowed) and not lowers_two(forbidden)
assert all(allowed[i,i]==forbidden[i,i]==0 for i in range(3))
# Dual homological motion starts in weight -4 and ends in weight -6.
gamma=s.Matrix([1,0,0]);dual_motion=-allowed.T*gamma
assert dual_motion==s.Matrix([0,0,-1])
assert weights[2]*(-1)==-6
primitive_tube_twist=-2
assert -6-primitive_tube_twist==-4

# The ambient quadric Gysin map has its one-dimensional Tate kernel.
G=s.Matrix([[1,1]])
assert G.rank()==1 and G*s.Matrix([1,-1])==s.zeros(1,1)

# Two vertices with r parallel edges contribute exactly r-1 graph cycles.
graphs=[]
for r in range(1,7):
    incidence=s.Matrix([[-1]*r,[1]*r])
    assert incidence.rank()==1 and len(incidence.nullspace())==r-1
    # Adding a leaf leaves the Betti number unchanged.
    leaf=incidence.col_join(s.zeros(1,r)).row_join(s.Matrix([0,-1,1]))
    assert leaf.cols-leaf.rank()==r-1
    graphs.append({'parallel_edges':r,'first_Betti_number':r-1})

# Universal family formula, retaining the two different weight pieces.
q=s.symbols('q',integer=True,positive=True)
family=[]
for N_group,e,k in [(24,12,48),(60,30,120)]:
    assert N_group==2*e
    assert s.cancel(2*(q*e)/N_group-1)==q-1
    genus=k*(k-2)//N_group
    assert k*(k-2)%N_group==0
    family.append({'projective_group_order':N_group,'projective_exponent':e,
                   'control_k':k,'quotient_boundary_genus':genus,
                   'larger_radical_genus_lower_bound':2*genus,
                   'Tate_boundary_lower_bound':2*k//N_group-1})

out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'rank_one_symmetric_skew_obstruction_exact':True,
 'alternating_parity_countercontrol_exact':True,
 'adjacent_weight_extension_fails_relative_drop_two_exact':True,
 'Tate_weight_extension_and_dual_motion_exact':True,
 'dual_Haar_motion_weights':{'ambient':-4,'complement_motion':-6,'primitive_cycle':-4},
 'quadric_Gysin_kernel_control_exact':True,
 'boundary_graph_incidence_controls':graphs,
 'safe_family_formula_controls':family,
 'all_degree_Tate_lower_bound_in_k_over_e':'q - 1',
 'universal_MHS_admissibility_and_topological_proofs_machine_verified':False,
 'all_exceptional_pure_fibres_semisimple_at_eigenvalue1':False,
 'invariant_Tate_space_zero_on_entire_pure_locus':False,
 'actual_Haar_variation_proved_zero_for_every_exceptional_phase':False,
 'abstract_weight_controls_realized_as_exceptional_phases':False,
 'all_coefficient_ring_topological_promotion_claimed':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('haar_tate_variation_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
