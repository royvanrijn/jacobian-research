"""Exact transcription controls for the written cycle/linking theorem."""
from pathlib import Path
import hashlib,json
import sympy as s

def linking_identity(B,T,d,l):
    assert T.T*B*T==B
    assert d.T*B*T==l.T*(s.eye(T.rows)-T)

# Nondegenerate intersection and semisimple eigenvalue 1: its motion vanishes.
B=s.diag(2,3);T=s.diag(1,-1);b=s.Matrix([5,7])
d=(T-s.eye(2))*b;l=B*b
linking_identity(B,T,d,l);assert d[0]==0
# Degenerate invariant intersection direction permits a nonzero translation.
Br=s.diag(0,2);Tr=s.diag(1,-1);dr=s.Matrix([7,-2]);lr=s.Matrix([0,2])
linking_identity(Br,Tr,dr,lr)
Ar=Tr.row_join(dr).col_join(s.Matrix([[0,0,1]]))
assert (Ar-s.eye(3))**2!=s.zeros(3)
# Use exact Bezout interpolation with multiplicity for the Jordan block.
q=s.symbols('q')
pol=s.rem((q+1)*s.invert(q+1,(q-1)**2),(q-1)**2*(q+1),q)
proj=s.zeros(3)
for (k,),a in s.Poly(pol,q).terms():proj+=a*Ar**k
assert proj**2==proj
gamma=s.Matrix([0,0,1]);assert (Ar-s.eye(3))*proj*gamma==s.Matrix([7,0,0])
# A split abstract extension can move the specified lift when fibre T is
# unipotent, even with nondegenerate symmetric intersection.
Bu=s.zeros(4);Bu[:2,2:]=s.eye(2);Bu[2:,:2]=s.eye(2)
Tu=s.eye(4);Tu[0,3]=1;Tu[1,2]=-1
bu=s.Matrix([0,0,0,1]);du=(Tu-s.eye(4))*bu;lu=Bu*bu
linking_identity(Bu,Tu,du,lu);assert Bu.det()!=0 and du!=s.zeros(4,1)
Au=Tu.row_join(du).col_join(s.Matrix([[0,0,0,0,1]]))
fixed=s.Matrix(list(-bu)+[1]);assert Au*fixed==fixed
assert Au*s.Matrix([0,0,0,0,1])!=s.Matrix([0,0,0,0,1])
# Period derivatives genuinely supply higher-pole original forms.
x,y=s.symbols('x y')
for k in range(1,7):
    assert s.diff(1/(x-y),x,k-1)==(-1)**(k-1)*s.factorial(k-1)/(x-y)**k
assert s.simplify(x/(x-y)-y/(x-y)-1)==0
# Genuine quadric control retained from the prior tame theorem: trivial
# fibre H1 monodromy alone does not rule out a complement logarithm.
C=s.log((x+2)/x)/2
assert s.simplify(s.diff(C,x)+1/(x*(x+2)))==0
out={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'linking_matrix_identity_controls_exact':True,
 'semisimple_nondegenerate_eigenvalue1_variation':[str(v) for v in d],
 'invariant_boundary_radical_translation':[str(v) for v in dr],
 'actual_eigenvalue1_projected_motion_in_radical_control':['7','0'],
 'split_extension_can_move_a_specified_lift':True,
 'higher_pole_derivative_and_original_form_spanning_controls_exact':True,
 'quadric_Cauchy_log_derivative_control_exact':True,
 'universal_Thom_linking_period_and_monodromy_proofs_machine_verified':False,
 'sufficient_fibre_boundary_hypotheses_verified_for_all_exceptional_phases':False,
 'exceptional_variation_cycle_proved_zero_for_all_pure_phases':False,
 'abstract_controls_claimed_as_exceptional_phase_realizations':False,
 'all_coefficient_ring_topological_criterion_claimed':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'new_degree_support_or_prefix_campaign':False}
Path(__file__).with_name('haar_cycle_linking_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
