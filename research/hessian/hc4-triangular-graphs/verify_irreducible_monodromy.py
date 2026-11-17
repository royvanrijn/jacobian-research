"""New target-frame and puncture-extension controls;run before registration."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,a,b=0):check(name,s.cancel(a-b)==0,True)
u,v,t,b,c=s.symbols('u v t b c')
g=v**5+3*v**2-v
P=u+g+b*v;Q=v
eq('NEW target shear inverse recovers original u',P-g-b*Q,u)
eq('NEW target shear preserves Jacobian',s.det(s.Matrix([[s.diff(P,u),s.diff(P,v)],[s.diff(Q,u),s.diff(Q,v)]])),1)
# A general target shear needs no source-map formula or source-coordinate change.
z=s.symbols('z')
for e in [1,2,5]:
 eq('NEW local Q pole differential order '+str(e),s.diff(z**(-e),z),-e*z**(-e-1))
 eq('NEW finite Q puncture differential order '+str(e),s.diff(c+z**e,z),e*z**(e-1))
solutions=[(genus,degree,punctures) for genus in range(4) for degree in range(1,8) for punctures in range(1,8) if 2*genus-2==-degree-punctures]
check('NEW all-poles degree bookkeeping control only',solutions==[(0,1,1)],True)
A=s.Matrix([[-1,0,1],[0,-1,0],[0,0,1]])
B=s.diag(-1,-1,1);I=s.eye(3)
F=(A-I).col_join(B-I)
check('NEW affine H fixed space zero',F.rank()==3,True)
check('NEW compact V fixed space zero',(B[:2,:2]-s.eye(2)).rank()==2,True)
check('NEW puncture W fixed space dimension one',A[2,2]==B[2,2]==1,True)
check('NEW compact images finite',(A[:2,:2]**2)==(B[:2,:2]**2)==s.eye(2),True)
J=s.Matrix([[0,1],[-1,0]])
check('NEW compact image symplectic',A[:2,:2].T*J*A[:2,:2]==J,True)
N=A*B-I
check('NEW full affine image nontrivial unipotent',N!=s.zeros(3) and N*N==s.zeros(3),True)
for n in [2,3,7]:
 check('NEW full affine product power '+str(n),(A*B)**n==I+n*N,True)
q1,q2=s.symbols('q1 q2');lift=s.Matrix([q1,q2,1])
check('NEW no equivariant W section',s.linsolve(list((A-I)*lift)+list((B-I)*lift),[q1,q2]) is s.EmptySet,True)
# A comparison split extension: same graded pieces,nonzero affine invariants.
split=B
check('NEW split comparison has one fixed class',3-(split-I).rank()==1,True)
for orbit_sizes in [[1,1],[2,3],[1,2,4]]:
 h=len(orbit_sizes);ss=sum(orbit_sizes);perm=s.zeros(ss);off=0
 for size in orbit_sizes:
  for j in range(size):perm[off+(j+1)%size,off+j]=1
  off+=size
 check('NEW permutation residues invariant rank '+str(orbit_sizes),ss-(perm-s.eye(ss)).col_join(s.ones(1,ss)).rank()==h-1,True)
context=json.loads((HERE/'sources/irreducible-monodromy-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
for p,h in context['new_primary_sha256'].items():check('new primary byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('ONE explicit published reduction,no unproved import',len(context['newly_imported_proved_external_statements'])==1 and context['newly_used_unproved_external_statements']==[])
positive={'published_actual_all_irreducible_target_frame_adopted','target_ring_full_normalization_and_trace_dual_preserved','fixed_cohomology_dimension_and_residue_injectivity_proved','quantified_puncture_extension_rank_proved','counterexample_requires_nonsplit_infinite_affine_monodromy','canonical_boundary_exclusions_apply_in_supplied_frame','counterexample_normalized_raw_quotient_torsion_free_positive_rank','newly_adopted_external_theorem','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'IRREDUCIBLE_TARGET_FRAMES_AND_PUNCTURE_EXTENSIONS.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-irreducible-monodromy.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
