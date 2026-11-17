"""New half-lattice/linking/cubic controls; no registered checker replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];started=time.perf_counter()
checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,x):check(name,s.simplify(s.expand(x))==0,True)
m=s.symbols('m',integer=True,positive=True);odd=2*m+1;even=2*m
zero('odd lower half exponent',s.floor((odd-1)/2)-m)
zero('odd upper half exponent',s.ceiling((odd-1)/2)-m)
zero('even lower half exponent',s.floor((even-1)/2)-(m-1))
zero('even upper half exponent',s.ceiling((even-1)/2)-m)
zero('odd trace dual half exponent',(odd-1)-m-m)
zero('even trace dual complementary exponent',(even-1)-(m-1)-m)
zero('odd residual layer vanishes',s.ceiling((odd-1)/2)-s.floor((odd-1)/2))
zero('even residual layer is one',s.ceiling((even-1)/2)-s.floor((even-1)/2)-1)
check('odd quadratic exponents stay below ramification',(odd-2*m).is_positive,True)
zero('universal cubic mixed exponent',m+m+1-odd)
zero('universal nonintegral multiplication coefficient exponent',-2*m-(1-odd))
check('odd half lattice nontrivial',m.is_positive,True)
zero('square discriminant exponent',2*m-(odd-1))
zero('odd self-duality parity',s.Mod(odd-1,2))
zero('even self-duality obstruction parity',s.Mod(even-1,2)-1)
t=s.symbols('t',nonzero=True)
def root_matrix(e):
    U=s.zeros(e)
    for j in range(e-1):U[j+1,j]=1
    U[0,e-1]=t
    return U
# These are independent finite matrix transcriptions, NOT a degree census.
# Universal statements use the written projective/divisorial proofs.
for e in [3,5]:
    n=(e-1)//2;U=root_matrix(e);powers={a:U**a for a in range(-2*n,2*n+1)}
    G=s.Matrix([[s.trace(powers[a+b]) for b in range(-n,n+1)] for a in range(-n,n+1)])
    target=e*s.zeros(e)
    for j in range(e):target[j,e-1-j]=e
    check('new half-lattice constant Gram e'+str(e),G==target,True)
    zero('new half-lattice unimodular determinant e'+str(e),G.det()-(-1)**(e*(e-1)//2)*e**e)
    cubic=s.trace((U**(-n))*(U**(-n))*(U**(-1)))
    zero('new mixed cubic has pole e'+str(e),cubic-e/t)
    # Reconstruct the (-n,-n) product from its cubic tensor and Gram inverse.
    cub=s.Matrix([s.trace((U**(-n))*(U**(-n))*(U**a)) for a in range(-n,n+1)])
    coeff=s.simplify(G.inv()*cub)
    expected=s.zeros(e,1);expected[n+1]=1/t
    check('new cubic reconstructs polar structure constant e'+str(e),coeff==expected,True)
    link=s.Matrix([[s.trace(U**(-a-b)) for b in range(1,e)] for a in range(1,e)])
    check('new odd linking is nondegenerate e'+str(e),s.simplify(link.det())!=0,True)
    check('new odd half subspace is isotropic e'+str(e),link[:n,:n]==s.zeros(n),True)
    check('new odd half subspace has half dimension e'+str(e),2*n==e-1,True)
U=root_matrix(6)
zero('new even residual self-pair',s.trace((U**(-3))**2)-6/t)
check('new even lower half isotropic',s.trace(U**(-2)*U**(-2))==0 and s.trace(U**(-1)*U**(-2))==0,True)
zero('new even lower half orthogonal to residual',s.trace(U**(-2)*U**(-3)))
x,y,z=s.symbols('x y z')
zero('cubic polarization coefficient',s.diff((x+y+z)**3,x,y,z)-6)
c=json.loads((H/'sources/trace-linking-half-lattice-context.json').read_text())
ledger=json.loads((R/'MATH_STATUS.json').read_text());known={a['id']:a for a in ledger['entries']}
for a in c['retained_structural_inputs']:check('prior '+a['claim'],known[a['claim']]['state']=='proved' and known[a['claim']]['canonical_source']==a['canonical_source'])
for p,h in c['retained_proof_and_source_sha256'].items():check('byte '+p,hashlib.sha256((R/p).read_bytes()).hexdigest()==h)
yes={'full_goal_active','previous_goal_turn_progress','global_perfect_trace_linking_proved','half_ramification_reflexive_lattices_base_projective_proved','canonical_residual_linking_supported_only_at_even_indices','odd_S_stable_selfdual_lattice_uniqueness_proved','nonzero_odd_metabolic_linking_controls_supplied','cubic_trace_multiplication_criterion_proved','conditional_actual_odd_three_moment_witness_proved','constant_Gram_Laurent_half_lattice_controls_supplied','external_review_deferred'}
for p,v in c['scope_flags'].items():check('scope '+p,v is (p in yes))
proof=H/'TRACE_LINKING_HALF_RAMIFICATION_AND_CUBIC_MULTIPLICATION.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=c['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-trace-linking-half-lattice.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
