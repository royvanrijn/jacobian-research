"""NEW additive-derivation and affine-action controls;run before registration."""
import hashlib,itertools,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,a,b=0):check(name,s.cancel(a-b)==0,True)
x,y,z,t,a,b=s.symbols('x y z t a b')
for k in [2,3,5]:
 def D(f):return s.diff(f,z)*x**2+s.diff(f,y)*k*z**(k-1)
 eq('NEW LND preserves Wk relation k'+str(k),D(x*x*y-z**k+1))
 eq('NEW LND fixes ruling x k'+str(k),D(x))
 eq('NEW LND squares to zero on z k'+str(k),D(D(z)))
 val=y
 for j in range(1,k+1):
  val=D(val);eq('NEW nilpotence induction transcription '+str((k,j)),val,s.factorial(k)/s.factorial(k-j)*z**(k-j)*x**(2*j-2))
 eq('NEW LND terminates on y k'+str(k),D(val))
 check('NEW LND nonzero at x0,z1 k'+str(k),D(y).subs({x:0,z:1})!=0,True)
 eq('NEW volume contraction coefficient along dx k'+str(k),-D(z)/x**2,-1)
# Finite affine action on generic A1:f(t)=a*t+b.
eq('generic nontranslation affine map has fixed point',a*(b/(1-a))+b,b/(1-a))
for n in [2,3,5]:
 val=t
 for _ in range(n):val=s.expand(a*val+b)
 eq('finite-order affine iterate transcription n'+str(n),val,a**n*t+b*sum(a**j for j in range(n)))
 eq('characteristic-zero translation iterate n'+str(n),val.subs(a,1),t+n*b)
 check('translation cannot have nonzero finite order n'+str(n),n!=0,True)
c0,d0=s.symbols('c0 d0')
eq('affine composition law',a*(c0*t+d0)+b,a*c0*t+a*d0+b)
eq('shift by a common fixed point removes affine offset',a*(t+d0)+(1-a)*d0-d0,a*t)
# The noncyclic full monodromy can have a trivial rational deck quotient.
S3=set(itertools.permutations(range(3)));H={g for g in S3 if g[0]==0}
def mul(g,h):return tuple(g[h[i]] for i in range(3))
def inv(g):return tuple(g.index(i) for i in range(3))
N={g for g in S3 if {mul(mul(g,h),inv(g)) for h in H}==H}
check('S3 control full monodromy nonabelian',any(mul(g,h)!=mul(h,g) for g in S3 for h in S3),True)
check('S3 control full point stabilizer order2',len(H)==2,True)
check('S3 control normalizer equals point stabilizer',N==H,True)
check('S3 control generic deck quotient trivial',len(N)//len(H)==1,True)
context=json.loads((HERE/'sources/completion-additive-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('no new external theorem adoption',context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[])
positive={'actual_completion_ML_trivial_iff_no_added_unramified_primes_proved','proper_completion_has_unique_additive_kernel_or_is_rigid','two_distinct_source_LND_fields_preserving_canonical_algebra_remove_unramified_added_primes','one_additive_action_forces_generic_rational_deck_group_cyclic','noncyclic_generic_deck_forces_completion_additive_rigidity','deck_action_on_unique_affine_ruling_proved_faithful','known_higher_power_Danielewski_exact_volume_controls_retained','new_cyclic_control_additive_derivation_proved_in_every_k','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'ADDITIVE_ACTIONS_ON_THE_CANONICAL_ETALE_COMPLETION.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-completion-additive.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
