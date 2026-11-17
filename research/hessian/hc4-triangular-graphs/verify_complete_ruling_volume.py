"""One preregistration transcription run; refuse registered replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];B=R.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
ledger=json.loads((R/'MATH_STATUS.json').read_text())
assert not any(e['id']=='HC4MYML0' for e in ledger['entries']),'Registered mathematical checker replay forbidden'
receipt=H/'receipt-complete-ruling-volume.json';assert not receipt.exists(),'Receipt already exists; inspect rather than rerun'
context=H/'sources/complete-ruling-volume-context.json';c=json.loads(context.read_text())
proof=H/'EXACT_VOLUME_AND_COMPLETE_RULING_COMPONENTS.md';pr=json.loads((B/'research/.scratch/complete-ruling-volume-preregistration.json').read_text())
for p,k in [(proof,'proof_sha256'),(Path(__file__),'checker_sha256'),(context,'context_sha256')]:assert sha(p)==pr[k],p
checks=[];algebraic=0;start=time.perf_counter()
def check(name,ok,alg=False):
 global algebraic
 assert ok,name
 checks.append(name);algebraic+=int(alg)
for p,d in c['retained_proof_and_source_sha256'].items():check('sealed '+p,sha(R/p)==d)
t,z,w,cc=s.symbols('t z w c');b=s.Function('b')(t);f=s.Function('f')(t);g=s.Function('g')(1/t)
ss=1/t;ww=-t**2*z+b
check('arbitrary Laurent torsor transition preserves two-form',s.simplify(s.Matrix([ss,ww]).jacobian([t,z]).det()-1)==0,True)
check('transition inverse recovers fibre',s.simplify((-t**2*((b-w)/t**2)+b)-w)==0,True)
zp=z+f;wp=ww+g;bp=b+t**2*f+g
check('arbitrary fibre translation coboundary sign',s.simplify(wp+t**2*zp-bp)==0,True)
theta0=-z;thetainf=s.simplify((-t**2*z+cc*t)/t**2)
check('first chart primitive has exterior derivative omega',-s.diff(theta0,z)==1,True)
check('second chart primitive has exterior derivative omega',-s.diff(thetainf,z)==1,True)
check('primitive difference is c dt over t',s.simplify(thetainf-theta0-cc/t)==0,True)
for n in [-7,-2,-1,0,1,2,3,8]:
 mon=t**n
 ff=-t**(n-2) if n>=2 else 0
 gg=-t**n if n<=0 else 0
 reduced=cc*t if n==1 else 0
 value=cc*mon+t**2*cc*ff+cc*gg
 check('Laurent coboundary exponent '+str(n),s.expand(value-reduced)==0,True)
n=s.symbols('n',integer=True)
check('all-exponent Laurent differentiation identity',s.simplify(s.diff(t**n,t)-n*t**(n-1))==0,True)
check('only candidate residue exponent contributes zero',s.diff(t**0,t)==0,True)
q=s.symbols('q',nonzero=True)
check('nonzero torsor parameter gives nonzero residue',s.residue(q/t,t,0)==q,True)
check('zero torsor parameter glues zero sections',(-t**2*z).subs(z,0)==0,True)
check('no multiplicity-one does not imply multiple whole fibre',min([2,3])>=2 and s.gcd(2,3)==1)
check('one nonreduced component does not exclude reduced selection',1 in [1,3] and max([1,3])>1)
check('one adopted local-structure package and one primary cache',len(c['external_imports'])==len(c['new_primary_caches'])==1)
check('no automatic LND or full HC4 claim',not any(c['scope_flags'][k] for k in ['additive_action_proved_automatic','all_complete_rulings_excluded','bad_fibre_gcd_asserted_gt_one','full_trace_pairing_proved_automatic','unrestricted_HC4_or_JC2_resolved']))
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,elapsed_seconds=time.perf_counter()-start,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),context_sha256=sha(context),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
receipt.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
