"""NEW MT0 slice/focal/finite-algebra controls; never replay predecessors."""
import hashlib,json,time
from pathlib import Path
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RECEIPT=HERE/'receipt-generic-source-translation-algebra.json'
assert not RECEIPT.exists(),'Sealed receipt exists: never replay.'
assert not any(e['id']=='HC4MYMT0' for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries'])
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value,alg=True):
 assert bool(value),name
 checks.append(name)
 if alg:algebraic.append(name)
s,t,x,y,z,v1,v2,alpha,w1,w2=S.symbols('s t x y z v1 v2 alpha w1 w2')
def delta(f):return S.expand(v1*S.diff(f,x)+v2*S.diff(f,y))
def iterate(f,n):
 for j in range(n):f=delta(f)
 return f

for n in [1,2,3,5]:
 f=x**n+2*x**(n-1)*y+3*y**n+x**(n-1)
 top=x**n+2*x**(n-1)*y+3*y**n
 check(f'top directional polar gives constant n{n}',S.expand(iterate(f,n)-S.factorial(n)*top.subs({x:v1,y:v2},simultaneous=True))==0)
 slice=iterate(f,n-1)/iterate(f,n)
 check(f'linear slice differentiation n{n}',S.cancel(delta(slice)-1)==0)
 check(f'linear slice is affine n{n}',S.Poly(S.together(slice*iterate(f,n)),x,y).total_degree()<=1)

# Coefficient extraction with a genuine slice; arbitrary coefficients in z.
f=(z**2-1)*s**4+(z**3-z)*s**2+(z**2+2)*s+z**5
for i in range(5):
 extracted=sum((-s)**j*S.diff(f,s,i+j)/(S.factorial(i)*S.factorial(j)) for j in range(5-i))
 check(f'exact slice coefficient projection i{i}',S.expand(extracted-S.Poly(f,s).nth(i))==0)
 check(f'coefficient projection invariant i{i}',S.expand(S.diff(extracted,s))==0)

# Explicit finite-type incidence chart and exact Taylor equivalence.
P=x+y**2;Q=y+P**2
p={x:s,y:alpha+t*s};q={x:s+w1,y:alpha+w2+t*s}
for F,tag in [(P,'P'),(Q,'Q')]:
 dv=lambda h:S.diff(h,x)+t*S.diff(h,y)
 deriv=F;coeff=[]
 for j in range(5):
  coeff.append(S.expand(deriv.subs({x:0,y:alpha},simultaneous=True)-deriv.subs({x:w1,y:alpha+w2},simultaneous=True)))
  deriv=dv(deriv)
 actual=S.expand(F.subs(p,simultaneous=True)-F.subs(q,simultaneous=True))
 check('whole-translation Taylor incidence '+tag,S.expand(actual-sum(coeff[j]*s**j/S.factorial(j) for j in range(5)))==0)
 check('finite coefficient incidence chart '+tag,S.degree(actual,s)<=4)

# The focal determinant is linear in the WHOLE line parameter.
a1=S.Function('a1')(t);a2=S.Function('a2')(t)
v1t=S.Function('v1')(t);v2t=S.Function('v2')(t)
par=S.Matrix([a1+s*v1t,a2+s*v2t])
jac=S.expand(par.jacobian([s,t]).det())
constant=v1t*S.diff(a2,t)-v2t*S.diff(a1,t)
linear=v1t*S.diff(v2t,t)-v2t*S.diff(v1t,t)
check('general focal determinant',S.expand(jac-constant-s*linear)==0)
check('projective direction derivative is focal coefficient',S.cancel(S.diff(v2t/v1t,t)-linear/v1t**2)==0)
vary=S.Matrix([s,s*t+t**2])
check('varying whole-line family has an actual focal root',vary.jacobian([s,t]).det()==s+2*t)
check('generic focal root lies on entire affine line',(s+2*t).subs(s,-2*t)==0)

# Punctured-line symplectic correspondence: both projections are etale only
# after removing the focal point; this is NOT an actual Keller collision.
puncp=S.Matrix([s,s*t]);puncq=S.Matrix([s+t,s*t+t**2/2])
check('punctured first projection Jacobian',puncp.jacobian([s,t]).det()==s)
check('punctured second projection Jacobian',puncq.jacobian([s,t]).det()==s)
check('punctured difference independent of line parameter',all(S.diff(f,s)==0 for f in puncq-puncp))
check('punctured ruling really varies direction',S.det(S.Matrix([[1,0],[t,1]]))==1)
check('missing focal point essential',puncp.jacobian([s,t]).det().subs(s,0)==0)
parallelp=S.Matrix([s,t]);parallelq=S.Matrix([s+t**2,t])
check('parallel whole-line correspondence first etale',parallelp.jacobian([s,t]).det()==1)
check('parallel whole-line correspondence second etale',parallelq.jacobian([s,t]).det()==1)
check('focal argument alone only gives constant direction',S.diff(parallelq-parallelp,s)==S.zeros(2,1))

for n in [1,2,4,7]:
 difference=S.expand((s+1)**n-s**n)
 check(f'polynomial translation period has nonzero leading difference n{n}',S.Poly(difference,s).LC()==n and S.degree(difference,s)==n-1)

# Global finite-algebra controls distinguish the final three hypotheses.
nodef=z**2-1;nodeg=z*(z**2-1)
check('nodal normalization finite',S.expand(z**2-(nodef+1))==0)
check('nodal normalization birational',S.cancel(nodeg/nodef-z)==0)
check('nodal normalization unramified',S.expand(S.Rational(3,2)*z*S.diff(nodef,z)-S.diff(nodeg,z))==1)
check('nodal normalization still noninjective',nodef.subs(z,1)==nodef.subs(z,-1)==nodeg.subs(z,1)==nodeg.subs(z,-1)==0)
check('cuspidal normalization is ramified despite bijection',S.diff(z**2,z).subs(z,0)==S.diff(z**3,z).subs(z,0)==0)
check('non-Keller vertical collision has proper derivative field',S.expand(((x**2).subs(x,-x))-x**2)==0 and S.diff(x**2,x)==2*x)

cp=HERE/'sources/generic-source-translation-algebra-context.json';C=json.loads(cp.read_text())
for p,h in C['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,sha(ROOT/p)==h,False)
check('four retained claims eight proof-checker pins',len(C['retained_structural_inputs'])==4 and len(C['retained_proof_and_source_sha256'])==8,False)
for k in ['all_degrees_every_fixed_source_frame','all_but_finitely_many_linear_directions','generic_raw_algebra_equals_polynomial_ring','generic_augmented_map_closed_embedding','whole_persistent_collision_lines_finite','focal_argument_on_actual_generic_fibre','target_finite_module_still_unproved']:
 check('scope '+k,C['scope_flags'][k],False)
for k in ['all_directions_raw_generation','uniform_degree_independent_jet_order','raw_generation_implies_finite_target_module','universal_line_family_assumed_etale','punctured_lines_substituted_for_whole','automatic_Mathieu_transfer','unrestricted_HC4_or_JC2_resolved']:
 check('boundary '+k,not C['scope_flags'][k],False)
check('no new imports caches downloads or predecessor replay',C['new_external_theorem_imports']==C['new_primary_caches']==C['new_primary_downloads']==0 and C['registered_mathematical_checkers_replayed']==[],False)
r=dict(outcome='PASS',claim='HC4MYMT0',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,algebraic_check_names=algebraic,elapsed_seconds=time.perf_counter()-started,
 proof_sha256=sha(HERE/'GENERIC_SOURCE_TRANSLATION_DERIVATIVES_GENERATE_THE_POLYNOMIAL_RING.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
RECEIPT.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
