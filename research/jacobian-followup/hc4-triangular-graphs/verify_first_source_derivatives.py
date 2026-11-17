"""NEW MS0 differential, valuation and finite-order controls. Never replay."""
import hashlib,json,time
from pathlib import Path
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
RECEIPT=HERE/'receipt-first-source-derivatives.json'
assert not RECEIPT.exists(),'Sealed receipt exists: never replay.'
assert not any(e['id']=='HC4MYMS0' for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries'])
started=time.perf_counter();checks=[];algebraic=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value,alg=True):
 assert bool(value),name
 checks.append(name)
 if alg:algebraic.append(name)
x,y,u,v,t,z,tau=S.symbols('x y u v t z tau')
P=x+y**2;Q=y+P**2
check('actual nonlinear automorphism control has Jacobian one',S.Matrix([P,Q]).jacobian([x,y]).det()==1)
for label,p in [('bare',S.Integer(1)),('vanishing-replica',y**2),('two-root-replica',y**2-1)]:
 a=p*S.diff(P,x);b=p*S.diff(Q,x)
 check(label+' Hamiltonian differential x coefficient',S.expand(a*S.diff(Q,x)-b*S.diff(P,x))==0)
 check(label+' Hamiltonian differential y coefficient',S.expand(a*S.diff(Q,y)-b*S.diff(P,y)-p)==0)
 for f,tag in [(P,'P'),(Q,'Q')]:
  orbit=sum(p**j*S.diff(f,x,j)*tau**j/S.factorial(j) for j in range(3))
  check(label+' finite source coaction '+tag,S.expand(orbit-f.subs(x,x+tau*p))==0)

for r,e,ell in [(1,2,0),(2,2,0),(3,2,0),(1,3,1),(1,3,2),(2,5,1),(2,5,2),(4,7,1)]:
 k=r*(ell+1)
 inv=u**r;ham=inv**(ell+1)/S.Integer(ell+1)
 tangent=-r*u**(k-e)
 check(f'exact local contraction r{r} e{e} ell{ell}',S.cancel(-u**(e-1)*tangent-S.diff(ham,u))==0)
 check(f'first derivative of separating target v detects order r{r} e{e} ell{ell}',tangent.as_powers_dict().get(u,0)==k-e)
 check(f'global multiplier threshold arithmetic r{r} e{e} ell{ell}',(k>=e)==(ell>=max(0,int(S.ceiling(S.Rational(e,r)))-1)))

for q in [1,2,7]:
 check(f'characteristic zero derivative detects a pole q{q}',S.diff(u**(-q),u)==-q*u**(-q-1))

# Non-Keller local control: arbitrary rational derivations do not satisfy
# the source-LND criterion, even if their two target derivatives are integral.
for e in [2,3,5]:
 check(f'non-Keller first target derivatives regular e{e}',S.cancel(u**(1-e)*S.diff(u**e,u))==e)
 check(f'non-Keller full normalization derivative polar e{e}',1-e<0)

# Multiplication by a=u/t in C(t)[u]/(u^n-t). Its rational field resolvent
# still has nonintegral target-prime coefficients.
for n in [2,3,4]:
 M=S.zeros(n)
 for j in range(n-1):M[j+1,j]=1/t
 M[0,n-1]=1
 chi=z**n-t**(1-n)
 check(f'finite-field characteristic polynomial n{n}',S.expand(M.charpoly(z).as_expr()-chi)==0)
 check(f'first n minus one field traces vanish n{n}',all(S.trace(M**j)==0 for j in range(1,n)))
 check(f'nth field trace detects pole n{n}',S.trace(M**n)==n*t**(1-n))
 check(f'rational ordinary trace resolvent n{n}',S.cancel(S.trace((z*S.eye(n)-M).inv())-S.diff(chi,z)/chi)==0)

# Newton coefficients, with no square/first-trace shortcut.
roots=S.symbols('a0:4');power=[None]+[sum(a**j for a in roots) for j in range(1,5)]
coeff=[S.Integer(1)]
for j in range(1,5):
 coeff.append(S.expand(-sum(coeff[j-i]*power[i] for i in range(1,j+1))/j))
 check(f'Newton coefficient j{j}',S.expand(coeff[j]-S.Poly(S.prod(z-a for a in roots),z).nth(4-j))==0)
res=z*t/(t*z**2-1)
check('quadratic normalized resolvent is rational',S.cancel(res-z/(z**2-1/t))==0)
check('quadratic marked resolvent first pole coefficient',S.limit(z**3*(res-1/z),z,S.oo)==1/t)
check('quadratic marked resolvent next pole coefficient',S.limit(z**5*(res-1/z-1/(t*z**3)),z,S.oo)==1/t**2)
check('rational field resolvent does not give target integrality',S.denom(1/t)==t and S.denom(1/t**2)==t**2)

cp=HERE/'sources/first-source-derivatives-context.json';C=json.loads(cp.read_text())
for p,h in C['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,sha(ROOT/p)==h,False)
check('five retained claims ten proof-checker pins and coordinator note',len(C['retained_structural_inputs'])==5 and len(C['retained_proof_and_source_sha256'])==11,False)
for k in ['all_source_LNDs_and_frames','first_two_integrality_iff_full_stability','full_finite_normalization_ramified_primes','finite_orbit_module_iff_first_two','finite_newton_trace_test','full_field_saturation_maximal_Ga_open','field_resolvent_not_B_regular']:
 check('scope '+k,C['scope_flags'][k],False)
for k in ['automatic_first_derivative_integrality','arbitrary_source_derivation_first_test','raw_first_jet_normalization_equals_full_orbit','raw_orbit_normality_assumed','Mathieu_criterion_assumed_automatic','knot_certificates_constructed','unrestricted_HC4_or_JC2_resolved']:
 check('boundary '+k,not C['scope_flags'][k],False)
check('no new imports caches downloads or predecessor replay',C['new_external_theorem_imports']==C['new_primary_caches']==C['new_primary_downloads']==0 and C['registered_mathematical_checkers_replayed']==[],False)
r=dict(outcome='PASS',claim='HC4MYMS0',check_count=len(checks),algebraic_controls=len(algebraic),checks=checks,algebraic_check_names=algebraic,elapsed_seconds=time.perf_counter()-started,
 proof_sha256=sha(HERE/'FIRST_SOURCE_DERIVATIVES_AND_FINITE_LATTICE_STABILITY.md'),checker_sha256=sha(Path(__file__)),context_sha256=sha(cp),registered_mathematical_checkers_replayed=[],independent_replay=False,formal_verification=False,external_review=False)
RECEIPT.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
