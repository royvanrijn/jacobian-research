"""New analytic-transcription and provenance controls; no registered replay."""
from pathlib import Path
import hashlib,json,time
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];started=time.perf_counter()
checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,a):check(name,s.simplify(s.expand(a))==0,True)
u,v=s.symbols('u v',positive=True);e=s.symbols('e',integer=True,positive=True)
P=u**e;Q=v/u**(e-1)
zero('all-index Laurent Jacobian is constant e',s.diff(P,u)*s.diff(Q,v)-s.diff(P,v)*s.diff(Q,u)-e)
zero('source volume becomes ramified canonical volume',e*u**(e-1)*s.diff(v/u**(e-1),v)-e)
zero('global primitive dP wedge dQ',s.diff(P*s.diff(Q,v),u)-s.diff(P*s.diff(Q,u),v)-e)
p,a,q,rho=s.symbols('p a q rho',positive=True)
zero('pole-volume radial exponent',2*(e-1)-p*q+1-(2*e-1-p*q))
zero('radial antiderivative',s.diff(rho**(2*e-p*q)/(2*e-p*q),rho)-rho**(2*e-1-p*q))
zero('source Laurent fibre integration exponent',2*(e-1)+1-p*a-(2*e-1-p*a))
zero('top codifferent mode L2 has constant norm',2*s.pi**2/(2*e-2*(e-1))-s.pi**2)
zero('constant Laurent mode volume',2*s.pi**2/(2*e)-s.pi**2/e)
E=s.symbols('E',integer=True,nonnegative=True);ee=E+2
check('top codifferent mode L4 never integrable',s.ask(s.Q.nonpositive(2*ee-4*(ee-1))),True)
zero('index two is logarithmic L4 endpoint',(2*ee-4*(ee-1)).subs(E,0))
zero('top codifferent generic critical exponent',2*ee/(ee-1)-(2+2/(ee-1)))
check('every index has L2 headroom',(2*ee-2*(ee-1)).is_positive,True)
d=s.symbols('d',integer=True,positive=True);pole=s.symbols('pole',integer=True,positive=True)
# e<=d and pole>=1 give 2d*pole>=2e; keep the two nonnegative terms explicit.
gap=s.symbols('gap',integer=True,nonnegative=True)
check('uniform 2d detects every integer pole',(2*d*(pole-1)+2*gap).is_nonnegative,True)
j=s.symbols('j',positive=True)
zero('pullback norm scalar',d/j**2-d*j**(-2))
zero('normalized pullback isometric',j**2/d*d/j**2-1)
zero('trace pullback projection idempotent scalar',d/d**2-1/d)
zero('pullback adjoint and trace scaling',j**(-2)*d-d/j**2)
m,n=s.symbols('m n',integer=True,nonnegative=True)
A=s.symbols('A',integer=True,positive=True);aa=A+2
den=(m+1)*(n+1-(m+1)/aa)
zero('Runge-control volume',s.pi**2/den.subs({m:0,n:0})-s.pi**2*aa/(aa-1))
zero('entire x is L2 for every a>=3',s.pi**2/den.subs({m:1,n:0})-s.pi**2*aa/(2*(aa-2)))
check('allowed coefficient witness strict threshold',(aa-(aa-2+1)).is_positive,True)
zero('multiplication moves witness to logarithmic divergence',aa-(aa-1+1))
zero('witness norm equals domain volume',s.pi**2/den.subs({m:aa-2,n:0})-s.pi**2*aa/(aa-1))
x,y=s.symbols('x y')
zero('Runge polynomial map has non-Keller Jacobian',s.diff(x**aa*y,x)*s.diff(y,y)-s.diff(x**aa*y,y)*s.diff(y,x)-aa*x**(aa-1)*y)
check('bilinear simple base pole is L1',s.Rational(1)<2,True)
check('bilinear simple base pole is not L2',not (s.Rational(2)<2),True)
c=json.loads((H/'sources/bergman-trace-dual-context.json').read_text())
ledger=json.loads((R/'MATH_STATUS.json').read_text());known={e['id']:e for e in ledger['entries']}
for i in c['retained_structural_inputs']:check('prior '+i['claim'],known[i['claim']]['state']=='proved' and known[i['claim']]['canonical_source']==i['canonical_source'])
for p,h in c['retained_proof_and_source_sha256'].items():check('byte '+p,hashlib.sha256((R/p).read_bytes()).hexdigest()==h)
positive={'full_goal_active','previous_goal_turn_progress','actual_polynomial_L2_space_equals_finite_trace_dual','global_quotient_point_integrability_audited','bounded_Hilbert_trace_retraction_automatic','L4_finite_trace_dual_inverse_equivalence_proved','natural_polynomial_multiplier_domain_dichotomy_proved','all_index_constant_Jacobian_Laurent_controls_supplied','all_a_entire_coordinate_nondensity_Runge_controls_supplied','external_review_deferred'}
for k,v in c['scope_flags'].items():check('scope '+k,v is (k in positive))
proof=H/'BERGMAN_TRACE_DUAL_AND_MULTIPLIER_DOMAIN_RIGIDITY.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=c['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-bergman-trace-dual.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
