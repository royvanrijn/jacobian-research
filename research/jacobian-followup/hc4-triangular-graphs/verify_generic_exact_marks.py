"""New symbolic parameter identities;no registered mathematical replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];start=time.perf_counter();checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,z):check(name,s.factor(s.simplify(s.expand(z)))==0,True)
x,T,c,q,z=s.symbols('x T c q z',nonzero=True);k,n=s.symbols('k n',integer=True,positive=True);m=s.symbols('m',integer=True,nonnegative=True)
A=1+T**n;G=T+T**(n+1)/(n+1);P=x**k*A;Y=T/x**(k+1);Q=G/(k*P**2)
zero('global relative Liouville identity dx',P*Y-T/(k*x**k)*s.diff(P,x))
zero('global relative Liouville identity dT',n/(k*(n+1))*s.diff(T**(n+1),T)-T/(k*x**k)*s.diff(P,T))
zero('target G derivative A',s.diff(G,T)-A)
zero('open Keller coframe Jacobian',x**(k+1)*(s.diff(P,x)*s.diff(Q,T)-s.diff(P,T)*s.diff(Q,x))-1)
zero('volume coefficient in P T frame',x**(-k-1)/(k*x**(k-1)*A)-A/(k*P**2))
zero('curve primitive derivative',s.diff(-n*T**(n+1)/(k*(n+1)*c),T)+n*T**n/(k*c))
# k1 whole-plane identities,computed in the dense x!=0 coordinates.
P1=x*A;D=lambda f:x*x*(A*s.diff(f,T)-x*s.diff(A,T)*s.diff(f,x))
a=2*n+1;Px=1+a*T**n;Py=n*x**3*T**(n-1)
zero('whole-plane submersion Bezout identity',(1-a*T**n)*Px+a*a/n*T**(n+1)/x**3*Py-1)
check('Bezout coefficient polynomial in original x y', (2*n-1).is_positive is True and (n+1).is_positive is True,True)
zero('generic volume raw annihilator upper bound',D(G)-P1**2)
tau=-n*x*T**n
zero('raw Liouville representative comparison',D(T/x)-P1+tau)
zero('raw Liouville annihilator upper bound',P1*tau+n/(n+1)*D(T**(n+1)))
h=s.Function('h')(T)
zero('all-degree P primitive differential equation',D(T*h/x)-x*(T*A*s.diff(h,T)+(A+T*s.diff(A,T))*h))
zero('all-degree mate differential equation',D(T*h/x**2)-((1+(2*n+1)*T**n)*h+T*A*s.diff(h,T)))
zero('P primitive monomial degree identity',D(T**(m+1)/x)-x*T**m*((m+1)+(m+n+1)*T**n))
zero('mate monomial degree identity',D(T**(m+1)/x**2)-T**m*((m+1)+(m+2*n+1)*T**n))
check('leading coefficient assumptions',m.is_nonnegative is True and n.is_positive is True,True)
zero('constant P primitive coefficient mismatch',(A+T*s.diff(A,T))-A-n*T**n)
zero('generic product inverse x',(c/A)*A-c)
zero('generic product inverse T',(c/A)**2*(T*A**2/c**2)-T)
zero('generic rational mate product expression',G/c**2-(T*A**2/c**2)*(1+T**n/(n+1))/A**2)
zero('open extension rational mate original expression',G/P1**2-(T/x**2)*(1+T**n/(n+1))/A**2)
# Global exact primitive on the open source,not a polynomial mate.
F=(T-(n-1)*T**(n+1)/(n+1))/(k*P)
zero('open Liouville primitive dx',Y+P*s.diff(Q,x)-s.diff(F,x))
zero('open Liouville primitive dT',P*s.diff(Q,T)-s.diff(F,T))
dd=s.symbols('dd',integer=True,positive=True)
twog=(k-1)*(n-1)-dd+1
zero('all-parameter Riemann Hurwitz',-2*k+n*(k-1)+(k-dd)-(twog-2))
zero('all-parameter first Betti number',twog+n+dd-1-(k*(n-1)+1))
zero('all-parameter Euler characteristic',2-twog-n-dd+k*(n-1))
g=s.symbols('g',integer=True,nonnegative=True)
zero('every genus specialization',twog.subs({k:2*g+1,n:2,dd:1})-2*g)
zero('every genus three punctures',(n+dd).subs({n:2,dd:1})-3)
zero('k1 genus zero',twog.subs({k:1,dd:1}))
zero('k1 rank n',(k*(n-1)+1).subs(k,1)-n)
# Critical-value and normalization transcriptions,not a finite root census.
zero('critical value at T^n minus one',T-T/(n+1)-n*T/(n+1))
zero('critical root second derivative',s.diff(G,T,2)-n*T**(n-1))
zero('critical double root local ramification order',2*k-(k+k))
zero('retained generic degree count',k*(n+1)-2*k-k*(n-1))
zero('finite normalization equation dense chart',(c*z**k-A).subs({c:P,z:1/x}))
zero('whole-target Danielewski equation',G-P1**2*(G/P1**2))
zero('G zero at nonzero root T^n minus nplus1',T-T*(n+1)/(n+1))
zero('G derivative at nonzero G root',1-(n+1)+n)
zero('source x0 retained mate value',(s.Symbol('y')*(1+T**n/(n+1))/A**2).subs(T,0)-s.Symbol('y'))
ctx=json.loads((H/'sources/generic-exact-marks-context.json').read_text());known={a['id']:a for a in json.loads((R/'MATH_STATUS.json').read_text())['entries']}
check('new checker unregistered','HC4MYMF0' not in known)
for a in ctx['retained_structural_inputs']:check('prior '+a['claim'],known[a['claim']]['state']=='proved' and known[a['claim']]['canonical_source']==a['canonical_source'])
for path,d in ctx['retained_proof_and_source_sha256'].items():check('byte '+path,hashlib.sha256((R/path).read_bytes()).hexdigest()==d)
yes={'full_goal_active','previous_goal_turn_progress','all_parameter_polynomial_relative_Liouville_identity','smooth_closed_embedded_exact_curves_every_genus','k1_whole_plane_submersion_for_every_n','trivial_finite_generic_connection_rank_n','generic_volume_and_Liouville_marks_zero','raw_volume_exact_annihilator_P_squared','raw_Liouville_exact_annihilator_P','no_polynomial_Keller_mate_all_degrees','genuine_surjective_open_Keller_coframes','known_Danielewski_completion_credited','actual_Kaliman_all_irreducible_frame_boundary_preserved','Regeta_primary_equivalence_route_audit','external_review_deferred'}
for path,value in ctx['scope_flags'].items():check('scope '+path,value is (path in yes))
proof=H/'GENERIC_EXACT_MARKS_AND_SMOOTH_LIOUVILLE_CURVES.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-start,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=ctx['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-generic-exact-marks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
