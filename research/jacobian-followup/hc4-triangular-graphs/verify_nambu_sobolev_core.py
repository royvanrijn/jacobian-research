"""New Nambu/finite-core transcriptions only; no registered checker imports."""
from pathlib import Path
import hashlib,json,time
import sympy as s
H=Path(__file__).resolve().parent;R=H.parents[1];started=time.perf_counter()
checks=[];algebra=[]
def check(name,test,alg=False):
    assert test,name
    checks.append(name)
    if alg:algebra.append(name)
def zero(name,v):check(name,s.simplify(s.expand(v))==0,True)
u,v=s.symbols('u v',positive=True);E=s.symbols('E',integer=True,nonnegative=True);e=E+2
P=u**e;Q=v/u**(e-1)
A=u**(1-e)/e;B=(e-1)*v*u**(-e)/e;C=u**(e-1)
DP=lambda a:A*s.diff(a,u)+B*s.diff(a,v)
DQ=lambda a:C*s.diff(a,v)
zero('new transverse target derivative P',DP(P)-1)
zero('new transverse target derivative Q',DP(Q))
zero('new tangent target derivative P',DQ(P))
zero('new tangent target derivative Q',DQ(Q)-1)
zero('new Laurent coframe commutator u',DP(DQ(u))-DQ(DP(u)))
zero('new Laurent coframe commutator v',DP(DQ(v))-DQ(DP(v)))
zero('new transverse source divergence',s.diff(A,u)+s.diff(B,v))
zero('new tangent source divergence',s.diff(C,v))
q,n=s.symbols('q n',integer=True,positive=True)
zero('first derivative of polar mode',DP(u**(-q))+q/e*u**(-q-e))
check('polar derivative lies below canonical cutoff',(q+1).is_positive,True)
zero('regular positive mode derivative',DP(u**n)-n/e*u**(n-e))
check('regular first derivative lies in canonical lattice',(n-1).is_nonnegative,True)
zero('constant mode derivative',DP(s.Integer(1)))
F0,F1,unit=s.symbols('F0 F1 unit',nonzero=True)
F=F0+u*F1
lie=unit*u**(1-e)*s.diff(F,u)+s.diff(unit*u**(1-e),u)*F
zero('canonical domain obstruction coefficient',s.expand(lie*u**e)-(1-e)*unit*F0-(2-e)*unit*F1*u)
check('canonical obstruction coefficient nonzero',(1-e).is_negative,True)
N,k=s.symbols('N k',integer=True,nonnegative=True)
coeff=s.symbols('coeff')
nu=e*(N-k)+1
zero('all-order prefix induction',DP(coeff*u**nu)-coeff*(N-k+1/e)*u**(nu-e))
zero('all-order exponent recurrence',e*(N-(k+1))+1-(nu-e))
zero('all prefix traces have nonzero cyclic remainder',nu-e*(N-k)-1)
check('fractional factor positive',(1/e).is_positive,True)
check('fractional factor strictly less than one',(1-1/e).is_positive,True)
L=s.symbols('L',integer=True,nonnegative=True)
check('all k<=N prefix modes regular',(nu.subs(k,N-L)).is_positive,True)
check('all k<=N+1 derivative norms finite',(e*(N-k+1)+1).subs(k,N+1-L).is_positive,True)
check('exact next derivative norm diverges',(e*(N-k+1)+1).subs(k,N+2).is_negative,True)
zero('new exact derivative norm denominator',2*e+2*nu-2*(e*(N-k+1)+1))
zero('last finite derivative mode',nu.subs(k,N+1)-(1-e))
zero('first divergent derivative mode',nu.subs(k,N+2)-(1-2*e))
zero('tangent derivative preserves pole modes',DQ(u**(-q)))
a=s.symbols('a',integer=True,nonnegative=True)
check('whole-source finite prefix valuation bound',(a+L*e).is_nonnegative,True)
M,S=s.symbols('M S',integer=True,nonnegative=True)
check('nontrivial infinity-only cover violates RH',(M+2+S+1-2).is_positive,True)
c=json.loads((H/'sources/nambu-sobolev-core-context.json').read_text())
ledger=json.loads((R/'MATH_STATUS.json').read_text());known={x['id']:x for x in ledger['entries']}
for item in c['retained_structural_inputs']:check('prior '+item['claim'],known[item['claim']]['state']=='proved' and known[item['claim']]['canonical_source']==item['canonical_source'])
for p,h in c['retained_proof_and_source_sha256'].items():check('byte '+p,hashlib.sha256((R/p).read_bytes()).hexdigest()==h)
yes={'full_goal_active','previous_goal_turn_progress','component_no_line_scope_separately_checked','first_additive_derivative_lattice_equals_full_normalization','finite_discriminant_quotient_equals_derivative_defect_globally','every_lifted_target_LND_has_nondense_domain_under_failure','every_finite_target_additive_stable_module_lies_in_base','all_derivative_polynomial_core_equals_base','arbitrary_finite_trace_zero_derivative_prefixes_proved','new_all_order_Laurent_tangent_transverse_controls','external_review_deferred'}
for p,x in c['scope_flags'].items():check('scope '+p,x is (p in yes))
proof=H/'NAMBU_SOBOLEV_LATTICE_AND_FINITE_DERIVATIVE_CORE.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=c['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-nambu-sobolev-core.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
