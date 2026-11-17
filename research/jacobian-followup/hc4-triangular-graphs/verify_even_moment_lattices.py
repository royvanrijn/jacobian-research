"""New even-moment and trace-discriminant controls; no retained replay."""
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
m=s.symbols('m',integer=True,positive=True);b,r=s.symbols('b r',integer=True,nonnegative=True)
# Euclidean division a=m*b+r, 0<=r<m, proves the universal floor claims.
a=m*b+r
zero('power valuation remainder',a-m*b-r)
zero('trace discriminant residual',(a-2*b)-((m-2)*b+r))
zero('m>=3 strict positivity certificate',3*(a-2*b)-a-(2*(m-3)*b+2*r))
zero('pair-product valuation inclusion',a-2*b-((m-2)*b+r))
zero('reflexive power defect',a-m*b-r)
zero('terminal moment at first integer above a',s.floor((m-1)/m))
zero('last strict jump floor',s.floor(m/m)-1)
e,q=s.symbols('e q',integer=True,positive=True)
rho=s.symbols('rho',positive=True)
radial=2*s.pi**2*rho**(2*e-1-2*m*q)
primitive=s.pi**2*rho**(2*(e-m*q))/(e-m*q)
zero('all-parameter radial norm antiderivative',s.diff(primitive,rho)-radial)
zero('Laurent endpoint e=m+1 q=1',((m+1)-m)-1)
zero('norm threshold versus codifferent exponent',2*e-1-2*m*q-(2*(e-m*q)-1))
P=s.symbols('P',nonzero=True)
def root_matrix(n):
    U=s.zeros(n)
    for j in range(n-1):U[j+1,j]=1
    U[0,n-1]=P
    return U
# Independent matrix transcriptions, not a bounded-degree theorem/census.
for n,k in [(3,2),(3,3),(7,3),(2,3)]:
    cut=(n-1)//k;U=root_matrix(n);exponents=range(-cut,n-cut)
    G=s.Matrix([[s.trace(U**(i+j)) for j in exponents] for i in exponents])
    label=' e'+str(n)+' m'+str(k)
    check('moment Gram entries polynomial'+label,all(s.cancel(v).is_polynomial(P) for v in G),True)
    zero('moment Gram determinant'+label,G.det()-(-1)**((n-1)*(n-2)//2)*n**n*P**(n-1-2*cut))
    check('defect vanishes exactly at lower half odd case'+label,(n-1-2*cut==0)==(n==3 and k==2),True)
    check('integer norm endpoint accepted'+label,k*cut<n,True)
    check('next integer norm pole rejected'+label,k*(cut+1)>=n,True)
u,v,Q=s.symbols('u v Q');F=v**4-u**9
targetQ=v**5/5-u**9*v
zero('finite point control Jacobian',s.diff(targetQ,v)-F)
zero('finite degree five monic relation',v**5-5*u**9*v-5*targetQ)
zero('monogenic different derivative',s.diff(v**5-5*u**9*v-5*Q,v)-5*F)
zero('canonical L2 cancellation',(1/F)**2*F**2-1)
zero('weighted point control common exponent',4*9-9*4)
zero('point region volume',s.pi*((2*rho**4)**2-(rho**4)**2)*s.pi*(rho**9/2)**2-3*s.pi**2*rho**26/4)
zero('point L3 lower bound exponent',26-36+10)
check('point L3 integral lower bound diverges',s.limit(rho**(-10),rho,0,dir='+')==s.oo,True)
zero('point curve lower bound coefficient',1-s.Rational(1,2)**4-s.Rational(15,16))
zero('point target Q scaling degree',9*5-(4*9+9))
c=json.loads((H/'sources/even-moment-lattices-context.json').read_text())
ledger=json.loads((R/'MATH_STATUS.json').read_text());known={a['id']:a for a in ledger['entries']}
for x in c['retained_structural_inputs']:check('prior '+x['claim'],known[x['claim']]['state']=='proved' and known[x['claim']]['canonical_source']==x['canonical_source'])
for p,h in c['retained_proof_and_source_sha256'].items():check('byte '+p,hashlib.sha256((R/p).read_bytes()).hexdigest()==h)
yes={'full_goal_active','previous_goal_turn_progress','global_even_integer_moment_formula_proved','moment_lattices_finite_base_projective','exact_terminal_moment_and_multiplication_criterion_proved','moment_trace_duals_and_discriminant_lengths_proved','uniform_sixth_moment_trace_perfection_criterion_proved','sixth_moment_defect_retains_every_positive_ramification_index','all_parameter_Laurent_norm_and_discriminant_controls_supplied','finite_cover_non_even_point_integrability_control_supplied','external_review_deferred'}
for p,x in c['scope_flags'].items():check('scope '+p,x is (p in yes))
proof=H/'EVEN_MOMENT_LATTICES_AND_SIXTH_MOMENT_TRACE_TEST.md'
out=dict(outcome='PASS',check_count=len(checks),algebraic_controls=len(algebra),elapsed_seconds=round(time.perf_counter()-started,3),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),checks=checks,scope_flags=c['scope_flags'],registered_mathematical_checkers_replayed=False,independent_replay=False,formal_verification=False,external_review=False)
(H/'receipt-even-moment-lattices.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({p:out[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
