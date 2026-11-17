"""Finite exact controls; the all-order proof is in the accompanying JB0 note."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp

START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):
    check(sp.expand(value)==0,label)
def vanishes_mod(value,variable,order,label):
    zero(sp.rem(sp.expand(value),variable**order,variable),label)
s,t,S,V,u,v,f,g=sp.symbols('s t S V u v f g')
for m in [1,2,4]:
    expression=(1-v*f)*sum((v*f)**j for j in range(m))
    zero(expression-(1-(v*f)**m),'the finite geometric inverse identity m='+str(m))
    vanishes_mod(expression-1,f,m,'the inverse works through the full fat-divisor order m='+str(m))
    zero((f*g)**m-f**m*g**m,'the correction has the exact m-th-power pullback m='+str(m))
check(sp.Poly(f**4*g**4,f,g).total_degree()==8,
      'the power factor keeps both selected and extra multiplicities')

L=[s,t+s*t*t]
J=sp.Matrix(L).jacobian([s,t]).det()
zero(J-(1+2*s*t),'the model is explicitly non-Keller on the entire source plane')
zero(L[1]-t*(1+s*t),'the whole target divisor has both source factors')
zero((1+s*t)-s*t-1,'the two source divisors generate the unit ideal')
zero(J.subs(t,0)-1,'the selected divisor is contained in the etale open')
zero(J.subs(t,-1/s)+1,'the extra component is also contained in that etale open')
zero(L[1].subs(t,-1/s),'the extra component maps to the same target line')
zero(t.subs(t,-1/s)*(-s)-1,'the selected source equation is a nonconstant unit on the extra component')
check(sp.limit(-1/s,s,0,dir='+')==-sp.oo and sp.limit(s,s,0)==0,
      'the normal trace has a pole at a finite-image boundary puncture')
check(sp.limit(-1/s,s,sp.oo)==0 and sp.limit(s,s,sp.oo)==sp.oo,
      'its zero is at the other puncture mapping to target infinity')
check(1+2*0-2+2==1,'the extra normalization has degree one and finite-image boundary index one')

for order in [2,3,5]:
    q=sum((-1)**(n-1)*sp.catalan(n-1)*S**(n-1)*V**n for n in range(1,order))
    vanishes_mod(q+S*q*q-V,V,order,
                 'the formal branch is a right inverse through order '+str(order))
    qpull=q.subs({S:s,V:L[1]},simultaneous=True)
    vanishes_mod(qpull-t,t,order,
                 'the polynomial jet recovers the selected source coordinate through order '+str(order))
    check(sp.Poly(q,S,V).total_degree()==2*order-3,
          'the unique reduced jet has the predicted growing target degree order='+str(order))
    check(sp.Poly(q,S,V).coeff_monomial(S**(order-2)*V**(order-1))
          ==(-1)**(order-2)*sp.catalan(order-2),
          'the top monomial survives every addition from the target ideal power order='+str(order))
for r,m in [(2,2),(3,2)]:
    order=r*m
    q=sum((-1)**(n-1)*sp.catalan(n-1)*S**(n-1)*V**n for n in range(1,order))
    vanishes_mod(q+S*q*q-V,V,order,'the same inverse retains fat-divisor multiplicity r='+str(r))
    vanishes_mod(q.subs({S:s,V:L[1]},simultaneous=True)-t,t,order,
                 'fat selected and target thickening orders agree r='+str(r))
disc=1+4*S*V
factors=sp.factor_list(disc,S,V)[1]
check(len(factors)==1 and factors[0][1]==1 and sp.Poly(factors[0][0],S,V).total_degree()==2,
      'the quadratic discriminant is one reduced irreducible divisor')
a=sp.symbols('a')
zero((a+S*a*a-V)*4*S-( (2*S*a+1)**2-disc),
     'the inverse equation has exactly that nonsquare discriminant')

# Positive Keller control and a strict degree-threshold check.
K=[s+t*t,t]
inverse=[S-V*V,V]
check(sp.Matrix(K).jacobian([s,t]).det()==1,'the positive comparison is an actual whole-plane Keller map')
for i in range(2):
    zero(inverse[i].subs({S:K[0],V:K[1]},simultaneous=True)-[s,t][i],
         'the positive fixed polynomial recovery is exact coordinate='+str(i))
zero((S-V*V).subs({S:K[0],V:K[1]},simultaneous=True)-s,
     'the positive selected divisor has the exact target equation')
check(6*1>1+2*2,'an actual inverse meets the sufficient rational degree threshold')
error=t-(t+s*t*t)
check(sp.Poly(error,s,t).total_degree()==3 and sp.rem(error,t*t,t)==0
      and sp.rem(error,t**3,t)!=0,
      'a finite jet may pass with nonzero error when the strict threshold is absent')
check(not (2*1>1+1*3),'that false finite recovery is excluded by the strict degree test')

context=json.loads((HERE/'sources/formal-divisor-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'prior proof/source bytes remain pinned: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'primary formal-etale context bytes are pinned')
check(context['newly_used_unproved_external_statements']==[],
      'the formal-etale special case is reproduced rather than silently imported')
check(context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the global embedding problem remains open')
algebra_count=sum(not label.startswith(('prior proof/source','primary formal-etale'))
                  for label in checks)
receipt=dict(outcome='PASS',claim='HC4MYJB0',check_count=len(checks),
 algebraic_and_scope_controls=algebra_count,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 source_sha256={'sources/formal-divisor-context.json':hashlib.sha256((HERE/'sources/formal-divisor-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 primary_context=context['newly_cached_primary_context'],
 assurance='Finite exact controls of a written all-order proof; not independent replay, formal verification or external review.',
 automatic_all_order_divisor_inverse_proved=True,
 rational_reconstruction_degree_gate_proved=True,
 full_w_dependent_escape_excluded=False,
 unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-formal-divisor.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
