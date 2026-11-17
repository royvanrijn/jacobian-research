"""Exact local controls and immutable input audit for the written classification.

This checker does not prove the source-layer/IC theorems or the global
classification. It checks their local algebra transcriptions and source locks.
Run once before registering;never replay registered mathematical checkers.
"""
import hashlib,json,time
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
started=time.monotonic();checks=[]
def check(ok,label):
 assert ok,label
 checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
u,a,b,c,lam,lam1,lam2,v,z=s.symbols('u a b c lam lam1 lam2 v z')

# A diagonal W-submodule need not be stable under the coefficient algebra.
w=s.Matrix([1,c]);T=s.diag(a,b)
check(s.expand(s.det(s.Matrix.hstack(w,T*w))-c*(b-a))==0,'a mixed same-image slope is invariant only if c*(b-a)=0')
check(s.det(s.Matrix.hstack(w,s.diag(1,0)*w))==-c,'one actual coefficient projector destroys every nonzero mixed slope')
check(s.diag(1,0)**2==s.diag(1,0),'the first unramified prime projector is idempotent')
check(s.diag(0,1)**2==s.diag(0,1),'the second unramified prime projector is idempotent')
check(s.diag(1,0)*s.diag(0,1)==s.zeros(2),'the distinct prime projectors are orthogonal')
check(s.diag(1,0)+s.diag(0,1)==s.eye(2),'prime projectors sum to the whole same-image coefficient space')

# A residue extension can have target character summands;its FIELD action
# does not preserve an arbitrary proper character line.
M=s.Matrix([[0,v],[1,0]])
check(M*M==v*s.eye(2),'quadratic residue-field multiplication obeys alpha squared=v')
check(s.det(z*s.eye(2)-M)==z*z-v,'the coefficient-field action has the irreducible quadratic generic polynomial')
check(s.Poly(z*z-v,z,domain=s.QQ.frac_field(v)).is_irreducible,'generic quadratic residue multiplication has no rational invariant line')
check(M*s.Matrix([0,1])==s.Matrix([v,0]),'multiplication exits the trace-zero character line,so character choice alone is insufficient')

# Normal derivative is essential:regular coefficients alone leave bounded
# pole-order submodules;in characteristic zero u^-1 generates all orders.
for k in [0,1,2,5,11]:
 check(s.diff(u**-1,u,k)==(-1)**k*s.factorial(k)*u**(-k-1),'principal-parts normal derivatives recover pole order'+str(k+1))
for m in [1,2,7]:
 check(s.simplify(u**(m-1)*u**(-m)-u**-1)==0,'a leading pole reduces to u^-1 at order'+str(m))
check(s.diff(u**-1,u)==-u**-2,'the pole-order-one coefficient submodule is NOT derivative stable')
check((-2)%2==0,'the characteristic-zero restriction matters:the derivative of u^-2 vanishes in characteristic2')

# Generic linear primitive and valuation exceptions are finite,not sampled
# evidence for actual Keller existence or an all-degree exclusion.
V=s.Matrix([[1,lam1],[1,lam2]])
check(s.det(V)==lam2-lam1,'two different linear parameters distinguish unequal embedding pairs')
check(s.Poly(a+lam*b,lam).degree()==1,'a nonzero leading-pole cancellation condition has at most one parameter')
check(s.solve(a+lam*b,lam)==[-a/b],'the unique potential cancellation parameter is explicit over the leading residue field')
check(s.Poly(a+lam*b,lam).diff().as_expr()==b,'two cancellation parameters force both leading residues zero')

# Local nonsingular-order identity,credited to Makar-Limanov Section2.
Pi,Pt,ai,ay=s.symbols('Pi Pt ai ay')
check(s.expand((-ai*(-Pi/Pt)+ay)-(ai*Pi+ay*Pt)/Pt)==0,'the hypersurface derivative identity recovers reciprocal P_t only under the global nonsingular-order premise')
algebra=len(checks)

context=json.loads((HERE/'sources/differential-subalgebra-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check(sha(ROOT/p)==h,'immutable retained input '+p)
for p,h in context['new_primary_sha256'].items():check(sha(HERE/p)==h,'inspected primary '+p)
check(context['newly_imported_proved_external_statements']==[],'no new proved external theorem is silently imported')
check(context['newly_used_unproved_external_statements']==[],'neither Makar-Limanov conjecture nor tensor-map surjectivity is imported')
check(context['finite_generation_assumed_for_input_subalgebra'] is False,'finite generation is a conclusion of the written classification')
check(context['last_HC4_source_to_regular_inverse_Gram_implication_proved'] is False,'the full HC4 last implication remains explicitly OPEN')
proof=HERE/'DIFFERENTIAL_SUBALGEBRA_CLASSIFICATION.md'
flags={k:v for k,v in context.items() if isinstance(v,bool)}
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-started,3),proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],**flags)
(HERE/'receipt-differential-subalgebras.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
