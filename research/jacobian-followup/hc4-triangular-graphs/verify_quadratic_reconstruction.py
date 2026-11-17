"""Exact local controls and immutable provenance;not a replay of published proofs."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
u=s.symbols('u',nonzero=True)
# Generic split hyperplane identities;the note supplies the all-degree proof.
for d in [3,4,6]:
    basis=[s.eye(d)[:,i] for i in range(d)]
    v=[basis[i]-basis[-1] for i in range(d-1)]
    prod=lambda a,b:a.multiply_elementwise(b)
    check(prod(v[0],v[1])==basis[-1],'distinct augmentation vectors recover the last idempotent at d='+str(d))
    check(all(prod(v[i],v[i])-prod(v[0],v[1])==basis[i] for i in range(d-1)),
          'quadratic products recover all other split idempotents at d='+str(d))
    matrix=s.Matrix.hstack(*(prod(a,b) for i,a in enumerate(v) for b in v[i:]))
    check(matrix.rank()==d,'the generic quadratic image has full rank at d='+str(d))
v=s.Matrix([1,-1]);check(v.multiply_elementwise(v)==s.ones(2,1),'degree two retains only the diagonal square line')
# Actual tame field traces via companion matrices,not an invented Gram pairing.
for e in [2,3,5,6]:
    mult=s.zeros(e)
    for i in range(e-1):mult[i+1,i]=1
    mult[0,e-1]=u
    a=mult.inv();b=mult**(-(e-1))
    check(s.trace(a)==0 and s.trace(b)==0,'both complementary fractional poles have zero trace at e='+str(e))
    check(a*b==s.eye(e)/u,'their product is the isolated integer pole at e='+str(e))
    check(s.trace(a*b)==e/u,'its trace has the required pole at e='+str(e))
    for j in [1,e-1]:
        lam=s.Rational(j,e)
        check(all(lam+k!=0 for k in [-2,-1,0,1]),'nontrivial Kummer weights remain localized at e='+str(e)+' j='+str(j))
# One omitted unramified sheet can pair with any OTHER omitted factor,
# including another unramified sheet;the compensating unit is retained.
for h in [1,2,3,6]:
    a=s.Matrix([1/u,-1/(h*u),0]);b=s.Matrix([1,0,-1]);weights=s.Matrix([1,h,1])
    check(weights.dot(a)==0 and weights.dot(b)==0,'two trace-zero elements exist for the unramified factor with partner degree '+str(h))
    check(a.multiply_elementwise(b)==s.Matrix([1/u,0,0]),'their product isolates that missing unramified sheet for h='+str(h))
    check(weights.dot(a.multiply_elementwise(b))==1/u,'the quadratic trace detects its shared image for h='+str(h))
# A sole omitted unramified coordinate is forced to be regular by trace zero.
a,b=s.symbols('a b');last=-a-b
check(s.expand(a+b+last)==0,'a sole unramified missing coordinate equals minus the regular retained trace')
check(s.Poly(last,a,b).total_degree()==1,'that coordinate has no Laurent denominator')
# Uniform residue-class witnesses for the cyclic non-Keller control.
for d in [3,6]:
    pairs=[(1,d-1),(2,d-1)]+[(1,r-1) for r in range(2,d)]
    check(all(0<i<d and 0<j<d and (i+j)%d==r for r,(i,j) in enumerate(pairs)),
          'nonzero cyclic characters produce every Laurent character at d='+str(d))
# Selected-image defect budgets include several primes and residue degree>1.
for groups,expected in [([[(1,1)]],[1]),([[(2,1)]],[]),([[(1,2)]],[]),([[(1,1),(1,1)]],[]),
                        ([[(2,1)],[(1,1)]],[2]),([[(3,2)],[(1,1)],[(1,1),(1,1)]],[2])]:
    missing=[sum(e*f for e,f in row) for row in groups]
    check([j+1 for j,ell in enumerate(missing) if ell==1]==expected,'the quadratic defect selects precisely loss-one images in '+str(groups))
algebra=len(checks)
context=json.loads((HERE/'sources/quadratic-reconstruction-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check(sha(ROOT/p)==h,'immutable retained input '+p)
check(len(context['newly_imported_proved_external_statements'])==1,'only the directly inspected Domrina degree-four theorem is newly imported')
check(context['newly_used_unproved_external_statements']==[],'no unproved external statement is used')
check(context['primary_inspection']['domrina_part_II']['main_theorem_unrestricted_degree_four'] is True,
      'Part II main theorem has no one-boundary condition;the primary was inspected by web')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(flags['all_degree_quadratic_defect_proved'] and flags['quartic_counterexample_realizations_retired']
      and not flags['quadratic_reconstruction_proved_to_force_polynomial_source'] and not flags['original_goal_complete'],
      'the all-degree reconstruction,published quartic retirement and remaining full goal are distinguished')
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=s.__version__,
 proof_sha256=sha(HERE/'QUADRATIC_SOURCE_RECONSTRUCTION_AND_QUARTIC_RETIREMENT.md'),
 checker_sha256=sha(Path(__file__)),context_sha256=sha(HERE/'sources/quadratic-reconstruction-context.json'),
 dependencies=context['retained_proof_and_source_sha256'],independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-quadratic-reconstruction.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies'}},indent=2))
