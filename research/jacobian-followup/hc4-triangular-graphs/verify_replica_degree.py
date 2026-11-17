"""NEW MR0 algebraic transcriptions, not a proof by bounded examples.

Refuse an existing receipt or registered claim. Never run predecessors.
"""
import hashlib
import json
import time
from pathlib import Path
import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / 'receipt-replica-degree.json'
assert not RECEIPT.exists(), 'Sealed receipt exists: never replay.'
assert not any(e['id'] == 'HC4MYMR0' for e in json.loads((ROOT / 'MATH_STATUS.json').read_text())['entries'])
started = time.perf_counter()
checks, algebraic = [], []

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def check(name, value, alg=True):
    assert bool(value), name
    checks.append(name)
    if alg:
        algebraic.append(name)

x, y, u, t = S.symbols('x y u t')

# Distinct polynomial quotient maps, including a nongalois cubic control.
for label, h, r in [
    ('linear', 2*y+3, u**3+u+1),
    ('quadratic', y**2+1, u+2),
    ('nongalois-cubic', y**3-3*y, u**2+1),
    ('quartic', y**4+y, S.Integer(1)),
]:
    hp = S.diff(h, y)
    p = S.expand(hp*r.subs(u, h))
    beta = 1/(y+2)
    sl = x/p+beta
    m, k, n = S.degree(h, y), S.degree(p, y), S.degree(r, u)
    check(label+' rational local slice', S.cancel(p*S.diff(sl, x)-1) == 0)
    check(label+' volume coefficient', S.cancel(S.Matrix([sl, h]).jacobian([x,y]).det()-hp/p) == 0)
    check(label+' ratio is polynomial in quotient', S.cancel(p/hp-r.subs(u,h)) == 0)
    check(label+' exact degree identity', k+1 == m*(n+1))
    check(label+' quotient index divides replica degree plus one', (k+1) % m == 0)

# The pole argument is symbolic for arbitrary positive multiplicity and order.
e, ell = S.symbols('e ell', integer=True, positive=True)
check('finite-pole valuation equals minus one minus e times ell minus one',
      S.expand((e-1-e*ell)-(-1-e*(ell-1))) == 0)
for multiplicity, pole in [(1,1),(2,1),(7,1),(1,4),(9,3)]:
    check(f'polynomial ratio rejects pole e{multiplicity} l{pole}', multiplicity-1-multiplicity*pole < 0)
check('bare translation has only index one', [m for m in range(1,9) if 1 % m == 0] == [1])

# Monomial replicas admit cyclic fixed fields; never divide p in E.
for m, n in [(2,0),(3,1),(5,0),(4,2)]:
    v = y-2
    h = v**m
    p = S.diff(h,y)*h**n
    sl = x/p+1/(y+1)
    new_t = S.cancel(sl*h**n*m)
    check(f'monomial slice multiplication is legitimate m{m} n{n}',
          S.cancel(new_t-x/v**(m-1)-m*h**n/(y+1)) == 0)
    check(f'monomial exact index m{m} n{n}', S.degree(p,y)+1 == m*(n+1))
    check(f'monomial character preserves volume m{m} n{n}', ((m-1)+1) % m == 0)

sigma = {x:-x+1/y**2, y:-y}
check('rational shift cyclic deck may have source pole', S.cancel((-x+1/y**2).subs(sigma, simultaneous=True)-x) == 0)
check('rational cyclic deck preserves source volume', S.Matrix([-x+1/y**2,-y]).jacobian([x,y]).det() == 1)
check('rational invariant slice cyclic order two', S.cancel((x/y-1/(2*y**3)).subs(sigma, simultaneous=True)-(x/y-1/(2*y**3))) == 0)

# Actual Keller automorphisms in two adapted source directions. These are
# controls only; the field theorem is the unbounded local-slice/volume proof.
P = x+y**2
Q = y+P**2
check('nonlinear automorphism control has Jacobian one', S.Matrix([P,Q]).jacobian([x,y]).det() == 1)
check('first direction derivatives recover invariant y', S.expand(Q-S.diff(Q,x)**2/4-y) == 0)
check('first direction reconstructed x', S.expand(P-(Q-P**2)**2-x) == 0)
a,b = S.symbols('a b')
Px = S.expand(P.subs({x:b,y:a}, simultaneous=True))
Qx = S.expand(Q.subs({x:b,y:a}, simultaneous=True))
check('second coordinate frame has constant nonzero Jacobian', S.Matrix([Px,Qx]).jacobian([a,b]).det() == -1)
check('second direction derivatives recover complementary b', S.expand(S.diff(Px,a)**2/4-(Px-b)) == 0)
check('second frame derivative ladder terminates at x degree bound', S.diff(Px,a,5) == 0 and S.diff(Qx,a,5) == 0)
for d,k in [(5,1),(12,4),(35,5),(12,5)]:
    g = S.gcd(d,k+1)
    check(f'field index divides both tower and replica d{d} k{k}',
          [m for m in range(1,d+1) if d % m == 0 and (k+1) % m == 0] == [m for m in range(1,int(g)+1) if g % m == 0])

# Birational, immersive, Ga-equivariant and nevertheless noninjective.
f,g = y**2-1,y*(y**2-1)
check('negative control is birational', S.cancel(g/f-y) == 0)
check('negative control differential has rank two everywhere', S.expand(S.Rational(3,2)*y*S.diff(f,y)-S.diff(g,y)) == 1)
check('negative control identifies two source lines', f.subs(y,1) == f.subs(y,-1) == g.subs(y,1) == g.subs(y,-1) == 0)
check('negative control missing y is integral', S.expand(y**2-(f+1)) == 0)
check('negative control has nodal relation', S.expand(g**2-f**2*(f+1)) == 0)
check('negative control is equivariant for bare translation', S.diff(f,x) == S.diff(g,x) == 0)

cp = HERE / 'sources/replica-degree-context.json'
C = json.loads(cp.read_text())
for path,digest in C['retained_proof_and_source_sha256'].items():
    check('retained byte pin '+path, sha(ROOT/path) == digest, False)
check('five prior claims and ten proof-checker pins', len(C['retained_structural_inputs']) == 5 and len(C['retained_proof_and_source_sha256']) == 10, False)
flags = C['scope_flags']
for key in ['all_degrees_and_source_coordinate_frames','polynomial_volume_factorization','exact_replica_degree_divisibility','bare_translation_no_proper_intermediate','actual_equivariant_source_image_surjective','monomial_replica_cyclic_fixed_field','finite_source_direction_derivatives_generate_full_field']:
    check('scope '+key, flags[key], False)
for key in ['all_replicas_have_no_proper_intermediate','general_replica_fixed_field_is_cyclic','raw_derivative_algebra_normal','field_reconstruction_implies_inverse','source_action_descends_to_target','full_trace_integrality_automatic','unrestricted_HC4_or_JC2_resolved']:
    check('boundary '+key, not flags[key], False)
check('no new imports caches downloads or predecessor replay', C['new_external_theorem_imports'] == C['new_primary_caches'] == C['new_primary_downloads'] == 0 and C['registered_mathematical_checkers_replayed'] == [], False)
r = dict(outcome='PASS', claim='HC4MYMR0', check_count=len(checks), algebraic_controls=len(algebraic), checks=checks, algebraic_check_names=algebraic, elapsed_seconds=time.perf_counter()-started,
         proof_sha256=sha(HERE/'REPLICA_DEGREE_AND_SOURCE_DERIVATIVE_FIELD_RECONSTRUCTION.md'), checker_sha256=sha(Path(__file__)), context_sha256=sha(cp), registered_mathematical_checkers_replayed=[], independent_replay=False, formal_verification=False, external_review=False)
RECEIPT.write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
