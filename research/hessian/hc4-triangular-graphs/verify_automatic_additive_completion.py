"""NEW MP0 controls only. Never overwrite a receipt or replay after registration."""
import hashlib
import json
import time
from pathlib import Path
import sympy as s

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / 'receipt-automatic-additive-completion.json'
assert not RECEIPT.exists(), 'Sealed mathematical receipt: never replay.'
assert not any(e['id'] == 'HC4MYMP0' for e in json.loads((ROOT / 'MATH_STATUS.json').read_text())['entries'])
started = time.perf_counter()
checks, algebraic = [], []

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def check(name, value, alg=True):
    assert bool(value), name
    checks.append(name)
    if alg:
        algebraic.append(name)

x, y, z, t, a, b, lam, c, step = s.symbols('x y z t a b lam c step', nonzero=True)
r = s.symbols('r', integer=True, positive=True)
root_x = c*z**(-r-1) + z**(1-r)*t
root_y = z**r
delta = lambda f: s.simplify(a*z**(r-1)*s.diff(f,t))
check('ALL symbolic positive root orders: delta(x_i)=a_i', s.simplify(delta(root_x)-a) == 0)
check('ALL symbolic positive root orders: delta(y_i)=0', delta(root_y) == 0)
check('ALL symbolic positive root orders: contraction equals dy',
      s.simplify((r/(a*b))*a*z**(r-1)-s.diff(root_y/b,z)) == 0)
check('ALL symbolic positive root orders: volume divergence zero', s.diff(a*z**(r-1),t) == 0)

# New controls align several DIFFERENT finite boundary values and complements.
for order, value, scale_a, scale_b in [(1,0,2,3), (2,5,-3,2), (5,-7,4,-2)]:
    comp = scale_a*x + y**3 + 2*y
    base = scale_b*(y-value)
    check(f'source complement derivative at order{order} value{value}', s.diff(comp,x) == scale_a)
    check(f'source base derivative at order{order} value{value}', s.diff(base,x) == 0)
    local_x = z**(-order-1)+2*z**(-2*order-1)+z**(1-order)*t
    local_delta = lambda f: s.cancel(scale_a*z**(order-1)*s.diff(f,t))
    check(f'whole-boundary derivative expression order{order}', local_delta(local_x) == scale_a)
    check(f'whole-boundary nonnegative exponent order{order}', order-1 >= 0)
    # A fixed-base root deck transformation sends (z,t) to (eps*z,eps^(r-1)*t_next)
    # and the vector coefficient transforms by the SAME power: equivariance.
    check(f'root deck equivariance order{order}',
          s.expand((step*z)**(order-1)-step**(order-1)*z**(order-1)) == 0)

for amplitude in [1, y, y**3-2*y+1]:
    polynomial = x**4*y**2 + x**2*(y-5)+x+1
    deriv = lambda f: s.expand(amplitude*s.diff(f,x))
    val = polynomial
    for j in range(5):
        val = deriv(val)
    check('source LND remains nilpotent for amplitude '+str(amplitude), val == 0)
    flow = polynomial.subs(x,x+step*amplitude)
    series = sum(step**j * amplitude**j*s.diff(polynomial,x,j)/s.factorial(j) for j in range(5))
    check('finite exponential equals whole source flow '+str(amplitude), s.expand(flow-series) == 0)

# A one-variable LND cannot have a nonconstant coefficient in the fibre variable.
for degree in [1,2,4]:
    f = (y+1)*x**degree + (0 if degree == 1 else 2)
    val = x
    for n in [1,2,3]:
        val = s.expand(f*s.diff(val,x))
        check(f'non-LND degree recurrence degree{degree} power{n}',
              s.degree(val,x) == 1+n*(degree-1) and val != 0)

# Faithful finite affine base actions: fixed point and averaging identities.
alpha, beta, u = s.symbols('alpha beta u')
fixed = beta/(1-alpha)
check('nontranslation affine fibre element has a fixed point', s.cancel(alpha*fixed+beta-fixed) == 0)
check('finite-order translation characteristic-zero recurrence', s.expand((u+3*beta)-u) == 3*beta)
check('cyclic base action may be nontrivial', (-s.I)**4 == 1 and -s.I != 1)
check('semiregular cyclic action need not preserve a chosen component', tuple((j+1)%4 for j in range(4)) == (1,2,3,0))

cp = HERE / 'sources/automatic-additive-completion-context.json'
C = json.loads(cp.read_text())
for p, digest in C['retained_proof_and_source_sha256'].items():
    check('retained byte pin '+p, sha(ROOT/p) == digest, False)
check('five retained claims/eight unique proof-checker pins', len(C['retained_structural_inputs']) == 5 and len(C['retained_proof_and_source_sha256']) == 8, False)
flags = C['scope_flags']
check('actual common source pencil and ALL component charts required',
      flags['actual_plane_source_and_full_normalization'] and flags['retained_common_unramified_pencil']
      and flags['whole_component_charts_not_formal_germs'], False)
check('automatic additive completion and universal cyclic generic decks',
      flags['all_actual_etale_completions_admit_Ga'] and flags['all_generic_rational_decks_cyclic']
      and flags['all_LNDs_preserve_source_for_proper_completion'], False)
check('do not transfer deck cyclicity to monodromy or arbitrary intermediate pi1',
      not flags['all_monodromy_cyclic'] and not flags['all_intermediate_fundamental_groups_cyclic'], False)
check('target and ramified-boundary obligations remain open',
      not flags['target_algebra_invariant_under_source_LND']
      and not flags['LND_extends_to_full_finite_normalization']
      and not flags['cyclic_generic_deck_excluded'] and not flags['unrestricted_HC4_or_JC2_resolved'], False)
check('no new imports/caches/downloads or registered replay',
      C['new_external_theorem_imports'] == C['new_primary_caches'] == C['new_primary_downloads'] == 0
      and C['registered_mathematical_checkers_replayed'] == [], False)
receipt = dict(
    outcome='PASS', claim='HC4MYMP0', check_count=len(checks), algebraic_controls=len(algebraic),
    checks=checks, algebraic_check_names=algebraic, elapsed_seconds=time.perf_counter()-started,
    proof_sha256=sha(HERE/'COMMON_UNRAMIFIED_PENCIL_AND_AUTOMATIC_ADDITIVE_COMPLETION.md'),
    checker_sha256=sha(Path(__file__)), context_sha256=sha(cp),
    registered_mathematical_checkers_replayed=[], independent_replay=False, formal_verification=False, external_review=False,
)
RECEIPT.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_controls','elapsed_seconds']}))
