"""Exact controls for the written trace theorem; never replays prior checkers."""
import hashlib
import json
import platform
import time
from pathlib import Path
import sympy as s

START = time.monotonic()
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
checks = []


def check(value, label):
    assert bool(value), label
    checks.append(label)


def zero(value, label):
    check(s.cancel(value) == 0, label)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


T, t = s.symbols('T t')
# Newton's identities for unrestricted symbolic coefficients, with the
# constant coefficient deliberately withheld from the supplied traces.
for degree in [1, 2, 3, 5, 7]:
    cs = s.symbols('c1:'+str(degree+1))
    M = s.zeros(degree)
    for j in range(degree-1):
        M[j+1, j] = 1
    for j in range(degree):
        M[j, degree-1] = -cs[degree-1-j]
    ps = [None]+[s.trace(M**k) for k in range(1, degree)]
    recovered = [s.Integer(1)]
    for k in range(1, degree):
        recovered.append(s.expand(-sum(recovered[k-i]*ps[i] for i in range(1, k+1))/k))
    check(recovered[1:] == list(cs[:-1]), 'Newton recovery omits only the constant coefficient at degree '+str(degree))
    f = T**degree+sum(cs[j]*T**(degree-1-j) for j in range(degree))
    check(s.Poly(f, T).all_coeffs()[-1] == cs[-1], 'evaluation at an integral sheet determines the last coefficient at degree '+str(degree))

# A split local root model with one integral root and q negative roots.
# These are companion multiplication matrices, not Keller maps.
for q in [1, 2, 3, 5, 8]:
    M = s.zeros(q)
    for j in range(q-1):
        M[j+1, j] = 1
    M[0, q-1] = t**(-q)
    check(M**q == t**(-q)*s.eye(q), 'negative root cluster has q-th powers t^-q at q='+str(q))
    check(all(s.trace(M**k) == 0 for k in range(1, q)), 'all earlier power sums cancel at q='+str(q))
    zero(s.trace(M**q)-q*t**(-q), 'the last allowed trace detects a genuine pole at q='+str(q))
    check(s.det(M) == (-1)**(q-1)*t**(-q), 'negative conjugate product is not integral at q='+str(q))

x, y, z = s.symbols('x y z')
u = 1+x*y
A = u**3*z+y**2*u*(4+3*x*y)
B = y+3*x*u**2*z+3*x*y**2*(4+3*x*y)
C = 2*x-3*x*x*y-x**3*z
zero(s.det(s.Matrix([A, B, C]).jacobian([x, y, z]))+2, 'the actual whole-source three-dimensional map has constant Jacobian minus two')
for point in [(0, 0, -s.Rational(1, 4)), (1, -s.Rational(3, 2), s.Rational(13, 2)), (-1, s.Rational(3, 2), s.Rational(13, 2))]:
    check(tuple(s.expand(f).subs(dict(zip((x, y, z), point))) for f in [A, B, C]) == (-s.Rational(1, 4), 0, 0), 'the known Keller map has the exact repeated image at '+str(point))
gamma = 2-3*x*y-x*x*z
w = u*gamma
zero(B*C-2*gamma-4*w+3*w*w, 'the weighted first target invariant has the stated identity')
zero(A*C*C-w*w-w*gamma+w**3, 'the weighted second target invariant has the stated identity')
zero(w**3-2*w*w+B*C*w-2*A*C*C, 'the source integral polynomial satisfies its global monic cubic')
P, Q, W, cc = s.symbols('P Q W cc')
G = (P-4*W+3*W**2)/2
xx = cc/G
yy = (W/G-1)/xx
zz = (5-3*W/G-G)/xx**2
f = W**3-2*W**2+P*W-2*Q
for output, expected, label in [(A, Q/cc**2, 'first'), (B, P/cc, 'second'), (C, cc, 'third')]:
    num = s.cancel(output.subs({x: xx, y: yy, z: zz}, simultaneous=True)-expected).as_numer_denom()[0]
    check(s.rem(num, f, W) == 0, 'generic cubic-root reconstruction recovers the '+label+' target coordinate')
check(s.Poly(f, Q).degree() == 1 and s.diff(f, Q) == -2, 'the generic hidden equation is linear with a unit coefficient in the independent target Q')

# Exact field traces in F[W]/(W^3-2W^2+P W-2Q).
M = s.Matrix([[0, 0, 2*Q], [1, 0, -P], [0, 1, 2]])
I = s.eye(3)
GM = (P*I-4*M+3*M*M)/2
U = (M*GM.inv()).applyfunc(s.cancel)
c = P**3-P**2-18*P*Q+27*Q**2+16*Q
b = P**2-12*Q
check(M**3-2*M*M+P*M-2*Q*I == s.zeros(3), 'the actual regular representation has its monic cubic')
check(s.trace(M) == 2 and s.trace(M*M) == 4-2*P, 'an integral source polynomial has both regular power traces in the noninvertible map')
check((c*U**3+b*U-4*Q*I).applyfunc(s.cancel) == s.zeros(3), 'the actual source polynomial u has the stated denominator-bearing cubic')
zero(s.trace(U), 'the nonintegral source polynomial has complete first-trace cancellation')
zero(s.trace(U*U)+2*b/c, 'its second trace has the exact rational pole formula')
check(c.subs({P: 0, Q: -s.Rational(16, 27)}) == 0 and b.subs({P: 0, Q: -s.Rational(16, 27)}) == s.Rational(64, 9), 'the second trace denominator has a zero at which its numerator stays nonzero')
check(s.cancel(s.trace(U*U)+b/c) != 0, 'a hostile missing-factor-two trace mutation fails')
gram = s.Matrix(3, 3, lambda i, j: s.trace(M**(i+j)))
check(s.det(gram) != 0, 'the separable cubic field trace pairing is nondegenerate')

# Differential trace pairing control on a quadratic function field.
a0, a1 = s.Function('a0')(t), s.Function('a1')(t)
N = s.Matrix([[0, t], [1, 0]])
AA = a0*s.eye(2)+a1*N
dAA = s.diff(a0, t)*s.eye(2)+(s.diff(a1, t)+a1/(2*t))*N
zero(s.diff(s.trace(AA), t)-s.trace(dAA), 'the extended field derivative commutes with trace')
zero(s.diff(s.trace(AA*N), t)-s.trace(dAA*N+AA*N/(2*t)), 'differential trace pairing has the full Leibniz term')
check(s.Matrix([[s.trace(s.eye(2)), s.trace(N)], [s.trace(N), s.trace(N*N)]]).det() == 4*t, 'trace duality detects both coefficients after field localization')
algebra = len(checks)
context = json.loads((HERE/'sources/whole-source-trace-context.json').read_text())
for path, digest in context['retained_proof_and_source_sha256'].items():
    check(sha(ROOT/path) == digest, 'immutable retained input '+path)
for entry in context['newly_cached_primary_context']:
    check(sha(HERE/'sources'/entry['retained_file']) == entry['sha256'], 'inspected primary context '+entry['retained_file'])
check(context['newly_imported_proved_external_statements'] == [] and context['newly_used_unproved_external_statements'] == [], 'all new implications are written proofs; no external splitting theorem or trace-regularity premise is imported')
flags = {k: v for k, v in context.items() if isinstance(v, bool)}
check(flags['whole_source_trace_regularity_proved_automatic'] is False and flags['unrestricted_HC4_or_JC2_resolved'] is False and flags['original_goal_complete'] is False, 'the global trace-regularity gap and full original goal remain open')
proof = HERE/'WHOLE_SOURCE_POWER_TRACES_AND_DIFFERENTIAL_SPLITTING.md'
paths = ['sources/whole-source-trace-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result = dict(outcome='PASS', check_count=len(checks), checks=checks, algebraic_and_scope_controls=algebra,
              elapsed_seconds=round(time.monotonic()-START, 3), python=platform.python_version(), sympy=s.__version__,
              proof_sha256=sha(proof), checker_sha256=sha(Path(__file__)), dependencies=context['retained_proof_and_source_sha256'],
              source_sha256={p: sha(HERE/p) for p in paths}, independent_replay=False, formal_verification=False, external_review=False, **flags)
(HERE/'receipt-whole-source-traces.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k: v for k, v in result.items() if k not in {'checks', 'dependencies', 'source_sha256'}}, indent=2))
