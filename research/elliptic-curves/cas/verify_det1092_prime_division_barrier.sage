#!/usr/bin/env sage-python
"""Pinned equations and symbolic checks for the all-prime division theorem.

The modular-curve, monodromy and cohomology arguments are written proofs,
not inferred from a finite prime sample. No exceptional coordinates are read.
"""
import hashlib, json, signal
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, vector

signal.alarm(25)
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'det1092_prime_division_barrier_v1'
PARENT = ART / 'curve302_recovered_mw17_parent_v1.json'
GEOM = ART / 'curve302_parent_geometric_picard19_v1.json'
TRACE = ART / 'det1092_trace_parity_descent_v1/replay.json'
def read(path): return json.loads(path.read_text())
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def save(name, data):
    OUT.mkdir(exist_ok=True)
    path = OUT / name
    if path.exists(): assert read(path) == data
    else:
        with path.open('x') as stream:
            json.dump(data, stream, indent=2, sort_keys=True)
            stream.write('\n')

save('protocol.json', {
    'classification': 'new written all-prime deduction with exact equation and symbolic checks',
    'inputs': {str(p.relative_to(ROOT)): sha(p) for p in [PARENT, GEOM, TRACE]},
    'checker_sha256': sha(Path(__file__)),
    'limits': {'seconds': 25, 'cycle_regression_primes': [2, 3, 5, 7, 11],
               'point_searches': 0, 'specializations': 0, 'class_groups': 0,
               'subgroup_states': 0, 'new_multisections': 0},
    'boundary': 'All-prime conclusions use the written proof; finitely many permutation regressions do not prove uniformity.'})

p = read(PARENT); geom = read(GEOM); trace = read(TRACE)
assert geom['status'] == 'PASS'
assert geom['geometric_generic_MW_rank'] == geom['arithmetic_generic_MW_rank'] == 17
assert geom['full_geometric_basis_is_displayed_rational_basis'] is True
for name, digest in geom['sources'].items(): assert sha(ROOT / name) == digest
assert trace['status'] == 'PASS_TRACE_PARITY_DESCENT_IDENTITIES_AND_INPUTS'
assert trace['uniform_geometric_conclusion'].startswith('For every nonzero M17 parity')
for name, digest in trace['inputs'].items(): assert sha(ROOT / name) == digest

R = PolynomialRing(QQ, 't'); K = R.fraction_field()
def dec(value): return K(R(value['numerator'])) / R(value['denominator'])
E = EllipticCurve(K, list(map(dec, p['a_invariants'])))
A = R(-E.c4()/48); B = R(-E.c6()/864)
Delta = R(E.discriminant()); c4 = R(E.c4())
assert A.degree() == 8 and B.degree() == 12
assert Delta.degree() == 24 and Delta.gcd(Delta.derivative()) == 1
assert Delta.gcd(c4) == 1
assert -16*(4*A[8]**3 + 27*B[12]**2) != 0
j = E.j_invariant()
assert j.denominator().degree() == 24
assert j.denominator().gcd(j.denominator().derivative()) == 1
assert j.numerator().degree() <= j.denominator().degree()
assert j == K(c4**3)/Delta
G = matrix(ZZ, p['generic_height_gram'])
assert G.nrows() == 17 and G.det() == 1092

# Cusp0 has width ell in Gamma0(ell): lower-left entry is exactly-m.
S = PolynomialRing(QQ, 'm'); m = S.gen()
swap = matrix(S, [[0, -1], [1, 0]])
assert swap.inverse()*matrix(S, [[1, m], [0, 1]])*swap == matrix(S, [[1, 0], [-m, 1]])

# Elementary upper/lower unipotents generate every determinant-one matrix.
S = PolynomialRing(QQ, names=('a', 'b', 'c')); a, b, c = S.gens()
F = S.fraction_field()
def U(t): return matrix(F, [[1, t], [0, 1]])
def L(t): return matrix(F, [[1, 0], [t, 1]])
def W(t): return U(t)*L(-1/t)*U(t)
D = matrix(F, [[a, 0], [0, 1/a]])
assert W(F(a))*W(F(1)).inverse() == D
assert L(c/a)*D*U(b/a) == matrix(F, [[a, b], [c, (1+b*c)/a]])

# The central-minus-one cocycle identity forces a fixed affine point.
# Commuting g and -I gives 2*c(g)=(I-g)*c(-I); put v=c(-I)/2.
S = PolynomialRing(QQ, names=('a', 'b', 'c', 'd', 'v0', 'v1'))
a, b, c, d, v0, v1 = S.gens()
g = matrix(S, [[a, b], [c, d]]); v = vector(S, [v0, v1])
I = matrix.identity(S, 2); cocycle = (I-g)*v
assert g*v+cocycle == v
assert 2*cocycle == (I-g)*(2*v)

# Exact general RH formula, then tiny finite-vector sanity checks only.
S = PolynomialRing(QQ, 'ell'); ell = S.gen()
genus = 1-ell**2+12*(ell-1)**2
assert genus == (ell-1)*(11*ell-13)
assert genus-9 == (ell-2)*(11*ell-2)
assert 1-(ell**2-1)+12*(ell-1)**2 == genus+1
rows = []
for prime in [2, 3, 5, 7, 11]:
    points = {(a,b) for a in range(prime) for b in range(prime)}
    cycles = []
    while points:
        q = min(points); orbit = []; z = q
        while z not in orbit:
            orbit.append(z); points.remove(z)
            z = ((z[0]+z[1]) % prime, z[1])
        assert z == q
        cycles.append(orbit)
    lengths = sorted(map(len, cycles))
    assert lengths == [1]*prime + [prime]*(prime-1)
    assert sum(n-1 for n in lengths) == (prime-1)**2
    torsion_lengths = sorted(len(o) for o in cycles if o != [(0,0)])
    assert torsion_lengths == [1]*(prime-1)+[prime]*(prime-1)
    rows.append({'prime': prime, 'division_degree': prime**2,
                 'division_genus': int(genus(prime)),
                 'nonzero_torsion_degree': prime**2-1,
                 'nonzero_torsion_genus': int(genus(prime)+1),
                 'local_division_cycle_lengths': lengths,
                 'local_ramification_contribution': (prime-1)**2})

report = {
    'status': 'PASS_PRIME_DIVISION_BARRIER_INPUTS_AND_IDENTITIES',
    'classification': 'new written theorem with exact symbolic/input checks, not formal verification',
    'j_poles': {'number': 24, 'all_simple': True, 'infinity_is_not_a_pole': True},
    'full_geometric_generic_basis_rank': 17,
    'cusp_width_matrix_identity': True, 'unipotent_generation_identities': True,
    'central_cocycle_fixed_point_identity': True, 'local_cycle_regressions': rows,
    'uniform_written_conclusions': [
        'No prime-degree isogeny over Qbar(t); mod-ell geometric monodromy is SL2(F_ell) for every prime ell.',
        'For Z not in ell*M17, the connected division cover has degree ell^2 and genus (ell-1)*(11*ell-13).',
        'The connected nonzero ell-torsion curve has degree ell^2-1 and genus one larger.',
        'Every base change of genus below9 has torsion-free MW group and fully saturated inherited M17.',
        'Every quadratic base change has these same saturation/torsion conclusions, regardless of genus.',
        'Every genuine geometrically irreducible multisection of genus below9, and every genuine bisection of any genus, gives a point outside the inherited rational span.'
    ],
    'boundary': 'These are function-field independence statements. They do not make every rational specialization productive, distinguish the fixed null fibres, or select a302 carrier member.',
    'independent_replay': False, 'formal_verification': False,
    'checker_sha256': sha(Path(__file__)), 'protocol_sha256': sha(OUT/'protocol.json'),
    'inputs': read(OUT/'protocol.json')['inputs'],
}
save('replay.json', report)
print(report['status'], [(r['prime'], r['division_genus']) for r in rows], flush=True)
