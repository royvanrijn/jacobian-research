#!/usr/bin/env sage-python
"""Exact input/stalk/residue checks for the written all-degree obstruction.

No cohomological proof is claimed as machine verified. No new points are
computed: the immutable independence and halving replays are reused.
"""
import hashlib
import json
import math
import signal
from itertools import product
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, matrix

signal.alarm(25)
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'det1092_unramified_kummer_v1'

def read(p):
    return json.loads(p.read_text())

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

protocol = read(OUT/'protocol.json')
for name, digest in protocol['inputs'].items():
    assert sha(ROOT/name) == digest
parent = read(ART/'curve302_recovered_mw17_parent_v1.json')
brauer = read(ART/'det1092_surface_brauer_triviality_v1/independent-replay.json')
assert brauer['status'] == 'PASS_ARITHMETIC_INPUTS_FOR_GLOBAL_BRAUER_TRIVIALITY'
for name, digest in brauer['inputs'].items():
    assert sha(ROOT/name) == digest
assert [(r['p'], r['Picard_index'], r['Brauer_order']) for r in brauer['reductions']] == [(149,1,1),(151,1,1)]

R = PolynomialRing(QQ, 't')
K = R.fraction_field()
def dec(row):
    return K(R(row['numerator']))/R(row['denominator'])

a1,a2,a3,a4,a6 = [dec(row) for row in parent['a_invariants']]
assert (a1,a2,a3) == (1,1,1)
G = matrix(ZZ, parent['generic_height_gram'])
assert G.nrows() == G.rank() == 17 and G.det() == 1092
T = PolynomialRing(K, 'X')
x = T.gen()
f = x**3+5*x**2+(16*a4+8)*x+64*a6+16
disc = R(f.discriminant())
assert disc.degree() == 24 and disc.gcd(disc.derivative()).degree() == 0
assert R(a4).degree() == 8 and R(a6).degree() == 12
finf = PolynomialRing(QQ, 'X')([64*R(a6)[12],16*R(a4)[8],0,1])
assert finf.discriminant() and finf[0]  # smooth infinity, all three roots nonzero

# The split norm sequence on the good and nodal geometric stalks.
triples = list(product(range(2), repeat=3))
even = [v for v in triples if sum(v)%2 == 0]
fixed = [v for v in triples if v[0] == v[1]]
nodal_kernel = [v for v in fixed if sum(v)%2 == 0]
assert len(even) == 4 and nodal_kernel == [(0,0,0),(1,1,0)]
assert {sum(v)%2 for v in fixed} == {0,1}
assert all(sum((b,b,b))%2 == b for b in (0,1))

construction = read(ART/'det1092_seed_norm_lift_v1/construction.json')
norm_replay = read(ART/'det1092_seed_norm_lift_v1/replay.json')
assert norm_replay['status'] == 'PASS_INDEPENDENT_NORM_LIFT_AND_EIGHT_NON_SELMER_OBSTRUCTIONS'
for row in (construction, norm_replay):
    for name, digest in row['inputs'].items():
        assert sha(ROOT/name) == digest
k = QQ(construction['Xstar'])
D = dec(construction['D'])
assert D == f(k) and D.denominator().degree() == 0
assert R(D).degree() == 12 and R(D).gcd(R(D).derivative()) == 1
assert R(D).gcd(disc) == 1
alpha = T([dec(row) for row in construction['alpha_coefficients']])
assert alpha == D*(k-x)
assert f.resultant(alpha) == D**4
assert QQ(D(0)).is_square()
# Simple zero of D over a good parameter: one simple x=k branch, two units.
local_orders = [1+v for v in (1,0,0)]
assert local_orders == [2,1,1]
assert -R(D).degree()-4 == -16  # each infinity branch; finf[0] != 0 above

controls = []
for row, checked in zip(construction['cases'], norm_replay['cases']):
    assert row['index'] == checked['index']
    tau = QQ(row['parameter'])
    value = QQ(f(k)(tau))
    assert value == QQ(row['D_value'])
    support = math.prod(int(v) for v in row['support_factors'])
    rem = abs(int(value.numerator()))
    for divisor in row['gcd_strips']:
        g0 = math.gcd(rem, support)
        assert g0 == int(divisor) and g0 > 1
        rem //= g0
    assert rem == int(row['coprime_remainder']) and math.gcd(rem,support) == 1
    if tau == 0:
        assert checked['status'] == 'PASS_RATIONAL_SEED_CLASS_OUTSIDE_GENERIC_SPAN'
        proof = checked['independence']
        m = matrix(GF(2), proof['matrix_rows'])
        assert m.rank() == proof['rank'] == 18 and m[:,:17].rank() == 17
        controls.append({'index':row['index'],'old_result':'rational primitive seed class'})
    else:
        assert math.isqrt(rem)**2 != rem
        assert checked['status'] == 'PASS_EXACT_GOOD_PRIME_LOCAL_OBSTRUCTION'
        controls.append({'index':row['index'],'old_result':'not Selmer: good-prime odd valuation'})
assert len(controls) == 9

split = read(ART/'det1092_split_descent_v1/independent-replay.json')
assert split['status'] == 'PASS_INDEPENDENT_HALVING_CYCLE_REPLAY'
for name, digest in split['inputs'].items():
    assert sha(ROOT/name) == digest
counts = {'primitive_productive_branches':0,'productive_antitraces_depth_one':0,
          'generic_recoveries_over_M16':0,'generic_recoveries_over_M17':0,'dependent_conic':0}
for row in split['cases']:
    fam, label = row['family'], row['label']
    if fam == 'dependent-conic':
        assert row['status'] == 'INHERITED_RATIONAL_SPAN' and row['steps'] in (2,3)
        counts['dependent_conic'] += 1
    elif label.endswith('core16'):
        assert row['status'] == 'NEW_INDEPENDENT_DIRECTION' and row['steps'] == 0
        counts['generic_recoveries_over_M16'] += 1
    elif label.endswith('full17'):
        assert row['status'] == 'INHERITED_RATIONAL_SPAN' and row['steps'] in (2,3)
        counts['generic_recoveries_over_M17'] += 1
    elif label.startswith('anti-trace'):
        assert row['status'] == 'NEW_INDEPENDENT_DIRECTION' and row['steps'] == 1
        counts['productive_antitraces_depth_one'] += 1
    else:
        assert label.startswith('branch')
        assert row['status'] == 'NEW_INDEPENDENT_DIRECTION' and row['steps'] == 0
        counts['primitive_productive_branches'] += 1
assert list(counts.values()) == [8,8,9,9,4] and split['case_count'] == 38

report = {
    'status':'PASS_HYPOTHESES_AND_RESIDUE_REGRESSION',
    'classification':'verified arithmetic application; written cohomological deduction',
    'generic_fibres':'24I1; smooth infinity',
    'spectral_genus_by_Riemann_Hurwitz':10,
    'norm_identity':'N(D*(k-theta))=D^4',
    'retrospective_residues':{'good_parameter_degree':12,'orders_above_each':[2,1,1],
        'geometric_ramification_points':24,'orders_at_infinity':[-16,-16,-16]},
    'control_decisions_reused':counts,'old_norm_lift_controls':controls,
    'new_points_computed':0,'cohomology_formally_verified':False,
    'proof_boundary':'The all-degree equality is proved in the canonical note using Kummer, the split cubic norm sheaf, Leray and the already proved Br(X)=Br(Q). This script checks inputs and local algebra, not those theorems.',
    'inputs':protocol['inputs'],'protocol_sha256':sha(OUT/'protocol.json'),
    'checker_sha256':sha(Path(__file__))}
dest = OUT/'replay.json'
if dest.exists():
    assert read(dest) == report
else:
    with dest.open('x') as stream:
        json.dump(report,stream,indent=2,sort_keys=True)
        stream.write('\n')
print(report['status'],counts,flush=True)
