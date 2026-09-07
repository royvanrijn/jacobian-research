#!/usr/bin/env python3
"""Bounded retrospective parent-relative block audit; never a rank search.

build: reduce the quotient metric once, freeze integral lifts, then audit.
check: replay the frozen lifts and exact arithmetic without a height/CAS run.
Additional parents can be supplied to build as a JSON list of {id, columns}.
Columns are an ambient-rank by parent-rank matrix in the public point basis.
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT/'elliptic-curves/cas'), str(ROOT/'elliptic-curves')]
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
import icarm_curve302 as curve
from latent_lattice.elliptic import EllipticCurve
from latent_lattice.integer import modular_rank
import retrospective as r
from verify_explicit_governing_octic import finite_polynomials

OUT = ROOT/'artifacts/generated-results/elliptic-curves'
INPUT = OUT/'curve302_parent_blocks_inputs_v1.json'
OUTPUT = OUT/'curve302_parent_blocks_v1.json'
CORE = OUT/'record_rank17_core_candidates_v1.json'
HEIGHT = OUT/'record_height_lattices_28_29_273_302_v1.json'
COMPONENT = OUT/'record_first17_subgroups_v1.json'
CODE = OUT/'icarm_curve302_point_cloud_v1.json'
LOCAL = ROOT/'artifacts/local/curve302-parent-block'
PRIME_BOUND = 1999


def matrix(a):
    return [[int(x) for x in row] for row in a.tolist()]


def gp_matrix(a):
    return '['+';'.join(','.join(map(str, row)) for row in a)+']'


def row302(path):
    return next(x for x in r.read(path)['curves'] if x['label'] == 'curve302')


def smith(a):
    d = smith_normal_form(a, domain=sp.ZZ)
    return [abs(int(d[i, i])) for i in range(min(d.shape)) if d[i, i]]


def prepare(extra=None):
    core = row302(CORE)
    parents = [{'id': 'additive_core17',
                'columns': matrix(sp.Matrix(core['saturated_basis_columns_in_public_point_coordinates']) *
                                  sp.Matrix(core['core_lll_transform_columns'])),
                'provenance': 'stored saturated additive core, with its stored numerical LLL basis; generic parent UNKNOWN'}]
    if extra:
        parents.extend(r.read(extra))
    assert len(parents) <= 8 and len({p['id'] for p in parents}) == len(parents)
    h = row302(HEIGHT)['height_gram']
    LOCAL.mkdir(parents=True, exist_ok=True)
    rows = []
    for parent in parents:
        c = sp.Matrix(parent['columns']); n, m = c.shape
        assert n == 31 and 1 <= n-m <= 20 and smith(c) == [1]*m
        # Exact Smith transforms give a completion for arbitrary primitive input.
        from sympy.polys.matrices import DomainMatrix
        from sympy.polys.matrices.normalforms import smith_normal_decomp
        d, left, right = smith_normal_decomp(DomainMatrix.from_Matrix(c).convert_to(sp.ZZ))
        left = left.to_Matrix(); right = right.to_Matrix(); d = d.to_Matrix()
        assert left*c*right == d
        complement = left.inv()[:, m:]
        initial = c.row_join(complement)
        assert abs(initial.det()) == 1
        program = f'''default(realprecision,100);
H={gp_matrix(h)}; C={gp_matrix(matrix(c))}; Q={gp_matrix(matrix(complement))};
G=C~*H*C; T=G^-1*C~*H*Q; S=Q~*H*Q-Q~*H*C*T;
V=qflllgram(S); R=round(T*V); W=Q*V-C*R;
print(V);print(R);print(W);print(vector(matsize(V)[1],i,(V~*S*V)[i,i]));
'''
        completed = subprocess.run(['gp', '-q'], input=program, text=True, capture_output=True, timeout=60, check=True)
        if completed.stderr.strip():
            raise RuntimeError(completed.stderr)
        lines = completed.stdout.strip().splitlines()
        assert len(lines) == 4
        def parse_mat(s):
            return [[int(x.strip()) for x in row.split(',')] for row in s.strip('[]').split(';')]
        v, shifts, lifts = map(parse_mat, lines[:3])
        v = sp.Matrix(v); shifts = sp.Matrix(shifts); lifts = sp.Matrix(lifts)
        assert abs(v.det()) == 1 and lifts == complement*v-c*shifts
        assert abs(c.row_join(lifts).det()) == 1
        rows.append({**parent, 'initial_complement_columns': matrix(complement),
                     'quotient_lll_columns': matrix(v), 'parent_shift_columns': matrix(shifts),
                     'residual_lift_columns': matrix(lifts),
                     'projected_height_diagonal_numerical': lines[3],
                     'reduction': '100-digit stored Gram; quotient Schur-complement LLL; rounded parent coordinates; not a closest-vector or exact-height claim'})
        (LOCAL/(parent['id']+'.gp')).write_text(program)
        (LOCAL/(parent['id']+'.out')).write_text(completed.stdout)
    return {'schema': 'curve302.parent-blocks-inputs.v1', 'parents': rows,
            'prime_bound': PRIME_BOUND, 'model': list(map(str, curve.short_coefficients())),
            'bindings': {str(p.relative_to(ROOT)): r.digest(p.read_bytes())
                         for p in (CORE, HEIGHT, COMPONENT, CODE, Path(curve.__file__))}}


def polynomial_degrees(coefficients, p):
    trim, sub, mul, power, gcd = finite_polynomials(p)
    h = trim([r.mod(x, p) for x in coefficients])
    if len(h) != 9 or len(gcd(h, trim([i*h[i] for i in range(1, 9)]))) != 1:
        return None
    xpk = [0, 1]; counts = {}; degrees = []
    for k in range(1, 9):
        xpk = power(xpk, p, h)
        total = len(gcd(sub(xpk, [0, 1]), h))-1
        old = sum(d*n for d, n in counts.items() if k % d == 0)
        assert (total-old) % k == 0
        counts[k] = (total-old)//k
        degrees.extend([k]*counts[k])
    assert sum(degrees) == 8
    return degrees


def radicals_at(points, a, b, p):
    trim, sub, mul, power, gcd = finite_polynomials(p)
    cubic = trim([r.mod(b, p), r.mod(a, p), 0, 1])
    size = p**3; odd = size-1; valuation = 0
    while odd % 2 == 0:
        odd //= 2; valuation += 1
    ns = next(i for i in range(2, p) if pow(i, (p-1)//2, p) == p-1)
    def sqrt(a):
        assert power(a, (size-1)//2, cubic) == [1]
        z = power([ns], odd, cubic); x = power(a, (odd+1)//2, cubic)
        t = power(a, odd, cubic); m = valuation
        while t != [1]:
            probe = t; i = 0
            while probe != [1]:
                probe = mul(probe, probe, cubic); i += 1
            assert i < m
            v = power(z, 2**(m-i-1), cubic)
            x = mul(x, v, cubic); z = mul(v, v, cubic); t = mul(t, z, cubic); m = i
        assert mul(x, x, cubic) == a
        return x
    roots = []
    for x, y in points:
        root = sqrt(trim([r.mod(x, p), -1]))
        if power(root, p*p+p+1, cubic) != [r.mod(y, p)]:
            root = trim([-v for v in root])
        assert power(root, p*p+p+1, cubic) == [r.mod(y, p)]
        roots.append([root, power(root, p, cubic), power(root, p*p, cubic)])
    norms = []
    for i, j in combinations(range(len(points)), 2):
        conjugates = [sub(roots[i][k], [-v for v in roots[j][k]]) for k in range(3)]
        norm = mul(mul(conjugates[0], conjugates[1], cubic), conjugates[2], cubic)
        assert len(norm) == 1 and norm[0]
        norms.append(norm[0])
    return norms


def audit_parent(parent, inputs, public_signatures):
    c = sp.Matrix(parent['columns']); w = sp.Matrix(parent['residual_lift_columns'])
    n, m = c.shape; q = n-m; combined = c.row_join(w)
    assert smith(c) == [1]*m and abs(combined.det()) == 1
    inverse = combined.inv(); projection = inverse[m:, :]
    assert projection*c == sp.zeros(q, m) and projection*w == sp.eye(q)
    assert all(x.q == 1 for x in inverse)
    initial = sp.Matrix(parent['initial_complement_columns'])
    v = sp.Matrix(parent['quotient_lll_columns']); shifts = sp.Matrix(parent['parent_shift_columns'])
    assert abs(v.det()) == 1 and w == initial*v-c*shifts
    assert abs(c.row_join(initial).det()) == 1
    ec = EllipticCurve(curve.short_coefficients())
    points = [ec.linear_combination(curve.SHORT_POINTS, list(map(int, w[:, i]))) for i in range(q)]
    assert all(p is not None and ec.is_on_curve(p) for p in points)
    signatures = []
    for j in range(n):
        code = 0
        for i in range(n):
            if int(combined[i, j]) % 2:
                code ^= public_signatures[i]
        signatures.append(code)
    assert r.rank(signatures) == n and r.rank(signatures[:m]) == m
    # All reduced lift points are independently checked in the same local code.
    model = inputs['model']; a, b = map(F, model[3:])
    blocks = [(p, r.roots_at(str(a), str(b), p)) for p in r.primes(1000)]
    blocks = [(p, roots) for p, roots in blocks if roots]
    assert [r.point_signature(model, p, blocks) for p in points] == signatures[m:]
    mod3 = row302(CODE)['finite_reduction_kummer_codes']['mod_3']
    local3 = [row for block in mod3['local_blocks'] for row in block['canonical_image_row_space_rref']]
    assert modular_rank(local3, 3) == 31
    transformed3 = matrix(sp.Matrix(local3)*combined)
    assert modular_rank(transformed3, 3) == n
    assert modular_rank([row[:m] for row in transformed3], 3) == m
    component = row302(COMPONENT)['bad_fibre_component_code']
    places = component['places']; moduli = component['ambient_component_group_moduli']
    local = sp.Matrix([p['g17_classes']+p['remaining_classes'] for p in places])
    diag = sp.diag(*moduli)
    component_smith = smith(diag.row_join(local*c))
    all_smith = smith(diag.row_join(local))
    # This verifies the full product, not just surjectivity at each place.
    component_data = {'primes': [p['prime'] for p in places], 'moduli': moduli,
                      'parent_image_columns': matrix(local*c), 'residual_image_columns': matrix(local*w),
                      'product_mod_parent_smith': component_smith, 'product_mod_displayed_smith': all_smith}
    pairs = list(combinations(range(q), 2)); octics = []; deltas = []
    for i, j in pairs:
        xp, yp = points[i]; xq, yq = points[j]; delta = xq-xp
        assert delta and yp and yq
        deltas.append(delta)
        octics.append([delta**6, F(0), -4*delta**3*(yq-yp), F(0),
                       6*delta**2*(xp+xq), F(0), -4*(yp+yq), F(0), F(1)])
    table = []; exclusions = []; squarefree = {}; pivot = {}; rank_primes = []
    for p in r.primes(PRIME_BOUND):
        if r.roots_at(str(a), str(b), p) != ():
            continue
        if any(F(t).denominator % p == 0 for point in points for t in point):
            exclusions.append({'prime': p, 'reason': 'lift denominator'}); continue
        if any(r.mod(points[j][0]-points[i][0], p) == 0 for i, j in pairs):
            exclusions.append({'prime': p, 'reason': 'pair x-coordinate collision'}); continue
        norms = radicals_at(points, a, b, p)
        bits = [int(pow(z, (p-1)//2, p) == p-1) for z in norms]
        # One exact squarefree reduction and cycle-type check per octic;
        # plus every pair at the first eligible prime. The rest is the norm formula.
        for k, h in enumerate(octics):
            if k not in squarefree:
                degrees = polynomial_degrees(h, p)
                if degrees:
                    assert degrees == ([1, 1, 3, 3] if bits[k] == 0 else [2, 6])
                    squarefree[k] = {'prime': p, 'degrees': degrees, 'psi': bits[k]}
        old_rank = len(pivot)
        value = r.pack(bits)
        while value:
            top = value.bit_length()-1
            if top not in pivot:
                pivot[top] = value; break
            value ^= pivot[top]
        if len(pivot) != old_rank:
            rank_primes.append(p)
        table.append({'prime': p, 'norms': norms, 'psi': bits})
    assert len(squarefree) == len(pairs)
    def square(value):
        return value > 0 and isqrt(value.numerator)**2 == value.numerator and isqrt(value.denominator)**2 == value.denominator
    scalar_matches = [list(pair) for pair in combinations(range(len(deltas)), 2)
                      if square(deltas[pair[0]]/deltas[pair[1]])]
    identical_bit_columns = [list(pair) for pair in combinations(range(len(pairs)), 2)
                             if all(row['psi'][pair[0]] == row['psi'][pair[1]] for row in table)]
    # Explicit commuting translations isolate each cup-product slot.
    # e_i=(1,0), f_j=(0,1) for an ordered i<j give exactly that commutator.
    commutators = [[int((i, j) == (u, v)) for u, v in pairs] for i, j in pairs]
    assert modular_rank(commutators, 2) == len(pairs)
    return {'id': parent['id'], 'parent_rank': m, 'displayed_rank': n, 'residual_free_rank': q,
            'generic_parent': 'UNKNOWN', 'full_E_Q_quotient': 'UNKNOWN',
            'parent_smith': smith(c), 'completion_determinant': int(combined.det()),
            'projection_rows': matrix(projection), 'full_inverse_rows': matrix(inverse),
            'residual_points_short_model': [[str(x), str(y)] for x, y in points],
            'finite_mod2_rank': n, 'parent_mod2_rank': m, 'residual_mod2_rank': q,
            'finite_mod3_rank': n, 'parent_mod3_rank': m, 'residual_mod3_rank': q,
            'finite_mod2_lift_signatures': signatures[m:],
            'bad_components': component_data,
            'relative_halving_degree': str(4**q),
            'residual_CT_rank': 0,
            'CT_reason': 'rational Kummer images annihilate the entire Selmer space; tautological retrospective vanishing, not an independent CT computation',
            'full_Selmer_dimension': 'UNKNOWN',
            'pair_count': len(pairs), 'central_degree_over_full_displayed_halving_field': str(2**len(pairs)),
            'governing_commutator_rank': len(pairs),
            'cochains': [{'indices_zero_based': list(pair), 'octic_ascending': list(map(str, h)),
                          'squarefree_cycle_witness': squarefree[k]}
                         for k, (pair, h) in enumerate(zip(pairs, octics))],
            'inert_prime_table': table, 'inert_prime_exclusions': exclusions,
            'sample_bit_matrix_rank': len(pivot), 'sample_rank_increase_primes': rank_primes,
            'sample_bit_matrix_kernel_dimension': len(pairs)-len(pivot),
            'identical_sample_bit_columns': identical_bit_columns,
            'scalar_norm_dictionary': {'definition': 'c_ij=x(Q_j)-x(Q_i), with i<j and the frozen reduced lifts',
                                       'scalars': list(map(str, deltas)),
                                       'equal_rational_squareclass_pairs_zero_based': scalar_matches,
                                       'scope': 'exact finite dictionary, not an invariant under parent translations or cochain gauge changes'},
            'sample_interpretation': 'cochain-gauge-dependent exact Frobenius values; no CT matrix or solubility frequencies',
            'low_genus_carrier': 'UNKNOWN; conditional one-point fixed-map explanation of all q directions requires genus at least q when M spans the full generic rational group',
            'shared_solubility_event': 'UNKNOWN; supplied rational witnesses do not identify a family-wide event'}


def compute(inputs):
    for path, digest in inputs['bindings'].items():
        assert r.digest((ROOT/path).read_bytes()) == digest
    assert inputs['prime_bound'] == PRIME_BOUND
    model = inputs['model']; assert model == list(map(str, curve.short_coefficients()))
    gal = r.galois(model); assert gal['galois_group'] == 'S3'
    a, b = model[3:]
    blocks = [(p, r.roots_at(a, b, p)) for p in r.primes(1000)]
    blocks = [(p, roots) for p, roots in blocks if roots]
    assert all(curve.on_short_curve(point) for point in curve.SHORT_POINTS)
    signatures = [r.point_signature(model, p, blocks) for p in curve.SHORT_POINTS]
    assert r.rank(signatures) == 31
    rows = []
    for parent in inputs['parents']:
        result = audit_parent(parent, inputs, signatures)
        rows.append(result)
        print(f"{parent['id']}: Z^{result['residual_free_rank']}, {result['pair_count']} cochains, finite bit rank {result['sample_bit_matrix_rank']}", flush=True)
    return {'schema': 'curve302.parent-blocks.v1', 'status': 'PASS',
            'bindings': {str(p.relative_to(ROOT)): r.digest(p.read_bytes()) for p in
                         (INPUT, Path(__file__), Path(r.__file__), Path(curve.__file__),
                          ROOT/'elliptic-curves/latent_lattice/elliptic.py',
                          ROOT/'elliptic-curves/latent_lattice/integer.py',
                          ROOT/'elliptic-curves/rank-jump/verify_explicit_governing_octic.py')},
            'two_division': gal, 'parents': rows,
            'boundary': 'Exact displayed-subgroup quotients and governing arithmetic. Parent identification, full Mordell-Weil quotient, low-genus carrier and new solubility mechanism remain UNKNOWN.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['build', 'check'])
    parser.add_argument('--additional-parents', type=Path)
    args = parser.parse_args()
    if hasattr(sys, 'set_int_max_str_digits'):
        sys.set_int_max_str_digits(100000)
    if args.mode == 'build':
        if INPUT.exists() or OUTPUT.exists():
            raise SystemExit('refusing to overwrite frozen outputs')
        r.write_new(INPUT, prepare(args.additional_parents))
        r.write_new(OUTPUT, compute(r.read(INPUT)))
    else:
        assert not args.additional_parents
        assert r.read(OUTPUT) == compute(r.read(INPUT))
    print('PASS: parent-relative quotient and governing replay')
