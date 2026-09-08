#!/usr/bin/env sage-python
"""Exact soluble/Sha comparison, including marked-cover presentation costs.

Build with --write; replay stored rational transformations with --check.
The only new CAS reduction is hyperellred on 63 pointed quartics. No point
search, class group, or new full pairing is part of this checker. Imported
CT arithmetic is replayed, not accepted from a matrix attestation.
"""
from __future__ import annotations

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import runpy
import sys

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ, ZZ, pari, vector
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / 'elliptic-curves/cas'
sys.path.insert(0, str(CAS))
from fixed_cubic_field_curve_family import field_product

RESULTS = ROOT / 'artifacts/generated-results/elliptic-curves'
PANEL = RESULTS / 'exceptional_soluble_selmer_panel_v1.json'
CT = RESULTS / 'fixed_cubic_u_minus1_cassels_tate_v1.json'
EVIDENCE = RESULTS / 'fixed_cubic_u_minus1_cassels_tate_evidence_v1.json.gz'
FAMILY = RESULTS / 'fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json'
PROBES = RESULTS / 'exceptional_selmer_feasibility_v1.json'
OUTPUT = RESULTS / 'exceptional_soluble_vs_sha_comparison_v1.json'
LOCAL = ROOT / 'artifacts/local/elliptic-curves/exceptional-soluble-vs-sha-v1'
IDS = (356, 385, 398, 400, 543)
R = PolynomialRing(QQ, 't')
t = R.gen()
COVER_MAP = ('pi_short=2*phi-(a,b); raw t=(m00*X+m01*Z)/(m10*X+m11*Z); '
             'raw w=2*k*Y/(m10*X+m11*Z)^2; '
             'phi_x=(t^2-a+w)/2; phi_y=t*(t^2-3*a+w)/2-b; '
             'short to input: x=x_short-b2/12, y=y_short-(a1*x+a3)/2')


def read(p):
    return json.loads(p.read_text())


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def invariants(g):
    e, d, c, b, a = [g[i] for i in range(5)]
    return 12*a*e-3*b*d+c*c, 72*a*c*e-27*a*d*d-27*b*b*e+9*b*c*d-2*c**3


def hom(g, x, z):
    return sum(g[i]*x**i*z**(4-i) for i in range(5))


def rational_cost(q):
    """Numeric slots: sign bit plus numerator and denominator bit lengths."""
    q = QQ(q)
    return 1 + max(1, abs(q.numerator()).nbits()) + q.denominator().nbits()


def slots_cost(values):
    return int(sum(rational_cost(q) for q in values))


def shape_cost(g):
    I, J = invariants(g)
    H = max(abs(c) for c in g)
    # Invariant under scalar multiplication of g, but NOT under GL2.
    ratio = H**6 / max(abs(I)**3, J**2)
    return {'coefficient_slots_bits': slots_cost([g[i] for i in range(5)]),
            'coefficient_max_numerator_denominator_bits': int(max(
                max(abs(c.numerator()).nbits(), c.denominator().nbits()) for c in g)),
            'scalar_normalized_height_sixth': str(ratio)}


def pairing_profile(rows, witnesses):
    """Certified-input linear algebra; callers bind entries to arithmetic.

    A null/missing entry is rejected. Completeness of this restricted square
    block never implies completeness of Sel_2 or a whole-curve rank bound.
    """
    n = len(rows)
    require(all(len(r) == n and all(type(x) is int and x in (0, 1) for x in r)
                for r in rows), 'matrix must be a complete binary square block')
    B = Matrix(GF(2), n, n, sum(rows, []))
    require(B == B.transpose() and all(B[i, i] == 0 for i in range(n)), 'not alternating')
    require(all(len(w) == n and all(type(x) is int and x in (0, 1) for x in w)
                for w in witnesses), 'invalid witness coordinates')
    Q = Matrix(GF(2), len(witnesses), n, sum(witnesses, []))
    require(Q*B == 0, 'a claimed rational class has nonzero CT pairing')
    q, rho = int(Q.rank()), int(B.rank())
    return {'tested_dimension': n, 'pairing_rank': rho,
            'restricted_radical_dimension': n-rho,
            'soluble_intersection_dimension_interval': [q, n-rho],
            'sha_image_dimension_lower_bound': rho,
            'provably_insoluble_class_count': 2**n-2**(n-rho),
            'provably_insoluble_fraction': str(QQ(1)-QQ(1)/2**rho),
            'unresolved_radical_mod_witness_span_dimension': n-rho-q,
            'full_residual_selmer_dimension': None, 'full_radical_dimension': None,
            'whole_curve_rank_upper_bound': None}


def raw_pointed(E, P):
    # Rational short model, with no factorization or choice of isomorphism.
    a = P[0]+E.b2()/12
    b = P[1]+(E.a1()*P[0]+E.a3())/2
    A, B = -E.c4()/48, -E.c6()/864
    require(b*b == a**3+A*a+B, 'short-model point transport')
    g = (t**4-6*a*t*t-8*b*t-3*a*a-4*A)/4
    require(invariants(g) == (E.c4()/16, E.c6()/32), 'pointed invariants')
    return g, a, b


def reduce_pointed(E, P, label):
    raw, a, b = raw_pointed(E, P)
    red = pari('(q)->{my(m);my(z=hyperellred(q,&m));[z,m]}')
    # Clearing by a square preserves the covering and permits integral input.
    z, m = red(raw*raw.denominator()**2)
    g = R(z[0])+R(z[1])**2/4
    I, J = invariants(raw)
    scale = QQ(I/invariants(g)[0]).nth_root(4)
    require(invariants(g)[1]*scale**6 == J, 'quartic reduction changed Jacobian')
    g *= scale**2
    M = Matrix(QQ, m[1])
    x, zz = M*vector([t, 1])
    pull = R(hom(raw, x, zz))
    k = QQ(pull.leading_coefficient()/g.leading_coefficient()).sqrt()
    require(k in QQ and pull == k*k*g, 'quartic reduction map')
    X, Z = M.inverse()*vector(QQ, [1, 0])
    Y = 1/(2*k)
    return {'label': label, 'input_point': [str(P[0]), str(P[1])],
            'short_point': [str(a), str(b)], 'quartic': [str(g[i]) for i in range(5)],
            'raw_parameter_from_reduced': [[str(v) for v in r] for r in M.rows()],
            'raw_quartic_y_scale': str(k),
            'witness_X_Z_Y': list(map(str, [X, Z, Y])),
            'cover_map': COVER_MAP}


def verify_pointed(rec, E):
    require(rec['cover_map'] == COVER_MAP, 'wrong marked covering map')
    P = E(*map(QQ, rec['input_point']))
    raw, a, b = raw_pointed(E, P)
    require(rec['short_point'] == list(map(str, [a, b])), 'short point changed')
    g = R(list(map(QQ, rec['quartic'])))
    M = Matrix(QQ, rec['raw_parameter_from_reduced'])
    k = QQ(rec['raw_quartic_y_scale'])
    require(M.nrows() == 2 and M.ncols() == 2 and M.det() != 0 and k != 0, 'singular transport')
    x, z = M*vector([t, 1])
    require(hom(raw, x, z) == k*k*g, 'wrong covering transformation')
    require(invariants(g) == invariants(raw), 'wrong quartic invariants')
    X, Z, Y = map(QQ, rec['witness_X_Z_Y'])
    require((X, Z) != (0, 0) and hom(g, X, Z) == Y*Y, 'invalid rational witness')
    require(M*vector([X, Z]) == vector([1, 0]) and k*Y == QQ(1)/2, 'wrong infinity transport')
    cost = shape_cost(g)
    cost['curve_model_slots_bits'] = slots_cost(E.a_invariants())
    cost['transport_and_point_slots_bits'] = slots_cost(list(M.list())+[k, a, b])
    cost['witness_slots_bits'] = slots_cost([X, Z, Y])
    cost['total_presentation_slots_bits'] = (cost['curve_model_slots_bits']+cost['coefficient_slots_bits']+
        cost['transport_and_point_slots_bits']+cost['witness_slots_bits'])
    return cost


def check_sources():
    panel = read(PANEL)
    pmod = runpy.run_path(str(CAS / 'certify_exceptional_soluble_selmer_panel.sage'))
    for path, value in panel['input_hashes'].items():
        require(digest(ROOT/path) == value, 'stale soluble-panel input: '+path)
    for case in pmod['inputs']():
        if case['curve_id'] in IDS:
            row = next(r for r in panel['curves'] if r['curve_id'] == case['curve_id'])
            require(pmod['build_case'](case, panel['prime_bound']) == row, 'soluble panel replay differs')
    ctmod = runpy.run_path(str(CAS / 'verify_fixed_cubic_cassels_tate.sage'))
    evidence, summary = json.loads(gzip.decompress(EVIDENCE.read_bytes())), read(CT)
    require(summary['evidence_sha256'] == digest(EVIDENCE), 'CT evidence hash')
    for path, value in summary['source_hashes'].items():
        require(digest(ROOT/path) == value, 'stale CT source: '+path)
    # No need to repeat point searches: they prove no obstruction used here.
    arithmetic = ctmod['verify']({k: v for k, v in evidence.items() if k != 'radical_search'})
    for key in ('matrix', 'pairing_rank', 'restricted_radical_dimension', 'verified_cover_count'):
        require(arithmetic[key] == summary['arithmetic'][key], 'CT replay differs')
    return panel, evidence, arithmetic


def verify_control_independence(evidence, family):
    """Good-prime characters independently distinguish eta and all 18 classes."""
    source = family['anchor']
    B, A = map(QQ, source['base_polynomial_ascending'][:2])
    f = t**3+A*t+B
    require(f.is_irreducible(), 'control must have no rational two-torsion')
    run = next(r for r in family['runs'] if r['parameter_u'] == '-1')
    covers = {r['anchor_mask']: r for r in evidence['covers']}
    betas = [list(map(QQ, covers[b['mask']]['beta'])) for b in run['W_u_basis']]
    for beta, b in zip(betas, run['W_u_basis']):
        expected = field_product([list(map(QQ, v)) for i, v in enumerate(
            source['known_kummer_basis_beta_power_coordinates']) if b['mask'] >> i & 1], A, B)
        require(list(expected) == beta, 'wrong inherited class')
    betas.append(list(map(QQ, [A+1, -1, 1])))
    from sage.all import prime_range
    rows, blocks = [[] for _ in betas], []
    for p in prime_range(3, 2001):
        if f.discriminant().numerator() % p == 0:
            continue
        if any(c.denominator() % p == 0 for b in betas for c in b):
            continue
        roots = sorted(f.change_ring(GF(p)).roots(multiplicities=False), key=int)
        if not roots:
            continue
        values = [[sum(GF(p)(b[i])*root**i for i in range(3)) for root in roots] for b in betas]
        if any(v == 0 for row in values for v in row):
            continue
        trial = [r+[int(not v.is_square()) for v in vs] for r, vs in zip(rows, values)]
        if Matrix(GF(2), trial).rank() == Matrix(GF(2), rows).rank():
            continue
        rows = trial
        blocks.append({'prime': int(p), 'roots': list(map(int, roots))})
        if Matrix(GF(2), rows).rank() == 19:
            break
    require(Matrix(GF(2), rows).rank() == 19, 'control independence not proved by bounded characters')
    return {'rows': rows, 'blocks': blocks, 'rank': 19}


def summary_for(panel, evidence, arithmetic, records):
    target_rows = []
    for cid in IDS:
        src = next(r for r in panel['curves'] if r['curve_id'] == cid)
        E = EllipticCurve(list(map(QQ, src['model'])))
        rs = records[str(cid)]
        # Bind each reduced cover to the actual selected, independent point.
        expected = [(c['label'], [str(QQ(c['cubic_point'][0])/4),
            str((QQ(c['cubic_point'][1])-E.a1()*QQ(c['cubic_point'][0])-4*E.a3())/8)])
                    for c in src['basis_covers']]
        require([(r['label'], r['input_point']) for r in rs] == expected, 'wrong target class binding')
        costs = [verify_pointed(r, E) for r in rs]
        q = len(rs)
        prof = pairing_profile([[0]*q for _ in range(q)], [[int(i == j) for j in range(q)] for i in range(q)])
        target_rows.append({'curve_id': cid, 'generic_rank': src['generic_rank'],
            'rank_lower_bound': src['generic_rank']+q, 'epsilon': 0,
            'basis_construction': 'KNOWN_TARGET_POINTS_RETROSPECTIVE',
            'pairing': prof, 'complementary_selmer_classes': None,
            'costs': costs,
            'dictionary_total_cost_thresholds': sorted(c['total_presentation_slots_bits'] for c in costs),
            'source_quadric_coefficient_bits': src['complexity']['basis_max_coefficient_bits']})
    family = read(FAMILY)
    run = next(r for r in family['runs'] if r['parameter_u'] == '-1')
    E = EllipticCurve(list(map(QQ, run['raw_curve_ainvariants'])))
    B, A = map(QQ, family['anchor']['base_polynomial_ascending'][:2])
    qr = records['control_point']
    require(qr['input_point'] == list(map(str, [A+1, A-B+1])), 'wrong control point')
    qcost = verify_pointed(qr, E)
    I, J = E.c4()/16, E.c6()/32
    inherited = []
    for i, b in enumerate(run['W_u_basis']):
        rec = next(r for r in evidence['covers'] if r['anchor_mask'] == b['mask'])
        g = R(list(map(QQ, rec['quartic'])))
        require(invariants(g) == (I, J), 'Sha/point quartic invariants differ')
        # Store a nonzero partner for every insoluble basis class.
        partner = next(j for j in range(18) if arithmetic['matrix'][i][j])
        inherited.append({'anchor_mask': b['mask'], 'nonzero_pairing_partner_index': partner,
                          'cost': shape_cost(g), 'point_or_sha': 'CERTIFIED_SHA_NONZERO'})
    mixed = [row+[0] for row in arithmetic['matrix']]+[[0]*19]
    return {'targets': target_rows,
        'control': {'curve': 'fixed cubic field u=-1', 'I': str(I), 'J': str(J),
            'reference_subgroup_rank': 1, 'reference_subgroup_generator': 'eta point Q',
            'reference_subgroup_epsilon': 0,
            'inherited_basis_construction': 'FIXED_FIELD_SQUARECLASSES_WITHOUT_POINT_WITNESSES_ON_THIS_CURVE',
            'inherited_subspace': pairing_profile(arithmetic['matrix'], []),
            'point_augmented_subspace': pairing_profile(mixed, [[0]*18+[1]]),
            'independence_certificate': verify_control_independence(evidence, family),
            'known_point_cost': qcost, 'inherited_basis': inherited,
            'all_insoluble_basis_equations_smaller_than_point_control': all(
                c['cost']['coefficient_max_numerator_denominator_bits'] <
                qcost['coefficient_max_numerator_denominator_bits'] for c in inherited),
            'same_curve_same_invariants_opposite_solubility': True}}


def verify_feasibility(panel):
    probes = read(PROBES)
    new = probes['new_probes']
    runner = runpy.run_path(str(CAS/'run_exceptional_selmer_feasibility.py'))
    require(new['source_sha256'] == digest(PANEL), 'wrong probe equation source')
    require(new['runner_sha256'] == digest(CAS/'run_exceptional_selmer_feasibility.py'), 'stale probe runner')
    require(new['worker_source'] == runner['WORKER'] and new['point_searches'] == 0, 'wrong probe worker')
    require([r['curve_id'] for r in new['records']] == [398, 400, 543], 'wrong probe targets')
    outcomes = []
    for r in new['records']:
        row = next(x for x in panel['curves'] if x['curve_id'] == r['curve_id'])
        require(r['input'] == {'model': row['model'], 'stack_bytes': 256_000_000}, 'probe received extra input')
        require(sha256(r['log_text'].encode()).hexdigest() == r['measurement']['log_sha256'], 'probe log hash')
        require(r['measurement']['outcome'] in ('strict_wall_timeout', 'strict_rss_limit'), 'unexpected completed probe')
        require(r['worker']['full_selmer_dimension'] is None and r['worker']['covers'] is None,
                'failed probe claims Selmer data')
        outcomes.append({'curve_id': r['curve_id'], 'outcome': r['measurement']['outcome'],
                         'last_stage': r['worker']['stages'][-1]})
    require([r['curve_id'] for r in probes['historical_probes']] == [356, 385], 'wrong historical probes')
    for r in probes['historical_probes']:
        require(sha256(r['source_text'].encode()).hexdigest() == r['source_sha256'], 'historical record hash')
        require(json.loads(r['source_text']) == r['record'], 'historical transcription changed')
        require(sha256(r['log_text'].encode()).hexdigest() == r['record']['log_sha256'], 'historical log hash')
        require(r['record']['status'] == 'strict_wall_timeout' and r['record']['checkpoint'] is None,
                'unexpected historical completion')
        outcomes.append({'curve_id': r['curve_id'], 'outcome': 'historical_1200_second_timeout',
                         'last_stage': 'class_group_relation_collection'})
    return sorted(outcomes, key=lambda r: r['curve_id'])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write', action='store_true')
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    require(args.write != args.check, 'choose exactly --write or --check')
    panel, evidence, arithmetic = check_sources()
    feasibility = verify_feasibility(panel)
    sources = [PANEL, CT, EVIDENCE, FAMILY, PROBES, Path(__file__),
               CAS/'run_exceptional_selmer_feasibility.py']
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in sources}
    if args.write:
        pari.setrand(1)
        records = {}
        for cid in IDS:
            row = next(r for r in panel['curves'] if r['curve_id'] == cid)
            E = EllipticCurve(list(map(QQ, row['model'])))
            rs = []
            for c in row['basis_covers']:
                X, Y = map(QQ, c['cubic_point'])
                P = E(X/4, (Y-E.a1()*X-4*E.a3())/8)
                rs.append(reduce_pointed(E, P, c['label']))
            records[str(cid)] = rs
            LOCAL.mkdir(parents=True, exist_ok=True)
            (LOCAL/'cover-checkpoint.json').write_text(json.dumps(records, indent=2)+'\n')
            print('REDUCED', cid, len(rs), flush=True)
        fam = read(FAMILY)
        run = next(r for r in fam['runs'] if r['parameter_u'] == '-1')
        E = EllipticCurve(list(map(QQ, run['raw_curve_ainvariants'])))
        B, A = map(QQ, fam['anchor']['base_polynomial_ascending'][:2])
        records['control_point'] = reduce_pointed(E, E(A+1, A-B+1), 'eta')
        result = {'schema': 'elliptic-curves.exceptional-soluble-vs-sha.v1',
            'status': 'PASS_EXACT_RESTRICTED_COMPARISON_FULL_COMPLEMENTS_UNKNOWN',
            'input_hashes': hashes, 'software': {'sage': version, 'pari': str(pari.version())},
            'records': records, 'comparison': summary_for(panel, evidence, arithmetic, records),
            'descent_feasibility': feasibility,
            'claim_boundary': ['No full target Selmer group or complementary Sha class is certified.',
                'Point-derived covers and their cost thresholds are retrospective.',
                'CT rank is basis-invariant on the declared subspace; a zero restricted block is not a solubility test.',
                'Numeric-slot costs include maps and witnesses; reduction is not a minimum over equivalent covers.',
                'The same-curve control rules out Jacobian quartic invariants as a class-solubility separator.']}
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    else:
        result = read(OUTPUT)
        require(result['input_hashes'] == hashes, 'stale comparison inputs')
        require(result['descent_feasibility'] == feasibility, 'stale feasibility summary')
        require(result['comparison'] == summary_for(panel, evidence, arithmetic, result['records']), 'stale comparison')
    print('EXCEPTIONAL_VS_SHA|PASS|target_covers=62|control_CT_rank=16|full_complements=UNKNOWN')


if __name__ == '__main__':
    main()
