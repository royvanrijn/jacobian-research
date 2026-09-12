#!/usr/bin/env sage-python
"""Frozen, matched degree-two ancestry atlas: 302 +14 versus 11952 +11.

No point search, local arithmetic census, or claim of global genus minimality.
Run under an outer 300-second timeout; completed target files are immutable.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import importlib.machinery
import importlib.util
from pathlib import Path

from sage.all import QQ, ZZ, EllipticCurve, PolynomialRing, prime_range, version
from memory_rank_certificate import checked_rank

GPATH = Path(__file__).with_name('rank_triangle_geometry.sage')
loader = importlib.machinery.SourceFileLoader('ancestry_geometry', str(GPATH))
spec = importlib.util.spec_from_loader(loader.name, loader)
g = importlib.util.module_from_spec(spec)
loader.exec_module(g)
R, F, t, ART, ROOT = g.R, g.F, g.t, g.ART, g.ROOT
OUT = ART/'rank_ancestry_principal28_v1'
INPUT = ART/'rank_accessibility_subsets_v1/inputs.json.gz'
PUBLIC = ART/'inventory188_public28_reproduction_v1.json'
LITERAL = ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
ATLAS = ART/'compact_six_r17_atlas_v1.json'


def source_paths():
    return [Path(__file__), GPATH, Path(__file__).with_name('verify_principal28_ancestry.sage'),
            Path(__file__).with_name('memory_rank_certificate.py'), INPUT, PUBLIC, LITERAL, ATLAS,
            ART/'curve302_recovered_mw17_parent_v1.json',
            ROOT/'artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-inputs/parents/11952.json',
            ART/'rank_triangle_v1/geometry.json'] + [
                ART/('x1092_class1_realization_'+s+'_v1.json') for s in ['trace', 'rr', 'equation']]


def protocol():
    return {'schema': 'elliptic-curves.principal28-ancestry-plan.v1',
            'sources': {str(p.relative_to(ROOT)): g.sha(p) for p in source_paths()},
            'runtime': version(), 'selection_mode': 'retrospective_known_high_rank',
            'targets': {'302': 14, '11952': 11},
            'parameters': {'302': '0', '11952': '110314/102227'},
            'atlas_per_target': 'vertical x; chords through +/-17 generic sections; one existing alternate-fibration fibre',
            'limits': {'targets': 25, 'covers_per_target': 36, 'wall_seconds_per_run': 300, 'jobs': 1},
            'complexity_order': ['degree_over_base', 'normalization_genus'],
            'cover_equivalence': 'over the fixed original base: monic branch polynomial equal and twist ratio a rational square',
            'global_minimum_genus': 'UNKNOWN', 'new_point_searches': False,
            'basis_invariant_rank_jump_decomposition': False}


def selected_data():
    old = next(r for r in g.read(INPUT)['curves'] if r['id'] == '302')
    pub = g.read(PUBLIC)
    assert pub['status'] == 'PASS' and pub['rank_lower_bound'] == 28
    assert pub['parameter'] == '110314/102227'
    assert pub['independent_column_indices'] == list(range(27)) + [53]
    assert pub['independent_points'][:17] == pub['points'][:17]
    return [dict(id='302', model=old['model'], basis_points=old['basis_points'],
                 basis_words_in_D=old['basis_words_in_D']),
            dict(id='11952', model=pub['curve'], basis_points=pub['independent_points'],
                 original_union_indices=pub['independent_column_indices'])]


def parent(cid):
    A, B, sections, t0 = g.parent(cid)
    if cid == '11952':
        t0 = QQ(110314)/102227
    return A, B, sections, t0


def transport(data):
    A, B, sections, t0 = parent(data['id'])
    Es = EllipticCurve(QQ, list(map(QQ, data['model'])))
    E0 = EllipticCurve(QQ, [g.at(A, t0), g.at(B, t0)])
    points = [Es(list(map(QQ, P))) for P in data['basis_points']]
    expected = [E0(g.at(x, t0), g.at(y, t0)) for x, y in sections]
    # Match every column in order, allowing the recorded generic sign convention.
    iso = next(i for i in Es.isomorphisms(E0)
               if all(i(P) == Q or i(P) == -Q for P, Q in zip(points[:17], expected)))
    signs = [1 if iso(P) == Q else -1 for P, Q in zip(points[:17], expected)]
    assert len(expected) == 17
    assert all(y*y == x*x*x+A*x+B for x, y in sections)
    return A, B, sections, t0, Es, points, iso, signs


def native_carrier_factory(A, B):
    """The retained direct map, with target-independent model transport."""
    literal = g.read(LITERAL)
    compact = next(r for r in g.read(ATLAS)['families'] if r['family'] == '11952')
    assert g.sha(LITERAL) == compact['source_sha256']
    a, b, c, d = map(QQ, compact['base_matrix_a_b_c_d'])
    scale = QQ(compact['total_scale_from_literal_source'])
    u = F((a*t+b)/(c*t+d))
    substitute = lambda rec: g.dec(rec)(u)
    gm = literal['genus_one_model']
    qpoly = PolynomialRing(F, 'z')([substitute(r) for r in gm['q_coefficients_in_t_low_to_high']])
    z = qpoly.parent().gen()
    tzero = substitute(gm['distinguished_point_from_old_zero']['t0'])
    vzero = substitute(gm['distinguished_point_from_old_zero']['W0'])
    gauge = substitute(literal['weierstrass_model']['gauge'])
    shift = qpoly(z+tzero)
    ee, dd, cc, bb, aa = [shift[i] for i in range(5)]
    assert ee == vzero*vzero
    a1 = dd/vzero
    a2 = cc-dd*dd/(4*vzero*vzero)
    a3 = 2*vzero*bb
    a4 = -4*vzero*vzero*aa
    a6 = a2*a4
    b2 = a1*a1+4*a2
    factor = F(3*gauge*(c*t+d)**2/scale)
    Egen = EllipticCurve(F, [a1, a2, a3, a4, a6])
    assert -Egen.c4()/48*factor**4 == A and -Egen.c6()/864*factor**6 == B

    def carrier(P, tstar):
        px, py = P
        at = lambda f: g.at(f, tstar)
        xg = px/at(factor)**2-at(b2)/12
        yg = py/at(factor)**3-(at(a1)*xg+at(a3))/2
        old_t = at(tzero)+(2*at(vzero)*(xg+at(cc))-at(dd)**2/(2*at(vzero)))/yg
        zs = old_t-at(tzero)
        W = (xg*zs**2-at(dd)*zs-2*at(vzero)**2)/(2*at(vzero))
        raw = F(qpoly(old_t))
        zc = F(old_t-tzero)
        x0 = (2*vzero*vzero+dd*zc)/zc**2
        x1 = 2*vzero/zc**2
        y0 = (4*vzero**3+2*vzero*dd*zc+(2*vzero*cc-dd*dd/(2*vzero))*zc**2)/zc**3
        y1 = 4*vzero*vzero/zc**3
        y0 = (y0+(a1*x0+a3)/2)*factor**3
        y1 = (y1+a1*x1/2)*factor**3
        x0 = (x0+b2/12)*factor**2
        x1 = x1*factor**2
        return raw, (x0, x1, y0, y1), W, {'family': 'old_R17_fibre', 'parameter': str(old_t)}
    return carrier


def record_cover(candidate, A, B, P, t0):
    raw, (x0, x1, y0, y1), v0, tag = candidate
    raw = F(raw)
    assert y0*y0+y1*y1*raw == x0**3+3*x0*x1*x1*raw+A*x0+B
    assert 2*y0*y1 == 3*x0*x0*x1+x1**3*raw+A*x1
    assert g.at(raw, t0) == v0*v0
    assert g.at(x0, t0)+g.at(x1, t0)*v0 == P[0]
    assert g.at(y0, t0)+g.at(y1, t0)*v0 == P[1]
    constant, q, s = g.squarefree_cover(raw)
    assert q.degree() > 0 and g.at(s, t0) != 0
    v = v0/g.at(s, t0)
    assert v*v == constant*q(t0)
    return {'tag': tag, 'degree_over_base': 2, 'normalization_genus': int((q.degree()-1)//2),
            'branch_monic_polynomial': list(map(str, q.list())), 'branch_at_infinity': bool(q.degree()%2),
            'constant_twist_representative': str(constant), 'square_multiplier': g.enc(s),
            'raw_radical': g.enc(raw), 'lift': [str(t0), str(v)],
            'maps': {k: g.enc(f) for k, f in zip(['x0', 'x1', 'y0', 'y1'], [x0, x1*s, y0, y1*s])},
            'rational_bisection_orbit_class': 'NOT_APPLICABLE_POSITIVE_GENUS' if q.degree() > 2 else 'UNKNOWN_ORBIT_LABEL',
            'geometric_orbit_identification': 'construction pencil only; no full NS orbit canonicalization'}


def equivalence_classes(rows):
    buckets = defaultdict(list)
    classes = []
    for row in rows:
        for k, cover in enumerate(row['covers']):
            key = row['fibre'], tuple(cover['branch_monic_polynomial'])
            old = next((c for c in buckets[key] if
                        (QQ(c['constant'])/QQ(cover['constant_twist_representative'])).is_square()), None)
            if old is None:
                old = {'fibre': row['fibre'], 'branch': list(key[1]),
                       'constant': cover['constant_twist_representative'], 'members': []}
                buckets[key].append(old)
                classes.append(old)
            old['members'].append([row['target'], k])
    for c in classes:
        c['targets'] = sorted({m[0] for m in c['members']})
    return classes


def run():
    g.save(OUT/'plan.json', protocol())
    roster = selected_data()
    g.save(OUT/'inputs.json', roster)
    allrows = []
    subgroup_checks = []
    for data in roster:
        cid = data['id']
        A, B, sections, t0, Es, points, iso, signs = transport(data)
        primes = [int(p) for p in prime_range(3, 998) if Es.discriminant().valuation(p) == 0
                  and all(a.denominator()%p for a in Es.a_invariants())]
        torsion_prime = next(p for p in primes if not Es.change_ring(__import__('sage.all', fromlist=['GF']).GF(p)).division_polynomial(2).roots())
        proof = checked_rank(list(map(Fraction, data['model'])),
                             [list(map(Fraction, P)) for P in data['basis_points']], primes, torsion_prime)
        g.save(OUT/f'rank-{cid}.json', proof)
        subgroup_checks.append({'fibre': cid, 'rank_lower_bound': len(points), 'generic_rank': 17,
                                'known_quotient_rank': len(points)-17, 'generic_column_signs': signs,
                                'fibre_isomorphism': list(map(str, iso.tuple())),
                                'minimum_multisection_degree': 2,
                                'degree_argument': 'A degree-one horizontal curve is a generic section. Generic rank17 and independence of M17 plus these extras exclude its specialization.'})
        native = native_carrier_factory(A, B) if cid == '11952' else None
        for j, original in enumerate(points[17:], 1):
            dest = OUT/'targets'/f'{cid}-E{j}.json'
            if dest.exists():
                allrows.append(g.read(dest))
                continue
            P = iso(original)
            px, py = P.xy()
            candidates = [(F(px**3+A*px+B), (F(px), F(0), F(0), F(1)), py, {'family': 'vertical_x'})]
            for i, (tx, ty) in enumerate(sections, 1):
                for sign in [-1, 1]:
                    sy = sign*ty
                    m = (py+g.at(sy, t0))/(px-g.at(tx, t0))
                    D = m**4-6*tx*m*m-8*sy*m-3*tx*tx-4*A
                    x0 = (m*m-tx)/2
                    candidates.append((D, (x0, F(QQ(1)/2), m*(x0-tx)-sy, F(m/2)),
                                       2*px-(m*m-g.at(tx, t0)),
                                       {'family': 'constant_slope', 'generic_index': i, 'sign': sign, 'slope': str(m)}))
            candidates.append(g.class1(A, B, (px, py), t0) if cid == '302' else native((px, py), t0))
            covers = [record_cover(c, A, B, P, t0) for c in candidates]
            row = {'fibre': cid, 'parameter': str(t0), 'target': j,
                   'source_point': list(map(str, original.xy())), 'parent_point': list(map(str, P.xy())),
                   'fibre_isomorphism': list(map(str, iso.tuple())), 'covers': covers,
                   'minimum_degree': 2, 'best_genus_in_atlas': min(c['normalization_genus'] for c in covers),
                   'global_minimum_genus': 'UNKNOWN', 'cohort': 'historical_external',
                   'selection_mode': 'retrospective_known_high_rank'}
            g.save(dest, row)
            allrows.append(row)
            print('ANCESTRY', cid, j, 'covers', len(covers), 'best genus', row['best_genus_in_atlas'], flush=True)
    classes = equivalence_classes(allrows)
    result = {'status': 'PASS_FIXED_14_VS_11_ATLAS', 'rows': allrows, 'subgroup_checks': subgroup_checks,
              'cover_equivalence_classes': classes, 'global_minimum_genus': 'UNKNOWN',
              'boundary': 'Existing alternate-fibration incidence, not a shared-cover rank decomposition or a proved cause of either jump.'}
    g.save(OUT/'geometry.json', result)
    for cid in ['302', '11952']:
        rows = [r for r in allrows if r['fibre'] == cid]
        print(cid, 'genera', dict(Counter(c['normalization_genus'] for r in rows for c in r['covers'])),
              'shared target covers', sum(len(c['targets']) > 1 for c in classes if c['fibre'] == cid), flush=True)


if __name__ == '__main__':
    run()
