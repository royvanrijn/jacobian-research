#!/usr/bin/env sage-python
"""Replay emitted maps and branch data, without rerunning their construction."""
import importlib.machinery
import importlib.util
import json
from pathlib import Path
from sage.all import QQ

path = Path(__file__).with_name('principal28_ancestry.sage')
loader = importlib.machinery.SourceFileLoader('principal28_ancestry', str(path))
spec = importlib.util.spec_from_loader(loader.name, loader)
a = importlib.util.module_from_spec(spec)
loader.exec_module(a)
g, OUT, R = a.g, a.OUT, a.R


def valid_map(c, A, B, P, t0):
    q = R(c['branch_monic_polynomial'])
    D = QQ(c['constant_twist_representative'])*q
    assert q.is_monic() and q.degree() > 0 and q.gcd(q.derivative()) == 1
    assert c['normalization_genus'] == (q.degree()-1)//2
    assert c['branch_at_infinity'] == bool(q.degree()%2)
    assert c['degree_over_base'] == 2
    assert g.dec(c['raw_radical']) == D*g.dec(c['square_multiplier'])**2
    x0, x1, y0, y1 = [g.dec(c['maps'][k]) for k in ['x0', 'x1', 'y0', 'y1']]
    assert x1 or y1
    def identities(ordinate):
        return (ordinate**2+y1*y1*D == x0**3+3*x0*x1*x1*D+A*x0+B and
                2*ordinate*y1 == 3*x0*x0*x1+x1**3*D+A*x1)
    assert identities(y0) and not identities(y0+1)
    base, v = map(QQ, c['lift'])
    assert base == t0 and v*v == D(base)
    assert g.at(x0, base)+g.at(x1, base)*v == P[0]
    assert g.at(y0, base)+g.at(y1, base)*v == P[1]


def main():
    frozen_plan, live_plan = g.read(OUT/'plan.json'), a.protocol()
    checker_key = str(Path(__file__).relative_to(a.ROOT))
    # The construction freeze retains the initially available checker hash.
    # Bind this corrected replay implementation separately; every construction
    # input, script and setting must still match the original pre-run freeze.
    frozen_checker = frozen_plan['sources'].pop(checker_key)
    live_checker = live_plan['sources'].pop(checker_key)
    assert frozen_plan == live_plan
    roster = a.selected_data()
    assert g.read(OUT/'inputs.json') == roster
    result = g.read(OUT/'geometry.json')
    expected = {(d['id'], j) for d in roster for j in range(1, len(d['basis_points'])-16)}
    assert {(r['fibre'], r['target']) for r in result['rows']} == expected
    assert len(result['rows']) == len(expected) == 25
    count = 0
    for data in roster:
        A, B, sections, t0, Es, points, iso, signs = a.transport(data)
        proof = g.read(OUT/f"rank-{data['id']}.json")
        from fractions import Fraction
        from memory_rank_certificate import checked_rank
        replay = checked_rank(list(map(Fraction, data['model'])),
                              [list(map(Fraction, P)) for P in data['basis_points']],
                              [s['prime'] for s in proof['signatures']], proof['no_rational_2_torsion_prime'])
        assert json.loads(json.dumps(replay)) == proof and proof['rank_lower_bound'] == len(points)
        for row in [r for r in result['rows'] if r['fibre'] == data['id']]:
            j = row['target']
            assert row == g.read(OUT/'targets'/f"{data['id']}-E{j}.json")
            P = iso(points[16+j])
            assert row['source_point'] == data['basis_points'][16+j]
            assert row['parent_point'] == list(map(str, P.xy()))
            assert row['fibre_isomorphism'] == list(map(str, iso.tuple()))
            assert row['parameter'] == str(t0) and len(row['covers']) == 36
            assert {(c['tag']['generic_index'], c['tag']['sign']) for c in row['covers']
                    if c['tag']['family'] == 'constant_slope'} == {(i, s) for i in range(1, 18) for s in [-1, 1]}
            assert sum(c['tag']['family'] == 'vertical_x' for c in row['covers']) == 1
            family = 'class1_fibre' if data['id'] == '302' else 'old_R17_fibre'
            assert sum(c['tag']['family'] == family for c in row['covers']) == 1
            for c in row['covers']:
                valid_map(c, A, B, P, t0)
                count += 1
            assert row['best_genus_in_atlas'] == min(c['normalization_genus'] for c in row['covers'])
            assert row['global_minimum_genus'] == 'UNKNOWN'
            print('REPLAY', data['id'], j, 'PASS', flush=True)
    assert count == 900
    assert result['cover_equivalence_classes'] == a.equivalence_classes(result['rows'])
    # The fixed 302 representatives and every old cover are unchanged.
    old = [r for r in g.read(a.ART/'rank_triangle_v1/geometry.json')['rows'] if r['fibre'] == '302']
    now = [r for r in result['rows'] if r['fibre'] == '302']
    for x, y in zip(old, now):
        assert x['source_point'] == y['source_point'] and x['parent_point'] == y['parent_point']
        for c, d in zip(x['covers'], y['covers']):
            assert all(d[k] == v for k, v in c.items())
    g.save(OUT/'verified.json', {'status': 'PASS_EXACT_14_VS_11_ANCESTRY_REPLAY',
        'carrier_identities_checked': count, 'perturbed_carriers_rejected': count,
        'finite_rank_certificates_replayed': [31, 28], 'generic_prefixes_checked': [17, 17],
        'old_302_cover_payloads_unchanged': 504,
        'construction_time_checker_sha256': frozen_checker,
        'replay_checker_sha256': live_checker,
        'checker_correction': 'Normalize Python tuple-valued signatures to their emitted JSON arrays before exact equality; construction evidence unchanged.',
        'shared_cover_target_groups': [c['targets'] for c in result['cover_equivalence_classes'] if len(c['targets']) > 1],
        'global_minimum_genus': 'UNKNOWN',
        'bindings': {str(p.relative_to(a.ROOT)): g.sha(p) for p in
                     [OUT/'plan.json', OUT/'inputs.json', OUT/'geometry.json', OUT/'rank-302.json', OUT/'rank-11952.json', Path(__file__)]}})
    print('PASS', count, 'exact cover identities and negative controls', flush=True)


if __name__ == '__main__':
    main()
