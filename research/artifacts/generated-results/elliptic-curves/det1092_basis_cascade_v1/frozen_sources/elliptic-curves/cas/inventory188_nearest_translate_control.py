#!/usr/bin/env python3
"""Frozen retrospective CVP-translate coverage control; no point enumeration."""
import argparse
from pathlib import Path
from fractions import Fraction as F
import json
from hashlib import sha256
from research_runtime.store import checkpoint
from search_observability import multiply

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
ART = ROOT/'artifacts/generated-results/elliptic-curves'
D = ROOT/'artifacts/local/elliptic-curves/inventory188-nearest-translates-v1'
PUBLIC = ART/'inventory188_public28_reproduction_v1.json'
PREVIOUS = ART/'inventory188_two_witness_chart_comparison_v1.json'
OLD = ART/'full11952_late64_r17_results_v1.json'
OWN = ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1/maps.json'
OUT = ART/'inventory188_nearest_translate_visibility_v1.json'


def read(path):
    return json.loads(path.read_text())


def hashed(path):
    return sha256(path.read_bytes()).hexdigest()


def raw_path():
    row = next(r for r in read(OLD)['curves'] if r['id'] == '11952-0959582')
    path = ROOT/row['discovery_witness']['path']
    assert hashed(path) == row['discovery_witness']['sha256']
    return path


def sources():
    paths = [PUBLIC, PREVIOUS, OLD, OWN, raw_path(), Path(__file__).resolve(),
             CAS/'run_inventory188_nearest_translates.sage',
             CAS/'verify_inventory188_nearest_translates.py',
             CAS/'prospective_half_lattice_v2.sage', CAS/'search_observability.py']
    return {str(p.relative_to(ROOT)): hashed(p) for p in paths}


def protocol():
    p = read(D/'protocol.json')
    assert p['sources'] == sources(), 'frozen control inputs changed'
    return p


def charts():
    answer = []
    for i, c in enumerate(read(raw_path())['charts']):
        s = c['search']; q = s['pointed_chart']
        den, scale, shift = map(F, (q['point_denominator_root'], q['curve_coordinate_scale'], q['shift_mod_denominator_squared']))
        M = multiply((den/scale, shift/(den*scale), 0, 1), multiply(
            tuple(map(F, q['unimodular_horizontal_matrix'])), tuple(map(F, s['horizontal_matrix']))))
        answer.append({'arm': 'generic17', 'index': i,
                       'base_point': [str(F(s['base_point']['x'])), str(F(s['base_point']['y']))],
                       'matrix': list(map(str, M))})
    for i, r in enumerate(read(OWN)['rows']):
        raw = list(map(F, r['raw_coefficients']))
        answer.append({'arm': 'own27', 'index': i,
                       'base_point': [str(-raw[2]/6), str(-raw[1]/8)], 'matrix': r['matrix']})
    assert len(answer) == 98
    return answer


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['freeze']); p.parse_args()
    if (D/'protocol.json').exists():
        raise FileExistsError('preserve fixed translate control')
    checkpoint(D/'protocol.json', {
        'schema': 'elliptic-curves.inventory188-nearest-translates-protocol.v1', 'sources': sources(),
        'public_indices': [26, 27], 'old_subgroup_rank': 27, 'charts_per_arm': 49,
        'height': 125000, 'seconds_per_stage': 300, 'rss_bytes': 1073741824,
        'rule': 'For each public witness separately, compute its full28-point canonical-height Gram at384 bits and round at10^6. Unimodular LLL reduces the old27 principal block. Project the negative witness into that27-space and call floating CVP with a fixed covering radius for coordinate rounding. Keep the original witness, the returned nearest translate, and its54 translates by plus or minus one reduced basis vector; deduplicate in this fixed order. Test both signs against all49 saved generic17 and49 saved own27 charts. No coordinate-directed optimization, extra shell or enumeration follows.',
        'scope': 'Retrospective oracle-assisted representative diagnostic, not prospective selection, point discovery, optimality or full quotient coverage. Existing maps and all completed searches remain unchanged.',
        'following_campaign': None})
    print('FROZEN two witnesses / at most112 representatives / at most21952 coordinates')
