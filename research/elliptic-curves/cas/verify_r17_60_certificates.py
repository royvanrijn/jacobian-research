#!/usr/bin/env python3
"""Verify exported native rank certificates without raw local search files.

This proves subgroup lower bounds, not the complete search-policy history.
Use audit_r17_60_panel.py for the latter. No point searches are performed.
"""
import argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path

import compact_atlas_specialization as atlas
import certify_compact_r17_candidates as cert
from r17_60_arithmetic import native_check
from select_r17_60_panel import ROOT, ART, FAMILIES
from v3_warm_support import read, require, sha


def verify(result_path):
    protocol_path = ART / 'r17_60_panel_protocol_v1.json'
    roster_path = ART / 'r17_60_panel_roster_v1.json'
    protocol, roster, result = read(protocol_path), read(roster_path), read(result_path)
    require(result['status'] == 'COMPLETE_BOUNDED_R17_60_PANEL' and
            result['completed'] == result['independently_verified'] == 60, 'complete sixty-case result required')
    require(result['protocol_sha256'] == sha(protocol_path), 'exported protocol differs')
    for name, digest in protocol['sources'].items():
        require(sha(ROOT / name) == digest, 'frozen source changed: ' + name)
    require(sha(atlas.ATLAS) == protocol['inputs'][str(atlas.ATLAS.relative_to(ROOT))], 'atlas changed')
    roster_hashes = [v for k, v in protocol['inputs'].items() if k.endswith('/r17-60-panel-v1/roster.json')]
    require(roster_hashes == [sha(roster_path)], 'exported roster differs')
    rows = roster['rows']
    require(len(rows) == 60 and len({r['id'] for r in rows}) == 60, 'sixty distinct roster entries required')
    require(Counter((r['family'], r['stratum']) for r in rows) ==
            Counter({(f, s): 5 for f in FAMILIES for s in ('low', 'high')}), 'stratification differs')
    by_id = {r['id']: r for r in rows}
    require(len(result['results']) == 60 and {r['id'] for r in result['results']} == set(by_id), 'result identities differ')
    models = []
    ranks = Counter()
    for record in result['results']:
        row = by_id[record['id']]
        require(record['status'] == 'PASS_INDEPENDENT_R17_60_CASE', 'unresolved case')
        for key in ('family', 'parameter', 'stratum'):
            require(record[key] == row[key], 'case identity differs')
        packet = record['packet']
        native_check(packet, row)
        require(packet['rank_lower_bound'] == record['rank_lower_bound'], 'reported rank differs')
        model = tuple(map(Fraction, packet['curve']))
        require(model == tuple(map(Fraction, row['model'])), 'selected equation differs')
        require(not any(cert.isomorphic(model, other) for other in models), 'duplicate rational isomorphism class')
        models.append(model)
        ranks[record['rank_lower_bound']] += 1
        print('R17_60_NATIVE_CERTIFICATE', record['id'], record['rank_lower_bound'], flush=True)
    print('R17_60_EXPORTED_CERTIFICATES_PASS|cases=60|charts=0|ranks=' + str(dict(sorted(ranks.items()))), flush=True)
    return dict(sorted(ranks.items()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result', type=Path, default=ART / 'r17_60_panel_results_v1.json')
    verify(parser.parse_args().result)
