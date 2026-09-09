#!/usr/bin/env python3
"""Index already-certified panel/cohort packets; no point or conductor search."""
import argparse
import copy
import json
from pathlib import Path
from v3_warm_support import read, sha, require, atomic
from research_curve_refresh import ROOT, MANIFEST

ART = ROOT / 'artifacts/generated-results/elliptic-curves'


def updated():
    result = copy.deepcopy(read(MANIFEST))
    entries = {r['id']: r for r in result['curves']}
    def add(identifier, family, parameter, rank, kind, source, claim, record=None):
        # Later certified foundry continuations can strengthen these same
        # curves. Replaying this older intake must never lower their bounds.
        if identifier in entries and entries[identifier]['rank_lower_bound'] > rank:
            return
        path = ART / source
        rel = str(path.relative_to(ROOT))
        entry = dict(id=identifier, family=family, parameter=parameter, rank_lower_bound=rank,
                     kind=kind, source=rel, claim=claim)
        if record is not None: entry['record'] = record
        entries[identifier] = entry
        result['sources'][rel] = sha(path)
    for family, parameter, rank in [('074d9','88/2551',20),('07ca9','-2475/2848',20),('08234','2570/2143',22)]:
        claim = 'EC-SECOND-FRESH6-' + ('08234-M22' if family=='08234' else '074D9-07CA9-M20') + '-20260909'
        add('r17-'+family+'-fresh-002',family,parameter,rank,'productive',
            f'fresh6_second_amplification_v1/{family}-complement/result.json',claim)
    low = read(ART / 'fresh6_lowheight_first_m18_v1/result.json')
    for row in low['records']:
        family, identifier = row['family'], row['id']
        source = f'fresh6_lowheight_amplification_v1/{identifier}/reconciled.json'
        kind, claim = 'reconciled', 'EC-LOWHEIGHT-FRESH6-AMPLIFICATION-20260909'
        record = None
        if family == '08234':
            source, kind, claim, record = 'fresh6_lowheight_first_m18_v1/result.json','cohort','EC-LOWHEIGHT-FRESH6-FIRST-M18-20260909',identifier
            packet = row['packet']
        else:
            if family == '103b2':
                source, kind, claim = 'fresh6_lowheight_amplification_v1/103b2-M22-complement/result.json','productive','EC-LOWHEIGHT-FRESH6-103B2-M25-20260909'
            elif family == '07ca9':
                source, claim = 'fresh6_lowheight_amplification_v1/07ca9-M20-complement/reconciled.json','EC-LOWHEIGHT-FRESH6-07CA9-M24-20260909'
            data = read(ART / source)
            packet = next(v for k,v in data['records'].items() if k.endswith('/terminal.json')) if kind=='productive' else data
        add('r17-'+identifier,family,row['parameter'],packet['rank_lower_bound'],kind,source,claim,record)
    for row in read(ART / 'r17_60_panel_results_v1.json')['results']:
        add('r17-panel-'+row['id'],row['family'],row['parameter'],row['rank_lower_bound'],'panel',
            'r17_60_panel_results_v1.json','EC-R17-SIXTY-SEED-COMPLEMENT-PANEL-20260909',row['id'])
    result['curves'] = list(entries.values())
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args()
    result = updated()
    if a.check:require(read(MANIFEST)==result,'ledger selection differs')
    else:atomic(MANIFEST,result)
    print('R17_PANEL_LEDGER_SELECTION_PASS',len(result['curves']))
