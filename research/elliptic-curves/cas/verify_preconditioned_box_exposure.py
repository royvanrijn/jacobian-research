#!/usr/bin/env python3
"""Seal a fixed exposure audit against the verified parent receipt chain."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path

from pointed_box_equivalence import box_key
import run_parent_seed_v3 as parent
from pointed_quartic_search import PointedQuarticSearch
import pari_pointed_backend as backend


def verify(source, audit):
    read, sha, require = parent.read, parent.sha, parent.require
    report = read(audit/'report.json')
    for category in ('inputs','sources'):
        require(all(sha(parent.ROOT/p)==h for p,h in report[category].items()), 'audit binding changed')
    protocol, terminal, verified = [read(source/n) for n in ('protocol.json','terminal.json','verified.json')]
    parent.guard(protocol)
    require(terminal['protocol_sha256']==sha(source/'protocol.json') and
            verified['status']=='PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY' and
            verified['terminal_sha256']==sha(source/'terminal.json') and
            verified['rank_lower_bound']==terminal['rank_lower_bound'] and
            verified['charts']==terminal['charts'], 'parent replay seal differs')
    require(len(terminal['stages'])==1, 'single parent epoch required')
    stage=terminal['stages'][0]
    epoch=source/'epoch-00'
    require(read(epoch/'stage.json')==stage and stage['before']==stage['after'], 'parent stage differs')
    selection=read(epoch/'landscape/selection.json')
    require(stage['selection_sha256']==sha(epoch/'landscape/selection.json'), 'parent selection differs')
    seed=read(source/'seed.json')
    require(seed['points']==selection['basis']==terminal['points'], 'parent basis differs')
    state=parent.certified_state(tuple(map(F,seed['curve'])),tuple(tuple(map(F,p)) for p in seed['points']),seed['proof'])
    previous=sha(epoch/'landscape/selection.json')
    boxes={}
    for i in range(stage['charts']):
        path=epoch/f'chart-{i:04d}.json'
        chart=read(path)
        require(chart['index']==i and chart['previous_sha256']==previous, 'parent receipt chain differs')
        previous=sha(path)
        mapping=chart['mapping']
        search=PointedQuarticSearch(state=state,centre={'coefficients':mapping['centre']['representative']},coordinate_policy=mapping['coordinate_policy'])
        backend.replay(search,mapping,chart['search'])
        if chart['search']['status']=='bounded_search_complete' and chart['search']['height_bound']==125000 and chart['search']['infinity_checked']:
            boxes.setdefault(tuple(mapping['centre']['point']),set()).add(box_key(mapping['matrix']))
    require(previous==stage['last_chart_sha256'], 'parent chain endpoint differs')
    require([row['centre_index'] for row in report['rows']]==list(range(16)), 'fixed roster differs')
    duplicates=[]
    for row in report['rows']:
        ci=row['centre_index']; wd=audit/f"map-{ci:04d}-{row['policy']}"
        for name in ('input','result','supervisor'):
            require(sha(wd/(name+'.json'))==row[name+'_sha256'], 'audit receipt differs')
        payload,result,sup=[read(wd/(n+'.json')) for n in ('input','result','supervisor')]
        require(payload==dict(curve=seed['curve'],points=seed['points'],centre=selection['centres'][ci],policy='preconditioned_full'), 'prospective map input differs')
        require(sup['outcome']=='completed' and sup['returncode']==0 and sup['worker_result_sha256']==sha(wd/'result.json') and sup['log_sha256']==sha(wd/'worker.log') and result['input_sha256']==sha(wd/'input.json'), 'map execution receipt differs')
        mapping=result['mapping']
        require(mapping['centre']==payload['centre'], 'prospective centre differs')
        search=PointedQuarticSearch(state=state,centre={'coefficients':payload['centre']['representative']},coordinate_policy=mapping['coordinate_policy'])
        backend.validate_map(search,mapping)
        duplicate=box_key(mapping['matrix']) in boxes.get(tuple(mapping['centre']['point']),set())
        require(duplicate==row['duplicate_completed_box'], 'box comparison differs')
        duplicates.append(duplicate)
    require(report['distinct_uncovered_boxes']==sum(not x for x in duplicates), 'exposure count differs')
    output=audit/'exposure-verified.json'
    require(not output.exists(), 'preserve prior verification')
    parent.checkpoint(output,dict(status='PASS_PARENT_CHAIN_AND_EXACT_BOX_EQUIVALENCE',
        report_sha256=sha(audit/'report.json'), source_sha256=sha(Path(__file__)),
        parent_terminal_sha256=sha(source/'terminal.json'), parent_charts_replayed=stage['charts'],
        duplicate_boxes=sum(duplicates), new_point_searches=0,
        claim_boundary='Fixed signed-permutation equivalence of completed height125000 projective boxes including infinity; not an upper bound or exclusion of untested coordinates.'))
    print('VERIFIED_EXPOSURE',len(duplicates),sum(duplicates))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent',type=Path,required=True)
    parser.add_argument('--audit',type=Path,required=True)
    args=parser.parse_args()
    verify(args.parent.resolve(),args.audit.resolve())
