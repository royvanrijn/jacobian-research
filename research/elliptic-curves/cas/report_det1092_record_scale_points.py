#!/usr/bin/env python3
"""Report certified gains and completed exposure, including unfinished curves."""
import argparse
import hashlib
import json
from pathlib import Path
import det1092_record_scale_points as c
from research_runtime.store import checkpoint

OUT = c.ART / 'det1092_record_scale_points_v1.json'


def report():
    p = c.campaign()
    raw = (c.BATCH / 'ledger.json').read_bytes()
    ledger = json.loads(raw)
    roster = p['rows']
    assert [r['id'] for r in ledger['rows']] == [r['id'] for r in roster[:len(ledger['rows'])]]
    rows, inputs, strata = [], {}, {}
    for chosen in roster:
        key = str(chosen['height_bin']) + ':' + chosen['stratum']
        group = strata.setdefault(key, dict(selected_curves=0, completed_curves=0,
            certified_gains=0, point_boxes_completed=0, certified_boxes=0,
            point_worker_seconds=0, completed_stage_seconds=0))
        group['selected_curves'] += 1
    for row, chosen in zip(ledger['rows'], roster):
        stages = list(row['stages'])
        rank, boxes, certified_boxes, waves, point_seconds = 17, 0, 0, [], 0
        for w in row['waves']:
            folder = c.BATCH / row['id'] / w['id']
            stages += w['stages']
            worker = next((s for s in w['stages'] if s['name'] == 'worker'), None)
            if worker is not None and worker['status'] == 'PASS':
                result = c.cert.read(folder / 'result.json')
                assert len(result['charts']) == 49
                assert all(x['search']['status'] == 'bounded_search_complete' for x in result['charts'])
                boxes += 49
                point_seconds += worker['supervision']['wall_seconds']
            if w['status'] != 'PASS':
                continue
            proof = c.cert.read(folder / 'certification-ledger.json')
            cp = c.cert.read(folder / 'certification-protocol.json')
            assert proof['status'] == 'PASS' and proof['rank_lower_bound'] == w['rank_lower_bound']
            assert proof['initial_rank'] == w['initial_rank'] == rank
            assert proof['completed_boxes'] == w['completed_boxes'] == 49
            assert all(v == proof['rank_lower_bound'] for v in proof['odd_modulus_ranks'].values())
            for name, digest in {**cp['sources'], **cp['inputs'], **proof['certificates']}.items():
                assert c.cert.hashed(c.ROOT / name) == digest
            for name in ['seed.json', 'protocol.json', 'maps.json', 'result.json',
                         'certification-ledger.json', 'certification-protocol.json']:
                path = folder / name
                inputs[str(path.relative_to(c.ROOT))] = c.cert.hashed(path)
            inputs.update(proof['certificates'])
            assert c.cert.hashed(c.ROOT / w['cloud_path']) == w['cloud_sha256']
            rank = proof['rank_lower_bound']
            certified_boxes += 49
            waves.append(dict(id=w['id'], initial_rank=w['initial_rank'],
                              certified_rank_lower_bound=rank, completed_boxes=49))
        if row['status'] == 'PASS':
            assert rank == row['rank_lower_bound'] and boxes == certified_boxes
            assert row['stop_reason'] in ['NO_CERTIFIED_GAIN', 'CERTIFIED_AT_LEAST32', 'FIVE_WAVE_LIMIT']
            if row['stop_reason'] == 'NO_CERTIFIED_GAIN':
                assert waves[-1]['initial_rank'] == rank
        seconds = sum(s['supervision']['wall_seconds'] for s in stages)
        item = dict(id=row['id'], parameter=chosen['parameter'], status=row['status'],
            height_bin=chosen['height_bin'], stratum=chosen['stratum'],
            certified_rank_lower_bound=rank if waves else None,
            certified_gains=rank-17, point_boxes_completed=boxes,
            certified_boxes=certified_boxes, point_worker_seconds=point_seconds,
            completed_stage_seconds=seconds, waves=waves,
            stop_reason=row.get('stop_reason'))
        rows.append(item)
        group = strata[str(row['height_bin'])+':'+row['stratum']]
        group['completed_curves'] += int(row['status'] == 'PASS')
        for key in ['certified_gains', 'point_boxes_completed', 'certified_boxes',
                    'point_worker_seconds', 'completed_stage_seconds']:
            group[key] += item[key]
    assert sum(r['certified_boxes'] for r in rows) == ledger['completed_boxes']
    if ledger['status'] == 'PASS':
        assert len(rows) == 48 and all(r['status'] == 'PASS' for r in rows)
    return dict(schema='elliptic-curves.det1092-record-scale-point-report.v1',
        status=ledger['status'], selected_curves=48, attempted_curves=len(rows),
        completed_curves=sum(r['status'] == 'PASS' for r in rows),
        certified_gains=sum(r['certified_gains'] for r in rows),
        certified_boxes=ledger['completed_boxes'], rows=rows, strata=strata,
        protocol_sha256=c.cert.hashed(c.BATCH / 'protocol.json'),
        ledger_sha256=hashlib.sha256(raw).hexdigest(), inputs=inputs,
        checker_sha256=c.cert.hashed(Path(__file__)),
        boundary='Only completed independent proofs count as certified gains. Completed workers can precede proof completion. Recorded stage seconds include finished failed/censored stages but exclude the live stage and must not be called total eventual cost. Rank lower bounds are not exact ranks. A running snapshot is not a terminal comparison. No external novelty claim.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--snapshot', action='store_true')
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    assert not (args.snapshot and args.check)
    result = report()
    if args.snapshot:
        checkpoint(c.BATCH / 'report-snapshot.json', result)
    else:
        assert result['status'] == 'PASS', 'Final report requires all48 completed curves'
        if args.check:
            assert c.cert.read(OUT) == result
        else:
            assert not OUT.exists()
            checkpoint(OUT, result)
    print(result['status'], result['completed_curves'], '/48 curves;',
          result['certified_boxes'], 'certified boxes;', result['certified_gains'], 'gains')
