#!/usr/bin/env python3
"""Calibration gate with exact public-span identities and bounded search history."""
import argparse
from pathlib import Path
import curve302_recovered_subgroup_followup as c
from research_runtime.store import checkpoint
OUT=c.ART/'curve302_recovery_calibration_v1.json'
def expected():
    p=c.campaign();ledger=c.cert.read(c.BATCH/'ledger.json');assert ledger['status']=='PASS'
    inputs={};waves=[];total=0
    paths=[c.BATCH/'protocol.json',c.BATCH/'ledger.json',c.ORIGINAL,c.CLOUD]
    for row in ledger['waves']:
        folder=c.BATCH/row['id'];proof=c.cert.read(folder/'certification-ledger.json');cp=c.cert.read(folder/'certification-protocol.json')
        assert proof['status']=='PASS' and proof['rank_lower_bound']==row['rank_lower_bound']
        assert all(v==row['rank_lower_bound'] for v in proof['odd_modulus_ranks'].values())
        assert c.cert.hashed(c.ROOT/row['cloud_path'])==row['cloud_sha256']
        for n,h in {**cp['inputs'],**cp['sources'],**proof['certificates']}.items():assert c.cert.hashed(c.ROOT/n)==h
        inputs.update(proof['certificates'])
        assert all(s['supervision']['outcome']=='completed' and s['supervision']['returncode']==0 for s in row['stages'])
        seconds=sum(s['supervision']['wall_seconds'] for s in row['stages']);total+=seconds
        for name in ['geometry','worker','replay']:
            path=folder/(name+'-data-access.json');reads=c.cert.read(path);paths.append(path)
            assert all((c.ROOT/r).is_relative_to(c.BATCH) or c.ROOT/r in (c.ORIGINAL,c.CLOUD) for r in reads)
            assert not any('public-span' in r or 'public-union' in r for r in reads)
        paths += [folder/n for n in ['seed.json','maps.json','result.json','protocol.json','certification-ledger.json','certification-protocol.json']]
        waves.append(dict(id=row['id'],initial_rank=row['initial_rank'],rank_lower_bound=row['rank_lower_bound'],completed_boxes=49,mask_floor=row['mask_floor'],stage_seconds=seconds))
    d=c.BATCH/'public-span';span=c.cert.read(d/'result.json');assert span['status']=='PASS' and span['recovered_exceptional_directions']==7
    paths += [d/n for n in ['protocol.json','input.json','result.json','build.supervisor.json','standalone/protocol.json','standalone/replay.supervisor.json']]
    for q in [d/'build.supervisor.json',d/'standalone/replay.supervisor.json']:
        s=c.cert.read(q);assert s['outcome']=='completed' and s['returncode']==0;total+=s['wall_seconds']
    for n,h in c.cert.read(d/'standalone/protocol.json')['files'].items():assert c.cert.hashed(d/'standalone'/n)==h
    u=c.BATCH/'public-union';union=c.cert.read(u/'ledger.json');assert union['status']=='PASS'
    paths += [u/n for n in ['protocol.json','input.json','ledger.json']]
    for n,h in union['files'].items():assert c.cert.hashed(c.ROOT/n)==h
    inputs.update(union['files']);total+=sum(s['supervision']['wall_seconds'] for s in union['stages'])
    inputs.update({str(q.relative_to(c.ROOT)):c.cert.hashed(q) for q in paths})
    return dict(schema='elliptic-curves.curve302-recovery-calibration.v1',status='PASS',
        starting_generic_rank=17,previous_initial_recovery_rank=19,rank_lower_bound=24,
        recovered_exceptional_directions=7,known_exceptional_dimension=14,significant_fraction_gate_met=True,
        gate_interpretation='The user requested a significant fraction after19-to24 recovery was under way. Seven of fourteen is accepted as half of the exceptional space; this is an operational decision, not a prospectively registered statistical success test.',
        waves=waves,stop_reason=ledger['stop_reason'],followup_boxes=ledger['completed_boxes'],
        total_calibration_boxes=49+ledger['completed_boxes'],public_union_lower_bound=union['rank_lower_bound'],
        public_union_points=union['point_count'],followup_and_postproof_seconds=total,
        prospective_policy='Initial49 generic-only charts, then at most4 updated-subgroup waves. Require certified gain to continue, exclude the preceding subgroup in centre parity, and certify every completed cloud. Same2048 SHA sample,384-bit rounded metric,49 largest computed norms, factor-free map,125000 height and10 seconds per chart. No public points or visibility coordinates in prospective selection.',
        inputs=inputs,source_sha256=c.cert.hashed(Path(__file__)),
        boundary='The original17-only control plus guarded recovered-only follow-ups reach24. All24 recovered basis points have exact rational identities in the public31 span; their quotient over generic17 has dimension7. Public data are used only after terminal point waves. Complete point-cloud union still supplies a lower bound only. No new curve/rank record or universal sensitivity theorem.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--check',action='store_true');a=ap.parse_args();r=expected()
    if a.check:assert c.cert.read(OUT)==r
    else:assert not OUT.exists();checkpoint(OUT,r)
    print('PASS calibration17->24,7/14; union',r['public_union_lower_bound'],flush=True)
