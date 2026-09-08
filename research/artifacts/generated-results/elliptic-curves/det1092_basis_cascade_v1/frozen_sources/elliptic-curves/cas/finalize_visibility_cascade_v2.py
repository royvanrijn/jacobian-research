#!/usr/bin/env python3
"""Package completed V2 calibration and independent replay without modifying V1."""
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
from research_runtime.store import checkpoint

ROOT = Path(__file__).resolve().parents[2]
D = ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v2.json'


def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p): return str(p.relative_to(ROOT))


def main():
    if OUT.exists(): raise FileExistsError('preserve final V2 result')
    terminal = read(D/'calibration302/terminal.json')
    replay = read(D/'replay.json')
    policy = read(D/'protocol.json')
    metric_replay = read(D/'metric-replay.json')
    assert metric_replay['protocol_sha256']==sha(D/'protocol.json')
    assert metric_replay['checker_sha256']==sha(Path(__file__).with_name('check_visibility_metric_v2.sage'))
    assert len(metric_replay['stages'])==len(terminal['stages'])
    assert all(sha(ROOT/s['selection'])==s['sha256'] for s in metric_replay['stages'])
    preserved = read(D/'v1-preservation.json')
    assert all(sha(ROOT/k)==h for k,h in preserved.items()), 'V1 changed'
    assert all(sha(ROOT/k)==h for k,h in policy['sources'].items())
    assert all(sha(ROOT/k)==h for k,h in policy['inputs'].items())
    assert replay['protocol_sha256']==sha(D/'protocol.json')
    assert replay['checker_sha256']==sha(Path(__file__).with_name('check_visibility_cascade_v2.sage'))
    assert len(replay['stages'])==len(terminal['stages'])
    assert replay['rank_lower_bound']==terminal['final_rank_lower_bound']
    statuses = {}
    charts = []
    bindings = 0
    gains = []
    for stage, checked in zip(terminal['stages'],replay['stages']):
        assert all(sha(ROOT/k)==h for k,h in checked['checkpoint_hashes'].items())
        assert (stage['before'],stage['after'],stage['charts'])==(checked['before'],checked['after'],checked['charts_replayed'])
        wd = D/'calibration302'/f"epoch-{stage['epoch']:02d}"
        selection = read(wd/'selection.json')
        anchor_parities = {tuple(int(x)%2 for x in a['anchor']['representative'][:17]) for a in selection['anchors']}
        assert len(anchor_parities)==2*policy['anchors_per_shell'], 'overlapping anchor fibres'
        points = [tuple(map(F,p)) for p in selection['basis']]
        seen = {(x,abs(y)) for x,y in points}
        prefix = []
        for j,path in enumerate(sorted(wd.glob('chart-*.json'))):
            chart = read(path)
            prefix.append(chart)
            for raw in chart['search']['finite_curve_points']:
                point = F(raw['x']),F(raw['y'])
                key = point[0],abs(point[1])
                if key not in seen:
                    seen.add(key);points.append(point)
            snapshot_path = wd/f'cloud-{j:03d}.json'
            snapshot = read(snapshot_path)
            audit = read(wd/f'mod2-{j:03d}.json')
            assert snapshot['charts']==prefix
            assert snapshot['final_state']['state']['reductions']['points']==selection['basis']
            assert audit['input_sha256']==sha(snapshot_path)
            assert audit['points']==[list(map(str,p)) for p in points]
            bindings += 1
            status = chart['search']['status']
            statuses[status] = statuses.get(status,0)+1
            charts.append({'path':rel(path),'sha256':sha(path)})
        if (wd/'modl.json').exists():
            extra = read(wd/'modl.json')
            assert extra['input_sha256']==sha(wd/stage['audit'])
            assert extra['points']==audit['points']
            assert all(sha(ROOT/k)==h for k,h in extra['sources'].items())
        if stage['after']>stage['before']:
            centre = prefix[-1]['centre']
            gains.append({'before':stage['before'],'after':stage['after'],'lane':centre['lane'],
                'base_shell':centre['shell'],'metric_norm':centre['metric_norm'],
                'coset_fingerprint':centre['fingerprint'],
                'extension_encoding_for_audit_only':centre.get('extension'),
                'outside_first_64_in_this_displayed_basis':centre.get('extension',0)>=64,
                'witness_chart':rel(path)})
    last = terminal['stages'][-1]
    wd = D/'calibration302'/f"epoch-{last['epoch']:02d}"
    certificate = read(wd/last['audit'])
    final_modl = read(wd/'modl.json') if (wd/'modl.json').exists() else None
    artifacts = [p for p in D.rglob('*') if p.is_file() and p.suffix in ('.json','.npz')]
    result = {'schema':'elliptic-curves.full-extension-cascade.v2','calibration_only':True,
        'status':'COMPLETE_INDEPENDENTLY_REPLAYED_BOUNDED_EXPERIMENT',
        'target_achieved':terminal['final_rank_lower_bound']>=31,
        'rank_lower_bound':terminal['final_rank_lower_bound'],
        'chain':[17]+[s['after'] for s in terminal['stages']],
        'protocol':policy,'protocol_sha256':sha(D/'protocol.json'),
        'replay':replay,'stages':terminal['stages'],'chart_status_counts':statuses,
        'metric_replay':metric_replay,
        'certificate_to_transcript_bindings_checked':bindings,
        'gain_witnesses':gains,
        'finalizer_sha256':sha(Path(__file__)),
        'runtime':read(D/'runtime.json'),
        'V1_preserved_files':len(preserved),'V1_preservation_manifest_sha256':sha(D/'v1-preservation.json'),
        'seed':read(D/'calibration302/seed.json'),'final_rank_certificate':certificate,
        'final_modl_certificate':final_modl,
        'final_modl_post_run_diagnostic':last['after']==last['before'],
        'charts':charts,'checkpoint_hashes':{rel(p):sha(p) for p in sorted(artifacts)},
        'claim_boundary':'One frozen calibration-only rule from generic MW17. Full extension masks were cheaply scored for every retained anchor; only a heuristic shortlist received exact rounded-metric CVP and finitely many charts were searched. This does not establish exact rank, full-policy basis invariance, absence of points beyond a stall, or a new rank record.'}
    checkpoint(OUT,result)
    print('FINALIZED',result['chain'],'charts',len(charts),'target achieved',result['target_achieved'],'sha256',sha(OUT))


if __name__=='__main__': main()
