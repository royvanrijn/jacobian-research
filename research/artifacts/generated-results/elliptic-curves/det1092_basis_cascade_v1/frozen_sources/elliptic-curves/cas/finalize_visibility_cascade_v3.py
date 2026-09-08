"""Bind frozen V3 schedules, transcripts, independent certificates and replays."""
import json
import hashlib
from pathlib import Path
from fractions import Fraction as F
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
D=ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v3.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(ROOT))
def package(start):
    out=D/('terminal-M30' if start==30 else 'replay-M17')
    terminal=read(out/'terminal.json')
    if start==17:
        terminal_reads=[k for k in terminal['read_paths'] if '/terminal-M30/' in k]
        assert terminal_reads==[rel(D/'terminal-M30/terminal.json')], 'M17 consumed first-trial witnesses'
    replay=read(D/f'replay-M{start}.json')
    metric=read(D/f'metric-replay-M{start}.json')
    assert replay['checker_sha256']==sha(Path(__file__).with_name('check_visibility_cascade_v3.sage'))
    assert metric['checker_sha256']==sha(Path(__file__).with_name('check_visibility_metric_v3.sage'))
    assert replay['protocol_sha256']==metric['protocol_sha256']==sha(D/'protocol.json')
    assert replay['rank_lower_bound']==terminal['final_rank_lower_bound']
    assert len(terminal['stages'])==len(replay['stages'])==len(metric['stages'])
    for s in metric['stages']:assert sha(ROOT/s['selection'])==s['sha256']
    statuses={};gains=[];bindings=0
    for stage,checked in zip(terminal['stages'],replay['stages']):
        assert all(sha(ROOT/k)==h for k,h in checked['checkpoint_hashes'].items())
        assert (stage['before'],stage['after'],stage['charts'])==(checked['before'],checked['after'],checked['charts_replayed'])
        wd=out/f"epoch-{stage['epoch']:02d}"
        selection=read(wd/'selection.json')
        assert len({tuple(int(v)%2 for v in a['anchor']['representative'][:17]) for a in selection['anchors']})==32
        points=[tuple(map(F,p)) for p in selection['basis']]
        seen={(x,abs(y)) for x,y in points};prefix=[]
        for j,path in enumerate(sorted(wd.glob('chart-*.json'))):
            chart=read(path);prefix.append(chart)
            for raw in chart['search']['finite_curve_points']:
                point=F(raw['x']),F(raw['y']);key=point[0],abs(point[1])
                if key not in seen:seen.add(key);points.append(point)
            snap=wd/f'cloud-{j:03d}.json';snapshot=read(snap);audit=read(wd/f'mod2-{j:03d}.json')
            assert snapshot['charts']==prefix
            assert snapshot['final_state']['state']['reductions']['points']==selection['basis']
            assert audit['input_sha256']==sha(snap)
            assert audit['points']==[list(map(str,q)) for q in points]
            bindings+=1
            status=chart['search']['status'];statuses[status]=statuses.get(status,0)+1
        if stage['after']>stage['before']:
            extra=read(wd/'modl.json')
            assert extra['input_sha256']==sha(wd/stage['audit']) and extra['points']==audit['points']
            gains.append({'before':stage['before'],'after':stage['after'],'chart':rel(path),'centre':prefix[-1]['centre']})
    last=terminal['stages'][-1];wd=out/f"epoch-{last['epoch']:02d}"
    return {'start':start,'rank_lower_bound':terminal['final_rank_lower_bound'],'chain':[start]+[s['after'] for s in terminal['stages']],
            'terminal':terminal,'replay':replay,'metric_replay':metric,'chart_status_counts':statuses,
            'certificate_to_transcript_bindings_checked':bindings,'gain_witnesses':gains,
            'seed':read(out/'seed.json'),'final_rank_certificate':read(wd/last['audit'])}
def main():
    if OUT.exists():raise FileExistsError('preserve completed V3')
    policy=read(D/'protocol.json')
    prior=read(ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2/protocol.json')
    unchanged=('anchors_per_shell','canonical_per_shell','height','seconds_per_chart','max_charts',
               'max_epochs','target_rank','prime_bound','exact_cvp_node_limit','gp_sha256','incremental','metric')
    assert all(policy[k]==prior[k] for k in unchanged)
    v2_result=read(ROOT/'artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v2.json')
    assert read(D/'fixed-M30.json')['points']==v2_result['final_rank_certificate']['independent_points']
    for block in ('sources','inputs'):assert all(sha(ROOT/k)==h for k,h in policy[block].items())
    preserved=read(D/'v2-preservation.json')
    assert all(sha(ROOT/k)==h for k,h in preserved.items())
    first=package(30);runs=[first]
    if first['rank_lower_bound']>=31:runs.append(package(17))
    result={'status':'COMPLETE_INDEPENDENTLY_REPLAYED_BOUNDED_CALIBRATION','schema':'visibility-cascade.v3',
            'protocol':policy,'protocol_sha256':sha(D/'protocol.json'),'runs':runs,'V2_preserved_files':len(preserved),
            'unchanged_V2_policy_fields_checked':list(unchanged),
            'finalizer_sha256':sha(Path(__file__)),'runtime':read(D/'runtime.json'),
            'claim_boundary':'Target-blind execution after retrospective calibration. Certified rank lower bounds only. Finite heuristic shortlist; no absence or exact-rank claim. Fixed-M30 starts with empty tested-chart history.'}
    calibration=ROOT/'artifacts/local/elliptic-curves/v2-terminal-retrospective-v1/v3-prefreeze-metrics.json'
    result['retrospective_prefreeze_calibration']=read(calibration)
    result['retrospective_prefreeze_calibration_sha256']=sha(calibration)
    result['fixed_M30_target_achieved']=first['rank_lower_bound']>=31
    result['M17_target_achieved']=len(runs)==2 and runs[1]['rank_lower_bound']>=31
    result['M17_dataflow_boundary']='Initial basis is generic MW17. The combined protocol hash-checks both seed files, but the M17 branch never parses the M30 seed or reads its discovery charts/certificates; only the terminal summary is read to enforce the success precondition.'
    checkpoint(OUT,result)
    print('PACKAGED V3',[(r['chain'],r['chart_status_counts']) for r in runs],sha(OUT))
if __name__=='__main__':main()
