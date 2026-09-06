#!/usr/bin/env python3
"""Render the audited class-target continuation and update its single status entry."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
D = ROOT/'artifacts/local/elliptic-curves/small-conductor-class-target-v1'
NOTE = ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_CLASS_TARGET_2026-09-06.md'
CHECKER = ROOT/'elliptic-curves/cas/pursue_small_conductor_class_target.sage'
ID = 'EC-SMALL-CONDUCTOR-CLASS-TARGET-20260906'


def main():
    paths = sorted(ART.glob('small_conductor_class_target_wave_*_v1.json'))
    if not paths:
        raise ValueError('no audited waves')
    rows = [json.loads(p.read_text()) for p in paths]
    for index,r in enumerate(rows,1):
        assert r['wave']==index and r['status']=='PASS'
        for name,h in r['sources'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    last = rows[-1]; dim = last['matrix']['quotient_dimension']
    lines = ['# MW16 at 3/17: reducing the class-group quotient toward 16', '',
        'Mathematical status: `MATH_STATUS.json`, entry `'+ID+'`.', '',
        'The current certified **GRH-conditional class-2-rank upper bound is '+str(dim)+'**.',
        'The unconditional curve-rank lower bound remains **22**. '+
        ('Matching bounds now prove **exact rank 22 under GRH**.' if dim<=16 else 'Exact rank remains unknown.'), '',
        'The [preceding descent study](SMALL_CONDUCTOR_DESCENT_SHORTCUT_2026-09-06.md) establishes',
        'the cubic field, exact local correction `rank <= dim Sel_2 <= g+7`, proved even',
        'Selmer parity, and the interval-certified generating cutoff 37,638 under GRH.',
        'Here `g=dim Cl(K)[2]`. A certified `g<=16` suffices for rank 22 under the same',
        'assumption. This criterion need not be attainable if the actual class-2-rank exceeds 16.', '',
        '## Authorized continuation', '',
        'The user explicitly set the goal of reducing this quotient to 16. Each wave fixes',
        'a finite target list, coefficient box, smoothness bound and single-worker resource',
        'limits before execution. Prime ideals are selected from free columns after exact',
        'supported row reduction, excluding unsuitable ramified or index-dividing primes.',
        'The Hessian-reduced index-prime lattices and previously searched boxes are replayed.',
        'Target checkpoints permit bounded resumes. A wave stops early only at dimension 16.', '',
        'Each saved witness is reconstructed and its principal ideal is factored exactly.',
        'The checker verifies the nonmonic norm identity including the fixed square factor.',
        'Relations may use primes up to the declared smoothness cutoff, never beyond the',
        'inherited complete 400,000 factor base. All coordinates above 37,638 are retained',
        'until exact elimination cancels them. The supported rank is independently checked',
        'as `rank(all relation rows) - rank(their outside projection)`.', '',
        'Only successful witnesses are required for the rank upper bound. Candidate counts',
        'and population digests describe the worker run; rejected values are not replayed',
        'and no exhaustive-search or smoothness-completeness assertion is made. Earlier',
        'full scalar sieve replays remain preserved as implementation calibration.', '',
        '## Audited waves', '',
        '| Wave | Targets completed | Box | Smooth bound | Candidate occurrences | Relation occurrences | Independent supported gain | Remaining dimension | Worker seconds |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r,path in zip(rows,paths):
        d = D/('wave_%03d'%r['wave']); p = json.loads((d/'protocol.json').read_text())
        seconds = sum(json.loads((ROOT/c['path']).read_text())['wall_seconds'] for c in r['chunks'])
        link = '../../'+str(path.relative_to(ROOT))
        lines.append('| ['+str(r['wave'])+']('+link+') | '+ ' | '.join(map(str,[len(r['chunks']),p['box'],p['smooth_bound'],r['candidate_occurrences'],r['relation_occurrences'],r['independent_supported_gain'],r['matrix']['quotient_dimension'],round(seconds,2)]))+' |')
    lines += ['', 'Starting dimension was **2,879**; the audited supported gain is **'+str(2879-dim)+'**.',
        '**'+str(max(0,dim-16))+' further independent supported rows** suffice for the target.',
        'The current coarse curve-rank upper bound is **'+str(last['conditional_on_grh_curve_rank_upper_bound'])+' under GRH**.',
        'A remaining dimension is an upper bound, not a proven number of actual class-group',
        'directions. It gives neither a new rational point nor an algebraic-rank parity claim.', '',
        '## Replay', '', 'The checker first replays the inherited proof chain, then every new principal',
        'relation, target selection, matrix transition and independent rank identity.', '',
        '```bash', 'sage -python elliptic-curves/cas/pursue_small_conductor_class_target.sage check --wave '+str(last['wave']), '```', '',
        'Protocols, successful witnesses, logs and supervisor outcomes are retained under',
        '`artifacts/local/elliptic-curves/small-conductor-class-target-v1/`.',
        'Every linked generated certificate pins its protocol, source and witness chunks.',
        'No new portable-archive replay is claimed until its separate report is available.', '']
    NOTE.write_text('\n'.join(lines))
    status = ROOT/'MATH_STATUS.json'; original = status.read_text(); data = json.loads(original)
    prior = next(r for r in data['entries'] if r['id']=='EC-SMALL-CONDUCTOR-SMALL-BASE-TARGETS-20260906')
    entry = json.loads(json.dumps(prior))
    entry.update({'id':ID,'title':'MW16 at 3/17: audited relation continuation reaches class quotient dimension '+str(dim),
        'scope':'The user-authorized goal is class quotient dimension at most16. '+str(last['wave'])+' bounded, checkpointed waves add '+str(2879-dim)+' independent supported relations to the baseline dimension2879. Every retained principal ideal and norm identity, deterministic free-column target selection, and supported matrix intersection is exactly replayed. All outside coordinates are retained until cancelled; rank(all)-rank(outside) verifies the intersection. Current GRH-conditional class-2-rank upper bound is '+str(dim)+' and the corresponding curve-rank upper bound is '+str(last['conditional_on_grh_curve_rank_upper_bound'])+'. '+('Matching the certified22-point lower bound gives exact rank22 under GRH.' if dim<=16 else 'Exact rank remains unknown; '+str(dim-16)+' further independent supported rows suffice for the target, whose attainability is not asserted.')+' Candidate counts and rejected-value digests are operational data, not an exhaustive miss certificate. New relations give no new point or unconditional upper bound. Local exact replay passes; a new isolated portable replay is not yet claimed.',
        'canonical_source':str(NOTE.relative_to(ROOT)), 'checker':str(CHECKER.relative_to(ROOT)),
        'dependencies':[prior['id']], 'independent_replay':False,
        'artifact_hash':'sha256:'+hashlib.sha256(CHECKER.read_bytes()).hexdigest(),
        'software_lock':[str(p.relative_to(ROOT)) for p in [CHECKER,Path(__file__).resolve(),NOTE,*paths]]})
    matches = [i for i,e in enumerate(data['entries']) if e['id']==ID]
    if matches: data['entries'][matches[0]]=entry
    else: data['entries'].append(entry)
    assert status.read_text()==original,'concurrent status update'
    status.write_text(json.dumps(data,indent=2)+'\n')
    print('Recorded',len(rows),'waves; dimension',dim)


if __name__=='__main__':
    main()
