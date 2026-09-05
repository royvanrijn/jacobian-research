#!/usr/bin/env python3
"""Independent Sage-free point certificates for new compact-atlas specializations."""
import argparse
import json
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec

ROOT=Path(__file__).resolve().parents[2]
PREVIOUS=('compact_r17_new_curves_v1.json','compact_r17_wide_new_curves_v1.json',
          'compact_r17_top64_interim_curves_v1.json','compact_r17_largest_gain_curve_v1.json')


def build(directory,output,minimum):
    if output.exists():raise FileExistsError('use a new immutable certificate path')
    families={r['family']:r for r in cert.read(spec.ATLAS)['families']}
    database=cert.read(cert.DATABASE);projection=[{'id':r['id'],'ainvs':r['ainvs']} for r in database['curves']]
    previous=[];previous_sources={}
    for name in PREVIOUS:
        path=ROOT/'artifacts/generated-results/elliptic-curves'/name;previous_sources[str(path.relative_to(ROOT))]=cert.hashed(path)
        previous.extend({'parameter':r['parameter'],'curve':r['curve']} for r in cert.read(path)['curves'])
    rows=[];seen=[]
    for path in sorted(directory.glob('*/candidate-*/result.json')):
        d=cert.read(path)
        if d['status'] not in ('COMPLETE_DECLARED_PILOT','TARGET_REACHED_PENDING_INDEPENDENT_REPLAY') or d.get('rank_lower_bound',0)<minimum:continue
        model=tuple(map(cert.F,d['curve']));points=tuple(tuple(map(cert.F,p)) for p in d['final_state']['state']['reductions']['points'])
        proof=cert.checked_rank(model,points);scale=spec.family_check(families[d['family']],d['parameter'],model,points)
        matches=[r['id'] for r in projection if cert.isomorphic(model,r['ainvs'])]
        own_matches=[r['parameter'] for r in previous if cert.isomorphic(model,r['curve'])]
        if matches or own_matches or any(cert.isomorphic(model,m) for m in seen):raise ArithmeticError('candidate is not new to the pinned comparison sets')
        if len(points)!=d['rank_lower_bound']:raise ArithmeticError('worker rank assertion changed')
        seen.append(model)
        rows.append({'family':d['family'],'parameter':d['parameter'],'curve':d['curve'],'points':[list(map(str,p)) for p in points],
            'rank_certificate':proof,'family_to_curve_scale_u':scale,'icarm_matches':matches,'previous_certificate_matches':own_matches,
            'discovery_witness':{'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path)},
            'completed_chart_records':len(d['charts']),'search_status':d['status']})
        print('CERTIFIED ATLAS',d['family'],d['parameter'],'rank >=',len(points),flush=True)
    paths=(Path(__file__).resolve(),Path(cert.__file__).resolve(),Path(spec.__file__).resolve(),spec.ATLAS,
           ROOT/'elliptic-curves/cas/mod2_reduction_independence.py',ROOT/'elliptic-curves/cas/elliptic_candidate_record.py')
    cert.write(output,{'schema':'elliptic-curves.compact-atlas-new-curves.v1','status':'PASS_EXACT_RANK_LOWER_BOUNDS','curves':rows,
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'protocol_sha256':cert.hashed(directory/'protocol.json'),'launch_sha256':cert.hashed(directory/'launch.json'),
        'novelty_snapshot':{'url':'https://elliptic-rank.icarm.cloud/database.json','date':'2026-09-05','raw_sha256':cert.hashed(cert.DATABASE),
            'equation_projection':projection,'curve_count':len(projection),'acknowledgement':'ICARM and NSF Grant DMS 2425401'},
        'previous_certificate_equations':previous,'previous_certificate_sources':previous_sources,
        'claim_boundary':'Exact lower bounds only; no exact ranks or conductor records. No match in the pinned public catalogue and fifteen earlier prospective certificates is not universal novelty.'})


def check(path):
    d=cert.read(path);families={r['family']:r for r in cert.read(spec.ATLAS)['families']}
    for name,h in d['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('certificate dependency changed: '+name)
    seen=[]
    for r in d['curves']:
        model=tuple(map(cert.F,r['curve']));points=tuple(tuple(map(cert.F,p)) for p in r['points']);old=r['rank_certificate']
        actual=cert.checked_rank(model,points,[s['prime'] for s in old['signatures']],old['no_rational_2_torsion_prime'])
        if json.dumps(actual,sort_keys=True)!=json.dumps(old,sort_keys=True):raise ArithmeticError('finite quotient certificate changed')
        if spec.family_check(families[r['family']],r['parameter'],model,points)!=r['family_to_curve_scale_u']:raise ArithmeticError('section transport changed')
        matches=[x['id'] for x in d['novelty_snapshot']['equation_projection'] if cert.isomorphic(model,x['ainvs'])]
        own=[x['parameter'] for x in d['previous_certificate_equations'] if cert.isomorphic(model,x['curve'])]
        if matches or own or matches!=r['icarm_matches'] or own!=r['previous_certificate_matches']:raise ArithmeticError('catalogue novelty comparison failed')
        if any(cert.isomorphic(model,x) for x in seen):raise ArithmeticError('duplicate exported curve')
        seen.append(model);print('REPLAYED ATLAS',r['family'],r['parameter'],'rank >=',len(points),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--directory',type=Path);g.add_argument('--check',type=Path)
    p.add_argument('--output',type=Path);p.add_argument('--minimum-rank',type=int,default=22);a=p.parse_args()
    check(a.check) if a.check else build(a.directory.resolve(),a.output,a.minimum_rank)
