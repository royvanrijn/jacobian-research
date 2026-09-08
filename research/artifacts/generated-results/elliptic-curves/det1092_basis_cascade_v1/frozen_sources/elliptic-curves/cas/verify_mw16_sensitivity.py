#!/usr/bin/env python3
"""Independent exact-map/group replay and deterministic sensitivity bundle.

The retained bundle includes the executed sources and all campaign/seed inputs.
Only verified rational points plus finite-reduction independence prove gains.
No bounded null result supplies an upper bound or a structural exclusion.
"""
import argparse
import ast
import base64
from fractions import Fraction as Q
import gzip
from hashlib import sha256
import json
from math import gcd,isqrt
from pathlib import Path
import sys
import subprocess

import mw16_sensitivity_backend as b
from alternate_quartic_covers import alternate_cover
from mod2_reduction_independence import (combined_mod2_rank,find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,mod2_reduction_signature)

ROOT=b.base.ROOT
FIXTURE='elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json'
PROSPECTIVE='artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_specializations_h300_v1.json'


def point(p): return Q(p['x']),Q(p['y'])
def canonical(p): return p[0],abs(p[1])
def digest(data): return sha256(data).hexdigest()


def check_group(model,generic,basis,classification,discoveries):
    if basis[:16]!=generic:
        raise ArithmeticError('this certificate requires the original generic prefix')
    if classification['status']!='PASS_BASIS_EQUALS_DISCOVERED_GROUP':
        raise ArithmeticError('unclassified discoveries cannot certify an exact recovered group')
    known={canonical(p) for p in generic}|{canonical(p) for p in discoveries}
    if any(canonical(p) not in known for p in basis):
        raise ArithmeticError('a claimed basis point was not recovered blindly')
    no_two=find_two_torsion_certificate_prime(model,prime_bound=1000)
    events=classification['events']
    if events and events[-1]['type']=='NEW_Q_INDEPENDENT_DIRECTION':
        saved=events[-1]['finite_reduction_certificate']['signatures']
        signatures=[]
        for row in saved:
            signature=mod2_reduction_signature(model,basis,row['prime'])
            if [list(v) for v in signature.rows]!=row['rows'] or signature.group_order!=row['group_order'] or signature.doubled_subgroup_order!=row['doubled_subgroup_order'] or signature.quotient_dimension!=row['quotient_dimension']:
                raise ArithmeticError('finite-reduction signature changed')
            signatures.append(signature)
    else:
        signatures=find_mod2_reduction_certificate(model,basis,prime_bound=1000)
    if combined_mod2_rank(signatures,len(basis))!=len(basis):
        raise ArithmeticError('claimed basis is not certified independent')
    accounted={canonical(p) for p in basis}
    for relation in classification['exact_integral_relations']:
        target=point(relation['point'])
        if b.base.linear_combination(model,basis,relation['coordinates'])!=target:
            raise ArithmeticError('integral relation failed exact group law')
        accounted.add(canonical(target))
    if any(canonical(p) not in accounted for p in discoveries):
        raise ArithmeticError('a returned point is missing from the classified group')
    return no_two


def check_chart(model,basis,record):
    if record.get('backend') == 'pointed_quartic_search_v1':
        from pointed_quartic_search import PointedQuarticSearch
        h=record['height_bound']
        if record['status']!='bounded_search_complete' or (record['denominator_start'],record['denominator_end'],record['completed_denominator'])!=(1,h,h):
            raise ArithmeticError('incomplete universal chart cannot certify a full control box')
        if len(record['representative'])!=len(basis) or sum((int(v)%2)<<i for i,v in enumerate(record['representative']))!=record['mask']:
            raise ArithmeticError('centre left its declared mod-two class')
        search=PointedQuarticSearch(model,basis,{'coefficients':record['representative']},record['specification'])
        return set(search.verify_record(record).curve_points)
    h=record['height_bound']
    if record['status']!='bounded_search_complete' or (record['denominator_start'],record['denominator_end'],record['completed_denominator'])!=(1,h,h):
        raise ArithmeticError('incomplete search cannot certify complete coverage')
    if record['integer_pairs_covered']!=h*(2*h+1) or not record['infinity_checked']:
        raise ArithmeticError('projective box census changed')
    if len(record['representative'])!=len(basis) or sum((int(v)%2)<<i for i,v in enumerate(record['representative']))!=record['mask']:
        raise ArithmeticError('centre left its declared mod-two class')
    centre=b.base.linear_combination(model,basis,record['representative'])
    if centre!=point(record['base_point']):
        raise ArithmeticError('centre is not the recorded blind basis combination')
    chart=b.base.make_chart(model,centre)
    if chart.record()!=record['pointed_chart']:
        raise ArithmeticError('pointed chart failed reconstruction')
    m=tuple(map(Q,record['horizontal_matrix'])); scale=Q(record['ordinate_scale'])
    if m[0]*m[3]==m[1]*m[2] or not scale:
        raise ArithmeticError('singular chart')
    f=tuple(map(int,record['coefficients']))
    if b.horizontal(chart,record['specification'])!=(m,f,scale):
        raise ArithmeticError('recorded coordinates differ from the declared sensitivity setting')
    aa,bb,cc,dd=b.multiply(chart.matrix,m)
    den,k,u=chart.denominator,chart.shift,chart.curve_scale
    raw=(Q(den*den*aa+k*cc,den)/u,Q(den*den*bb+k*dd,den)/u,cc,dd)
    expected=b.base.binary_transform(alternate_cover(model,centre).coefficients,raw)
    if tuple(v*u**4*scale**2/den**2 for v in expected)!=f:
        raise ArithmeticError('five-coefficient identity against original quartic failed')
    mapped=set(); finite=0; infinity=[]
    a,bb,c,d=m
    for n,q,r in (tuple(map(int,p)) for p in record['primitive_square_hits']):
        if q==0:
            if n!=1: raise ArithmeticError('noncanonical infinity')
            infinity.append((n,q,r))
        else:
            if not -h<=n<=h or not 1<=q<=h or gcd(n,q)!=1:
                raise ArithmeticError('hit outside primitive rational box')
            finite+=1
        if r<0 or r*r!=sum(v*n**i*q**(4-i) for i,v in enumerate(f)):
            raise ArithmeticError('square hit failed')
        for signed in {r,-r}:
            p=chart.map_point(a*n+bb*q,c*n+d*q,Q(signed)/scale)
            if p is not None:
                if not b.base.point_on_short_curve(model,p): raise ArithmeticError('point is off curve')
                mapped.add(p)
    expected_infinity=[(1,0,isqrt(f[4]))] if f[4]>=0 and isqrt(f[4])**2==f[4] else []
    if infinity!=expected_infinity or finite!=record['nonnegative_square_hits']:
        raise ArithmeticError('hit census failed')
    if mapped!=set(map(point,record['finite_curve_points'])):
        raise ArithmeticError('mapped point list changed')
    return mapped


def verify(payload,documents):
    def read(name): return json.loads(documents[name]['text'])
    fixture=read(FIXTURE)
    parents={p['parent_id']:p for p in fixture['parents']}
    mode=payload['declared_budget']['mode']
    if mode=='prospective':
        candidates=read(PROSPECTIVE)['candidates']
        if len(payload['results'])!=104 or {r['parent_id'] for r in payload['results']}!={c['candidate_id'] for c in candidates}:
            raise ArithmeticError('complete prospective certificate requires all 104 frozen fibres')
        parents={c['candidate_id']:{'target_short_model':c['raw_short_model'],
            'specialized_generic_points':c['raw_generic_points']} for c in candidates}
        gates=[read(n) for n in payload['inputs'] if n.endswith('.json') and n in documents and read(n).get('schema')=='elliptic-curves.mw16-sensitivity-prospective-gate.v1']
        if len(gates)!=1: raise ArithmeticError('missing independently replayed prospective gate')
        gate=gates[0]
        inputs=[read(n) for n in gate['inputs'] if n.endswith('.json') and n in documents]
        initial=[v for v in inputs if v.get('declared_budget',{}).get('mode')=='initial']
        adaptive=[v for v in inputs if v.get('declared_budget',{}).get('mode')=='adaptive']
        if len(initial)!=1 or len(adaptive)!=1: raise ArithmeticError('ambiguous calibrated pipeline')
        from build_mw16_sensitivity_gate import build
        expected=build(initial[0],adaptive[0],{'status':'PASS_EXACT_REPLAY'})
        if any(gate[k]!=expected[k] for k in expected): raise ArithmeticError('prospective gate differs from verified calibration')
        calibrated={(s['centre'],s['specification'],s['height']) for r in initial[0]['results'] for s in r['settings']}
        actual={(s['centre'],s['specification'],s['height']) for r in payload['results'] for s in r['settings']}
        if actual!=calibrated: raise ArithmeticError('prospective initial policy differs from its calibration')
    summaries=[]
    if payload['status']!='COMPLETE': raise ArithmeticError('campaign is incomplete')
    for result in payload['results']:
        parent=parents[result['parent_id']]
        model=tuple(map(Q,parent['target_short_model']))
        generic=tuple(map(point,parent['specialized_generic_points']))
        basis=tuple(map(point,result['search_basis']))
        prior=set()
        if mode=='adaptive':
            # Seed membership is checked from its frozen, independently replayed
            # initial campaign elsewhere in the same bundle.
            candidates=[read(n) for n in payload['inputs'] if n.endswith('.json') and n in documents and read(n).get('declared_budget',{}).get('mode')=='initial']
            seeds=[r for doc in candidates for r in doc.get('results',[]) if r.get('parent_id')==result['parent_id'] and 'current_basis' in r]
            if len(seeds)!=1: raise ArithmeticError('adaptive seed ambiguous or missing')
            prior=set(map(point,seeds[0]['current_basis']))
            if any(p not in prior for p in basis): raise ArithmeticError('adaptive centre uses a point outside the blind seed')
            trigger=payload['declared_budget'].get('adaptive_trigger')
            if trigger:
                gain=seeds[0]['exact_quotient_rank_recovered']
                active=trigger['initial_gain_minimum']<=gain<trigger['initial_gain_exclusive_maximum']
                if bool(result['settings'])!=active:
                    raise ArithmeticError('adaptive work differs from the frozen gain trigger')
        all_points=set(prior); settings=[]
        for setting in result['settings']:
            discovered=set()
            for record in setting['charts']:
                discovered.update(check_chart(model,basis,record))
            local_basis=tuple(map(point,setting['current_basis']))
            prime=check_group(model,generic,local_basis,setting['classification'],discovered)
            gain=len(local_basis)-16
            if gain!=setting['exact_quotient_rank_recovered']: raise ArithmeticError('setting score changed')
            all_points.update(discovered)
            settings.append({'centre':setting['centre'],'specification':setting['specification'],
                'height':setting['height'],'chart_count':len(setting['charts']),
                'exact_quotient_rank_recovered':gain,'no_two_torsion_prime':prime})
        final=tuple(map(point,result['current_basis']))
        prime=check_group(model,generic,final,result['classification'],all_points)
        if len(final)-16!=result['exact_quotient_rank_recovered']: raise ArithmeticError('union score changed')
        summaries.append({'parent_id':result['parent_id'],'curve_id':result['curve_id'],
            'exact_quotient_rank_recovered':len(final)-16,'settings':settings,'no_two_torsion_prime':prime})
        print(f'VERIFY_SENSITIVITY|parent={result["parent_id"]}|gain={len(final)-16}|QQ=PASS',flush=True)
    return {'mode':mode,'results':summaries,'total_quotient_rank':sum(r['exact_quotient_rank_recovered'] for r in summaries),
        'chart_count':sum(s['chart_count'] for r in summaries for s in r['settings']),
        'positive_candidate_count':sum(r['exact_quotient_rank_recovered']>0 for r in summaries)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--campaign',type=Path,action='append')
    parser.add_argument('--bundle',type=Path,required=True)
    parser.add_argument('--summary',type=Path,required=True)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.check:
        bundle=json.loads(gzip.decompress(args.bundle.read_bytes()))
    else:
        if not args.campaign: parser.error('--campaign required for packaging')
        files={}; roots=[]
        def include(path,recursive=False):
            name=str(path.resolve().relative_to(ROOT))
            if name in files: return
            data=path.read_bytes()
            files[name]={'sha256':digest(data),**({'base64':base64.b64encode(data).decode()} if name.endswith('.gz') else {'text':data.decode()})}
            if recursive:
                payload=json.loads(data)
                for n,expected in payload.get('inputs',{}).items():
                    source=ROOT/n
                    if digest(source.read_bytes())!=expected: raise ArithmeticError('changed executed input: '+n)
                    include(source,n.endswith('.json') and '/mw16-sensitivity/' in n)
        for path in args.campaign:
            roots.append(str(path.resolve().relative_to(ROOT))); include(path,True)
        include(Path(__file__))
        include(Path(__file__).with_name('replay_mw16_sensitivity.sh'))
        include(ROOT/'elliptic-curves/tests/test_mw16_sensitivity_backend.py')
        include(Path(__file__).with_name('build_mw16_sensitivity_gate.py'))
        include(Path(__file__).with_name('freeze_mw16_sensitivity_setting.py'))
        # Pin local transitive imports of the relation engine as well as the
        # sieve. Reading source does not execute any historical target loader.
        visited=set()
        while True:
            pending=[n for n in files if n.endswith(('.py','.sage')) and n not in visited]
            if not pending: break
            for name in pending:
                visited.add(name)
                tree=ast.parse(files[name]['text'])
                modules=[]
                for node in ast.walk(tree):
                    if isinstance(node,ast.Import): modules.extend(v.name for v in node.names)
                    elif isinstance(node,ast.ImportFrom) and node.module: modules.append(node.module)
                for module in modules:
                    for folder in (ROOT/'elliptic-curves/cas',ROOT/'elliptic-curves'):
                        path=folder/(module.replace('.','/')+'.py')
                        if path.is_file(): include(path)
        pari=subprocess.run(['gp','--version'],text=True,capture_output=True,check=True)
        bundle={'schema':'elliptic-curves.mw16-sensitivity-evidence.v1','campaigns':roots,'files':files,
            'software_at_verification':{'python':sys.version,'pari':(pari.stdout+pari.stderr).strip()}}
    for n,entry in bundle['files'].items():
        data=base64.b64decode(entry['base64']) if 'base64' in entry else entry['text'].encode()
        if digest(data)!=entry['sha256']: raise ArithmeticError('bundle content checksum changed')
        if n.endswith(('.py','.sage','.cpp','.sh')) and digest((ROOT/n).read_bytes())!=entry['sha256']:
            raise ArithmeticError('executed source differs from working tree: '+n)
    results=[]
    campaigns={}
    for name,entry in bundle['files'].items():
        if not name.endswith('.json'): continue
        data=json.loads(entry['text'])
        if data.get('schema')=='elliptic-curves.mw16-sensitivity.v1':
            campaigns[name]=data
    order={'initial':0,'adaptive':1,'prospective':2}
    for name,payload in sorted(campaigns.items(),key=lambda item:(order[item[1]['declared_budget']['mode']],item[0])):
        results.append({'campaign':name,**verify(payload,bundle['files'])})
    summary={'schema':'elliptic-curves.mw16-sensitivity-summary.v1','status':'PASS_EXACT_REPLAY',
        'campaigns':results,'software_at_verification':bundle['software_at_verification'],
        'claim_boundary':'Finite-coordinate searches certify recovered subgroups only; no bounded miss gives a rank upper bound.'}
    if not args.check:
        args.bundle.parent.mkdir(parents=True,exist_ok=True)
        args.bundle.write_bytes(gzip.compress((json.dumps(bundle,sort_keys=True)+'\n').encode(),mtime=0))
    summary['bundle_sha256']=digest(args.bundle.read_bytes())
    if args.check:
        if summary!=json.loads(args.summary.read_text()): raise ArithmeticError('summary changed')
    else:
        args.summary.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print('SENSITIVITY_CERTIFICATE|PASS_EXACT_REPLAY',flush=True)


if __name__=='__main__': main()
