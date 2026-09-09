#!/usr/bin/env python3
"""Six equation-only choices from retained scores, excluding measured/queued curves."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import certify_compact_r17_candidates as cert

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/'artifacts/local/elliptic-curves'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
POOLS=[LOCAL/n/'result.json' for n in ('r17-retained-extended-primes-v1',
    'r17-discarded-shards-extended-v1','r17-retention512-extended-v1')]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(folder,replay=False):
    if replay:
        protocol=json.loads((folder/'protocol.json').read_text())
        paths=[ROOT/p for p in protocol['inputs']]
    else:
        folder.mkdir(exist_ok=False)
        paths=sorted(set(POOLS+list(ART.glob('*results*.json'))+
            list(LOCAL.glob('*pari*/protocol.json'))+[ART/'compact_atlas_endpoints_v2.json']))
        protocol=dict(schema='fresh6-retained-r17-selection.v1',per_family=1,
            inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},
            sources={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),Path(cert.__file__)]},
            rule='One per family in the6144 saved H4096 rows, descending combined_selection_units, '
                 'descending combined_good, denominator, signed numerator. Exclude exact Q-isomorphs '
                 'of equations in frozen result ledgers/previous-equation lists/endpoints and '
                 'addresses or equations in frozen point-run rosters. No measured rank, point, '
                 'or validation-band score enters ordering; no new public catalogue is loaded.',
            scope='Finite six-address selection only, no point search or predicted rank. '
                  'Previously unattempted means absent from these named snapshots, not every possible external search.')
        (folder/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for category in ('inputs','sources'):
        if any(sha(ROOT/p)!=h for p,h in protocol[category].items()):
            raise ArithmeticError('frozen selection binding changed')
    data={p:json.loads(p.read_text()) for p in paths}
    pool=[];addresses={}
    for p in POOLS:
        for raw in data[p]['rows']:
            row={k:v for k,v in raw.items() if not k.startswith('validation_')}
            key=(row['family'],str(F(row['parameter'])))
            if key in addresses:raise ArithmeticError('duplicate score address')
            if max(abs(row['numerator']),row['denominator'])>4096:raise ArithmeticError('height outside fixed population')
            row['source_score_file']=str(p.relative_to(ROOT))
            addresses[key]=row;pool.append(row)
    if len(pool)!=6144:raise ArithmeticError('fixed score population differs')
    equations={};scheduled=set();singular_exclusions=[]
    def add(model,label):
        if not isinstance(model,list) or len(model)!=5:return
        q=tuple(map(F,model))
        if not cert.weierstrass_invariants(q)['discriminant']:
            singular_exclusions.append(label);return
        equations.setdefault(q,label)
    for path,d in data.items():
        if path in POOLS or not isinstance(d,dict):continue
        label=str(path.relative_to(ROOT))
        for name in ('previous_equations','curves','rows','roster','cases'):
            rows=d.get(name,[])
            if not isinstance(rows,list):continue
            for row in rows:
                if not isinstance(row,dict):continue
                add(row.get('curve',row.get('model')),label)
                if path.name=='protocol.json' and 'family' in row and 'parameter' in row:
                    try:key=(row['family'],str(F(row['parameter'])))
                    except (ValueError,ZeroDivisionError):continue
                    scheduled.add(key)
                    if key in addresses:add(addresses[key]['model'],label)
        add(d.get('curve'),label)
    def invariant(model):
        inv=cert.weierstrass_invariants(tuple(map(F,model)))
        if not inv['discriminant']:raise ArithmeticError('singular selection equation')
        return inv['c4']**3/inv['discriminant']
    by_j={}
    for model,label in equations.items():
        by_j.setdefault(invariant(model),[]).append((model,label))
    chosen=[];skips=[];available={}
    for family in sorted({r['family'] for r in pool}):
        ordered=sorted((r for r in pool if r['family']==family),key=lambda r:
            (-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
        accepted=[]
        for position,row in enumerate(ordered,1):
            model=tuple(map(F,row['model']));key=(family,str(F(row['parameter'])))
            matches=[label for q,label in by_j.get(invariant(model),[]) if cert.isomorphic(model,q)]
            aliases=[q['id'] for q in chosen if cert.isomorphic(model,tuple(map(F,q['model'])))]
            if key in scheduled or matches or aliases:
                skips.append(dict(family=family,parameter=row['parameter'],score_order=position,
                    scheduled_address=key in scheduled,matches=matches,selected_aliases=aliases))
                continue
            accepted.append((position,row))
        available[family]=len(accepted)
        if not accepted:raise ArithmeticError('no fresh address in family')
        position,row=accepted[0]
        chosen.append({**row,'id':family+'-fresh-001','score_order':position})
    if len(chosen)!=6:raise ArithmeticError('six-family roster required')
    result=dict(status='PASS_FROZEN_FRESH6_SELECTION',point_searches=0,
        protocol_sha256=sha(folder/'protocol.json'),pool_rows=len(pool),
        unique_exclusion_equations=len(equations),scheduled_addresses=len(scheduled),
        eligible_addresses_by_family=available,rows=chosen,skips=skips,
        singular_nonelliptic_rows_skipped=len(singular_exclusions),
        claim_boundary=protocol['scope'])
    if any(sha(ROOT/p)!=h for p,h in protocol['inputs'].items()):raise ArithmeticError('inputs changed during selection')
    target=folder/'result.json'
    if replay:
        if json.loads(target.read_text())!=result:raise ArithmeticError('selection replay differs')
        (folder/'verified.json').write_text(json.dumps(dict(status='PASS_EXACT_SELECTION_REPLAY',
            result_sha256=sha(target),source_sha256=sha(Path(__file__))),indent=2)+'\n')
    else:target.write_text(json.dumps(result,indent=2)+'\n')
    print('FRESH6',[(r['family'],r['parameter'],r['score_order']) for r in chosen],available,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,required=True)
    parser.add_argument('--replay',action='store_true')
    args=parser.parse_args();run(args.folder.resolve(),args.replay)
