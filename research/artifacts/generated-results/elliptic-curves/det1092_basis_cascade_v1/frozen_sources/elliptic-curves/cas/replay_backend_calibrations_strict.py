#!/usr/bin/env python3
"""Strengthen calibration replay with exact rosters and raw-program bindings."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
import compare_native11952_search_engines as native
import compare_new_rank26_search_engines_v2 as fresh
import check_pari_height_boundaries as boundary
ROOT=Path(__file__).resolve().parents[2];OUTPUT=ROOT/'artifacts/generated-results/elliptic-curves/strict_backend_calibration_replay_v1.json'


def roster(rows,count):
    if [r['chart'] for r in rows]!=list(range(count)):raise ArithmeticError('missing, repeated or reordered chart index')


def expected():
    inputs={};results=[]
    for module,count in ((native,49),(fresh,43)):
        data=cert.read(module.OUTPUT);protocol=cert.read(module.D/'protocol.json');roster(data['rows'],count)
        if data['protocol_hash']!=digest(protocol):raise ArithmeticError('calibration protocol differs')
        module.check();maps=cert.read(module.MAPS);source=cert.read(module.INPUT)
        hits=0
        for row in data['rows']:
            i=row['chart'];m=maps['rows'][i];raw_path=module.D/f'chart-{i:02}.json';raw=cert.read(raw_path);poly=lambda a:'+'.join(f'({v})*x^{j}' for j,v in enumerate(a));program='C=['+poly(m['reduced_P'])+','+poly(m['reduced_Q'])+'];gettime();R=hyperellratpoints(C,100000);print("MS|",gettime());for(i=1,#R,print("X|",R[i][1]));quit\n'
            xs={cert.F(line[2:]) for line in raw['stdout'].splitlines() if line.startswith('X|')};old={cert.F(int(n),int(d)) for n,d,_ in source['charts'][i]['search']['primitive_square_hits'] if int(d)}
            if raw['program']!=program or raw['returncode']!=0 or '***' in raw['stderr'] or row['pari_hit_count']!=len(xs) or row['gmp_hit_count']!=len(old) or row['pari_only'] or row['gmp_only'] or xs!=old:raise ArithmeticError('raw calibration program or counts differ')
            hits+=len(xs);inputs[str(raw_path.relative_to(ROOT))]=cert.hashed(raw_path)
        inputs.update({str(p.relative_to(ROOT)):cert.hashed(p) for p in (module.OUTPUT,module.INPUT,module.MAPS,module.D/'protocol.json')});results.append({'charts':count,'affine_square_coordinates':hits,'status':'PASS'})
    data=cert.read(boundary.OUTPUT)
    wanted=[(list(c),h) for c in ((-56,0,0,0,1),(-1,3,0,0,81)) for h in (2,3)]
    if data['status']!='PASS' or [(r['coefficients'],r['height']) for r in data['rows']]!=wanted:raise ArithmeticError('boundary roster differs')
    for name,h in data['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('boundary source differs')
    for r in data['rows']:
        poly='+'.join(f'({a})*x^{i}' for i,a in enumerate(r['coefficients']));program=f'R=hyperellratpoints({poly},{r["height"]});for(i=1,#R,print("X|",R[i][1]));quit\n'
        if r['program']!=program or r['returncode']!=0 or '***' in r['stderr']:raise ArithmeticError('boundary program differs')
    boundary.main(True);inputs[str(boundary.OUTPUT.relative_to(ROOT))]=cert.hashed(boundary.OUTPUT)
    return {'schema':'elliptic-curves.strict-backend-calibration-replay.v1','status':'PASS','calibrations':results,'boundary_boxes':4,'inputs':inputs,'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(native.__file__).resolve(),Path(fresh.__file__).resolve(),Path(boundary.__file__).resolve())},'claim_boundary':'Checks exact complete chart rosters and raw program/count bindings in addition to prior square-set replays. Previous certificates and failed attempts remain immutable; no new search or universal engine equivalence claim.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();data=expected()
    if a.check:
        if cert.read(OUTPUT)!=data:raise ArithmeticError('strict replay output differs')
    else:
        if OUTPUT.exists():raise FileExistsError('preserve strict replay')
        checkpoint(OUTPUT,data)
    print('STRICT BACKEND CALIBRATIONS PASS')
