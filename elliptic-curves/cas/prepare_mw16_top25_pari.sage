#!/usr/bin/env sage-python
"""Re-reduce the same exact generic16 centres; no new metric or subgroup choice."""
import sys,argparse
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import mw16_top25_pari_followup as follow
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
mapper=SourceFileLoader('mw16_coordinate_mapper',str(CAS/'prepare_fresh_r17_pari_batch.sage')).load_module()
def main(index):
    p=follow.protocol();row=p['rows'][index];seed=cert.read(ROOT/row['seed_path']);out=follow.D/row['id']/'maps.json'
    if out.exists():raise FileExistsError('preserve MW16 coordinate maps')
    model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['generic_points']);data={'status':'RUNNING','protocol_hash':digest(p),'seed_sha256':row['seed_sha256'],'centres':seed['centres'],'rows':[]};checkpoint(out,data);mapper.pari.allocatemem(256000000,silent=True)
    if len(points)!=16:raise ArithmeticError('generic16 centre model required')
    for c in seed['centres']:
        if len(c['representative'])!=16 or any((c['representative'][j]-(c['mask']>>j))%2 for j in range(16)):raise ArithmeticError('old generic parity differs')
        data['rows'].append(mapper.mapping(model,points,c));checkpoint(out,data)
    data['status']='COMPLETE_DECLARED_MAPS';checkpoint(out,data);print('FROZEN SAME43 MW16 CENTRES WITH PARI COORDINATES',row['id'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,required=True);a=p.parse_args();main(a.index)
