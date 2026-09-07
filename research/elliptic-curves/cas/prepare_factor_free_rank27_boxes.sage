#!/usr/bin/env sage-python
"""Exact maps and coordinate-box witnesses on a fixed retained centre roster."""
import sys
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import factor_free_rank27_box_audit as batch
from research_runtime.store import checkpoint,digest
from research_runtime.projective_box_change import classify
from pointed_quartic_search import PointedQuarticSearch,point_record
from half_lattice_pointed_sieve import linear_combination
import pari_pointed_backend as backend

def main(index):
    p=batch.protocol();r=p['rows'][index];old=batch.cert.read(ROOT/r['folder']/'maps.json');seed=batch.cert.read(ROOT/r['folder']/'seed.json');out=batch.D/r['id']/'maps.json'
    if out.exists():raise FileExistsError('preserve exact maps')
    model=tuple(map(batch.cert.F,seed['curve']));points=tuple(tuple(map(batch.cert.F,P)) for P in seed['points'])
    mapper=SourceFileLoader('factor_free_mapper',str(CAS/'factor_free_pari_mapping.sage')).load_module();mapper.pari.allocatemem(256000000,silent=True)
    data={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(out,data)
    for i,original in enumerate(old['rows']):
        m=mapper.mapping(model,points,original['centre']);C=linear_combination(model,points,m['centre']['representative'])
        for mapping in (original,m):
            search=PointedQuarticSearch(curve=model,subgroup=[],centre={'point':point_record(C)},coordinate_policy=mapping['coordinate_policy']);backend.validate_map(search,mapping)
        box=classify(original['matrix'],m['matrix'],p['height']);data['rows'].append({'index':i,'mapping':m,'box_comparison':box});checkpoint(out,data)
    data['status']='PASS';checkpoint(out,data);print('PASS49 exact map pairs',flush=True)

if __name__=='__main__':main(int(sys.argv[1]))
