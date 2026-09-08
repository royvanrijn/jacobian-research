#!/usr/bin/env sage-python
"""Two fixed geometry feasibility fixtures; no point enumeration."""
import sys,time
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from half_lattice_pointed_sieve import linear_combination
from pointed_quartic_search import PointedQuarticSearch,point_record
from pari_pointed_backend import validate_map
mapper=SourceFileLoader('factor_free_map',str(CAS/'factor_free_pari_mapping.sage')).load_module()
D=ROOT/'artifacts/local/elliptic-curves/factor-free-map-probe-v1'

def main(index):
    p=cert.read(D/'protocol.json')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['bindings'].items()):raise ArithmeticError('probe binding changed')
    row=p['rows'][index];seed=cert.read(ROOT/row['seed']);old=cert.read(ROOT/row['maps']);centre=old['centres'][row['chart']]
    model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points'])
    out=D/(row['id']+'.json')
    if out.exists():raise FileExistsError('preserve fixed probe')
    started=time.monotonic();data={'status':'RUNNING','protocol_hash':digest(p),'phases':[]}
    def progress(stage):
        data['phases'].append({'stage':stage,'seconds':time.monotonic()-started});checkpoint(out,data);print(stage,flush=True)
    mapper.pari.allocatemem(256000000,silent=True)
    result=mapper.mapping(model,points,centre,progress)
    C=linear_combination(model,points,centre['representative'])
    search=PointedQuarticSearch(curve=model,subgroup=[],centre={'point':point_record(C)},coordinate_policy=result['coordinate_policy'])
    validate_map(search,result)
    data.update(status='PASS',mapping=result,seconds=time.monotonic()-started,
                reduced_coefficient_bits=max(abs(cert.F(v).numerator).bit_length() for v in result['reduced_P']))
    checkpoint(out,data);print('PASS',data['seconds'],data['reduced_coefficient_bits'],flush=True)

if __name__=='__main__':main(int(sys.argv[1]))
