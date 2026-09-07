#!/usr/bin/env sage-python
"""Retrospective single-chart known28 control, without exceptional worker points."""
import sys,time
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from pointed_quartic_search import PointedQuarticSearch,point_record
from half_lattice_pointed_sieve import linear_combination
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest
D=ROOT/'artifacts/local/elliptic-curves/factor-free-known28-control-v1'
OLD=ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/factor_free_known28_control_v1.json'

def main():
    if OUT.exists() or (D/'protocol.json').exists():raise FileExistsError('preserve control')
    paths=[Path(__file__).resolve(),CAS/'factor_free_pari_mapping.sage',CAS/'memory_rank_certificate.py',OLD/'seed.json',OLD/'maps.json']
    p={'sources':{**backend.sources(),**{str(x.relative_to(ROOT)):cert.hashed(x) for x in paths}},
       'chart_index_zero_based':4,'height':125000,'seconds':10,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),
       'scope':'One retrospectively chosen successful historical chart, using only the old27 seed and its frozen centre. No exceptional coordinate input;10-second new factor-free coordinate box. A positive control, not blind selection or new rank discovery.'}
    checkpoint(D/'protocol.json',p)
    seed=cert.read(OLD/'seed.json');old=cert.read(OLD/'maps.json');model=tuple(map(cert.F,seed['curve']));points=tuple(tuple(map(cert.F,P)) for P in seed['points'])
    mapper=SourceFileLoader('factor_free_map',str(CAS/'factor_free_pari_mapping.sage')).load_module();mapper.pari.allocatemem(256000000,silent=True)
    started=time.monotonic();m=mapper.mapping(model,points,old['centres'][4]);checkpoint(D/'map.json',m)
    C=linear_combination(model,points,m['centre']['representative'])
    search=PointedQuarticSearch(curve=model,subgroup=[],centre={'point':point_record(C)},coordinate_policy=m['coordinate_policy'])
    r,found=backend.execute(search,m,125000,10,p['gp_sha256']);checkpoint(D/'search.json',r)
    if backend.replay(search,m,r)!=found:raise ArithmeticError('exact transcript replay differs')
    cloud=list(points);seen={(x,abs(y)) for x,y in points}
    for P in found:
        if (P[0],abs(P[1])) not in seen:seen.add((P[0],abs(P[1])));cloud.append(P)
    proof=seed['rank_certificate'];rank=checked_rank(model,cloud,[a['prime'] for a in proof['signatures']],proof['no_rational_2_torsion_prime'])
    checkpoint(OUT,{'status':'PASS' if r['status']=='bounded_search_complete' and rank['rank_lower_bound']>=28 else 'CONTROL_NOT_RECOVERED',
        'protocol_hash':digest(p),'sources':p['sources'],'inputs':{str(x.relative_to(ROOT)):cert.hashed(x) for x in (D/'protocol.json',D/'map.json',D/'search.json')},
        'curve':seed['curve'],'points':[list(map(str,P)) for P in cloud],'rank_certificate':rank,
        'rank_lower_bound':rank['rank_lower_bound'],'completed_boxes':int(r['status']=='bounded_search_complete'),
        'seconds':time.monotonic()-started,'boundary':p['scope']})
    print(cert.read(OUT)['status'],rank['rank_lower_bound'],len(cloud),flush=True)
if __name__=='__main__':main()
