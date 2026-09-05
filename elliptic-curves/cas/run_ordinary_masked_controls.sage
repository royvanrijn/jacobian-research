#!/usr/bin/env sage-python
"""Frozen exposure controls on the 31 already searched ordinary fibres.

Withhold section zero; select twelve centres from a fixed 256-mask geometry
sample of the remaining subgroup. No exceptional point or withheld point is
read by the search worker. All charts use the existing H=100000, 20s limits.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys

from sage.all import matrix, ZZ, QQ, pari
from fpylll import GSO, IntegerMatrix, Enumeration

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
sys.set_int_max_str_digits(0)
from search_observability import masked_control, point_visibility
from pointed_quartic_search import PointedQuarticSearch, sources
from research_runtime.deep_centres import exact_coset_minimum
from research_runtime.store import digest, checkpoint
from research_runtime.supervisor import capture, Limits

ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/ordinary-masked-controls-v1'
POP=ART/'fibre_height_population_v1.json.gz'
OLD=ART/'fibre_height_search_v1.json.gz'


def read(path):
    raw=path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix=='.gz' else raw)


def geometry(gram):
    g=matrix(QQ,gram); u=matrix(ZZ,pari(g).qflllgram()).transpose()
    assert abs(u.det())==1
    reduced=u*g*u.transpose(); n=g.nrows()
    integral=(reduced*reduced.denominator()).change_ring(ZZ)
    rows=[list(map(int,r)) for r in integral.rows()]
    gso=GSO.Mat(IntegerMatrix.from_matrix(rows),gram=True,float_type='dd',update=True)
    mu=[[gso.get_mu(i,j) if i>j else float(i==j) for j in range(n)] for i in range(n)]
    masks=[]; i=0
    while len(masks)<256:
        mask=int(digest(['ordinary-mask-v1',n,i]),16)%(1<<n); i+=1
        if mask and mask not in masks:masks.append(mask)
    candidates=[]
    for mask in masks:
        residue=[(mask>>i)&1 for i in range(n)]
        target=tuple(-sum(residue[i]*mu[i][j] for i in range(n))/2 for j in range(n))
        radius=sum(abs(v) for r in rows for v in r)/4+1
        _, coords=Enumeration(gso).enumerate(0,n,radius,0,target=target)[0]
        rep=[residue[i]+2*round(coords[i]) for i in range(n)]
        value, exact, nodes=exact_coset_minimum(rows,mask,rep,node_budget=200000)
        word=(matrix(ZZ,1,n,exact)*u).row(0)
        candidates.append(dict(mask=mask,minimum_norm=value,nodes=nodes,
            reduced_representative=list(exact),coefficients=list(map(int,word))))
    chosen=sorted(candidates,key=lambda r:(-r['minimum_norm'],r['mask']))[:12]
    return dict(sample=candidates,selected=chosen,reduced_gram=rows,
        change_of_basis=[list(map(int,r)) for r in u.rows()],
        description='Deepest twelve of 256 predetermined masks, not the complete deepest stratum.')


def freeze():
    path=LOCAL/'protocol.json'
    if path.exists():raise FileExistsError(path)
    pop,old=read(POP),read(OLD)
    tmpl=read(ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json')
    gram16=next(r['generic_height_gram'] for r in tmpl['presentations'] if r['presentation_id']=='a1-presentation-01')
    gram18=next(r['height_gram'] for r in read(ART/'mw18_generic_height_geometry_v1.json')['covers'] if r['cover_label']=='07ca9-orbit-08c1e')
    gr={'mw16':gram16,'mw18':gram18}; geometries={}; cases={}
    for name,family in pop['families'].items():
        geometries[name]=geometry([r[1:] for r in gr[name][1:]])
        for row in family['rows']:
            if row['id'] not in old['results']:continue
            certificate=old['results'][row['id']]['generic_independence']
            if certificate['status']!='CERTIFIED_INDEPENDENT' or certificate['rank']!=len(row['subgroup']):
                raise ArithmeticError('missing original independence certificate')
            blind,oracle=masked_control(row['search_model'],row['subgroup'],gr[name],[0])
            # Search file contains only retained points and geometry; the oracle
            # and full subgroup certificate are in a different file.
            blind.update(centres=[{'coefficients':r['coefficients']} for r in geometries[name]['selected']])
            checkpoint(LOCAL/'blind'/f"{row['id']}.json",blind)
            oracle.update(original_independence=certificate,original_subgroup=row['subgroup'])
            checkpoint(LOCAL/'oracles'/f"{row['id']}.json",oracle)
            cases[row['id']]=dict(family=name,blind_sha256=digest(blind),oracle_sha256=digest(oracle))
    assert len(cases)==31
    protocol=dict(schema='elliptic-curves.ordinary-masked-controls.v1',cases=cases,geometry=geometries,
        original_population_sha256=sha256(POP.read_bytes()).hexdigest(),
        original_search_sha256=sha256(OLD.read_bytes()).hexdigest(),
        withhold_zero_based=[0],coordinate_policy='metric:16',height=100000,seconds=20,
        charts_per_curve=12,workers=4,cell_wall_seconds=900,rss_bytes=2147483648,
        source_hashes={**sources(),str(Path(__file__).relative_to(ROOT)):sha256(Path(__file__).read_bytes()).hexdigest(),
            'elliptic-curves/cas/search_observability.py':sha256((ROOT/'elliptic-curves/cas/search_observability.py').read_bytes()).hexdigest()},
        endpoint='WITHHELD_KNOWN_DIRECTIONS_NOT_NEW_RANK',
        claim_boundary='31 distinct curves selected by the previous height/Nagao pilot. '
        'Changed reference subgroup; this does not calibrate all exceptional directions or an unbiased population incidence.')
    protocol['protocol_hash']=digest(protocol); checkpoint(path,protocol)
    print('FROZEN',len(cases),'curves',len(cases)*12,'charts',protocol['protocol_hash'],flush=True)


def protocol():
    p=read(LOCAL/'protocol.json'); h=p.pop('protocol_hash')
    if digest(p)!=h:raise ArithmeticError('protocol drift')
    p['protocol_hash']=h
    for path,expected in p['source_hashes'].items():
        if sha256((ROOT/path).read_bytes()).hexdigest()!=expected:raise ArithmeticError('source drift: '+path)
    return p


def cell(identifier,verify=False):
    p=protocol(); blind=read(LOCAL/'blind'/f'{identifier}.json')
    if digest(blind)!=p['cases'][identifier]['blind_sha256']:raise ArithmeticError('blind input drift')
    path=LOCAL/'cells'/f'{identifier}.json'
    result=read(path) if path.exists() else dict(id=identifier,protocol_hash=p['protocol_hash'],charts=[])
    if result['protocol_hash']!=p['protocol_hash']:raise ArithmeticError('cell protocol drift')
    for i,centre in enumerate(blind['centres']):
        search=PointedQuarticSearch(curve=blind['curve'],subgroup=blind['points'],centre=centre,coordinate_policy=p['coordinate_policy'])
        if i<len(result['charts']):
            search.verify_record(result['charts'][i]);continue
        if verify:raise ArithmeticError('missing chart')
        result['charts'].append(search.search(p['height'],p['seconds'],checkpoint_dir=LOCAL/'charts').record)
        checkpoint(path,result)
        print('MASKED',identifier,i+1,result['charts'][-1]['status'],flush=True)
    result['status']='COMPLETE' if all(c['status']=='bounded_search_complete' for c in result['charts']) else 'INCOMPLETE'
    if not verify:checkpoint(path,result)


def run(verify=False):
    p=protocol()
    def job(identifier):
        args=['sage','-python',str(Path(__file__).resolve()),'cell','--id',identifier]
        if verify:args.append('--verify')
        return identifier,capture(args,limits=Limits(wall_seconds=p['cell_wall_seconds'],rss_bytes=p['rss_bytes']),
            log_path=LOCAL/('replay-logs' if verify else 'logs')/f'{identifier}.log',cwd=ROOT)
    with ThreadPoolExecutor(max_workers=p['workers']) as pool:
        futures=[pool.submit(job,i) for i in sorted(p['cases'])]
        for f in as_completed(futures):
            identifier,_=f.result();print('CELL_COMPLETE',identifier,flush=True)


def audit():
    p=protocol(); results={}
    for identifier,meta in p['cases'].items():
        c=read(LOCAL/'cells'/f'{identifier}.json');oracle=read(LOCAL/'oracles'/f'{identifier}.json')
        if digest(oracle)!=meta['oracle_sha256']:raise ArithmeticError('oracle drift')
        if c.get('status') not in ('COMPLETE','INCOMPLETE') or len(c['charts'])!=12:raise ArithmeticError('open unfinished search')
        target=oracle['withheld_points'][0]; variants=[target,dict(x=target['x'],y=str(-Fraction(target['y'])))]
        visibility=[dict(chart=i,sign=1 if sign==0 else -1,**point_visibility(r,v))
            for i,r in enumerate(c['charts']) for sign,v in enumerate(variants)]
        found={digest(q):q for ch in c['charts'] for q in ch['finite_curve_points']}
        results[identifier]=dict(family=meta['family'],status=c['status'],charts=c['charts'],
            withheld_point=target,visibility=visibility,finite_point_count=len(found),
            direct_signed_recovery=any(v['status']=='VISIBLE_AND_RECORDED' for v in visibility),
            minimum_signed_representative_height=min(v['minimum_affine_height'] for v in visibility if v['minimum_affine_height'] is not None),
            claim_boundary='Direct signed oracle recovery is exact. Other recovered representatives need a group-law relation before being labelled withheld-direction recovery.')
        print('MASK_AUDIT',identifier,len(found),results[identifier]['direct_signed_recovery'],flush=True)
    output=dict(schema='elliptic-curves.ordinary-masked-controls-results.v1',protocol=p,results=results,
        status='COMPLETE' if all(c['status']=='COMPLETE' for c in results.values()) else 'INCOMPLETE')
    raw=(json.dumps(output,sort_keys=True,indent=2)+'\n').encode()
    (ART/'ordinary_masked_controls_v1.json.gz').write_bytes(gzip.compress(raw,mtime=0))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('freeze','run','cell','audit'));parser.add_argument('--id');parser.add_argument('--verify',action='store_true')
    a=parser.parse_args()
    if a.mode=='freeze':freeze()
    elif a.mode=='run':run(a.verify)
    elif a.mode=='cell':cell(a.id,a.verify)
    else:audit()
