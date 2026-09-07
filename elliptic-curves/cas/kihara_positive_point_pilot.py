#!/usr/bin/env python3
"""Fixed six-fibre, torsion-compatible pointed exposure; ranks checked afterward."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from pointed_quartic_search import PointedQuarticSearch,point_record
from half_lattice_pointed_sieve import linear_combination
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';BATCH=ROOT/'artifacts/local/elliptic-curves/kihara-positive-point-pilot-v1';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def sources():
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in ['kihara_positive_point_pilot.py','prepare_kihara_positive_point_maps.sage','prospective_half_lattice_v3.sage','factor_free_pari_mapping.sage']}}
def freeze():
    assert not (BATCH/'protocol.json').exists();bundle=ART/'kihara_positive_fibre_seed_bundle_v1.json';replay=ART/'kihara_positive_fibre_seed_replay_v1.json';data=cert.read(bundle);assert cert.read(replay)['status']=='PASS'
    rows=[];inputs={str(f.relative_to(ROOT)):cert.hashed(f) for f in [bundle,replay]}
    for r in data['rows']:
        old=r['original'];points=[old['points'][i] for i in r['basis_indices']];rank=len(points);row={'id':old['id'],'family':old['family'],'parameter':old['parameter'],'parent_ratio':old['parent_ratio'],'initial_rank':rank};rows.append(row)
        seed={**row,'curve':old['curve'],'points':points,'known_points':old['points'],'source_basis_indices':r['basis_indices'],'scale_from_parent':old['scale_from_parent']};path=BATCH/row['id']/'seed.json';assert not path.exists();checkpoint(path,seed);inputs[str(path.relative_to(ROOT))]=cert.hashed(path)
    checkpoint(BATCH/'protocol.json',{'schema':'kihara-positive-point-pilot.v1','rows':rows,'sources':sources(),'inputs':inputs,'sample_domain':'kihara-positive-point-pilot-v1','charts':49,'height':125000,'seconds_per_chart':10,'rank_stop':False,'maximum_boxes':294,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'geometry_wall_seconds':120,'worker_wall_seconds':600,'replay_wall_seconds':120,'rss_bytes':2147483648,'maximum_workers':1,
        'gate':'Three new Q-distinct parents and six fixed fibres have independently certified displayed specialized spans6,7,4,9,12,12. Rational2-torsion occurs on the first parent; use empty rank state with explicit rational centres, and certify all returned points after exposure at odd primes. No adaptive rank gate or claim is used to schedule charts.',
        'selection':'Same six fixed parent/base pairs as the intake. No score, public target, validation-prime input, replacement or outcome-driven selection. All six map files precede every point search.',
        'centre_policy':'All nonzero parity masks when fewer than2048, otherwise2048 fixed SHA256 masks.384-bit canonical heights,10^6 rounding, unimodular LLL, numerical CVP with exact parity/norm transport. Choose49 largest computed-norm representatives after normalizing sign. If fewer than49 masks exist (the rank4 seed), fill deterministically from the same representatives plus signed twice-basis translations, excluding equal words up to sign. This adds chart representatives, not independent parity classes or rank directions. Exact coordinate maps precede search.',
        'rank_policy':'No production MWState admission: rational2-torsion makes its no-two-torsion gate inapplicable. The pointed backend uses an empty valid state and an exactly checked explicit centre. Retain original known points and every returned point; independently certify odd-prime full-cloud lower bounds afterward. Exact seed-span upper bounds distinguish discovered gains from old point reuse.',
        'failure_policy':'Retain failures and completed checkpoints. No automatic refill, larger boxes, second wave or parameter scan. Any map failure prevents all point attempts.','following_campaign':None})
    print('FROZEN6 FIBRES;294 BOXES;RANKS',[r['initial_rank'] for r in rows],flush=True)
def protocol():
    p=cert.read(BATCH/'protocol.json');assert p['sources']==sources()
    for n,h in p['inputs'].items():assert cert.hashed(ROOT/n)==h
    return p
def configure(i):
    global ROW,D,SEED
    ROW=protocol()['rows'][i];D=BATCH/ROW['id'];SEED=D/'seed.json'
def masks(p):
    rank=ROW['initial_rank'];n=(1<<rank)-1
    if n<=2048:return list(range(1,n+1))
    result=[];i=0
    while len(result)<2048:
        m=int(digest([p['sample_domain'],i]),16)%(1<<rank);i+=1
        if m and m not in result:result.append(m)
    return result
def centres(sample,g):
    order=sorted(sample,key=lambda r:(-r['metric_norm'],r['parity']));result=[];seen=set();rank=len(g)
    def admit(r,offset):
        w=[a+b for a,b in zip(r['representative'],offset)];nz=next((a for a in w if a),0)
        if not nz:return
        sign=1 if nz>0 else -1;w=tuple(sign*a for a in w)
        if w in seen:return
        seen.add(w);norm=sum(w[i]*g[i][j]*w[j] for i in range(rank) for j in range(rank))
        result.append({'parity':r['parity'],'representative':list(w),'metric_norm':norm,'origin_representative':r['representative'],'offset':offset,'sign':sign})
    for r in order:
        admit(r,[0]*rank)
        if len(result)==49:return result
    for axis in range(rank):
        for sign in (-1,1):
            for r in order:
                offset=[0]*rank;offset[axis]=2*sign;admit(r,offset)
                if len(result)==49:return result
    raise ArithmeticError('fixed translated representatives do not supply49 centres')
def search_for(seed,m):
    model=tuple(map(cert.F,seed['curve']));points=[tuple(map(cert.F,P)) for P in seed['points']];centre=linear_combination(model,points,m['centre']['representative']);assert centre is not None
    return PointedQuarticSearch(curve=model,subgroup=[],centre={'point':point_record(centre)},coordinate_policy=m['coordinate_policy'])
def worker():
    p=protocol();seed=cert.read(SEED);maps=cert.read(D/'maps.json');assert maps['status']=='COMPLETE_DECLARED_MAPS' and maps['protocol_hash']==digest(p) and len(maps['rows'])==49
    out=D/'result.json';assert not out.exists();cloud=[];seen=set()
    def admit(P):
        x,y=map(cert.F,P);key=(x,abs(y))
        if key not in seen:seen.add(key);cloud.append([str(x),str(y)])
    for P in seed['points']+seed['known_points']:admit(P)
    data={**seed,'status':'RUNNING','protocol_hash':digest(p),'maps_sha256':cert.hashed(D/'maps.json'),'charts':[],'retained_points':cloud};checkpoint(out,data)
    for i,m in enumerate(maps['rows']):
        search=search_for(seed,m);record,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256'])
        for P in points:admit(P)
        data['charts'].append({'index':i,'centre':m['centre'],'search':record});checkpoint(out,data);print(ROW['id'],i+1,record['status'],'retained',len(cloud),flush=True)
    data['status']='COMPLETE_DECLARED_POINT_ATTEMPT';checkpoint(out,data)
def replay():
    p=protocol();seed=cert.read(SEED);maps=cert.read(D/'maps.json');data=cert.read(D/'result.json');assert data['protocol_hash']==digest(p) and data['maps_sha256']==cert.hashed(D/'maps.json') and len(data['charts'])==49
    cloud=[];seen=set()
    def admit(P):
        x,y=map(cert.F,P);key=(x,abs(y))
        if key not in seen:seen.add(key);cloud.append([str(x),str(y)])
    for P in seed['points']+seed['known_points']:admit(P)
    for i,(m,r) in enumerate(zip(maps['rows'],data['charts'])):
        assert r['index']==i and r['centre']==m['centre'];record=r['search'];assert record['height_bound']==p['height'] and record['timeout_seconds']==p['seconds_per_chart'] and record['gp_binary_sha256']==p['gp_sha256']
        for P in backend.replay(search_for(seed,m),m,record):admit(P)
    assert cloud==data['retained_points'];print('PASS',ROW['id'],'49 exact charts and full retained cloud',len(cloud),flush=True)
def launch():
    p=protocol();out=BATCH/'ledger.json';assert not out.exists();ledger={'status':'RUNNING_MAPS','maps':[],'rows':[]};checkpoint(out,ledger)
    for i,row in enumerate(p['rows']):
        configure(i);s=run([SAGE,str(CAS/'prepare_kihara_positive_point_maps.sage'),'--index',str(i)],limits=Limits(p['geometry_wall_seconds'],p['rss_bytes']),log_path=D/'maps.log',checkpoint_path=D/'maps.supervisor.json',cwd=ROOT);ledger['maps'].append({'id':row['id'],'supervision':s});checkpoint(out,ledger);assert s['outcome']=='completed' and s['returncode']==0
    ledger['status']='RUNNING_POINTS';checkpoint(out,ledger)
    for i,row in enumerate(p['rows']):
        configure(i);entry={'id':row['id'],'stages':[]};ledger['rows'].append(entry)
        for stage in ('worker','replay'):
            s=run([sys.executable,str(Path(__file__).resolve()),stage,'--index',str(i)],limits=Limits(p[stage+'_wall_seconds'],p['rss_bytes']),log_path=D/(stage+'.log'),checkpoint_path=D/(stage+'.supervisor.json'),cwd=ROOT);entry['stages'].append({'stage':stage,'supervision':s});checkpoint(out,ledger);print(row['id'],stage,s['outcome'],s['wall_seconds'],flush=True);assert s['outcome']=='completed' and s['returncode']==0
        entry['status']='PASS';checkpoint(out,ledger)
    ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','launch','worker','replay']);p.add_argument('--index',type=int);a=p.parse_args()
    if a.stage in ('worker','replay'):configure(a.index)
    globals()[a.stage]()
