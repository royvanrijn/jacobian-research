#!/usr/bin/env python3
"""Post-batch union of retained points on exactly matched prior curve equations."""
import argparse,sys
from pathlib import Path
from math import isqrt
import certify_compact_r17_candidates as cert
import certify_extended20_mw16_results as cohort
from replay_retention24_geometry import cloud_check
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=cohort.ROOT;ART=cohort.ART;D=ROOT/'artifacts/local/elliptic-curves/extended20-mw16-point-unions-v1';OUT=ART/'extended20_mw16_union_results_v1.json'

def points_at(address):
    path=ROOT/address
    if path.is_file():
        d=cert.read(path);return [(path,d)]
    name,family,t=address.split(':',2);path=ART/name;d=cert.read(path);rows=[r for r in d['curves'] if r['family']==family and r['parameter']==t]
    if len(rows)!=1:raise ArithmeticError('prior certificate address ambiguous')
    r=rows[0];records=[(path,r)]
    witness=r.get('discovery_witness',{})
    if isinstance(witness,dict) and 'path' in witness:
        p=ROOT/witness['path']
        if cert.hashed(p)!=witness['sha256']:raise ArithmeticError('old raw witness binding differs')
        records.append((p,cert.read(p)))
    return records

def transport(model,target,P):
    inv=cert.weierstrass_invariants(model);out=cert.weierstrass_invariants(target)
    if not inv['c4'] or not inv['c6']:raise ArithmeticError('exceptional j outside finite union transport')
    s2=(out['c6']/inv['c6'])/(out['c4']/inv['c4'])
    if s2<=0:raise ArithmeticError('nonrational positive scale')
    n,d=isqrt(s2.numerator),isqrt(s2.denominator)
    if n*n!=s2.numerator or d*d!=s2.denominator:raise ArithmeticError('nonsquare isomorphism scale')
    s=cert.F(n,d)
    if out['c4']!=s**4*inv['c4'] or out['c6']!=s**6*inv['c6'] or target[:3]!=(0,0,0):raise ArithmeticError('short isomorphism invariants differ')
    x,y=P;Q=(s*s*(x+inv['b2']/12),s**3*(y+(model[0]*x+model[2])/2))
    if not cert.is_on_weierstrass_curve(model,P) or not cert.is_on_weierstrass_curve(target,Q):raise ArithmeticError('point transport failed')
    return Q

def input_for(row):
    current=ROOT/row['discovery_witness']['path'];data=cert.read(current);target=tuple(map(cert.F,row['curve']));records=[(current,data)]
    # Catalogue points never enter this union. Known rediscoveries retain only their own search cloud.
    if not row['icarm_matches']:
        for address in row['previous_matches']:records+=points_at(address)
    unique={};raw=[]
    for path,d in records:
        model=tuple(map(cert.F,d['curve']))
        if not cert.isomorphic(model,target):raise ArithmeticError('prior matched model differs')
        points=d.get('points',d.get('final_state',{}).get('state',{}).get('reductions',{}).get('points',[]))
        if not points:raise ArithmeticError('prior witness contains no point list')
        for P in points:
            x,y=transport(model,target,tuple(map(cert.F,P)));raw.append({'x':str(x),'y':str(y)})
        for c in d.get('charts',[]):
            for P in c['search']['finite_curve_points']:
                x,y=transport(model,target,(cert.F(P['x']),cert.F(P['y'])));raw.append({'x':str(x),'y':str(y)})
        unique[str(path.relative_to(ROOT))]=cert.hashed(path)
    return {'status':'POINT_ONLY_EXACT_ISOMORPHISM_UNION_NOT_AN_ADMISSION_HISTORY','family':row['family'],'parameter':row['parameter'],'curve':row['curve'],'final_state':data['final_state'],'charts':[{'search':{'finite_curve_points':raw}}],'sources':{str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve()),**unique},'source_inputs':[{'path':p,'sha256':h} for p,h in sorted(unique.items())],'claim_boundary':'Exact Q-isomorphism transport of current and earlier retained point witnesses only, after terminal cohort verification. No catalogue points, point search, rank upper bound or absence claim.'}

def run_all():
    if OUT.exists() or (D/'ledger.json').exists():raise FileExistsError('preserve point-union results')
    source=cert.read(cohort.OUT)
    if cert.read(cohort.D/'verification-ledger.json')['status']!='PASS' or len(source['curves'])!=20:raise ArithmeticError('terminal verified20 cohort required')
    rows=[];ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    for row in source['curves']:
        folder=D/row['id'];ip=folder/'input.json';payload=input_for(row);checkpoint(ip,payload);cloud=ART/('extended20_mw16_'+row['id'].replace('-','_')+'_union_mod2_v1.json')
        for label,args in [('build',['--input',str(ip),'--input-sha256',cert.hashed(ip),'--output',str(cloud),'--prime-bound','997']),('check',['--check',str(cloud)])]:
            s=run([sys.executable,str(cohort.batch.CAS/'audit_recorded_point_mod2_rank_v3.py'),*args],limits=Limits(180,1610612736),log_path=folder/(label+'.log'),checkpoint_path=folder/(label+'.supervisor.json'),cwd=ROOT)
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('union cloud failed/censored')
        c=cert.read(cloud);r={**row,'points':c['independent_points'],'rank_certificate':c['rank_certificate'],'rank_lower_bound':c['rank_lower_bound'],'initial_cohort_rank_lower_bound':row['rank_lower_bound'],'point_union_input':{'path':str(ip.relative_to(ROOT)),'sha256':cert.hashed(ip)},'complete_cloud_certificate':{'path':str(cloud.relative_to(ROOT)),'sha256':cert.hashed(cloud)}};rows.append(r);ledger['rows'].append({'id':row['id'],'status':'PASS','source_records':len(payload['source_inputs']),'before':row['rank_lower_bound'],'after':c['rank_lower_bound'],'retained_points':len(c['points'])});checkpoint(D/'ledger.json',ledger);print('MW16 POINT UNION',row['id'],row['rank_lower_bound'],'->',c['rank_lower_bound'],'from',len(payload['source_inputs']),'records',flush=True)
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    result={**source,'schema':'elliptic-curves.extended20-mw16-retained-point-unions.v1','sources':{**source['sources'],str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve())},'curves':rows,'initial_cohort_certificate':{'path':str(cohort.OUT.relative_to(ROOT)),'sha256':cert.hashed(cohort.OUT)},'union_ledger_sha256':cert.hashed(D/'ledger.json'),'claim_boundary':'Exact point lower bounds from current plus previously retained points transported through checked Q-isomorphisms. Fixed20 initial-only exposure is recorded separately. Catalogue points excluded. Pinned593-catalogue absence is not universal novelty; no exact rank, conductor or upper bound.'};checkpoint(OUT,result)

def check():
    d=cert.read(OUT);initial=cert.read(cohort.OUT);families={f['fibration_id']:f for f in cert.read(cohort.spec.ATLAS)['families']}
    if d['initial_cohort_certificate']['sha256']!=cert.hashed(cohort.OUT) or d['union_ledger_sha256']!=cert.hashed(D/'ledger.json') or len(d['curves'])!=20 or d['sources']!={**initial['sources'],str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve())}:raise ArithmeticError('union source binding differs')
    for source,row in zip(initial['curves'],d['curves']):
        ip=ROOT/row['point_union_input']['path'];payload=input_for(source);cloud=cert.read(ROOT/row['complete_cloud_certificate']['path'])
        if cert.read(ip)!=payload or cert.hashed(ip)!=row['point_union_input']['sha256'] or cloud['input_sha256']!=cert.hashed(ip) or cert.hashed(ROOT/row['complete_cloud_certificate']['path'])!=row['complete_cloud_certificate']['sha256']:raise ArithmeticError('union provenance differs')
        if row['points']!=cloud['independent_points'] or row['rank_certificate']!=cloud['rank_certificate'] or row['initial_cohort_rank_lower_bound']!=source['rank_lower_bound']:raise ArithmeticError('union proof extraction differs')
        cloud_check(payload,payload['charts'],ROOT/row['complete_cloud_certificate']['path'],ip)
        cohort.verify(row,families,d['catalogue']['equations'],d['previous_equations'])
    print('REPLAYED20 MW16 RETAINED-POINT UNIONS AND EXACT PROOFS',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();check() if a.check else run_all()
