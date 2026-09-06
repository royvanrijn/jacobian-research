#!/usr/bin/env python3
"""Exact geometry and union-cloud audits of the fixed coordinate/height pilots."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import mw16_top25_pari_followup as mw
import new27_million_height_pilot as million
from replay_retention24_geometry import geometry,cloud_check,tuples
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=mw.ROOT;ART=mw.ART;CAS=mw.CAS

def inputs(kind,index):
    module=mw if kind=='mw16' else million;p=module.protocol();row=p['rows'][index];folder=module.D/row['id'];path=folder/'result.json';data=cert.read(path)
    if cert.read(module.D/'ledger.json')['status']!='PASS':raise ArithmeticError('terminal history replay required')
    paths=[ROOT/row['seed_path']]
    if kind=='million':paths.append(ROOT/'artifacts/local/elliptic-curves/retention-rank27-11952-113-adaptive-v1/result.json')
    paths.append(path);records=[cert.read(q) for q in paths]
    if any(d['curve']!=data['curve'] or d['parameter']!=data['parameter'] or d['family']!=data['family'] for d in records):raise ArithmeticError('union model differs')
    mp=folder/'maps.json' if kind=='mw16' else ROOT/row['maps_path'];maps=cert.read(mp);initial=tuples(data['initial_state']['state']['reductions']['points'])
    geometry(data,maps,initial,{**p,'maps_path':mp})
    sources={str(q.relative_to(ROOT)):cert.hashed(q) for q in [Path(__file__).resolve(),CAS/'replay_retention24_geometry.py',module.D/'protocol.json',module.D/'ledger.json',folder/'replay.supervisor.json',mp,*paths]}
    charts=[r for d in records for r in d['charts']]
    return p,row,folder,data,paths,charts,sources

def main(kind,index,check):
    p,row,folder,data,paths,charts,sources=inputs(kind,index);stem=('mw16_top25_' if kind=='mw16' else 'new27_million_')+row['id'].replace('-','_');out=folder/'all-retained-point-cloud-only.json';mod2=ART/(stem+'_all_retained_mod2_v1.json');modl=ART/(stem+'_all_retained_modl_v1.json');summary=ART/(stem+'_coverage_v1.json')
    payload={'status':'POINT_ONLY_CONCATENATION_NOT_AN_ADMISSION_TRANSCRIPT','family':data['family'],'parameter':data['parameter'],'curve':data['curve'],'final_state':data['final_state'],'charts':[{'search':{'finite_curve_points':r['search']['finite_curve_points']}} for r in charts],'sources':sources,'source_inputs':[{'path':str(q.relative_to(ROOT)),'sha256':cert.hashed(q)} for q in paths],'scope':'Union of point witnesses only. Original histories have separate replays; this is not a synthetic admission transcript.'}
    if check:
        if cert.read(out)!=payload:raise ArithmeticError('union inputs differ')
        cloud_check(data,charts,mod2,out)
        old=cert.read(summary)
    else:
        if out.exists() or summary.exists():raise FileExistsError('preserve followup cloud audit')
        checkpoint(out,payload)
        stages=[('mod2-build',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--input',str(out),'--input-sha256',cert.hashed(out),'--output',str(mod2),'--prime-bound','997'],180),('mod2-check',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(mod2)],180),('modl-build',[str(CAS/'audit_retained_cloud_modl.py'),'--input',str(mod2),'--output',str(modl)],240),('modl-check',[str(CAS/'audit_retained_cloud_modl.py'),'--check',str(modl)],240)]
        ledger={'status':'RUNNING','rows':[]};checkpoint(folder/'cloud-verification-ledger.json',ledger)
        for name,args,seconds in stages:
            s=run([sys.executable,*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=ROOT);ok=s['outcome']=='completed' and s['returncode']==0;ledger['rows'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(folder/'cloud-verification-ledger.json',ledger)
            if not ok:raise ArithmeticError('cloud audit failed/censored')
        ledger['status']='PASS';checkpoint(folder/'cloud-verification-ledger.json',ledger);cloud_check(data,charts,mod2,out)
    cloud=cert.read(mod2);odd=cert.read(modl)
    result={'schema':'elliptic-curves.fixed-coordinate-height-followup-coverage.v1','status':'PASS','sources':sources,'kind':kind,'id':row['id'],'family':data['family'],'parameter':data['parameter'],'height':p['height'],'seconds_per_chart':p['seconds_per_chart'],'attempted_charts':len(data['charts']),'completed_boxes':sum(r['search']['status']=='bounded_search_complete' for r in data['charts']),'union_chart_inputs':len(charts),'retained_point_count':len(cloud['points']),'mod2_lower_bound':cloud['rank_lower_bound'],'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in odd['audits']},'certificates':{str(q.relative_to(ROOT)):cert.hashed(q) for q in [out,mod2,modl,folder/'cloud-verification-ledger.json']},'claim_boundary':'Exact rational geometry and point-union provenance replayed. Separate finite-quotient checks certify lower bounds only. Recorded box coverage trusts pinned PARI execution. No exact rank, saturation, upper bound or absence of further points.'}
    if check:
        if old!=result:raise ArithmeticError('coverage summary differs')
    else:checkpoint(summary,result)
    print('AUDITED',kind,row['id'],result['completed_boxes'],'boxes',len(cloud['points']),'points',result['mod2_lower_bound'],result['odd_modulus_lower_bounds'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--kind',choices=['mw16','million'],required=True);p.add_argument('--index',type=int,default=0);p.add_argument('--check',action='store_true');a=p.parse_args();main(a.kind,a.index,a.check)
