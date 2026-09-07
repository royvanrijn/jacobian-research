#!/usr/bin/env python3
"""Finite mod2/3/5 proofs for old and new retained points on three new27 curves."""
import sys
from pathlib import Path
import new27_height125_followup as follow
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=follow.ROOT;D=follow.D;ART=follow.ART;CAS=follow.CAS
OLD=('paired_rank27_all_retained_mod2_v1.json','paired_second27_all_retained_mod2_v1.json','next24_rank27_all_retained_mod2_v1.json')
def main():
    p=follow.protocol();ledger=cert.read(D/'ledger.json')
    if ledger['status']!='PASS':raise ArithmeticError('terminal exact histories required')
    out=ART/'new27_height125_followup_coverage_v1.json'
    if out.exists():raise FileExistsError('preserve height125 coverage')
    rows=[]
    for index,row in enumerate(p['rows']):
        folder=D/row['id'];path=folder/'result.json';data=cert.read(path);old_path=ART/OLD[index];old=cert.read(old_path);cloudinput=folder/'all-retained-point-cloud-only.json'
        if cloudinput.exists():raise FileExistsError('preserve point-only input')
        if old['curve']!=data['curve'] or old['rank_lower_bound']!=27:raise ArithmeticError('old certified cloud differs')
        sources={str(q.relative_to(ROOT)):cert.hashed(q) for q in (Path(__file__).resolve(),path,old_path,D/'ledger.json',D/'protocol.json')}
        checkpoint(cloudinput,{'status':'POINT_ONLY_UNION_NOT_AN_ADMISSION_TRANSCRIPT','family':data['family'],'parameter':data['parameter'],'curve':data['curve'],'final_state':data['final_state'],'charts':[{'search':{'finite_curve_points':[{'x':x,'y':y} for x,y in old['points']]}}]+[{'search':{'finite_curve_points':c['search']['finite_curve_points']}} for c in data['charts']],'sources':sources,'scope':'Prior exact cloud and new raw points only; both chronological histories remain separate and replayed.'})
        mod2=ART/('new27_height125_'+row['id'].replace('-','_')+'_mod2_v1.json');modl=mod2.with_name(mod2.name.replace('_mod2_','_modl_'));stages=[('mod2-build','audit_recorded_point_mod2_rank_v3.py',['--input',str(cloudinput),'--input-sha256',cert.hashed(cloudinput),'--output',str(mod2),'--prime-bound','997'],120),('mod2-check','audit_recorded_point_mod2_rank_v3.py',['--check',str(mod2)],120),('modl-build','audit_retained_cloud_modl.py',['--input',str(mod2),'--output',str(modl)],180),('modl-check','audit_retained_cloud_modl.py',['--check',str(modl)],180)];records=[]
        for name,script,args,seconds in stages:
            s=run([sys.executable,str(CAS/script),*args],limits=Limits(seconds,1610612736),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=ROOT);records.append({'name':name,'status':'PASS' if s['outcome']=='completed' and s['returncode']==0 else 'FAILED_OR_CENSORED','supervision':s});checkpoint(folder/'cloud-verification.json',{'status':'RUNNING','rows':records})
            if records[-1]['status']!='PASS':raise ArithmeticError('height follow-up cloud stage failed/censored')
        checkpoint(folder/'cloud-verification.json',{'status':'PASS','rows':records});cloud=cert.read(mod2);odd=cert.read(modl);r={'id':row['id'],'sources':sources,'declared_charts':row['charts'],'completed_boxes':sum(c['search']['status']=='bounded_search_complete' for c in data['charts']),'point_count':len(cloud['points']),'mod2_lower_bound':cloud['rank_lower_bound'],'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in odd['audits']},'certificates':{str(q.relative_to(ROOT)):cert.hashed(q) for q in (mod2,modl)}};rows.append(r);checkpoint(out,{'status':'RUNNING','rows':rows});print('NEW27 HEIGHT125 FULL CLOUD',r['id'],r['point_count'],r['mod2_lower_bound'],r['odd_modulus_lower_bounds'],flush=True)
    checkpoint(out,{'schema':'elliptic-curves.new27-height125-followup-coverage.v1','status':'PASS','height':125000,'rows':rows,'protocol_sha256':cert.hashed(D/'protocol.json'),'claim_boundary':'Full retained-point finite lower bounds on three previously certified new27 curves after a fixed141-map height increase. No exact rank, saturation, point absence or universal novelty.'})
if __name__=='__main__':main()
