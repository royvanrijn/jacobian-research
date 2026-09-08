#!/usr/bin/env python3
"""Independent full-cloud proofs for each terminal prospective specialized-parity trial."""
import sys
from pathlib import Path
import retained_native19_trial_v3 as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;ART=batch.ART;CAS=batch.CAS

def main():
    out=batch.BATCH/'verification-ledger.json'
    if out.exists():raise FileExistsError('preserve fixed full11952 follow-up cloud audit')
    p=batch.protocol();ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger)
    for index,row in enumerate(p['rows']):
        batch.configure(index);D=batch.D;input=D/'result.json';data=cert.read(input)
        for label in ('worker','replay'):
            s=cert.read(D/(label+'.supervisor.json'))
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('terminal worker and history replay required')
        mod2=ART/('retained_native19_trial_v3_'+row['id'].split('-')[-1]+'_mod2_v1.json');modl=mod2.with_name(mod2.name.replace('_mod2_','_modl_'))
        jobs=[('mod2-build',['audit_recorded_point_mod2_rank_v3.py','--input',str(input),'--input-sha256',cert.hashed(input),'--output',str(mod2),'--prime-bound','997'],180),('mod2-check',['audit_recorded_point_mod2_rank_v3.py','--check',str(mod2)],180),('modl-build',['audit_retained_cloud_modl.py','--input',str(mod2),'--output',str(modl)],300),('modl-check',['audit_retained_cloud_modl.py','--check',str(modl)],300)]
        for label,args,seconds in jobs:
            s=run([sys.executable,str(CAS/args[0]),*args[1:]],limits=Limits(seconds,1610612736),log_path=D/(label+'.log'),checkpoint_path=D/(label+'.supervisor.json'),cwd=ROOT)
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('full-cloud proof failed or censored')
        cloud=cert.read(mod2);ell=cert.read(modl)
        ledger['rows'].append({'id':row['id'],'input':str(input.relative_to(ROOT)),'input_sha256':cert.hashed(input),'status':'PASS','attempted_charts':len(data['charts']),'completed_boxes':sum(c['search']['status']=='bounded_search_complete' for c in data['charts']),'rank_lower_bound':cloud['rank_lower_bound'],'retained_points':len(cloud['points']),'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in ell['audits']},'mod2_certificate':str(mod2.relative_to(ROOT)),'mod2_sha256':cert.hashed(mod2),'modl_certificate':str(modl.relative_to(ROOT)),'modl_sha256':cert.hashed(modl)})
        checkpoint(out,ledger);print('AUDITED SPECIALIZED FULL11952',row['id'],len(data['charts']),cloud['rank_lower_bound'],len(cloud['points']),flush=True)
    ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':main()
