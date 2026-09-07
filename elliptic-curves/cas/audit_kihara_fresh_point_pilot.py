#!/usr/bin/env python3
"""Full retained-cloud mod2/3/5 proofs after the fixed parent calibration."""
import sys
from pathlib import Path
import kihara_fresh_point_pilot as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits

def main():
    out=batch.BATCH/'verification-ledger.json';p=batch.protocol()
    if out.exists():raise FileExistsError('preserve full-cloud verification ledger')
    if cert.read(batch.BATCH/'ledger.json')['status']!='PASS':raise ArithmeticError('complete search and history replays required')
    ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger)
    for row in p['rows']:
        d=batch.BATCH/row['id'];source=d/'result.json'
        mod2=batch.ART/('kihara_fresh_point_pilot_'+row['id'].replace('-','_')+'_mod2_v1.json')
        modl=mod2.with_name(mod2.name.replace('_mod2_','_modl_'))
        jobs=[('mod2-build',['audit_recorded_point_mod2_rank_v3.py','--input',str(source),'--input-sha256',cert.hashed(source),'--output',str(mod2),'--prime-bound','997'],120),
          ('mod2-check',['audit_recorded_point_mod2_rank_v3.py','--check',str(mod2)],120),
          ('modl-build',['audit_retained_cloud_modl.py','--input',str(mod2),'--output',str(modl)],180),
          ('modl-check',['audit_retained_cloud_modl.py','--check',str(modl)],120)]
        stages=[]
        for label,args,seconds in jobs:
            s=run([sys.executable,str(batch.CAS/args[0]),*args[1:]],limits=Limits(seconds,1610612736),
              log_path=d/(label+'.log'),checkpoint_path=d/(label+'.supervisor.json'),cwd=batch.ROOT)
            stages.append({'name':label,'supervision':s})
            if s['outcome']!='completed' or s['returncode']!=0:
                ledger.update(status='FAILED_OR_CENSORED',failed_row=row['id'],failed_stages=stages);checkpoint(out,ledger)
                raise ArithmeticError('retained-cloud verification failed or censored')
        cloud=cert.read(mod2);odd=cert.read(modl)
        ledger['rows'].append({'id':row['id'],'status':'PASS','rank_lower_bound':cloud['rank_lower_bound'],
          'retained_points':len(cloud['points']),'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in odd['audits']},
          'mod2_certificate':str(mod2.relative_to(batch.ROOT)),'mod2_sha256':cert.hashed(mod2),
          'modl_certificate':str(modl.relative_to(batch.ROOT)),'modl_sha256':cert.hashed(modl),'stages':stages})
        checkpoint(out,ledger);print(row['id'],'VERIFIED',len(cloud['points']),'points; rank >=',cloud['rank_lower_bound'],flush=True)
    ledger['status']='PASS';checkpoint(out,ledger)

if __name__=='__main__':main()
