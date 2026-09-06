#!/usr/bin/env python3
"""Complete-cloud certificates and coverage of the known28 specialized-parity control."""
from pathlib import Path
import sys,argparse
import certify_compact_r17_candidates as cert
import native28_specialized_parity_control as follow
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=follow.ROOT;ART=follow.ART;CAS=follow.CAS

def main():
    D=follow.D
    output=D/'all-retained-point-cloud-only.json'
    if output.exists():raise FileExistsError('preserve point-only cloud input')
    p=follow.protocol();check=cert.read(D/'replay.supervisor.json');terminal=cert.read(D/'worker.supervisor.json');data=cert.read(D/'result.json');seed=cert.read(follow.SEED)
    if check['outcome']!='completed' or check['returncode']!=0 or terminal['outcome']!='completed' or terminal['returncode']!=0:raise ArithmeticError('adaptive terminal and exact history replay required')
    rows=seed['charts']+data['charts'];sources={str(path.relative_to(ROOT)):cert.hashed(path) for path in (Path(__file__).resolve(),follow.SEED,D/'result.json',D/'protocol.json',D/'maps.json',D/'replay.supervisor.json',D/'worker.supervisor.json')}
    checkpoint(output,{'status':'POINT_ONLY_CONCATENATION_NOT_AN_ADMISSION_TRANSCRIPT','family':data['family'],'parameter':data['parameter'],'curve':data['curve'],'final_state':data['final_state'],'charts':[{'search':{'finite_curve_points':r['search']['finite_curve_points']}} for r in rows],'sources':sources,'source_inputs':[{'path':str(path.relative_to(ROOT)),'sha256':cert.hashed(path)} for path in (follow.SEED,D/'result.json')],'scope':'Only raw point witnesses are concatenated; initial and adaptive histories each have their own exact replayer. This is not a synthetic chronological transcript.'})
    mod2=ART/('native28_specialized_parity_all_retained_mod2_v1.json');modl=ART/('native28_specialized_parity_all_retained_modl_v1.json');stages=[('mod2-build',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--input',str(output),'--input-sha256',cert.hashed(output),'--output',str(mod2),'--prime-bound','997'],120),('mod2-check',[str(CAS/'audit_recorded_point_mod2_rank_v3.py'),'--check',str(mod2)],120),('modl-build',[str(CAS/'audit_retained_cloud_modl.py'),'--input',str(mod2),'--output',str(modl)],180),('modl-check',[str(CAS/'audit_retained_cloud_modl.py'),'--check',str(modl)],180)];ledger={'status':'RUNNING','rows':[]};checkpoint(D/'cloud-verification-ledger.json',ledger)
    for name,args,seconds in stages:
        s=run([sys.executable,*args],limits=Limits(seconds,1073741824),log_path=D/(name+'.log'),checkpoint_path=D/(name+'.supervisor.json'),cwd=ROOT);ledger['rows'].append({'name':name,'status':'PASS' if s['outcome']=='completed' and s['returncode']==0 else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'cloud-verification-ledger.json',ledger)
        if ledger['rows'][-1]['status']!='PASS':raise ArithmeticError('complete cloud stage failed/censored')
    ledger['status']='PASS';checkpoint(D/'cloud-verification-ledger.json',ledger);cloud=cert.read(mod2);odd=cert.read(modl)
    checkpoint(ART/('native28_specialized_parity_adaptive_coverage_v1.json'),{'schema':'elliptic-curves.native28-specialized-parity-coverage.v1','sources':sources,'initial_charts':len(seed['charts']),'adaptive_declared_charts':p['charts'],'adaptive_attempted_charts':len(data['charts']),'adaptive_completed_boxes':sum(r['search']['status']=='bounded_search_complete' for r in data['charts']),'height':p['height'],'retained_point_count':len(cloud['points']),'mod2_lower_bound':cloud['rank_lower_bound'],'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in odd['audits']},'point_certificates':{str(path.relative_to(ROOT)):cert.hashed(path) for path in (mod2,modl)},'cloud_verification_ledger_sha256':cert.hashed(D/'cloud-verification-ledger.json'),'claim_boundary':'Recorded completed box coverage trusts the pinned PARI executions. Every attempted adaptive history replays; a target stop can leave part of the49-chart roster unsearched. Full retained-cloud finite quotient certificates are lower bounds only; no exact rank, saturation or point-absence conclusion.'})
    print('AUDITED KNOWN28 SPECIALIZED CONTROL',len(cloud['points']),'points; mod2',cloud['rank_lower_bound'],'odd',[(a['modulus'],a['finite_column_rank']) for a in odd['audits']],flush=True)
if __name__=='__main__':
    main()
