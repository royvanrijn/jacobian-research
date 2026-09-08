#!/usr/bin/env python3
"""Bounded whole-cloud and independent geometric/rank checks of one frozen run."""
import argparse,sys,shutil
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def main(d,protocol,prefix):
 out=d/'certification-ledger.json';assert not out.exists()
 source=d/'result.json';mod2=ART/(prefix+'_mod2_v1.json');modl=ART/(prefix+'_modl_v1.json')
 names=['audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py','verify_factor_free_exposure.sage','verify_factor_free_rank.sage','certify_factor_free_exposure.py']
 checkpoint(d/'certification-protocol.json',dict(inputs={str(x.relative_to(ROOT)):cert.hashed(x) for x in [source,protocol,d/'maps.json',d/'seed.json']},sources={str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names},maximum_workers=1,rss_bytes=2147483648,seconds_per_stage=180))
 jobs=[('mod2-build',[sys.executable,str(CAS/names[0]),'--input',str(source),'--input-sha256',cert.hashed(source),'--output',str(mod2),'--prime-bound','997']),
 ('mod2-check',[sys.executable,str(CAS/names[0]),'--check',str(mod2)]),
 ('modl-build',[sys.executable,str(CAS/names[1]),'--input',str(mod2),'--output',str(modl)]),
 ('modl-check',[sys.executable,str(CAS/names[1]),'--check',str(modl)]),
 ('geometry',[SAGE,str(CAS/names[2]),'--run',str(d),'--protocol',str(protocol),'--cloud',str(mod2)])]
 ledger={'status':'RUNNING','stages':[]};checkpoint(out,ledger)
 def stage(name,cmd,cwd):
  s=run(cmd,limits=Limits(180,2147483648),log_path=d/(name+'-certificate.log'),checkpoint_path=d/(name+'-certificate.supervisor.json'),cwd=cwd)
  ok=s['outcome']=='completed' and s['returncode']==0
  ledger['stages'].append(dict(name=name,status='PASS' if ok else 'FAILED_OR_CENSORED',supervision=s));checkpoint(out,ledger);print(d.name,name,s['outcome'],s['returncode'],s['wall_seconds'],flush=True)
  if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed certificate stage')
 for name,cmd in jobs:stage(name,cmd,ROOT)
 fresh=d/'independent-rank';fresh.mkdir(exist_ok=False)
 for path in [CAS/names[3],mod2]:shutil.copy2(path,fresh/path.name)
 checkpoint(fresh/'protocol.json',dict(files={x.name:cert.hashed(x) for x in [CAS/names[3],mod2]},seconds=180,rss_bytes=2147483648))
 stage('independent-rank',[SAGE,str(fresh/names[3]),'--input',str(fresh/mod2.name)],fresh)
 cloud=cert.read(mod2);odd=cert.read(modl);result=cert.read(source);initial=cert.read(d/'seed.json')['points']
 assert cloud['points'][:len(initial)]==initial
 assert all(r['search']['status']=='bounded_search_complete' for r in result['charts']) and len(result['charts'])==49
 ledger.update(status='PASS',initial_rank=len(initial),rank_lower_bound=cloud['rank_lower_bound'],discovered_rank_gain=cloud['rank_lower_bound']-len(initial),completed_boxes=49,point_count=len(cloud['points']),odd_modulus_ranks={str(a['modulus']):a['finite_column_rank'] for a in odd['audits']},certificates={str(x.relative_to(ROOT)):cert.hashed(x) for x in [mod2,modl]})
 checkpoint(out,ledger);print('CERTIFIED',d.name,'rank >=',cloud['rank_lower_bound'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--protocol',type=Path,required=True);p.add_argument('--prefix',required=True);a=p.parse_args();main(a.run.resolve(),a.protocol.resolve(),a.prefix)
