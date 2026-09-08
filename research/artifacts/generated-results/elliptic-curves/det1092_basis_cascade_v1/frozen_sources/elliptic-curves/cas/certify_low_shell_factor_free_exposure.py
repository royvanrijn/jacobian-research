#!/usr/bin/env python3
"""Independent rank and low-shell geometry certificates for one frozen wave."""
import argparse, shutil, sys
from pathlib import Path

import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'


def main(directory,protocol,prefix):
    out=directory/'certification-ledger.json'
    if out.exists():raise FileExistsError('preserve low-shell certificate ledger')
    source=directory/'result.json';mod2=ART/(prefix+'_mod2_v1.json');modl=ART/(prefix+'_modl_v1.json')
    names=['audit_recorded_point_mod2_rank_v3.py','audit_retained_cloud_modl.py','verify_low_shell_factor_free_exposure.sage','verify_factor_free_rank.sage','certify_low_shell_factor_free_exposure.py']
    checkpoint(directory/'certification-protocol.json',{'inputs':{str(path.relative_to(ROOT)):cert.hashed(path) for path in (source,protocol,directory/'maps.json',directory/'seed.json')},'sources':{str((CAS/name).relative_to(ROOT)):cert.hashed(CAS/name) for name in names},'seconds_per_stage':180,'rss_bytes':2147483648})
    jobs=[('mod2-build',[sys.executable,str(CAS/names[0]),'--input',str(source),'--input-sha256',cert.hashed(source),'--output',str(mod2),'--prime-bound','997']),('mod2-check',[sys.executable,str(CAS/names[0]),'--check',str(mod2)]),('modl-build',[sys.executable,str(CAS/names[1]),'--input',str(mod2),'--output',str(modl)]),('modl-check',[sys.executable,str(CAS/names[1]),'--check',str(modl)]),('geometry',[SAGE,str(CAS/names[2]),'--run',str(directory),'--protocol',str(protocol),'--cloud',str(mod2)])]
    ledger={'status':'RUNNING','stages':[]};checkpoint(out,ledger)
    def stage(name,command,cwd):
        supervised=run(command,limits=Limits(180,2147483648),log_path=directory/(name+'-certificate.log'),checkpoint_path=directory/(name+'-certificate.supervisor.json'),cwd=cwd);ok=supervised['outcome']=='completed' and supervised['returncode']==0
        ledger['stages'].append({'name':name,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':supervised});checkpoint(out,ledger)
        if not ok:ledger['status']='FAILED_OR_CENSORED';checkpoint(out,ledger);raise ArithmeticError('preserve failed low-shell certificate stage')
    for name,command in jobs:stage(name,command,ROOT)
    fresh=directory/'independent-rank';fresh.mkdir()
    for path in (CAS/names[3],mod2):shutil.copy2(path,fresh/path.name)
    checkpoint(fresh/'protocol.json',{'files':{path.name:cert.hashed(path) for path in (CAS/names[3],mod2)},'seconds':180,'rss_bytes':2147483648})
    stage('independent-rank',[SAGE,str(fresh/names[3]),'--input',str(fresh/mod2.name)],fresh)
    cloud,odd,result=cert.read(mod2),cert.read(modl),cert.read(source);initial=cert.read(directory/'seed.json')['points']
    if cloud['points'][:len(initial)]!=initial or any(row['search']['status']!='bounded_search_complete' for row in result['charts']) or len(result['charts'])!=49:raise ArithmeticError('low-shell completed exposure differs')
    ledger.update(status='PASS',initial_rank=len(initial),rank_lower_bound=cloud['rank_lower_bound'],discovered_rank_gain=cloud['rank_lower_bound']-len(initial),completed_boxes=49,point_count=len(cloud['points']),odd_modulus_ranks={str(row['modulus']):row['finite_column_rank'] for row in odd['audits']},certificates={str(path.relative_to(ROOT)):cert.hashed(path) for path in (mod2,modl)})
    checkpoint(out,ledger);print('CERTIFIED LOW-SHELL',directory.name,'rank >=',cloud['rank_lower_bound'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',type=Path,required=True);parser.add_argument('--protocol',type=Path,required=True);parser.add_argument('--prefix',required=True);args=parser.parse_args();main(args.run.resolve(),args.protocol.resolve(),args.prefix)
