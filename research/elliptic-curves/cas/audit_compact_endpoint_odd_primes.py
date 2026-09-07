#!/usr/bin/env python3
"""Check all transported endpoint sections modulo3 and5 before diagnosing specialization loss."""
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/compact-endpoint-odd-primes-v1';INPUT=ART/'compact_atlas_endpoints_v2.json'
def main():
    if (D/'protocol.json').exists():raise FileExistsError('preserve endpoint odd-prime protocol')
    data=cert.read(INPUT);rows=[r for r in data['rows'] if r['status']=='CERTIFIED_SPECIALIZED_SUBGROUP'];gate=ROOT/'artifacts/local/elliptic-curves/compact-atlas-endpoints-v2/check.supervisor.json';g=cert.read(gate)
    if data['status']!='COMPLETE_DECLARED_ENDPOINT_AUDIT' or len(rows)!=21 or g['outcome']!='completed' or g['returncode']!=0:raise ArithmeticError('complete endpoint replay required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.compact-endpoint-odd-primes.v1','sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),INPUT,gate,CAS/'audit_endpoint_section_cloud_modl.py')},'rows':[(r['family'],r['endpoint']) for r in rows],'prime_bound':997,'moduli':[3,5],'seconds_per_stage':180,'rss_bytes':1073741824,'maximum_workers':1,'gate':'The exact endpoint transport gives finite mod2 lower bounds11..17 from generic sections. A lower mod2 image may reflect nonsaturation rather than rational dependence. Check every finite transported generic point modulo3 and5 before interpreting a specialization loss.','boundaries':'No new rational point search, candidate score, whole-curve rank upper bound, saturation or exact generic specialization rank.'})
    ledger={'status':'RUNNING','rows':[]};checkpoint(D/'ledger.json',ledger)
    for r in rows:
        name=r['family']+'-'+r['endpoint'];folder=D/name;ip=folder/'generic-points.json';op=ART/('endpoint_sections_'+name.replace('-','_')+'_modl_v1.json')
        checkpoint(ip,{**r,'parameter':r['endpoint'],'points':r['generic_points'],'endpoint_proof_sha256':cert.hashed(INPUT),'input_kind':'EXACT_GENERIC_SECTION_VALUES_NOT_A_SEARCH_TRANSCRIPT','scope':'Includes all finite generic section values. The prior mod2 certificate supplies only its lower bound, without rational dependence or saturation claims.'})
        for label,args in [('build',['--input',str(ip),'--output',str(op)]),('check',['--check',str(op)])]:
            s=run(['/usr/bin/python3',str(CAS/'audit_endpoint_section_cloud_modl.py'),*args],limits=Limits(180,1073741824),log_path=folder/(label+'.log'),checkpoint_path=folder/(label+'.supervisor.json'),cwd=ROOT)
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('endpoint odd-prime audit failed/censored')
        proof=cert.read(op);ledger['rows'].append({'family':r['family'],'endpoint':r['endpoint'],'input':str(ip.relative_to(ROOT)),'input_sha256':cert.hashed(ip),'output':str(op.relative_to(ROOT)),'output_sha256':cert.hashed(op),'status':'PASS','mod2_lower_bound':r['rank_lower_bound'],'odd_modulus_lower_bounds':{str(a['modulus']):a['finite_column_rank'] for a in proof['audits']}});checkpoint(D/'ledger.json',ledger);print('ENDPOINT ODD PRIMES',name,ledger['rows'][-1]['mod2_lower_bound'],ledger['rows'][-1]['odd_modulus_lower_bounds'],flush=True)
    ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
if __name__=='__main__':main()
