#!/usr/bin/env python3
"""Six frozen masked generic-direction controls; worker inputs omit the oracle."""
import argparse,sys
from pathlib import Path
import certify_compact_r17_candidates as cert
import pari_pointed_backend as backend
from pointed_quartic_search import PointedQuarticSearch
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import run,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'higher-r17-masked-controls-v1';PARENT=LOCAL/'higher24-r17-pari-v1';SAGE='/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python'
def sources():
    names=['higher_r17_masked_controls.py','prepare_higher_r17_masked_controls.sage','run_ordinary_masked_controls.sage','prepare_product24_r17_pari_batch.sage','prospective_half_lattice_v2.sage','search_observability.py','research_runtime/deep_centres.py']
    return {**backend.sources(),**{str((CAS/n).relative_to(ROOT)):cert.hashed(CAS/n) for n in names}}
def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve masked protocol')
    parent=cert.read(PARENT/'protocol.json'); rows=[]
    if cert.read(ART/'product22_comparison_portable_replay_v1.json')['status']!='PASS':raise ArithmeticError('paired audit and visibility countercheck must complete first')
    for family in sorted({r['family'] for r in parent['rows']}):
        r=min((r for r in parent['rows'] if r['family']==family),key=lambda r:r['retained_index']);mp=PARENT/r['id']/'maps.json';cp=ART/('higher24_r17_'+r['id'].replace('-','_')+'_mod2_v1.json')
        rows.append({'id':r['id'],'family':family,'parameter':r['parameter'],'maps_path':str(mp.relative_to(ROOT)),'maps_sha256':cert.hashed(mp),'certificate_path':str(cp.relative_to(ROOT)),'certificate_sha256':cert.hashed(cp)})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher-r17-masked-controls.v1','sources':sources(),'rows':rows,'parent_protocol_sha256':cert.hashed(PARENT/'protocol.json'),'withhold_zero_based':[0],'generic_dimension':17,'retained_dimension':16,'sample_masks':256,'charts_per_curve':12,'height':125000,'seconds_per_chart':10,'gp_sha256':cert.hashed(Path('/usr/bin/gp')),'workers':1,'map_wall_seconds':600,'worker_wall_seconds':300,'replay_wall_seconds':120,'rss_bytes':1610612736,'selection':'One already exposed higher-cohort address per family: smallest retained index, independent of its point result. No replacement or outcome filter. Withhold original generic section zero uniformly. The sixteen-point principal specialized metric block alone selects twelve deepest of256 fixed hashed parity masks using the existing exact sampled-coset geometry. Neither withheld coordinates nor exceptional points enter centre selection or the worker.','gate':'The completed paired score experiment gives weak detection, while exact direct generic-representative visibility is not a valid failure diagnostic. The earlier ordinary-fibre masked controls motivate a separately bounded control on the actual larger-coefficient R17 fibres. This diagnoses recovery of known withheld directions under changed subgroup geometry; it does not measure the incidence of exceptional points.','endpoint':'WITHHELD_KNOWN_DIRECTIONS_NOT_NEW_RANK','boundaries':'Six curves and72 boxes only; no new parameter scan or public target. Freeze precedes masking, geometry and search. Open oracles only after all point attempts and geometry replays finish. Require exact rational group relations with nonzero withheld coefficient for a recovered-known-direction claim; unresolved relations stay UNKNOWN. Recovery and failure do not establish sensitivity to all exceptional directions, exact rank, saturation or upper bounds.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen masked sources changed')
    return p
def cell(index,check=False):
    p=protocol();row=p['rows'][index];folder=D/row['id'];mp=folder/'maps.json';maps=cert.read(mp);blind=cert.read(folder/'blind.json');binding=cert.read(D/'prepared.json')
    if binding['protocol_hash']!=digest(p) or binding['maps_hashes'][row['id']]!=cert.hashed(mp) or maps['blind_sha256']!=cert.hashed(folder/'blind.json'):raise ArithmeticError('blind inputs changed')
    if len(blind['points'])!=16 or len(maps['rows'])!=12:raise ArithmeticError('fixed masked geometry differs')
    out=folder/'result.json'
    if not check and out.exists():raise FileExistsError('preserve masked point attempt')
    data=cert.read(out) if check else {'status':'RUNNING','protocol_hash':digest(p),'maps_sha256':cert.hashed(mp),'family':row['family'],'parameter':row['parameter'],'curve':blind['curve'],'retained_points':blind['points'],'charts':[]}
    if data['protocol_hash']!=digest(p) or data['maps_sha256']!=cert.hashed(mp) or data['curve']!=blind['curve'] or data['retained_points']!=blind['points']:raise ArithmeticError('masked search binding differs')
    for i,m in enumerate(maps['rows']):
        search=PointedQuarticSearch(curve=blind['curve'],subgroup=blind['points'],centre={'coefficients':m['centre']['representative']},coordinate_policy=m['coordinate_policy'])
        if check:
            r=data['charts'][i]
            if r['index']!=i or r['search']['height_bound']!=p['height'] or r['search']['timeout_seconds']!=p['seconds_per_chart'] or r['search']['gp_binary_sha256']!=p['gp_sha256']:raise ArithmeticError('masked budget differs')
            backend.replay(search,m,r['search'])
        else:
            r,points=backend.execute(search,m,p['height'],p['seconds_per_chart'],p['gp_sha256']);data['charts'].append({'index':i,'search':r});checkpoint(out,data);print('MASKED HIGHER',row['id'],i+1,r['status'],flush=True)
    if check:
        if len(data['charts'])!=12 or data['status']!='COMPLETE_FIXED_MASKED_ATTEMPT':raise ArithmeticError('all72 fixed attempts required')
    else:data['status']='COMPLETE_FIXED_MASKED_ATTEMPT';checkpoint(out,data)
def launch():
    p=protocol();out=D/'ledger.json'
    if out.exists():raise FileExistsError('preserve masked ledger')
    ledger={'status':'RUNNING','rows':[]};checkpoint(out,ledger)
    for i,row in enumerate(p['rows']):
        stages=[]
        for name,budget in [('worker',p['worker_wall_seconds']),('replay',p['replay_wall_seconds'])]:
            folder=D/row['id'];s=run([sys.executable,str(Path(__file__).resolve()),name,'--index',str(i)],limits=Limits(budget,p['rss_bytes']),log_path=folder/(name+'.log'),checkpoint_path=folder/(name+'.supervisor.json'),cwd=ROOT);stages.append({'name':name,'supervision':s})
            if s['outcome']!='completed' or s['returncode']!=0:raise ArithmeticError('masked stage failed/censored; retain checkpoint')
        ledger['rows'].append({'id':row['id'],'status':'PASS','stages':stages,'result_sha256':cert.hashed(D/row['id']/'result.json')});checkpoint(out,ledger)
    ledger['status']='PASS';checkpoint(out,ledger)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','launch','worker','replay']);p.add_argument('--index',type=int);a=p.parse_args();cell(a.index,a.stage=='replay') if a.stage in ('worker','replay') else globals()[a.stage]()
