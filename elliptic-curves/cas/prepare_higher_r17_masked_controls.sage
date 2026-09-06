#!/usr/bin/env sage-python
"""Split the oracle before computing centres from only the retained metric block."""
import sys
from decimal import Decimal
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import higher_r17_masked_controls as run
import certify_compact_r17_candidates as cert
from search_observability import masked_control
from research_runtime.store import checkpoint,digest
from memory_rank_certificate import checked_rank
sample=SourceFileLoader('masked_sample',str(CAS/'run_ordinary_masked_controls.sage')).load_module()
mapping=SourceFileLoader('masked_map',str(CAS/'prepare_product24_r17_pari_batch.sage')).load_module()
def main():
    p=run.protocol();out=run.D/'prepared.json'
    if out.exists():raise FileExistsError('preserve masked preparation')
    hashes={}
    for row in p['rows']:
        mp=ROOT/row['maps_path'];cp=ROOT/row['certificate_path']
        if cert.hashed(mp)!=row['maps_sha256'] or cert.hashed(cp)!=row['certificate_sha256']:raise ArithmeticError('source geometry/proof changed')
        old=cert.read(mp);proof=cert.read(cp);points=old['generic_points'];model=tuple(map(cert.F,old['curve']));c=proof['rank_certificate']
        generic_proof=checked_rank(model,[tuple(map(cert.F,P)) for P in points],[r['prime'] for r in c['signatures']],c['no_rational_2_torsion_prime'])
        gram=[[int((Decimal(v)*1000000).to_integral_value()) for v in r] for r in old['metric_gram']]
        blind,oracle=masked_control(old['curve'],points,gram,[0]);folder=run.D/row['id'];checkpoint(folder/'blind.json',blind)
        oracle.update(original_generic_points=points,original_independence=generic_proof);checkpoint(folder/'oracle.json',oracle)
        # Only the sixteen-point input enters geometry below.
        g=sample.geometry(blind['metric_gram']);kept=[(cert.F(P['x']),cert.F(P['y'])) for P in blind['points']]
        centres=[{'mask':r['mask'],'representative':r['coefficients'],'sample_minimum_norm':r['minimum_norm']} for r in g['selected']]
        maps=[mapping.mapping(model,kept,c) for c in centres]
        checkpoint(folder/'maps.json',{'status':'COMPLETE_FIXED_MASKED_MAPS','blind_sha256':cert.hashed(folder/'blind.json'),'geometry':g,'rows':maps});hashes[row['id']]=cert.hashed(folder/'maps.json');print('MASKED HIGHER PREPARED',row['id'],len(maps),flush=True)
    checkpoint(out,{'status':'PASS','protocol_hash':digest(p),'maps_hashes':hashes,'oracle_hashes':{r['id']:cert.hashed(run.D/r['id']/'oracle.json') for r in p['rows']}})
if __name__=='__main__':main()
