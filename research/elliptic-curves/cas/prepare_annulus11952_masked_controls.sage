#!/usr/bin/env sage-python
"""Split the oracle before computing centres from only the retained metric block."""
import sys
from decimal import Decimal
from pathlib import Path
from importlib.machinery import SourceFileLoader
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import annulus11952_masked_controls as run
import certify_compact_r17_candidates as cert
from search_observability import masked_control
from research_runtime.store import checkpoint,digest
from memory_rank_certificate import checked_rank
from research_runtime.search_state import raw_state
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
sample=SourceFileLoader('masked_sample',str(CAS/'run_ordinary_masked_controls.sage')).load_module()
mapping=SourceFileLoader('masked_map',str(CAS/'prepare_product24_r17_pari_batch.sage')).load_module()
def main():
    p=run.protocol();out=run.D/'prepared.json'
    if out.exists():raise FileExistsError('preserve masked preparation')
    hashes={}
    for row in p['rows']:
        mp=ROOT/row['maps_path']
        if cert.hashed(mp)!=row['maps_sha256']:raise ArithmeticError('source geometry changed')
        old=cert.read(mp);points=old['generic_points'];model=tuple(map(cert.F,old['curve']));generic=tuple(tuple(map(cert.F,P)) for P in points)
        state=raw_state(model,generic,cache=ReductionCache(MemoryFactStore()),prime_bound=997)
        if state.rank!=17 or state.basis!=generic:raise ArithmeticError('all17 generic points must independently certify before masking')
        generic_proof=checked_rank(model,generic,state.reductions.primes,state.no_two_torsion_prime)
        gram=[[int((Decimal(v)*1000000).to_integral_value()) for v in r] for r in old['metric_gram']]
        blind,oracle=masked_control(old['curve'],points,gram,[0]);folder=run.D/row['id'];checkpoint(folder/'blind.json',blind)
        oracle.update(original_generic_points=points,original_independence=generic_proof);checkpoint(folder/'oracle.json',oracle)
        # Only the sixteen-point input enters geometry below.
        g=sample.geometry(blind['metric_gram']);kept=[(cert.F(P['x']),cert.F(P['y'])) for P in blind['points']]
        centres=[{'mask':r['mask'],'representative':r['coefficients'],'sample_minimum_norm':r['minimum_norm']} for r in g['selected']]
        maps=[mapping.mapping(model,kept,c) for c in centres]
        checkpoint(folder/'maps.json',{'status':'COMPLETE_FIXED_MASKED_MAPS','blind_sha256':cert.hashed(folder/'blind.json'),'geometry':g,'rows':maps});hashes[row['id']]=cert.hashed(folder/'maps.json');print('MASKED NEW11952 PREPARED',row['id'],len(maps),flush=True)
    checkpoint(out,{'status':'PASS','protocol_hash':digest(p),'maps_hashes':hashes,'oracle_hashes':{r['id']:cert.hashed(run.D/r['id']/'oracle.json') for r in p['rows']}})
if __name__=='__main__':main()
