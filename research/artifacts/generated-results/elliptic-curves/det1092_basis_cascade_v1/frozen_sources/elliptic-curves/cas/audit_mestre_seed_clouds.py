#!/usr/bin/env python3
"""Check all supplied divisor images, beyond the selected own11 seed."""
import argparse
import json
from dataclasses import asdict
from pathlib import Path
import mestre_parent_calibration as batch
import audit_retained_cloud_modl as odd
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
OUT=batch.ART/'mestre_parent_calibration_seed_clouds_v1.json'

def expected():
    p=batch.protocol();rows=[];bindings={str(Path(__file__).resolve().relative_to(batch.ROOT)):cert.hashed(Path(__file__))}
    bindings.update(odd.sources())
    for row in p['rows']:
        source=batch.BATCH/row['id']/'seed.json';seed=cert.read(source)
        bindings[str(source.relative_to(batch.ROOT))]=cert.hashed(source)
        model=tuple(map(cert.F,seed['curve']));points=[tuple(map(cert.F,P)) for P in seed['divisor_cloud']]
        if any(not cert.is_on_weierstrass_curve(model,P) for P in points):raise ArithmeticError('seed image off curve')
        audits=[]
        for ell in (3,5):
            torsion=odd.ml.find_no_rational_l_torsion_prime(model,modulus=ell);basis={};signatures=[]
            for prime in odd.ml._primes_up_to(997):
                if prime in (2,ell):continue
                try:sig=odd.finite.signature(model,points,prime,ell)
                except ValueError:continue
                before=len(basis)
                for r in sig.rows:odd.insert(basis,r,ell)
                if len(basis)>before:signatures.append(asdict(sig))
            audits.append({'modulus':ell,'finite_rank':len(basis),'no_rational_ell_torsion_prime':torsion,
              'independent_column_indices':sorted(basis),'signatures':signatures})
        rows.append({'id':row['id'],'curve':seed['curve'],'points':seed['divisor_cloud'],'selected_seed_rank':11,'audits':audits})
        print(row['id'],'ALL14 INPUT IMAGES',[(a['modulus'],a['finite_rank']) for a in audits],flush=True)
    return {'schema':'elliptic-curves.mestre-parent-calibration-seed-clouds.v1','status':'PASS',
      'rows':rows,'sources':bindings,'prime_bound':997,
      'scope':'Exact finite3/5 lower bounds on all14 originally supplied divisor/covariant images. A finite rank deficit is not a dependence or upper-bound certificate. Search gains are measured from the certified selected11 seed, not asserted specialization jumps above the unknown exact generic rank.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=expected()
    if a.check:assert json.dumps(result,sort_keys=True)==json.dumps(cert.read(OUT),sort_keys=True)
    else:
        if OUT.exists():raise FileExistsError('preserve full input-image audit')
        checkpoint(OUT,result)
