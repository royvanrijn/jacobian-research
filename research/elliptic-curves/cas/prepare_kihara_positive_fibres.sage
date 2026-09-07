#!/usr/bin/env sage-python
"""Fixed six-fibre seed intake before positive-parent point exposure."""
import json,hashlib,sys
from pathlib import Path
from dataclasses import asdict
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,lcm,prime_range
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from audit_recorded_point_mod2_rank_v3 import signature,insert
from memory_rank_certificate import checked_rank
D=ROOT/'artifacts/local/elliptic-curves/kihara-positive-fibre-intake-v1'
def write(p,r):
    with p.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
def main():
    p=json.loads((D/'protocol.json').read_text())
    for n,h in p['sources'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h
    parents=json.loads((ROOT/p['parent_input']).read_text())['rows'];R=PolynomialRing(QQ,'T');K=R.fraction_field();rows=[]
    for index,row in enumerate(p['rows']):
        parent=parents[row['parent_index']];t=QQ(row['fibre_T']);a=R(parent['raw_A'])(t);b=R(parent['raw_B'])(t);den=lcm(a.denominator(),b.denominator());scale=QQ(den)
        ai=ZZ(a*scale**4);bi=ZZ(b*scale**6)
        for ell in prime_range(2,998):
            power=min(ai.valuation(ell)//4,bi.valuation(ell)//6)
            if power:ai//=ell**(4*power);bi//=ell**(6*power);scale/=ell**power
        E=EllipticCurve(QQ,[ai,bi]);assert E.discriminant()!=0
        all_points=[E([scale**2*K(P[0])(t),scale**3*K(P[1])(t),K(P[2])(t)]) for P in parent['sections']]
        selected=[all_points[i] for i in parent['independent_indices']];assert all(P for P in selected)
        model=tuple(cert.F(str(c)) for c in E.a_invariants());cloud=[tuple(cert.F(str(c)) for c in P.xy()) for P in selected]
        cache=ReductionCache(MemoryFactStore());pivots={};signatures=[];torsion=None
        for prime in prime_range(3,998):
            prime=int(prime)
            if torsion is None and cert.short_curve_has_no_rational_2_torsion_modular_certificate(model,prime):torsion=prime
            try:sig=signature(cache,model,cloud,prime)
            except ValueError:continue
            before=len(pivots)
            for bits in sig.rows:insert(pivots,bits)
            if len(pivots)>before:signatures.append(asdict(sig))
        expected=parent['generic_section_span_rank'];ok=len(pivots)==expected and torsion is not None
        proof=checked_rank(model,cloud,[s['prime'] for s in signatures],torsion) if ok else None
        result={'id':row['id'],'family':'kihara-positive-'+parent['ratio'],'parameter':str(t),'parent_ratio':parent['ratio'],'fibre_T':str(t),'parent_index':row['parent_index'],'status':'PASS' if ok else 'UNRESOLVED_SEED_GATE','curve':list(map(str,model)),'points':[[str(c) for c in P] for P in cloud],'generic_points':[[str(c) for c in P] for P in cloud],'generic_indices':parent['independent_indices'],'scale_from_parent':str(scale),'expected_seed_rank':expected,'finite_mod2_rank':len(pivots),'no_rational_2_torsion_prime':torsion,'rank_certificate':proof,'rank_lower_bound':expected if ok else None,'model_coefficient_bits':max(abs(c.numerator).bit_length() for c in model),'scope':'Fixed fibre specialization of independently certified generic section basis. Model scaling removes only fourth/sixth powers at primes<=997; no global minimality claim. Full specialized subgroup requires the independent finite quotient gate; no rank gain inferred from generic ranks.'}
        write(D/('seed'+str(index)+'.json'),result);rows.append(result)
        print(row['id'],result['status'],'mod2',len(pivots),'expected',expected,'torsion witness',torsion,'bits',result['model_coefficient_bits'],flush=True)
    write(D/'result.json',{'status':'PASS' if all(r['status']=='PASS' for r in rows) else 'UNRESOLVED_SEED_GATE','rows':rows,'sources':p['sources']})
if __name__=='__main__':main()
