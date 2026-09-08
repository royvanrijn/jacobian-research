#!/usr/bin/env sage-python
"""Exact torsion and odd-prime diagnostics on the three failed seed gates."""
import json,hashlib,sys
from pathlib import Path
from dataclasses import asdict
from sage.all import QQ,PolynomialRing,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import audit_compact_r17_ambiguous as finite
from audit_retained_cloud_modl import insert
D=ROOT/'artifacts/local/elliptic-curves/kihara-positive-seed-gaps-v1'
def write(p,r):
    with p.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
def main():
    protocol=cert.read(D/'protocol.json')
    for p,h in protocol['sources'].items():assert cert.hashed(ROOT/p)==h
    out=[]
    for i in protocol['seed_indices']:
        seed=cert.read(ROOT/'artifacts/local/elliptic-curves/kihara-positive-fibre-intake-v1'/('seed'+str(i)+'.json'))
        model=tuple(map(cert.F,seed['curve']));points=[tuple(map(cert.F,p)) for p in seed['points']];R=PolynomialRing(QQ,'x');x=R.gen();torsion=(x**3+QQ(model[3])*x+QQ(model[4])).roots();E=EllipticCurve(QQ,list(map(QQ,model)))
        row={'id':seed['id'],'two_torsion_points':[[str(r),'0'] for r,n in torsion],'audits':[]}
        assert all(2*E([r,0])==E(0) for r,n in torsion)
        for ell in (3,5):
            tp=finite.ml.find_no_rational_l_torsion_prime(model,modulus=ell);basis={};sigs=[]
            for p in finite.ml._primes_up_to(997):
                if p in (2,ell):continue
                try:sig=finite.signature(model,points,p,ell)
                except ValueError:continue
                before=len(basis)
                for r in sig.rows:insert(basis,r,ell)
                if len(basis)>before:sigs.append(asdict(sig))
            row['audits'].append({'modulus':ell,'rank_lower_bound':len(basis),'independent_indices':sorted(basis),'no_rational_ell_torsion_prime':tp,'signatures':sigs})
            print(seed['id'],'mod',ell,'rank',len(basis),'2-torsion points',len(torsion),flush=True)
        write(D/('seed'+str(i)+'.json'),row);out.append(row)
    write(D/'result.json',{'status':'PASS','rows':out,'sources':protocol['sources'],'scope':'Fixed three failed mod2 intake gates. Exact rational2-torsion and complete bounded mod3/5 column ranks with torsion exclusion. Deficient finite ranks do not prove rank loss or nonsaturation; no point search or replacement.'})
if __name__=='__main__':main()
