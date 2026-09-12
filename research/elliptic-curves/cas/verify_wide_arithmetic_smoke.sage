#!/usr/bin/env sage-python
"""Recheck completed smoke rows, independently of the profiling worker."""
import argparse,hashlib,json
from pathlib import Path
from sage.all import EllipticCurve,PolynomialRing,QQ,ZZ,pari

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path('research/artifacts/local/elliptic-curves/wide-arithmetic-profile-v1'));args=ap.parse_args();out=args.output.resolve()
    pop=read(out/'population.json')['curves'];records=[];bindings={}
    for row in pop[:3]:
        key=row['curve_key'];inp=out/'inputs'/f'{key}.json'
        base_path=out/'base'/f'{key}.json';local_path=out/'local'/f'{key}.json'
        b=read(base_path);l=read(local_path)
        for p in [inp,base_path,local_path]:bindings[str(p)]=sha(p)
        assert b['status']=='PASS' and b['input_sha256']==sha(inp)==l['input_sha256']
        original=EllipticCurve(QQ,row['ainvs']);E=EllipticCurve(QQ,list(map(QQ,b['minimal_ainvs'])))
        iso=original.isomorphism_to(E);u=iso.tuple()[0]
        assert original.discriminant()==u**12*E.discriminant()
        assert E.c4()**3-E.c6()**2==1728*E.discriminant()
        assert str(E.discriminant())==b['minimal_discriminant']
        R=PolynomialRing(QQ,'x');x=R.gen();f=R(list(map(QQ,b['two_division_cubic'])))
        expected=4*x**3+E.b2()*x**2+2*E.b4()*x+E.b6()
        assert f.monic()==expected.monic() and str(f.discriminant())==b['two_division_cubic_discriminant']
        assert f.is_irreducible() and b['rational_2torsion_rank']==0
        result={'curve_key':key,'base':'PASS','local':l['status']}
        if l['status']=='PASS':
            factors=[(ZZ(p),e) for p,e in l['minimal_discriminant_factorization']]
            product=ZZ(1)
            for p,e in factors:
                assert p.is_prime(proof=True);product*=p**e
            assert product==abs(E.discriminant())
            primes=[p for p,e in factors];pari.addprimes(primes)
            # Integral monic 2-division polynomial, unlike the worker's
            # primitive x-coordinate polynomial. X=4x gives the same field.
            F=x**3+E.b2()*x*x+8*E.b4()*x+16*E.b6()
            nf=pari.nfinit([pari(F),[2]+primes]);assert pari.nfcertify(nf)==[]
            assert str(nf.disc())==l['field_discriminant']
            assert list(map(int,nf.nf_get_sign()))==l['field_signature']
            offset=1 if E.discriminant()<0 else 2;N=ZZ(1)
            for r in l['local_data']:
                p=ZZ(r['p']);ld=pari.elllocalred(pari(E),p)
                assert E.discriminant().valuation(p)==r['v_delta'] and int(ld[0])==r['conductor_exponent']
                n=len(pari.idealprimedec(nf,p))
                offset += int(r['v_delta']%2==0) if int(ld[0])==1 else n-1
                if int(ld[0])!=1:assert n==r['cubic_prime_count']
                N*=p**int(ld[0])
            assert offset==l['bk_local_term'] and str(N)==l['conductor']
            assert int(pari.ellrootno(pari(E)))==l['root_number']
            result['bk_local_term']=offset
        else:
            assert l['status']=='UNKNOWN_TIMEOUT',l
        records.append(result)
    assert any(r['local']=='PASS' for r in records),'No complete local smoke witness'
    assert not list((out/'class').glob('*.json'))
    result={'status':'PASS_SMOKE_REPLAY','records':records,'bindings':bindings,'checker_sha256':sha(Path(__file__)),
            'boundary':'BASE transport and cubic identities; completed LOCAL rows replayed through monic PARI field and local reductions. Shared CAS, not an independent full descent. Timeouts remain UNKNOWN; no CLASS probes.'}
    p=out/'SMOKE_VALIDATION.json';raw=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if p.exists():assert p.read_text()==raw
    else:p.write_text(raw)
    print(json.dumps(result['records']),flush=True)

if __name__=='__main__':main()
