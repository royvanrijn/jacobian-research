#!/usr/bin/env sage-python
"""Exact prime certificates and two local conductor implementations, no factoring."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from sage.all import EllipticCurve, QQ, ZZ, pari
from sage.version import version

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert

WORK = ROOT/'artifacts/local/elliptic-curves/submitted627-630-conductors-v1'
OUT = ROOT/'artifacts/generated-results/elliptic-curves/submitted627_630_conductors_v1'


def kodaira(k):
    known = {1:'I0',2:'II',3:'III',4:'IV',-1:'I0*',-2:'II*',-3:'III*',-4:'IV*'}
    return known[k] if k in known else 'I'+str(abs(k)-4)+('*' if k<0 else '')


def prove(row, check):
    identifier = row['icarm_id']
    output = OUT/f'curve{identifier}.json'
    if check:
        saved = cert.read(output)
        factors = saved['discriminant_factorization']
    else:
        state_path = WORK/f'ecm_{identifier}/state.json'
        if not state_path.exists() or cert.read(state_path)['status'] != 'PROBABLE_PRIME_FACTORIZATION':
            print(identifier,'UNKNOWN: unresolved composite factors',flush=True)
            return
        counts = Counter(int(p) for p in cert.read(state_path)['factors'])
        for local in row['known_local_data']:
            counts[int(local['prime'])] += local['displayed_discriminant_valuation']
        factors = [[str(p), e] for p,e in sorted(counts.items()) if e]
    product = ZZ(1)
    certificates = {}
    for text,e in factors:
        p = ZZ(text)
        if p<2 or e<1:
            raise ArithmeticError('invalid prime/exponent')
        if check:
            primality = pari(saved['prime_certificates'][text])
        else:
            cache = WORK/f'prime_{p}.gp'
            primality = pari(cache.read_text()) if cache.exists() else pari(p).primecert()
            if not cache.exists():
                cache.write_text(str(primality)+'\n')
        root = ZZ(primality) if primality.type() == 't_INT' else ZZ(primality[0][0])
        if root != p or not pari.primecertisvalid(primality):
            raise ArithmeticError('invalid exact primality certificate')
        certificates[text] = str(primality)
        product *= p**e
    if product != abs(ZZ(row['discriminant'])):
        raise ArithmeticError('incomplete discriminant factorization')
    model = [QQ(a) for a in row['ainvs']]
    E = EllipticCurve(QQ,model)
    ep = pari.ellinit(model)
    if E.discriminant() != ZZ(row['discriminant']) or ZZ(ep.disc()) != E.discriminant():
        raise ArithmeticError('independent discriminants differ')
    conductor = ZZ(1)
    local_rows = []
    for text,e in factors:
        p = ZZ(text)
        local = E.local_data(p,algorithm='generic',proof=True)
        raw = ep.elllocalred(p)
        f = int(local.conductor_valuation())
        k = str(local.kodaira_symbol())
        if f != int(raw[0]) or k != kodaira(int(raw[1])):
            raise ArithmeticError('Sage/PARI local conductor disagreement')
        if QQ(raw[2][0]) != 1 or int(local.discriminant_valuation()) != e:
            raise ArithmeticError('submitted equation not globally minimal')
        if f<=0:
            raise ArithmeticError('discriminant prime is not bad on submitted minimal model')
        conductor *= p**f
        local_rows.append({'prime':text,'minimal_discriminant_valuation':e,
            'conductor_exponent':f,'kodaira':k,'pari_local_reduction':str(raw)})
    result = {'schema':'elliptic-curves.submitted-exact-conductor.v1','status':'PASS',
        'icarm_id':identifier,'inventory_id':row['inventory_id'],'ainvs':row['ainvs'],
        'discriminant':row['discriminant'],'discriminant_factorization':factors,
        'prime_certificates':certificates,'bad_primes':[p for p,e in factors],
        'conductor':str(conductor),'conductor_digits':len(str(conductor)),
        'global_minimality':'PROVED_AT_EVERY_DISCRIMINANT_PRIME',
        'local_data':local_rows,'sage_version':version,'pari_version':str(pari.version()),
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in
                   (Path(__file__).resolve(),WORK/'protocol.json',WORK/f'curve{identifier}.json')},
        'argument':'Exact ECPP certificates (or deterministic small-prime certificates) '
            'prove every factor prime and their product reconstructs the full discriminant. '
            'Independent Sage generic Tate and PARI elllocalred agree at every bad prime. '
            'The integral equation is minimal at all discriminant primes; all other primes '
            'are good. Thus the displayed list is complete and N=product(p^f_p) is exact.',
        'claim_boundary':'Exact bad-prime support and conductor; no rank change, universal '
            'record claim or external submission. Replay performs no factor search.'}
    if check:
        if saved != result:
            raise ArithmeticError('exact conductor replay differs')
    else:
        if output.exists() and cert.read(output) != result:
            raise FileExistsError('preserve exact conductor certificate')
        if not output.exists():
            cert.write(output,result)
        (OUT/f'curve{identifier}_bad_primes.txt').write_text(', '.join(result['bad_primes'])+'\n')
        (OUT/f'curve{identifier}_primes_payload.json').write_text(json.dumps({'primes':result['bad_primes']},indent=2)+'\n')
    print(identifier,'EXACT CONDUCTOR',result['conductor'],'BAD PRIMES',len(factors),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--id',type=int)
    args=parser.parse_args()
    for row in cert.read(WORK/'protocol.json')['roster']:
        if args.id is None or row['icarm_id']==args.id:
            prove(row,args.check)
