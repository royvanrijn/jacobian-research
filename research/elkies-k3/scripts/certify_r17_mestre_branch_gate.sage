#!/usr/bin/env sage-python
"""Small finite-polynomial witnesses for the unrestricted auxiliary-function gate."""
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import PolynomialRing, QQ, GF, prime_range

ROOT=Path(__file__).resolve().parents[2]
PATH=ROOT/'artifacts/generated-results/elkies-k3-r17-mestre-shared-twist-v1'


def main():
    resource.setrlimit(resource.RLIMIT_CPU,(20,25));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time();packet=json.loads((PATH/'input.json').read_text());result=json.loads((PATH/'result.json').read_text())
    R=PolynomialRing(QQ,'t');A,B=R(packet['A']),R(packet['B']);D=R(result['D']);Delta=4*A**3+27*B*B
    primes=list(prime_range(5,1000));rows=[]
    for name,f in [('A',A),('B',B)]:
        irreducible=None;constant_field=None
        for p in primes:
            if f.denominator()%p==0:continue
            S=PolynomialRing(GF(p),'t');g=S(f)
            if g.degree()!=f.degree() or g.gcd(g.derivative()).degree()!=0:continue
            if irreducible is None and g.is_irreducible():irreducible=int(p)
            congruence=(p%12==11) if name=='A' else (p%4==3)
            if constant_field is None and congruence:
                roots=g.roots(multiplicities=False)
                if roots:constant_field={'prime':int(p),'root':int(roots[0])}
            if irreducible and constant_field:break
        rows.append({'polynomial':name,'degree':int(f.degree()),'irreducible_reduction_prime':irreducible,
                     'degree_one_prime':constant_field,'status':'CERTIFICATE_READY' if irreducible and constant_field else 'UNKNOWN'})
    separability=None
    for p in primes:
        S=PolynomialRing(GF(p),'t');d,e=S(D),S(Delta)
        if d.degree()!=44 or e.degree()!=24:continue
        if d.gcd(d.derivative()).degree()==e.gcd(e.derivative()).degree()==d.gcd(e).degree()==0:
            separability=int(p);break
    out={'schema':'r17-mestre-branch-gate-v1','input_sha256':sha256((PATH/'input.json').read_bytes()).hexdigest(),
         'prime_bound':1000,'root_field_certificates':rows,'smooth_disjoint_branch_prime':separability,
         'cpu_seconds':time.process_time()-started,
         'proposed_general_conclusion':'For every rational u(t), the Mestre twist has odd valuation at the degree8 A divisor and degree12 B divisor, hence genus at least9. This requires the written valuation proof and independent finite-polynomial replay.'}
    with (PATH/'branch-gate.json').open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(out,sort_keys=True),flush=True)


if __name__=='__main__':main()
