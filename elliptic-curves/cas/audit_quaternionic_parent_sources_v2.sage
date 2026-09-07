#!/usr/bin/env sage-python
"""Finite Lin--Yang polarization/source comparison at D=546,314,206."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from sage.all import QQ,ZZ,QuaternionAlgebra,divisors,kronecker,prod
from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives

def build():
    records=[]
    for D in [546,314,206]:
        twists=[dict(m=int(m),algebra_discriminant=int(QuaternionAlgebra(QQ,-D,m).discriminant()))
                for m in divisors(D)]
        admitted=[r['m'] for r in twists if r['algebra_discriminant']==D]
        forms=BinaryQF_reduced_representatives(-16*ZZ(D),primitive_only=False)
        polarizations=set()
        for f in forms:
            a,b,c=map(int,f)
            if any((a*x*x+b*x*y+c*y*y)%4 not in [0,1] for x in range(4) for y in range(4)):continue
            if QuaternionAlgebra(QQ,-D,a).discriminant()!=D:continue
            polarizations.add((a,abs(b),c))
        h=len(BinaryQF_reduced_representatives(-4*ZZ(D),primitive_only=True))
        k=len(ZZ(D).prime_divisors());predicted=QQ(h)/2**k+QQ(len(admitted))/4
        assert len(polarizations)==predicted
        primes=ZZ(D).prime_divisors()
        e2=prod(1-kronecker(-4,p) for p in primes)
        e3=prod(1-kronecker(-3,p) for p in primes)
        genus=1+QQ(prod(p-1 for p in primes))/12-QQ(e2)/4-QQ(e3)/3
        # Full Fricke fixed points have discriminant -4D for these even D.
        quotient=1+(2*genus-2-h)/4
        assert quotient in ZZ
        records.append(dict(D=D,twisting_tests=twists,twisting_divisors=admitted,
                            class_number_minus_4D=h,polarization_forms=sorted(polarizations),
                            component_count=int(predicted),norm_one_genus=int(genus),
                            full_Fricke_genus=int(quotient)))
    return dict(schema='quaternionic-parent-source-audit.v1',status='PASS',records=records,
                external_inputs='Lin-Yang Theorem1/3, Proposition13; quaternionic genus and Fricke fixed-point formulas',
                boundary='Polarization moduli only. No rootless-frame enumeration, K3 coefficient search, or point search.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=build();r['checker_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(r))
