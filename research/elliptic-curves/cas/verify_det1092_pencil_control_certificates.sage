#!/usr/bin/env sage-python
"""Independent Frobenius-gcd replay of the13 root-free control certificates."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,gcd,lcm
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_pencil_multiples_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
data=json.loads((OUT/'controls.json').read_text())
for path,expected in data['inputs'].items():assert sha(ROOT/path)==expected
R=PolynomialRing(QQ,'z');passed=0;unknown=0
for row in data['cases']:
    if row['status']=='UNKNOWN':unknown+=1;continue
    assert row['status']=='NO_RATIONAL_PREIMAGE'
    cover=json.loads((OUT/('multiple%d.json'%row['n'])).read_text())
    fmap=cover['t_of_z'];deg=cover['parameter_map_degree']
    f=R(fmap['numerator'])-QQ(row['parameter'])*R(fmap['denominator'])
    f*=lcm(c.denominator() for c in f.list())
    f=R(f/gcd([ZZ(c) for c in f.list()]))
    p=row['certificate']['prime'];Rp=PolynomialRing(GF(p),'z');z=Rp.gen();fp=Rp(f)
    assert fp.degree()==deg
    assert list(map(int,fp.list()))==row['certificate']['primitive_polynomial_mod_p']
    # Roots in Fp are exactly gcd(f,z^p-z); infinity is checked separately.
    assert fp.gcd(pow(z,p,fp)-z)==1
    assert int(fp[deg])==row['certificate']['infinity_value']!=0
    assert [int(fp(i)) for i in range(p)]==row['certificate']['affine_values']
    passed+=1
assert (passed,unknown)==(13,5)
result={'classification':'independent verified modular application',
        'status':'PASS_CONTROL_ROOT_FREE_CERTIFICATES','exclusions':passed,'unknown':unknown,
        'checker_sha256':sha(Path(__file__)),'controls_sha256':sha(OUT/'controls.json'),
        'boundary':'UNKNOWN cases remain UNKNOWN; no additional primes or rational root searches.'}
path=OUT/'control-replay.json'
if path.exists():assert json.loads(path.read_text())==result
else:
    with path.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
print(json.dumps(result,sort_keys=True))
