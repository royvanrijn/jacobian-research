#!/usr/bin/env python3
"""One exact halving layer of seven prescribed subset sums on a fixed small curve."""
import argparse
from pathlib import Path
import retrospective as r

PROTOCOL=Path(__file__).with_name('SHARED_VALUE_HALVING_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_shared_value_halving_v1.json'


def compute():
    from sage.all import QQ,GF,AA,EllipticCurve,PolynomialRing
    spec=r.read(PROTOCOL);E=EllipticCurve(QQ,spec['curve_ainvariants'])
    points=[E(P) for P in spec['points']];halves=[];records=[]
    for mask in range(1,8):
        value=E(0)
        for i,P in enumerate(points):
            if mask>>i&1:value+=P
        found=value.division_points(2)
        assert all(2*H==value for H in found)
        records.append({'mask':mask,'sum':list(map(str,value)),
                        'halves':[list(map(str,H)) for H in found]})
        halves.extend(found)
    allpoints=points+halves;signatures=[0]*len(allpoints)
    R=PolynomialRing(QQ,'z');z=R.gen();f=z**3+7*z*z-10*z+9
    assert R.change_ring(GF(2))(f).is_irreducible()
    roots=f.roots(AA,multiplicities=False);offset=len(roots)
    for i,P in enumerate(allpoints):signatures[i]=r.pack(int(P[0]<a) for a in roots)
    fingerprints=[]
    for p in r.primes(spec['limits']['largest_fingerprint_prime']):
        if f.discriminant()%p==0:continue
        fp=PolynomialRing(GF(p),'z')(f);roots=sorted(int(t) for t in fp.roots(multiplicities=False))
        if len(roots)!=3:continue
        local=[]
        for i,P in enumerate(allpoints):
            bits=[]
            for root in roots:
                if QQ(P[0]).denominator()%p==0:bits.append(0);continue
                value=(r.mod(r.F(str(P[0])),p)-root)%p
                if not value:value=int(fp.derivative()(root))
                bits.append(int(pow(value,(p-1)//2,p)==p-1))
            sig=r.pack(bits);local.append(sig);signatures[i]|=sig<<offset
        fingerprints.append({'prime':p,'roots':roots,'signatures':local});offset+=3
    rank=r.rank(signatures)
    assert rank<=3 # every added point has its double in the original three-generator subgroup
    return {'schema':'rank-jump.shared-value-halving.v1',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL)},
        'subset_divisions':records,'all_points':[list(map(str,P)) for P in allpoints],
        'joint_fingerprints':signatures,'finite_fingerprints':fingerprints,
        'original_fingerprint_rank':r.rank(signatures[:3]),'enlarged_fingerprint_rank':rank,
        'original_subgroup_rank':3 if rank==3 else 'UNKNOWN',
        'boundary':'Every added point doubles into the original subgroup, so rational ranks agree. No full saturation or full curve rank is claimed.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args()
    data=compute()
    if args.mode=='check':assert r.read(OUTPUT)==data;print('PASS exact subgroup halving')
    else:r.write_new(OUTPUT,data);print(data['original_subgroup_rank'])
