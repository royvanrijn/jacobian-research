#!/usr/bin/env sage-python
"""Exact local arithmetic for the rational CM locus of X(546)/w546.

Completeness uses the class-number-one classification, quadratic genus theory,
and Gonzalez--Rotger Proposition5.6/Theorem5.12/Corollary5.14. Coordinates of the
CM points are not recovered. No point on the elliptic model is labelled non-CM.
"""
import argparse,json
from pathlib import Path
from hashlib import sha256
from sage.all import ZZ, QQ, GF, EllipticCurve, QuadraticField, kronecker, gcd
from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives

def build():
    # Complete external classification of negative quadratic order discriminants.
    ds=[-3,-4,-7,-8,-11,-12,-16,-19,-27,-28,-43,-67,-163]
    records=[]
    for d in ds:
        forms=BinaryQF_reduced_representatives(ZZ(d),primitive_only=True)
        assert len(forms)==1
        field_d=int(QuadraticField(d).discriminant())
        symbols=[int(kronecker(field_d,p)) for p in [2,3,7,13]]
        # These 13 orders have no conductor divisible by 2,3,7,13 unless their
        # underlying field already splits at one of those primes. All excluded
        # rows below therefore have a field-splitting obstruction, independently
        # of the possible additional optimal-embedding conductor obstruction.
        split=[p for p,s in zip([2,3,7,13],symbols) if s==1]
        records.append(dict(discriminant=d,field_discriminant=field_d,
                            primitive_reduced_forms=[list(map(int,f)) for f in forms],
                            kronecker_at_2_3_7_13=symbols,split_primes=split,
                            survives=not split))
    surviving=[r for r in records if r['survives']]
    assert [r['discriminant'] for r in surviving]==[-67,-163]
    assert all(r['kronecker_at_2_3_7_13']==[-1]*4 for r in surviving)
    E=EllipticCurve(QQ,[1,0,1,-137,380]);P=E(-9,-26)
    counts=[int(E.change_ring(GF(p)).cardinality()) for p in [5,11,17]]
    assert counts==[8,16,20] and gcd(counts)==4 and 4*P!=E(0)
    points=[n*P for n in range(17)]
    assert len(set(points))==17
    return dict(schema='det1092.rational-cm-locus.v1',status='EXACT_CM_DISCRIMINANTS_WITH_EXTERNAL_COMPLETENESS',
        class_number_one_order_rows=records,
        no_inert_prime_case=dict(required_ramified_primes=[2,3,7,13],
                                 genus_theory_class_number_divisor=8,
                                 maximum_class_number_for_quadratic_top_point=2,
                                 excluded=True),
        surviving_discriminants=[-67,-163],geometric_CM_points_upstairs_per_discriminant=16,
        rational_CM_points_on_full_Fricke_quotient_per_discriminant=8,total_rational_CM_points=16,
        seventeen_distinct_rational_points=[list(map(str,Q)) for Q in points],
        at_least_one_of_seventeen_is_non_CM=True,individual_non_CM_indices=[],
        boundary='CM coordinates and an explicit moduli identification remain unknown. Existence in a finite list does not identify a non-CM point or satisfy the individual-point construction gate.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=build();result['checker_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({'status':result['status'],'CM_discriminants':result['surviving_discriminants'],
                      'rational_CM_point_count':16,'individual_non_CM_indices':[]}))
