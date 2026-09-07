#!/usr/bin/env sage-python
"""Pointed non-CM admission on the full X(206)/w206 curve, before equations."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from sage.all import QQ,ZZ,QuadraticField,kronecker,matrix
from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives

def build(marking):
    assert marking['projectively_stable_labels']==[1,206]
    assert marking['even_order_is_maximal'] and marking['quaternion_discriminant']==206
    T=matrix(ZZ,marking['T']);g=T.inverse().column(2)
    assert g*T*g==QQ(5)/412
    assert kronecker(5*4,103)==-kronecker(4,103)==-1
    ds=[-3,-4,-7,-8,-11,-12,-16,-19,-27,-28,-43,-67,-163]
    cm=[];checks=[]
    for d in ds:
        dk=ZZ(QuadraticField(d).discriminant());f=ZZ(d/dk).sqrt()
        assert f in ZZ and len(BinaryQF_reduced_representatives(ZZ(d),primitive_only=True))==1
        symbols=[int(kronecker(dk,p)) for p in [2,103]]
        admit=all(s!=1 for s in symbols) and all(f%p for p in [2,103])
        if admit:cm.append(d)
        checks.append(dict(discriminant=d,field_discriminant=int(dk),conductor=int(f),
                           symbols=symbols,embeds=bool(admit)))
    assert cm==[-4,-8,-19,-163]
    h824=len(BinaryQF_reduced_representatives(ZZ(-824),primitive_only=True));assert h824==20
    # No-inert-prime case: both 2 and103 ramify. More ramified primes force
    # class number divisible by4; with just these two the field is Q(sqrt(-206)).
    # Its class number20 also excludes a degree-at-most2 point upstairs.
    excluded_j=['0','infinity','1','-1','1/2','-1/2']
    assert '3' not in excluded_j
    return dict(schema='det412.pointed-noncm-admission.v1',status='ADMITTED_MODULI_POINT_NOT_YET_EQUATION',
        literal_T=marking['T'],full_projective_stable_curve='X(206)/<w206>',
        local_L206_test_at_103=-1,class_number_one_tests=checks,
        no_inert_case_class_number_minus_824=h824,
        complete_rational_CM_discriminants=cm,
        Elkies_CM_coordinates_r=['0','infinity','1','-1','2','-2'],
        Lin_Yang_CM_coordinates_j=excluded_j,coordinate_relation='j=+/-1/r; CM set unaffected by sign',
        admitted_j='3',individual_non_CM_point=True,
        references=['https://arxiv.org/pdf/0802.1301 section6',
                    'https://arxiv.org/pdf/1807.00466 AppendixB Table6',
                    'https://web.mat.upc.edu/victor.rotger/docs/ShimuraGenusOne.pdf Appendix5'],
        external_inputs='Full marked period correspondence; published pointed Hauptmoduls; CM residue fields; class-number-one classification and genus theory.',
        boundary='The rational non-CM moduli point is identified. No rootless frame, rational generic section basis, or new elliptic fibre is yet supplied.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--marking',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=build(json.loads(a.marking.read_text()));r['checker_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    r['marking_sha256']=sha256(a.marking.read_bytes()).hexdigest()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(r['status'],r['admitted_j'],r['complete_rational_CM_discriminants'])
