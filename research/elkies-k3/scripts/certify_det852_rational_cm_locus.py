from sage.all import *
from sage.quadratic_forms.binary_qf import BinaryQF_reduced_representatives
D = ZZ(6)
N = ZZ(71)
M = ZZ(426)

# Complete Watkins lists of imaginary quadratic orders of class number one
# and two, represented by (conductor, fundamental discriminant), as pinned in
# Padurariu--Saia GenusAtMost2 commit 6cc368fe37aa67187783118f18d149b2b1fd6230.
CLASS_NUMBER_ONE_ORDERS = [
    (1, -3), (2, -3), (3, -3), (1, -4), (2, -4),
    (1, -7), (2, -7), (1, -8), (1, -11), (1, -19),
    (1, -43), (1, -67), (1, -163),
]
CLASS_NUMBER_TWO_ORDERS = [
    (4, -3), (5, -3), (7, -3), (3, -4), (4, -4),
    (5, -4), (4, -7), (2, -8), (3, -8), (3, -11),
    (1, -15), (2, -15), (1, -20), (1, -24), (1, -35),
    (1, -40), (1, -51), (1, -52), (1, -88), (1, -91),
    (1, -115), (1, -123), (1, -148), (1, -187),
    (1, -232), (1, -235), (1, -267), (1, -403),
    (1, -427),
]


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order_record(class_number, conductor, field_discriminant):
    conductor = ZZ(conductor)
    field_discriminant = ZZ(field_discriminant)
    order_discriminant = conductor**2*field_discriminant
    assert len(BinaryQF_reduced_representatives(
        order_discriminant, primitive_only=True
    )) == class_number

    def eichler_symbol(prime):
        return ZZ(1) if conductor % prime == 0 else ZZ(
            kronecker(field_discriminant, prime)
        )

    local_factors = {
        "2": int(1-eichler_symbol(2)),
        "3": int(1-eichler_symbol(3)),
        "71": int(1+eichler_symbol(71)),
    }
    top_count = ZZ(class_number)*prod(local_factors.values())

    d_r = prod(
        prime for prime in (2, 3)
        if conductor % prime and kronecker(field_discriminant, prime) == -1
    )
    n_r = prod(
        prime for prime in (71,)
        if conductor % prime == 0
        or kronecker(field_discriminant, prime) == 1
    )
    nstar_r = prod(
        prime for prime in (71,)
        if conductor % prime and kronecker(field_discriminant, prime) == 1
    )
    m_r = gcd(M, abs(order_discriminant)//gcd(N, conductor))
    quotient = M//m_r
    d_r_nstar_r = ZZ(d_r*nstar_r)

    # Corollary 5.14: for h=1 the criterion below is exact.  For h=2 a
    # rational field can occur only in the two-involution case.  No nonempty
    # class-number-two row reaches that case, so no Artin-class ambiguity is
    # left in this instance.
    if class_number == 1:
        rational_image = bool(
            top_count
            and (d_r_nstar_r == 1 or quotient == d_r_nstar_r)
        )
    else:
        rational_image = bool(
            top_count and d_r_nstar_r == 1 and quotient == 1
        )

    marked_count = 0
    if rational_image:
        # The fixed locus of w_426 has order discriminant -1704.  The two
        # surviving orders are different, so w_426 acts freely on each locus.
        assert order_discriminant != -1704
        assert top_count % 2 == 0
        marked_count = top_count//2

    return {
        "class_number": int(class_number),
        "conductor": int(conductor),
        "field_discriminant": int(field_discriminant),
        "order_discriminant": int(order_discriminant),
        "local_embedding_factors_p2_p3_p71": local_factors,
        "top_curve_cm_points": int(top_count),
        "D_R": int(d_r),
        "N_R": int(n_r),
        "N_star_R": int(nstar_r),
        "m_R": int(m_r),
        "m_over_m_R": int(quotient),
        "rational_image_on_w426_quotient": rational_image,
        "marked_curve_rational_cm_points": int(marked_count),
    }


import json,hashlib,argparse
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUTPUT=ROOT/'artifacts/generated-results/elkies-k3-det852-cm-gate-v1/certificate.json'
records=[order_record(h,f,d) for h,orders in [(1,CLASS_NUMBER_ONE_ORDERS),(2,CLASS_NUMBER_TWO_ORDERS)] for f,d in orders]
survivors=[r for r in records if r['rational_image_on_w426_quotient']]
assert [r['order_discriminant'] for r in survivors]==[-67,-163]
assert [r['marked_curve_rational_cm_points'] for r in survivors]==[4,4]
assert not any(r['class_number']==2 and r['top_curve_cm_points'] and r['D_R']*r['N_star_R']==1 and r['m_over_m_R']==1 for r in records)
E=EllipticCurve(QQ,[1,1,0,-286,1780]);P=E(7,-17);points=[n*P for n in range(9)]
assert len(set(points))==9
inputs=[ROOT/'artifacts/generated-results/elkies-k3-det852-marking-gate-v1/certificate.json',ROOT/'elkies-k3/scripts/certify_det1236_rational_cm_locus.sage']
out={'schema':'elkies-k3.det852-rational-CM-locus.v1','status':'PASS_COMPLETE_EIGHT_RATIONAL_CM_POINTS_COORDINATES_UNRESOLVED','curve':'X_0^6(71)/<w426>','orders_checked':records,'rational_CM_orders':[-67,-163],'rational_CM_point_count':8,'elliptic_test_points':[[str(c) for c in Q] for Q in points],'finite_test_set_conclusion':'Under any Q-isomorphism of the full marked curve with426b1, at least one of these nine distinct points is non-CM. No particular member is identified here.','theorem_inputs':['Complete Watkins class-number-one/two imaginary quadratic order lists','Gonzalez-Rotger Theorem5.8 and Corollary5.14','Optimal embedding count for squarefree Eichler level coprime to quaternion discriminant','Previously certified full marked curve and infinite-order point'],'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in inputs],'boundary':'Complete rational CM count and discriminants, not CM coordinates or a particular non-CM witness. No new surface/divisor descent or rootless/MW17 claim.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
if args.check:assert json.loads(OUTPUT.read_text())==out
else:OUTPUT.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS determinant852: exactly eight rational CM points, discriminants-67 and-163; nine-point test set contains a non-CM point but no particular member is labelled.')
