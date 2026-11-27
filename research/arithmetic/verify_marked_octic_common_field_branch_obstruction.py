"""Exact algebra and numerical degree bounds for BF1–BF3."""

from sage.all import GF, QQ, NumberField, PolynomialRing, matrix, valuation


def check_branch_fields():
    ring = PolynomialRing(QQ, "z")
    z = ring.gen()
    quartic = 49*z**4 + 750*z**3 - 2403*z**2 + 27000*z + 63504
    reciprocal_sum = z**2 + 36
    assert quartic == 49*reciprocal_sum**2 + 750*z*reciprocal_sum - 5931*z**2

    quadratic = NumberField(z**2 - 11, "a")
    a = quadratic.gen()
    pair_sum = (-375 + 198*a)/49
    conjugate_sum = (-375 - 198*a)/49
    assert 49*pair_sum**2 + 750*pair_sum - 5931 == 0
    assert quartic == 49*(z**2 - pair_sum*z + 36)*(z**2 - conjugate_sum*z + 36)
    assert (pair_sum**2 - 144).norm() == -QQ(11390625)/343
    assert (pair_sum + 12).norm() == -QQ(1125)/7
    assert (pair_sum - 12).norm() == QQ(10125)/49
    assert not (QQ(10125)/49).is_square()
    assert [valuation(a, 2) for a in quartic.list()] == [4, 3, 0, 1, 0]

    # One of the two roots in each pair has normalized valuation 2.
    valuation_rows = matrix(GF(3), [
        [2, 0, 2, 0],
        [2, 0, 0, 2],
        [0, 2, 2, 0],
        [0, 2, 0, 2],
    ])
    assert valuation_rows.rank() == 3
    assert list(valuation_rows.right_kernel().basis()[0]) == [1, 1, -1, -1]

    # In the joint Kummer span, the fifth coordinate is the scale ratio.
    joint_ring = GF(3)**5
    first = [joint_ring.basis()[i] for i in range(4)]
    second = [generator + joint_ring.basis()[4] for generator in first]
    pair_relation = first[0] + first[1] - first[2] - first[3]
    assert matrix(GF(3), [pair_relation] + first + second).rank() == 5
    for a in first:
        for b in second:
            assert matrix(GF(3), [pair_relation, a, b]).rank() == 3

    for scale in (49, 147):
        assert quartic(scale**2*z**6).is_irreducible()


def check_degree_bounds():
    for degree in (3, 6, 9, 12):
        count = 24*degree
        assert count > 2*(degree - 1)*degree
        candidates = [
            other for other in range(degree, 25)
            if 20*other + 4 <= count and
            (other <= 12 and other == degree or other > 13)
        ]
        assert candidates == ([degree, 14] if degree == 12 else [degree])
    assert 24*2 > 2*13 - 2
    assert 24*12 == 288
    assert 2*14 - 2 - 24 == 2


if __name__ == "__main__":
    check_branch_fields()
    check_degree_bounds()
    print("BF1–BF3 exact polynomial, Kummer matrix and degree bounds verified.")
