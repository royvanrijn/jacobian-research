"""Check the trace identities and local data for the beta scaling families."""
from pathlib import Path
from runpy import run_path

from sage.all import GF, QQ, ZZ, Matrix, Qp

existing = run_path(str(Path(__file__).with_name("verify_marked_octic_mirror_class_rank_four.py")))
F, t, beta = (existing[name] for name in ("F", "t", "beta"))
marks = existing["points"] + [beta]
trace = 49*t*t/9 + QQ(125)/3 + 196/(t*t)
assert trace == 2*beta[1]
assert t**8*F(6/t) == 6**4*F

possible = [
    2**i*3**j*7**k
    for i in (0, 1) for j in (0, 1) for k in (0, 1, 2)
    if (2**i * 7**k) % 9 in (4, 5)
]
assert set(possible) == {49, 147, 14, 42}
for c, partner in [(49, 42), (147, 14)]:
    assert c*partner == 6*7**3

for c in possible:
    for p in (2, 3, 5, 7, 11):
        exponent = 5
        if ((p in (2, 7) and c % 2) or (p == 3 and c % 3 == 0)):
            exponent = -5
        parameter = c*QQ(p)**(3*exponent)
        completion = Qp(p, prec=160)
        radical = completion(F(parameter)).sqrt()
        for coefficient, rational in marks:
            for sign in (1, -1):
                value = completion(rational(parameter)) + sign*completion(coefficient(parameter))*radical
                assert value.precision_relative() >= 3
                assert value.valuation() % 3 == 0
                if p == 3:
                    assert ZZ(value.unit_part().lift()) % 9 in (1, 8)

for c, data in [
    (49, [(13, 4, 2, 3, [1, 2, 1, 0]), (31, 1, 8, 25, [1, 1, 2, 0]),
          (37, 2, 13, 26, [1, 1, 1, 0]), (43, 5, 11, 36, [0, 2, 0, 1])]),
    (147, [(13, 1, 1, 3, [0, 2, 2, 1]), (19, 4, 2, 7, [2, 2, 1, 2]),
           (31, 4, 12, 25, [2, 1, 2, 2]), (37, 1, 10, 26, [2, 1, 2, 1])]),
]:
    rows = []
    for p, u, radical, root_of_unity, row in data:
        field = GF(p)
        parameter = field(c*u**3)
        radical = field(radical)
        root_of_unity = field(root_of_unity)
        assert root_of_unity.multiplicative_order() == 3
        assert radical != 0 and F.change_ring(field)(parameter) == radical**2

        def evaluate(f):
            return f.numerator().change_ring(field)(parameter)/f.denominator().change_ring(field)(parameter)

        for (coefficient, rational), exponent in zip(marks, row):
            value = evaluate(rational) + radical*evaluate(coefficient)
            assert value != 0 and value**((p-1)//3) == root_of_unity**exponent
        rows.append(row)
    assert Matrix(GF(3), rows).rank() == 4

print("Verified the beta trace and reciprocal identities, four local scaling slices, and two independence matrices.")
