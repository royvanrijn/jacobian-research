#!/usr/bin/env python3
"""Independent standard-library audit of both C15 reductions.

The generators use SymPy.  This implementation deliberately does not: sparse
polynomials are dictionaries and every scalar operation uses Fraction.  It
starts from C01, regenerates the 22 stable extensions and the 95-dimensional
form, then independently polarizes it and verifies the compact GZ matrices.
"""

from fractions import Fraction as Q
from itertools import permutations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
H_FILE = ROOT / "results" / "cubic_homogeneous_counterexample.json"
A_FILE = ROOT / "results" / "cubic_linear_counterexample.json"
Poly = dict[tuple[int, ...], Q]


def clean(p: Poly) -> Poly:
    return {e: c for e, c in p.items() if c}


def add(*ps: Poly) -> Poly:
    out: Poly = {}
    for p in ps:
        for e, c in p.items():
            out[e] = out.get(e, Q(0)) + c
    return clean(out)


def scale(p: Poly, c: Q) -> Poly:
    return clean({e: c * a for e, a in p.items()})


def mul(p: Poly, q: Poly) -> Poly:
    out: Poly = {}
    for e, a in p.items():
        for f, b in q.items():
            g = tuple(x + y for x, y in zip(e, f))
            out[g] = out.get(g, Q(0)) + a * b
    return clean(out)


def power(p: Poly, n: int) -> Poly:
    out = {tuple(0 for _ in next(iter(p))): Q(1)}
    for _ in range(n):
        out = mul(out, p)
    return out


def var(n: int, i: int) -> Poly:
    e = [0] * n
    e[i] = 1
    return {tuple(e): Q(1)}


def const(n: int, value) -> Poly:
    return {tuple([0] * n): Q(value)} if value else {}


def pad(p: Poly, amount=2) -> Poly:
    return {e + (0,) * amount: c for e, c in p.items()}


def mono(n: int, e, coefficient=1) -> Poly:
    return {tuple(e) + (0,) * (n - len(e)): Q(coefficient)}


def derivative(p: Poly, i: int) -> Poly:
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            out[tuple(f[:i] + [f[i] - 1] + f[i+1:])] = c * e[i]
    return out


def determinant3(m):
    out = {}
    for perm in permutations(range(3)):
        inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
        term = const(len(next(iter(m[0][0]), (0, 0, 0))), -1 if inversions % 2 else 1)
        for i in range(3):
            term = mul(term, m[i][perm[i]])
        out = add(out, term)
    return out


def evaluate(p: Poly, point) -> Q:
    return sum(c * product(x**a for x, a in zip(point, e)) for e, c in p.items())


def product(values):
    answer = Q(1)
    for value in values:
        answer *= value
    return answer


def split_two(e):
    first, second, left = [0] * len(e), list(e), 2
    for i in range(len(e)):
        take = min(second[i], left)
        first[i], second[i], left = take, second[i] - take, left - take
        if not left:
            break
    return tuple(first), tuple(second)


def dense(support, n):
    e = [0] * n
    for i, a in support:
        e[i] = a
    return tuple(e)


def regenerate_h():
    n = 3
    x, y, z = (var(n, i) for i in range(3))
    u = add(const(n, 1), mul(x, y))
    f = [
        scale(add(scale(x, 2), scale(mul(power(x, 2), y), -3), scale(mul(power(x, 3), z), -1)), Q(1, 2)),
        add(y, scale(mul(mul(x, power(u, 2)), z), 3), scale(mul(mul(x, power(y, 2)), add(const(n, 4), scale(mul(x, y), 3))), 3)),
        add(mul(power(u, 3), z), mul(mul(power(y, 2), u), add(const(n, 4), scale(mul(x, y), 3)))),
    ]
    jac = [[derivative(f[i], j) for j in range(3)] for i in range(3)]
    assert determinant3(jac) == {(0, 0, 0): Q(1)}
    points = [[Q(0), Q(0), -Q(1,4)], [Q(1), -Q(3,2), Q(13,2)], [-Q(1), Q(3,2), Q(13,2)]]

    steps = 0
    while max(sum(e) for p in f for e in p) > 3:
        degree = max(sum(e) for p in f for e in p)
        chosen = None
        for component, p in enumerate(f):
            candidates = [(e, c) for e, c in p.items() if sum(e) == degree]
            if candidates:
                e, coefficient = max(candidates)
                chosen = component, e, coefficient
                break
        component, e, coefficient = chosen
        a_e, b_e = split_two(e)
        for point in points:
            av = coefficient * product(value**exponent for value, exponent in zip(point, a_e))
            bv = product(value**exponent for value, exponent in zip(point, b_e))
            point.extend([-av, -bv])
        n += 2
        f = [pad(p) for p in f]
        a, b = mono(n, a_e, coefficient), mono(n, b_e)
        Y, Z = var(n, n-2), var(n, n-1)
        f[component] = add(f[component], scale(mul(add(Y, a), add(Z, b)), -1))
        f.extend([add(Y, a), add(Z, b)])
        steps += 1
    assert (steps, n) == (22, 47)
    assert max(sum(e) for p in f for e in p) == 3

    N, t = 2*n + 1, 2*n
    H = [dict() for _ in range(N)]
    cubic = []
    for i, p in enumerate(f):
        qpart = {e: c for e, c in p.items() if sum(e) == 2}
        cpart = {e: c for e, c in p.items() if sum(e) == 3}
        cubic.append(cpart)
        e = [0] * N
        e[n+i], e[t] = 1, 2
        H[i][tuple(e)] = -Q(1)
        for old, c in qpart.items():
            e = list(old) + [0] * (n+1)
            e[t] = 1
            H[i][tuple(e)] = c
        for old, c in cpart.items():
            H[n+i][old + (0,) * (n+1)] = c

    stored = json.loads(H_FILE.read_text())
    expected = [{dense(term["monomial"], N): Q(term["coefficient"]) for term in component} for component in stored["H"]]
    assert H == expected
    lifted = []
    for point in points:
        cp = [evaluate(component, point) for component in cubic]
        lifted.append(point + [-value for value in cp] + [Q(1)])
    assert lifted == [[Q(value) for value in point] for point in stored["collision_points"]]
    return H


def canonical(form, coefficient):
    first = next(value for value in form if value)
    return (tuple(-x for x in form), -coefficient) if first < 0 else (tuple(form), coefficient)


def polarize(H):
    cubes = {}
    for output, component in enumerate(H):
        for e, coefficient in component.items():
            indices = [i for i, a in enumerate(e) if a]
            target = -coefficient
            pieces = []
            if len(indices) == 1:
                pieces = [({indices[0]: 1}, Q(1))]
            elif len(indices) == 2:
                a = next(i for i in indices if e[i] == 1)
                b = next(i for i in indices if e[i] == 2)
                pieces = [({a:1,b:1},Q(1,6)), ({a:1,b:-1},Q(1,6)), ({a:1},-Q(1,3))]
            else:
                a,b,c = indices
                pieces = [({a:1,b:1,c:1},Q(1,24)), ({a:1,b:-1,c:-1},Q(1,24)),
                          ({a:1,b:1,c:-1},-Q(1,24)), ({a:1,b:-1,c:1},-Q(1,24))]
            for entries, factor in pieces:
                form = [0] * len(H)
                for i, value in entries.items(): form[i] = value
                form, value = canonical(form, target * factor)
                column = cubes.setdefault(form, {})
                column[output] = column.get(output, Q(0)) + value
                if not column[output]: del column[output]
    return {f:c for f,c in cubes.items() if c}


def rank(rows):
    a = [dict(row) for row in rows]
    pivots, r = [], 0
    columns = 1 + max((j for row in a for j in row), default=-1)
    for c in range(columns):
        pivot = next((i for i in range(r, len(a)) if a[i].get(c, 0)), None)
        if pivot is None: continue
        a[r], a[pivot] = a[pivot], a[r]
        value = a[r][c]
        a[r] = {j:v/value for j,v in a[r].items()}
        for i in range(len(a)):
            if i != r and a[i].get(c, 0):
                factor = a[i][c]
                a[i] = clean({j: a[i].get(j,Q(0))-factor*v for j,v in a[r].items()} | {
                    j:v for j,v in a[i].items() if j not in a[r]
                })
        pivots.append(c); r += 1
        if r == len(a): break
    return r


def rows(data):
    return [{int(i):Q(v) for i,v in row} for row in data]


def audit_linear(H):
    data = json.loads(A_FILE.read_text())
    cubes, forms = polarize(H), None
    forms = sorted(cubes)
    assert len(forms) == 415
    B = rows(data["pairing"]["B_columns"])
    D = rows(data["pairing"]["D_rows"])
    C = rows(data["pairing"]["C_rows"])
    A = rows(data["A_rows"])
    assert B[:415] == [cubes[f] for f in forms]
    assert D[:415] == [{i:Q(v) for i,v in enumerate(f) if v} for f in forms]
    b0_rows = [{j:B[j][i] for j in range(415) if i in B[j]} for i in range(95)]
    assert rank(b0_rows) == 59
    for i in range(95):
        for j in range(95):
            assert sum(B[k].get(i,0)*C[k].get(j,0) for k in range(451)) == (i == j)
    b_rows = [{j:B[j][i] for j in range(451) if i in B[j]} for i in range(95)]
    for i in range(451):
        expected = {}
        for k,d in D[i].items():
            for j,b in b_rows[k].items(): expected[j] = expected.get(j,Q(0))+d*b
        assert clean(expected) == A[i]
        for j in range(95):
            assert sum(A[i].get(k, 0)*C[k].get(j, 0) for k in A[i]) == D[i].get(j, 0)
    target = [Q(v) for v in data["common_image"]]
    points = [[Q(v) for v in p] for p in data["collision_points"]]
    for point in points:
        linear = [sum(v*point[j] for j,v in row.items()) for row in A]
        assert [x-y**3 for x,y in zip(point,linear)] == target
    assert len({tuple(p) for p in points}) == 3


def main():
    H = regenerate_h()
    print("PASS (stdlib): regenerated the 95D stable cubic-homogeneous form from C01")
    audit_linear(H)
    print("PASS (stdlib): regenerated 415 cube forms and rank(B0)=59")
    print("PASS (stdlib): verified the 451D GZ pairing and its exact collision")


if __name__ == "__main__":
    main()
