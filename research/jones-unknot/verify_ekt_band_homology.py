#!/usr/bin/env python3
"""Check selected EKT homology formulas against full integral Fox presentations."""

from collections import Counter
from itertools import product

from sympy import Matrix, Symbol, ZZ, expand
from sympy.matrices.normalforms import smith_normal_form


class FoxDiagram:
    def __init__(self):
        self.arc_count = 0
        self.relations = []

    def arc(self):
        index = self.arc_count
        self.arc_count += 1
        return index

    def join(self, first, second):
        if first != second:
            self.relations.append({first: 1, second: -1})

    def twist(self, first, second, sign):
        over, incoming_under = (first, second) if sign > 0 else (second, first)
        outgoing_under = self.arc()
        relation = Counter([incoming_under, outgoing_under])
        relation[over] -= 2
        self.relations.append(relation)
        return (outgoing_under, over) if sign > 0 else (over, outgoing_under)

    def integer(self, value):
        nw, sw = self.arc(), self.arc()
        ne, se = nw, sw
        for _ in range(abs(value)):
            ne, se = self.twist(ne, se, value)
        return nw, ne, sw, se

    def infinity(self):
        west, east = self.arc(), self.arc()
        return west, east, west, east

    def horizontal_sum(self, left, right):
        self.join(left[1], right[0])
        self.join(left[3], right[2])
        return left[0], right[1], left[2], right[3]

    def dot(self, left, right):
        return self.horizontal_sum(reflect(left), right)

    def omega(self, tangle, sign):
        tangle = self.horizontal_sum(tangle, self.integer(2 * sign))
        tangle = self.dot(tangle, self.integer(sign))
        return self.dot(tangle, self.integer(2 * sign))

    def hopf_clasp(self, first, second):
        strands = [first[1], first[3], second[1], second[3]]
        for index in [1, 0, 2, 1] * 2:
            strands[index : index + 2] = self.twist(*strands[index : index + 2], -1)
        for target, source in zip((first[0], first[2], second[0], second[2]), strands):
            self.join(target, source)

    def cover_homology(self):
        matrix = Matrix([[r.get(index, 0) for index in range(self.arc_count)]
                         for r in self.relations])
        diagonal = smith_normal_form(matrix, domain=ZZ).diagonal()
        free_rank = self.arc_count - sum(d != 0 for d in diagonal) - 1
        torsion = tuple(abs(int(d)) for d in diagonal if abs(d) > 1)
        return free_rank, torsion


def reflect(tangle):
    nw, ne, sw, se = tangle
    return nw, sw, ne, se


def iterate_slope(a, b, count, sign=1):
    for _ in range(count):
        a, b = 3 * a - 8 * sign * b, -sign * a + 3 * b
    return a, b


def check_ll2():
    for count in range(6):
        diagram = FoxDiagram()
        left = diagram.horizontal_sum(diagram.infinity(), diagram.integer(-2))
        right = diagram.horizontal_sum(diagram.infinity(), diagram.integer(2))
        for _ in range(count):
            left = diagram.omega(left, 1)
            right = diagram.omega(right, -1)
        diagram.hopf_clasp(left, right)
        a, b = iterate_slope(1, 0, count)
        order = a * a + 4 * b * b
        expected = (1, (order,) if order > 1 else ())
        assert diagram.cover_homology() == expected, ("LL2", count)


def check_ll1():
    for m, n, iterations in product(range(4), range(-3, 4), range(-2, 4)):
        diagram = FoxDiagram()
        left = diagram.dot(diagram.integer(n), diagram.integer(1))
        for _ in range(m):
            left = diagram.dot(left, diagram.integer(1))
            left = diagram.dot(left, reflect(diagram.integer(2)))
        left = diagram.dot(left, diagram.integer(-1))
        right = diagram.horizontal_sum(diagram.infinity(), diagram.integer(2))
        sign = 1 if iterations >= 0 else -1
        for _ in range(abs(iterations)):
            left = diagram.omega(left, sign)
            right = diagram.omega(right, -sign)
        diagram.hopf_clasp(left, right)

        x, y = n + 1, -n
        for _ in range(m):
            x, y = 3 * x - y, -2 * x + 2 * y
        assert x + y == 1
        a, c = iterate_slope(1, 0, abs(iterations), sign=-sign)
        order = abs(1 + 12 * c * (c - a * x))
        assert (order > 1) == (iterations != 0)
        expected = (1, (order,) if order > 1 else ())
        assert diagram.cover_homology() == expected, ("LL1", m, n, iterations)
        if iterations == 1:
            assert order == abs(12 * (3 * n + 2) * 4**m - 1)


def check_clasp_bracket():
    z = Symbol("z", nonzero=True)  # z = A**2
    delta = -z - 1 / z
    h00 = -z**7 - z**3 - 2 * z - 2 / z - z**-3 - z**-7
    assert expand(h00 - (delta**7 - 7 * delta**5 + 15 * delta**3 - 8 * delta)) == 0


if __name__ == "__main__":
    controls = ((0, 0, (1, (2, 2))), (0, None, (2, ())), (None, None, (1, ())))
    for first, second, expected in controls:
        diagram = FoxDiagram()
        left = diagram.infinity() if first is None else diagram.integer(first)
        right = diagram.infinity() if second is None else diagram.integer(second)
        diagram.hopf_clasp(left, right)
        assert diagram.cover_homology() == expected
    check_ll2()
    check_ll1()
    check_clasp_bracket()
    print("Smith-form agreement for 3 controls, 6 LL2 members, and 168 LL1 orbit members")
