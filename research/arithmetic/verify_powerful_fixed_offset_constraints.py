#!/usr/bin/env python3
"""Generate and verify rational certificates for powerful-number box constraints.

Generation uses SciPy to propose exact fractions. Verification uses only the
standard library and checks every branch and inequality in rational arithmetic.
"""
import argparse
import json
import math
import sys
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def combine(a, b, c, d):
    return tuple(c * x + d * y for x, y in zip(a, b, strict=True))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b, strict=True))


class Model:
    """Coordinates, a pooled tail per side, three maxima, and the cost D."""

    def __init__(self, u, v, r, cutoff, saving):
        require(0 <= v <= u and 1 <= cutoff <= r and r >= 5,
                "Require r >= 5, u >= v >= 0 and 1 <= cutoff <= r.")
        self.offsets = (u, v, 0)
        self.r = r
        self.cutoff = cutoff
        self.width = 3 * (cutoff + 1) + 4
        self.cost = self.row([(self.width - 1, 1)])
        self.threshold = u + v + saving

    def index(self, side, coord):
        return side * (self.cutoff + 1) + coord

    def row(self, entries):
        values = [0] * self.width
        for col, value in entries:
            values[col] += value
        return tuple(values)

    def base(self):
        n = self.cutoff
        eq = [(self.row((self.index(side, j), 1) for j in range(n + 1)), 1)
              for side in range(3)]
        ub = []
        for side in range(3):
            mx = self.width - 4 + side
            # Pool the tail at its lowest cost and cap only explicit entries.
            ub.extend((self.row([(self.index(side, j), 1), (mx, -1)]), 0)
                      for j in range(n))
        ub.append((self.row((j, 1) for j in range(self.width - 4, self.width - 1)), 1))
        for left, right in combinations(range(3), 2):
            terms = [(self.index(side, j), self.offsets[side] + j)
                     for side in (left, right) for j in range(n + 1)]
            ub.append((self.row(terms + [(self.width - 1, -1)]), 0))
        return eq, ub

    def pairs(self, point=None):
        groups = []
        for side, offset in enumerate(self.offsets):
            indices = [j for j in range(self.cutoff)
                       if point is None or point[self.index(side, j)]]
            groups.append([(group, math.gcd(*(self.r + offset + j for j in group)))
                           for k in range(1, len(indices) + 1)
                           for group in combinations(indices, k)])
        for i, k in combinations(range(3), 2):
            for (a, g), (b, h) in product(groups[i], groups[k]):
                if math.gcd(g, h) >= 2:
                    yield i, a, k, b

    def group_rows(self, pair):
        i, a, k, b = pair
        require(isinstance(i, int) and isinstance(k, int) and 0 <= i < k < 3
                and all(s and list(s) == sorted(set(s))
                        and all(isinstance(j, int) and 0 <= j < self.cutoff for j in s)
                        for s in (a, b)), "Invalid coordinate groups.")
        g = math.gcd(*(self.r + self.offsets[side] + j for side, inds in
                       ((i, a), (k, b)) for j in inds))
        require(g >= 2, "Graph pair has no common divisor.")
        left = self.row((self.index(i, j), 1) for j in a)
        right = self.row((self.index(k, j), 1) for j in b)
        # Singleton pairs have a height margin in competitive boxes. Groups
        # may require the complementary large-height branch.
        height = Q(g, 2) if len(a) + len(b) > 2 and g <= 4 else None
        return left, right, height

    def branches(self, pair):
        left, right, height = self.group_rows(pair)
        rows = [([(left, 0)], []), ([(right, 0)], []),
                ([], [(combine(left, right, 2, 1), 1),
                       (combine(left, right, 1, 2), 1)])]
        if height is not None:
            rows.append(([], [(combine(left, right, -1, -1), -height)]))
        return rows


def violation(a, b, height):
    return min(a, b, max(2*a+b, a+2*b)-1, 1 if height is None else height-a-b)


def verify_witness(model, point):
    eq, ub = model.base()
    require(len(point) == model.width and all(x >= 0 for x in point)
            and all(dot(a, point) == b for a, b in eq)
            and all(dot(a, point) <= b for a, b in ub),
            "Witness violates the probability, maximum, or pair-cost constraints.")
    require(dot(model.cost, point) < model.threshold, "Witness misses the target.")
    require(all(point[model.index(i, model.cutoff)] == 0 for i in range(3)),
            "A witness with pooled tail mass cannot certify the full model.")
    # Removing zero-mass coordinates can only raise the gcd, strengthening
    # the corresponding graph condition at the same masses.
    for pair in model.pairs(point):
        left, right, height = model.group_rows(pair)
        require(violation(dot(left, point), dot(right, point), height) <= 0,
                f"Witness violates grouped graph {pair}.")


def verify_leaf(model, eq, ub, leaf):
    kind = leaf["kind"]
    require(kind in ("dual", "infeasible"), "Unknown leaf type.")
    a, b = tuple(map(Q, leaf["eq"])), tuple(map(Q, leaf["ub"]))
    require(len(a) == len(eq) and len(b) == len(ub) and all(x <= 0 for x in b),
            "Wrong number or signs of LP multipliers.")
    rows, weights = eq + ub, a + b
    for col, bound in enumerate(model.cost):
        lhs = sum(w * row[col] for w, (row, _) in zip(weights, rows, strict=True))
        require(lhs <= (bound if kind == "dual" else 0),
                f"Invalid LP dual inequality in column {col}.")
    value = dot(weights, [rhs for row, rhs in rows])
    require(value >= model.threshold if kind == "dual" else value > 0,
            "LP multipliers do not certify the claimed lower bound.")


def verify_tree(model, node, eq, ub):
    if node["kind"] != "split":
        verify_leaf(model, eq, ub, node)
        return
    for (de, du), child in zip(model.branches(node["pair"]), node["children"], strict=True):
        verify_tree(model, child, eq + de, ub + du)


def verify(data):
    model = Model(data["u"], data["v"], data["r"], data["cutoff"], Q(data["saving"]))
    if data["kind"] == "witness":
        verify_witness(model, tuple(map(Q, data["point"])))
    else:
        require(data["kind"] == "proof", "Unknown certificate type.")
        verify_tree(model, data["tree"], *model.base())


def generate(model):
    import numpy as np
    from scipy.optimize import linprog

    def exact(values):
        return tuple(Q(float(x)).limit_denominator(1_000_000) for x in values)

    def infeasible(eq, ub):
        rows = eq + ub
        coeff = np.array([row for row, rhs in rows], dtype=float).T
        sol = linprog([-float(b) for a, b in rows], A_ub=coeff, b_ub=np.zeros(model.width),
                      bounds=[(-1, 1)] * len(eq) + [(-1, 0)] * len(ub), method="highs")
        require(sol.success, "Cannot obtain infeasibility multipliers.")
        x = exact(sol.x)
        return {"kind": "infeasible", "eq": list(map(str, x[:len(eq)])),
                "ub": list(map(str, x[len(eq):]))}

    graphs = []
    for pair in model.pairs():
        left, right, height = model.group_rows(pair)
        graphs.append((pair, [i for i, n in enumerate(left) if n],
                       [i for i, n in enumerate(right) if n], height))
    witness = None

    def search(eq, ub, used):
        nonlocal witness
        if witness is not None:
            return None
        sol = linprog(model.cost, A_eq=[a for a, b in eq], b_eq=[float(b) for a, b in eq],
                      A_ub=[a for a, b in ub], b_ub=[float(b) for a, b in ub], method="highs")
        if sol.status == 2:
            node = infeasible(eq, ub)
            verify_leaf(model, eq, ub, node)
            return node
        require(sol.success, sol.message)
        if sol.fun >= float(model.threshold) - 1e-9:
            node = {"kind": "dual", "eq": list(map(str, exact(sol.eqlin.marginals))),
                    "ub": list(map(str, exact(sol.ineqlin.marginals)))}
            verify_leaf(model, eq, ub, node)
            return node

        def score(item):
            _, la, lb, h = item
            a, b = sum(sol.x[i] for i in la), sum(sol.x[i] for i in lb)
            return violation(a, b, None if h is None else float(h))

        chosen = max((item for item in graphs if item[0] not in used),
                     key=score, default=None)
        if chosen is None or score(chosen) <= 1e-9:
            witness = exact(sol.x)
            verify_witness(model, witness)
            return None
        pair = chosen[0]
        children = [search(eq + e, ub + u, used | {pair})
                    for e, u in model.branches(pair)]
        return {"kind": "split", "pair": pair, "children": children}

    tree = search(*model.base(), frozenset())
    return ({"kind": "witness", "point": list(map(str, witness))} if witness is not None
            else {"kind": "proof", "tree": tree})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    gen = commands.add_parser("generate", help="Propose and verify a rational certificate.")
    for name in ("u", "v", "r"):
        gen.add_argument(name, type=int)
    gen.add_argument("--cutoff", type=int, default=6)
    gen.add_argument("--saving", type=Q, default=Q(2))
    val = commands.add_parser("verify", help="Verify a JSON certificate without SciPy.")
    val.add_argument("file", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "verify":
            data = json.loads(args.file.read_text())
            verify(data)
            print("Verified", data["kind"])
        else:
            model = Model(args.u, args.v, args.r, args.cutoff, args.saving)
            data = {"u": args.u, "v": args.v, "r": args.r, "cutoff": args.cutoff,
                    "saving": str(args.saving), **generate(model)}
            verify(data)
            print(json.dumps(data, indent=2))
    except (ImportError, KeyError, TypeError, ValueError) as error:
        sys.exit(str(error))


if __name__ == "__main__":
    main()
