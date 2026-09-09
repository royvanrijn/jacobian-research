"""Exact V3 parity CVP with integer-scaled LDL, for newly frozen runs only.

The rational decomposition is computed once by the reference implementation.
Clearing its fixed denominators makes every branch calculation integer exact.
Traversal, ties, node accounting and the closed ellipsoid are unchanged.
"""
from math import isqrt, lcm

from visibility_lattice_v2 import ExactParity as ReferenceParity


class IntegerExactParity(ReferenceParity):
    def __init__(self, gram):
        super().__init__(gram)
        self.shift_denominators = [lcm(*(self.mu[j][i].denominator
                                       for j in range(i + 1, self.n)))
                                   for i in range(self.n)]
        self.shift_coefficients = [
            [int(self.mu[j][i] * self.shift_denominators[i])
             for j in range(i + 1, self.n)] for i in range(self.n)]
        weights = [self.d[i] / self.shift_denominators[i]**2 for i in range(self.n)]
        self.norm_scale = lcm(*(x.denominator for x in weights))
        self.integer_weights = [int(x * self.norm_scale) for x in weights]
        if any(x <= 0 for x in self.integer_weights):
            raise ArithmeticError('non-positive integer LDL weight')

    def solve(self, residue, seed, node_limit=2000000):
        p = [int(x) % 2 for x in residue]
        if len(p) != self.n or len(seed) != self.n:
            raise ArithmeticError('wrong vector dimension')
        if [int(x) % 2 for x in seed] != p:
            raise ArithmeticError('wrong seed parity')
        best = self.norm(seed)
        minima = {tuple(map(int, seed))}
        w = [0] * self.n
        nodes = 0
        scale = self.norm_scale

        def visit(i, used):
            nonlocal best, minima, nodes
            nodes += 1
            if nodes > node_limit:
                raise RuntimeError('exact CVP node budget exhausted')
            if i < 0:
                value = self.norm(w)
                if used != value * scale:
                    raise ArithmeticError('integer LDL norm identity failed')
                if value < best:
                    best, minima = value, {tuple(w)}
                elif value == best:
                    minima.add(tuple(w))
                return
            remaining = best * scale - used
            if remaining < 0:
                return
            b = self.shift_denominators[i]
            a = sum(c * w[j] for j, c in enumerate(self.shift_coefficients[i], i + 1))
            weight = self.integer_weights[i]
            rad = isqrt(remaining // weight)
            lo, hi = -((rad + a) // b), (rad - a) // b
            lo += (p[i] - lo) % 2
            for x in sorted(range(lo, hi + 1, 2), key=lambda x: abs(b * x + a)):
                w[i] = x
                visit(i - 1, used + weight * (b * x + a)**2)

        visit(self.n - 1, 0)
        return {'norm': best, 'minima': sorted(minima), 'nodes': nodes,
                'certificate': 'complete exact rational LDL ellipsoid enumeration'}
