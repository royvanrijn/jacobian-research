"""Independently verify a parity minimum by exact fixed-radius enumeration.

The supplied minimum sets the radius, never a branch exclusion by itself.
An independent rational Schur elimination supplies every lower bound. The
closed ellipsoid is exhausted, so a smaller vector or an omitted tie fails.
No numerical height, optimizer traversal or optimizer node count is trusted.
"""
from math import isqrt
from operator import index

from gmpy2 import mpq


def integer(value):
    if isinstance(value, bool):
        raise ValueError('boolean is not an integer certificate field')
    return index(value)


class ParityMinimumVerifier:
    def __init__(self, gram):
        self.gram = tuple(tuple(integer(x) for x in row) for row in gram)
        n = self.dimension = len(self.gram)
        if not n or any(len(row) != n for row in self.gram):
            raise ValueError('a nonempty square Gram matrix is required')
        if any(self.gram[i][j] != self.gram[j][i] for i in range(n) for j in range(n)):
            raise ArithmeticError('asymmetric metric')
        residual = [[mpq(v) for v in row] for row in self.gram]
        self.upper = [[mpq(int(i == j)) for j in range(n)] for i in range(n)]
        self.diagonal = []
        # Symmetric Gaussian elimination, independent of the optimizer's
        # lower-triangular coefficient recurrence and denominator scaling.
        for k in range(n):
            pivot = residual[k][k]
            if pivot <= 0:
                raise ArithmeticError('metric is not positive definite')
            self.diagonal.append(pivot)
            for j in range(k+1, n):
                self.upper[k][j] = residual[k][j]/pivot
            for i in range(k+1, n):
                for j in range(i, n):
                    value = residual[i][j]-residual[k][i]*residual[k][j]/pivot
                    residual[i][j] = residual[j][i] = value
        # Check the exact factorization rather than relying on its algorithm.
        for i in range(n):
            for j in range(i, n):
                value = sum((self.diagonal[k]*self.upper[k][i]*self.upper[k][j]
                             for k in range(min(i, j)+1)), mpq(0))
                if value != self.gram[i][j]:
                    raise ArithmeticError('exact Schur factorization differs')

    def norm(self, vector):
        v = tuple(map(integer, vector)); n = self.dimension
        if len(v) != n:
            raise ValueError('vector dimension differs')
        return sum(v[i]*self.gram[i][j]*v[j] for i in range(n) for j in range(n))

    def verify(self, parity, radius, minima, node_limit=2000000):
        n = self.dimension; p = tuple(integer(x) for x in parity)
        radius = integer(radius); node_limit = integer(node_limit)
        if len(p) != n or any(x not in (0, 1) for x in p) or radius < 0 or node_limit < 1:
            raise ValueError('invalid parity, radius or node limit')
        claimed = [tuple(map(integer, v)) for v in minima]
        if not claimed or len(set(claimed)) != len(claimed):
            raise ArithmeticError('empty or duplicate minimum list')
        for vector in claimed:
            if len(vector) != n or tuple(v % 2 for v in vector) != p or self.norm(vector) != radius:
                raise ArithmeticError('claimed point is outside the parity sphere')
        vector = [0]*n; found = set(); nodes = 0

        def visit(i, used):
            nonlocal nodes
            nodes += 1
            if nodes > node_limit:
                raise RuntimeError('verification node limit; no minimum certified')
            if i < 0:
                exact = self.norm(vector)
                if used != exact:
                    raise ArithmeticError('leaf quadratic identity differs')
                if exact < radius:
                    raise ArithmeticError('a strictly shorter parity vector exists')
                if exact != radius:
                    raise ArithmeticError('leaf lies outside the closed ellipsoid')
                found.add(tuple(vector))
                return
            remaining = radius-used
            if remaining < 0:
                return
            shift = sum((self.upper[i][j]*vector[j] for j in range(i+1, n)), mpq(0))
            a, b = int(shift.numerator), int(shift.denominator)
            bound = remaining*b*b/self.diagonal[i]
            root = isqrt(int(bound.numerator)//int(bound.denominator))
            lower = -((root+a)//b)
            upper = (root-a)//b
            lower += (p[i]-lower) % 2
            # Fixed radius and increasing coordinates: no optimizer sorting,
            # incumbent shrinking or supplied discovery path is reused.
            for value in range(lower, upper+1, 2):
                vector[i] = value
                visit(i-1, used+self.diagonal[i]*(value+shift)**2)

        visit(n-1, mpq(0))
        if found != set(claimed):
            raise ArithmeticError('the complete minimum set differs')
        return {'status': 'PASS_EXACT_PARITY_MINIMUM', 'norm': radius,
                'minima': [list(v) for v in sorted(found)], 'verification_nodes': nodes,
                'boundary': 'Complete closed parity ellipsoid in the supplied exact positive metric. No canonical-height or optimizer-node-count assertion.'}
