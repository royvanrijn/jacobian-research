"""Target-free parity reduction: vectorized Babai and rational exhaustive CVP."""
from fractions import Fraction as Q
from math import isqrt
import numpy as np


class ExactParity:
    def __init__(self, gram):
        self.g = [[int(x) for x in row] for row in gram]
        self.n = len(self.g)
        n = self.n
        self.mu = [[Q(0) for _ in range(n)] for _ in range(n)]
        self.d = []
        for i in range(n):
            for j in range(i):
                self.mu[i][j] = (self.g[i][j] - sum(self.mu[i][k]*self.mu[j][k]*self.d[k] for k in range(j))) / self.d[j]
            self.d.append(Q(self.g[i][i]) - sum(self.mu[i][k]**2*self.d[k] for k in range(i)))
            if self.d[-1] <= 0:
                raise ArithmeticError('non-positive exact metric')
        self.float_mu = np.array(self.mu, dtype=float)

    def norm(self, w):
        return sum(int(w[i])*self.g[i][j]*int(w[j]) for i in range(self.n) for j in range(self.n))

    def babai(self, residues):
        w = np.asarray(residues, dtype=np.int64).copy()
        for i in range(self.n-1, -1, -1):
            shift = w[:, i+1:] @ self.float_mu[i+1:, i]
            w[:, i] += 2*np.rint((-shift-w[:, i])/2).astype(np.int64)
        gg = np.array(self.g, dtype=np.int64)
        if self.n**2*int(np.abs(w).max())**2*int(np.abs(gg).max()) >= 2**62:
            raise ArithmeticError('Babai norm integer overflow')
        norms = np.einsum('ij,jk,ik->i', w, gg, w, optimize=True)
        return w, norms

    def solve(self, residue, seed, node_limit=2000000):
        """Enumerate the complete closed ellipsoid; return all minimum words.

        The initial radius is an exactly evaluated feasible vector. Every branch
        exclusion uses rational LDL arithmetic and an integer square root.
        A node-limit failure is not a CVP certificate.
        """
        p = [int(x) % 2 for x in residue]
        if [int(x) % 2 for x in seed] != p:
            raise ArithmeticError('wrong seed parity')
        best = self.norm(seed)
        minima = {tuple(map(int, seed))}
        w = [0]*self.n
        nodes = 0

        def visit(i, used):
            nonlocal best, minima, nodes
            nodes += 1
            if nodes > node_limit:
                raise RuntimeError('exact CVP node budget exhausted')
            if i < 0:
                value = self.norm(w)
                if value < best:
                    best, minima = value, {tuple(w)}
                elif value == best:
                    minima.add(tuple(w))
                return
            remaining = Q(best)-used
            if remaining < 0:
                return
            shift = sum((self.mu[j][i]*w[j] for j in range(i+1, self.n)), Q(0))
            a, b = shift.numerator, shift.denominator
            radius = remaining/self.d[i]*b*b
            rad = isqrt(radius.numerator//radius.denominator)
            lo = -((rad+a)//b)
            hi = (rad-a)//b
            lo += (p[i]-lo) % 2
            values = range(lo, hi+1, 2)
            for x in sorted(values, key=lambda x: abs(Q(x)+shift)):
                w[i] = x
                visit(i-1, used+self.d[i]*(Q(x)+shift)**2)

        visit(self.n-1, Q(0))
        return {'norm': best, 'minima': sorted(minima), 'nodes': nodes,
                'certificate': 'complete exact rational LDL ellipsoid enumeration'}
