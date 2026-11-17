"""Complete SO(q) cyclic-frame quotient and exact moment circuits.

The normalized model fixes D=1/8, hence det(G)=-1/64 for q=2zw+t^2.
No Gram entry, eigenvalue or discriminant is separately inverted.
"""
from collections import Counter
from functools import lru_cache
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path

import sympy as S

p, r, s0, s1, s2, D, c0 = S.symbols('p r s0 s1 s2 D c0')
TRIPLES = tuple((i, j, k) for i in range(3) for j in range(i, 3) for k in range(j, 3))
U = {indices: S.Symbol('u'+''.join(map(str, indices))) for indices in TRIPLES}
VARIABLES = (p, r, s0, s1, s2, *U.values())
G = S.Matrix([[s0, s1, s2], [s1, s2, p*s1+r*s0],
              [s2, p*s1+r*s0, p*s2+r*s1]])
C = S.Matrix([[0, 0, r], [1, 0, p], [0, 1, 0]])
DET = S.expand(G.det())


def tensor(i, j, k):
    return U[tuple(sorted((i, j, k)))]


def trace_relations():
    return tuple(S.expand(sum(G[i, j]*tensor(i, j, k) for i in range(3) for j in range(3)))
                 for k in range(3))


def normalized_relations():
    return (S.expand(DET+S.Rational(1, 64)), *trace_relations())


def phase_polynomial():
    """The phase in dual frame coordinates xi_i=<B^i v,x>."""
    xi = S.Matrix(S.symbols('xi0:3'))
    quadratic = -64*C*G.adjugate()
    cubic = sum(S.Rational(6, math.prod(math.factorial(v) for v in Counter(indices).values()))
                *U[indices]*math.prod(xi[j] for j in indices) for indices in TRIPLES)
    return S.expand(xi[0]+(xi.T*quadratic*xi)[0]+cubic), tuple(xi)


def rodrigues_point(value=S.Rational(1, 4)):
    point = {variable: S.Integer(0) for variable in VARIABLES}
    point[s2] = value
    point[U[(2, 2, 2)]] = -128*value
    return point


def frame_data(v, B, metric):
    E = S.Matrix.hstack(v, B*v, B*B*v)
    return E, S.simplify(E.det()), S.simplify(E.T*metric*E)


def krylov_stratum(v, B, metric):
    frame, determinant, gram = frame_data(v, B, metric)
    rank, gram_rank = frame.rank(), gram.rank()
    allowed = {(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (3, 3)}
    assert (rank, gram_rank) in allowed
    return rank, gram_rank


def invariants_from_phase(v, B, cubic, coordinates):
    """Convert an original phase to complete cyclic-frame invariants.

    This does not rescale D or impose a spectral/Gram-entry condition.
    """
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    assert B.T*metric == metric*B and S.trace(B) == 0
    E, determinant, gram = frame_data(v, B, metric)
    if determinant == 0:
        raise ValueError('the full cyclic-frame adapter requires D nonzero')
    values = {p: S.trace(B*B)/2, r: B.det(),
              s0: (v.T*metric*v)[0], s1: (v.T*metric*B*v)[0],
              s2: (v.T*metric*B*B*v)[0]}
    xi = S.Matrix(S.symbols('xi0:3'))
    physical = (E.T*metric).inv()*xi
    transformed = S.Poly(S.expand(cubic.subs(dict(zip(coordinates, physical)), simultaneous=True)), *xi)
    for indices in TRIPLES:
        counts = tuple(indices.count(i) for i in range(3))
        multiplicity = S.Rational(6, math.prod(math.factorial(v) for v in counts))
        values[U[indices]] = transformed.coeff_monomial(counts)/multiplicity
    assert S.simplify(DET.subs(values)+determinant**2) == 0
    assert all(S.simplify(f.subs(values)) == 0 for f in trace_relations())
    return values, determinant


class Circuit:
    """Hash-consed arithmetic DAG over Q in the fifteen quotient coordinates."""
    def __init__(self, node_limit=1200000):
        self.nodes = []
        self.index = {}
        self.node_limit = node_limit
        self.zero = self.constant(0)
        self.one = self.constant(1)
        self.variables = {variable: self.node(('v', str(variable))) for variable in VARIABLES}

    def node(self, data):
        if data in self.index:
            return self.index[data]
        if len(self.nodes) >= self.node_limit:
            raise RuntimeError('declared circuit node limit reached; retain the last prefix checkpoint')
        number = len(self.nodes)
        self.index[data] = number
        self.nodes.append(data)
        return number

    def constant(self, value):
        value = Fraction(value)
        return self.node(('q', value.numerator, value.denominator))

    def add(self, left, right):
        if left == self.zero:
            return right
        if right == self.zero:
            return left
        return self.node(('+', min(left, right), max(left, right)))

    def multiply(self, left, right):
        if left == self.zero or right == self.zero:
            return self.zero
        if left == self.one:
            return right
        if right == self.one:
            return left
        return self.node(('*', min(left, right), max(left, right)))

    def summation(self, items):
        result = self.zero
        for item in items:
            result = self.add(result, item)
        return result

    def from_sympy(self, expr):
        expr = S.sympify(expr)
        if expr.is_Rational:
            return self.constant(Fraction(int(expr.p), int(expr.q)))
        if expr in self.variables:
            return self.variables[expr]
        if expr.is_Add:
            return self.summation(self.from_sympy(arg) for arg in expr.args)
        if expr.is_Mul:
            result = self.one
            for arg in expr.args:
                result = self.multiply(result, self.from_sympy(arg))
            return result
        if expr.is_Pow and expr.exp.is_Integer and expr.exp >= 0:
            result = self.one
            base = self.from_sympy(expr.base)
            for exponent in range(int(expr.exp)):
                result = self.multiply(result, base)
            return result
        raise ValueError('unsupported polynomial expression '+str(expr))

    def evaluate(self, values, outputs):
        result = []
        for node in self.nodes:
            if node[0] == 'q':
                value = Fraction(node[1], node[2])
            elif node[0] == 'v':
                raw = S.Rational(values[S.Symbol(node[1])])
                value = Fraction(int(raw.p), int(raw.q))
            elif node[0] == '+':
                value = result[node[1]]+result[node[2]]
            else:
                value = result[node[1]]*result[node[2]]
            result.append(value)
        return [result[index] for index in outputs]

    def save(self, path, outputs):
        payload = {'schema': 1, 'coefficient_field': 'QQ',
                   'variables': list(map(str, VARIABLES)),
                   'nodes': self.nodes, 'moment_outputs': outputs,
                   'relations': list(map(str, normalized_relations())),
                   'scope': 'Full normalized D=1/8 quotient. Global dimension/length not asserted.'}
        # A fixed gzip timestamp keeps the exact circuit artifact reproducible.
        with Path(path).open('wb') as stream:
            with gzip.GzipFile(fileobj=stream, mode='wb', mtime=0, filename='') as compressed:
                compressed.write(json.dumps(payload, separators=(',', ':')).encode())


def moment_circuit(prefix=12, progress=None, node_limit=1200000):
    graph = Circuit(node_limit)
    phase, xi = phase_polynomial()
    phase_terms = [(powers, graph.from_sympy(coefficient))
                   for powers, coefficient in S.Poly(phase, *xi).terms()]
    covariance = [[graph.from_sympy(G[i, j]) for j in range(3)] for i in range(3)]

    @lru_cache(maxsize=None)
    def wick(powers):
        total = sum(powers)
        if total == 0:
            return graph.one
        if total % 2:
            return graph.zero
        i = next(i for i, exponent in enumerate(powers) if exponent)
        remaining = list(powers)
        remaining[i] -= 1
        terms = []
        for j in range(3):
            if remaining[j]:
                exponent = remaining[j]
                reduced = list(remaining)
                reduced[j] -= 1
                value = graph.multiply(covariance[i][j], wick(tuple(reduced)))
                value = graph.multiply(graph.constant(exponent), value)
                terms.append(value)
        return graph.summation(terms)

    power_terms = {(0, 0, 0): graph.one}
    outputs = []
    for n in range(1, prefix+1):
        updated = {}
        for exponent, coefficient in power_terms.items():
            for powers, factor in phase_terms:
                target = tuple(i+j for i, j in zip(exponent, powers))
                term = graph.multiply(coefficient, factor)
                updated[target] = graph.add(updated.get(target, graph.zero), term)
        power_terms = updated
        summands = []
        for powers, coefficient in power_terms.items():
            if sum(powers) % 2 == 0:
                term = graph.multiply(coefficient, wick(powers))
                term = graph.multiply(graph.constant(Fraction(1, math.prod(range(sum(powers)+1, 0, -2)))), term)
                summands.append(term)
        outputs.append(graph.summation(summands))
        if progress:
            progress(n, len(graph.nodes), len(power_terms))
    return graph, outputs


def hilbert_residue():
    x = S.Symbol('x')
    integrand = x*x/((1+x*x)**3*(1+4*x*x)**2*(1+9*x*x))
    contour = S.simplify(S.I*sum(S.residue(integrand, x, S.I/k) for k in [1, 2, 3])/2)
    complete = [S.Integer(1)]+[S.Integer(0)]*10
    for value in [2, 2, 2, 4, 4, 6]:
        for degree in range(1, 11):
            complete[degree] += value*complete[degree-1]

    def schur(partition):
        return S.det(S.Matrix([[complete[partition[i]-i+j] if partition[i]-i+j >= 0 else 0
                                for j in range(6)] for i in range(6)]))

    published = S.factor(2*schur([3, 3, 3, 2, 1, 0])/schur([5, 4, 3, 2, 1, 0]))
    return contour, published


def quadratic_moments(prefix=6):
    """Exact full H3=0 regression in the five cyclic-frame base invariants."""
    tau = S.Symbol('tau')
    numerator = s0+2*s1*tau+(4*s2-4*p*s0)*tau*tau
    deviation = 4*p*tau*tau+8*r*tau**3
    moments = []
    for n in range(1, prefix+1):
        value = S.Integer(0)
        for k in range(n//2+1):
            target = n-2*k
            denominator_series = sum(S.rf(S.Rational(1, 2)+k, j)*deviation**j/S.factorial(j)
                                     for j in range(target//2+1))
            coefficient = S.expand(numerator**k*denominator_series).coeff(tau, target)
            value += S.factorial(n)*coefficient/(2**k*S.factorial(k)*S.factorial2(2*n-2*k+1))
        moments.append(S.expand(value))
    return tuple(moments)
