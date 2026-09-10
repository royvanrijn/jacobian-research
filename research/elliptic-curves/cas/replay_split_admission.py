"""Independent admission replay: Fraction group law and Sage finite groups.

No producer imports, artifact reads, point search, or new finite places.
This checks supplied witnesses, including exact complete quartic factorizations
for nonhalving. Shared CAS polynomial arithmetic is used only over Q.
"""
from fractions import Fraction as F
from sage.all import QQ, ZZ, GF, EllipticCurve, PolynomialRing, matrix, vector


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def point(row):
    return None if row is None else tuple(map(F, row))


def add(A, P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P; u, v = Q
    if x == u and y == -v:
        return None
    slope = (3*x*x+A)/(2*y) if P == Q else (v-y)/(u-x)
    z = slope*slope-x-u
    return z, slope*(x-z)-y


def multiply(A, n, P):
    n = int(n)
    if n < 0:
        P = None if P is None else (P[0], -P[1])
        n = -n
    result = None
    for bit in bin(n)[2:]:
        result = add(A, result, result)
        if bit == '1':
            result = add(A, result, P)
    return result


def word_point(A, basis, word):
    require(len(word) == len(basis), 'word length')
    result = None
    for n, P in reversed(list(zip(word, basis))):
        result = add(A, result, multiply(A, int(n), P))
    return result


class ReplayFrame:
    def __init__(self, frame):
        self.A, self.B = map(F, frame['curve'][3:])
        require(list(map(F, frame['curve'][:3])) == [0, 0, 0], 'short model')
        self.basis = [self.checked_point(row) for row in frame['basis']]
        self.groups = []
        primes = [row['prime'] for row in frame['records']]
        require(len(primes) == len(set(primes)), 'duplicate place')
        torsion = False
        for record in frame['records']:
            p = record['prime']
            require(p > 2 and ZZ(p).is_prime(proof=True), 'prime')
            require(self.A.denominator % p and self.B.denominator % p, 'integral coefficients')
            ep = EllipticCurve(GF(p), [QQ(self.A), QQ(self.B)])
            require(ep.discriminant(), 'good reduction')
            # Preserve the sealed affine order, while using an independent
            # group implementation to enumerate doubling and every coset.
            points = [ep(0)] + sorted((P for P in ep.points() if not P.is_zero()),
                                     key=lambda P: (int(P[0]), int(P[1])))
            key = lambda P: tuple(map(int, P))
            doubles = {key(2*P): 2*P for P in points}
            labels = {k: 0 for k in doubles}
            reps = [ep(0)]
            dim = 0
            for P in points:
                if key(P) in labels:
                    continue
                new = []
                for i, Q in enumerate(reps):
                    V = P+Q
                    new.append(V)
                    for D in doubles.values():
                        k = key(V+D)
                        require(k not in labels, 'disjoint quotient cosets')
                        labels[k] = i | (1 << dim)
                reps += new
                dim += 1
            require(len(labels) == len(points) == len(doubles)*2**dim, 'complete finite quotient')
            require(len(points) == record['order'] and dim == record['dimension'], 'sealed quotient metadata')
            if p == frame['no_two_torsion_prime']:
                require(len(points) % 2 == 1, 'odd-order reduction')
                torsion = True
            self.groups.append((p, ep, labels, dim))
        require(torsion, 'no rational two-torsion witness')
        codes = [self.code(P) for P in self.basis]
        self.M = matrix(GF(2), codes).transpose()
        require(self.M.rank() == len(self.basis), 'inherited injection')
        require([list(map(int, row)) for row in self.M.rows()] == frame['finite_rows'], 'sealed finite rows')
        offset = 0
        for record in frame['records']:
            dim = record['dimension']
            require(record['codes'] == [sum(int(c[offset+j]) << j for j in range(dim)) for c in codes],
                    'sealed local codes')
            offset += dim
        self.frame = frame

    def checked_point(self, row):
        P = point(row)
        require(P is None or P[1]**2 == P[0]**3+self.A*P[0]+self.B, 'rational curve membership')
        return P

    def code(self, P):
        out = []
        for p, ep, labels, dim in self.groups:
            if P is None or any(v.denominator % p == 0 for v in P):
                red = ep(0)
            else:
                red = ep([v.numerator*pow(v.denominator, -1, p) % p for v in P])
            label = labels[tuple(map(int, red))]
            out.extend((label >> j) & 1 for j in range(dim))
        return vector(GF(2), out)

    def word(self, w):
        return word_point(self.A, self.basis, w)

    def quartic(self, target, certificate, *, nonhalving=False):
        if target is None:
            require(not nonhalving and certificate['kind'] == 'ZERO_TARGET_NO_RATIONAL_TWO_TORSION',
                    'zero halving target')
            return
        R = PolynomialRing(QQ, 'x'); x = R.gen()
        A, B, a = map(QQ, (self.A, self.B, target[0]))
        polynomial = x**4-4*a*x**3-2*A*x*x-(8*B+4*A*a)*x+A*A-4*B*a
        require(certificate['kind'] == 'EXACT_DUPLICATION_QUARTIC'
                and R(certificate['polynomial']) == polynomial, 'duplication quartic')
        product = R(certificate['unit'])
        for factor in certificate['factors']:
            g = R(factor['coefficients']); e = factor['exponent']
            require(e > 0 and g.degree() > 0, 'factor degree and multiplicity')
            product *= g**e
            if nonhalving:
                require(g.is_irreducible(), 'complete rational factorization')
                if g.degree() == 1:
                    xx = -g[0]/g[1]; yy2 = xx**3+A*xx+B
                    if yy2.is_square():
                        for yy in {yy2.sqrt(), -yy2.sqrt()}:
                            P = (F(str(xx)), F(str(yy)))
                            require(multiply(self.A, 2, P) != target, 'rational half exists')
        require(product == polynomial, 'factorization product')

    def replay(self, result):
        original = self.checked_point(result['original']); current = original
        r = len(self.basis); W = [0]*r
        states, words = [current], [W[:]]
        for n, h in enumerate(result['history']):
            require(n == h['step'] and current == self.checked_point(h['point']), 'chain step')
            require(original == add(self.A, multiply(self.A, 2**n, current), self.word(W)), 'chain invariant')
            bits = h['parity_word']
            require(len(bits) == r and all(type(b) is int and b in (0, 1) for b in bits), 'binary parity')
            c = self.code(current)
            require(list(map(int, c)) == h['code'] and self.M*vector(GF(2), bits) == c, 'unique compatible parity')
            T = self.word(bits)
            target = add(self.A, current, multiply(self.A, -1, T))
            require(target == self.checked_point(h['target']), 'exact subtraction')
            self.quartic(target, h['halving'], nonhalving='next' not in h)
            if 'next' not in h:
                require(n == len(result['history'])-1 and result['reason'] == 'GLOBAL_NONHALVING_ESCAPE',
                        'nonhalving terminal only')
                break
            current = self.checked_point(h['next'])
            require(multiply(self.A, 2, current) == target, 'exact rational doubling')
            W = [w+2**n*b for w, b in zip(W, bits)]
            states.append(current); words.append(W[:])
        steps = result['steps']
        require(steps == len(states)-1 and steps <= result['max_successful_halves'], 'successful halving count')
        require(len(set(states[:-1])) == len(states[:-1]) and None not in states[:-1], 'first terminal reached')
        require(original == add(self.A, multiply(self.A, 2**steps, current), self.word(W)), 'final chain invariant')
        status = result['status']
        if status == 'INHERITED_RATIONAL_SPAN':
            if result['reason'] == 'EXACT_HALVING_CYCLE':
                a = result['cycle_start']
                require(0 <= a < steps and states[a] == current, 'exact repeated point')
                multiplier = 2**(steps-a)-1
                relation = [2**(steps-a)*u-v for u, v in zip(words[a], W)]
            else:
                require(result['reason'] == 'ZERO_REACHED' and current is None, 'zero terminal')
                multiplier, relation = 1, W
            require(str(multiplier) == result['relation_multiplier']
                    and list(map(str, relation)) == result['relation_word'], 'cycle-derived relation')
            require(multiply(self.A, multiplier, original) == self.word(relation), 'independent exact relation')
        elif status == 'NEW_INDEPENDENT_DIRECTION':
            require(current == self.checked_point(result['terminal'])
                    and W == list(map(int, result['chain_word'])), 'independence transport')
            if result['reason'] == 'FINITE_FOOTPRINT_ESCAPE':
                c = self.code(current)
                require(list(map(int, c)) == result['terminal_code'], 'terminal finite column')
                require(self.M.augment(matrix(GF(2), len(c), 1, list(c))).rank() == r+1, 'finite escape')
            else:
                require(result['reason'] == 'GLOBAL_NONHALVING_ESCAPE'
                        and 'next' not in result['history'][-1], 'global nonhalving escape')
        else:
            require(status == 'UNKNOWN' and result['reason'] == 'HALVING_STEP_CAP'
                    and steps == result['max_successful_halves'], 'UNKNOWN proof bound')
            require(current is not None and current not in states[:-1], 'no closed terminal at cap')
            pending = result['pending']; bits = pending['parity_word']
            require(pending['point'] == result['original'] if steps == 0 else
                    self.checked_point(pending['point']) == current, 'pending point')
            require(len(bits) == r and set(bits) <= {0, 1}
                    and list(map(int, self.code(current))) == pending['code']
                    and self.M*vector(GF(2), bits) == self.code(current), 'pending unique parity')
            require(self.checked_point(pending['target']) == add(self.A, current, multiply(self.A, -1, self.word(bits))),
                    'pending subtraction')
        return dict(status=status, steps=steps, finite_rank=r, finite_places=len(self.groups),
                    finite_rows=self.M.nrows(), exact_halves_checked=steps,
                    reason=result['reason'])
