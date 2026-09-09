"""Exact rational-span classification from a fixed injective mod-2 footprint.

The unbounded algorithm terminates by Mordell--Weil and height contraction.
The implementation is bounded and returns UNKNOWN at its step cap. No rational
point search, numerical heights, full rank, class group or Selmer oracle.
"""
from pathlib import Path
import runpy
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, prime_range, gcd, lcm

KERNEL = Path(__file__).with_name('verify_det1092_funnel_small_conic_seed.sage')
finite_group = runpy.run_path(str(KERNEL))['finite_group']


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def point_record(P):
    return None if P.is_zero() else list(map(str, P.xy()))


def decode_point(E, row):
    return E(0) if row is None else E(list(map(QQ, row)))


def reduce_point(P, p):
    if P.is_zero():
        return None
    x, y = P.xy()
    z = lcm(x.denominator(), y.denominator())
    coords = [ZZ(x*z), ZZ(y*z), ZZ(z)]
    g = gcd(coords)
    X, Y, Z = [int(v/g % p) for v in coords]
    if Z == 0:
        require(X == 0 and Y != 0, 'invalid projective reduction')
        return None
    return X*pow(Z,-1,p) % p, Y*pow(Z,-1,p) % p


def build_frame(E, basis, prime_cap=1009):
    """Uses only the curve and inherited basis; stops before any candidate read."""
    require(list(E.a_invariants()[:3]) == [0,0,0], 'short model required')
    r = len(basis)
    rows, records, exposures = [], [], []
    torsion_prime = None
    for p0 in prime_range(3, prime_cap+1):
        p = int(p0)
        if any(a.denominator() % p == 0 for a in E.a_invariants()):
            exposures.append(dict(prime=p, status='COEFFICIENT_DENOMINATOR'))
            continue
        if E.discriminant() % p == 0:
            exposures.append(dict(prime=p, status='BAD_REDUCTION'))
            continue
        A, B = [int(a % p) for a in E.a_invariants()[3:]]
        points, doubles, labels, reps, dim = finite_group(p, 0, A, B)
        codes = [labels[reduce_point(P,p)] for P in basis]
        new = [[(v >> j)&1 for v in codes] for j in range(dim)]
        rows.extend(new)
        if len(points) % 2 == 1 and torsion_prime is None:
            torsion_prime = p
        records.append(dict(prime=p, order=len(points), dimension=dim, codes=codes))
        rank = matrix(GF(2), len(rows), r, sum(rows, [])).rank()
        exposures.append(dict(prime=p, status='GOOD', quotient_dimension=dim, inherited_rank=int(rank)))
        if rank == r and torsion_prime is not None:
            break
    M = matrix(GF(2), len(rows), r, sum(rows, []))
    require(M.rank() == r and torsion_prime is not None, 'UNKNOWN_GENERIC_FOOTPRINT_CAP')
    return dict(status='PASS_GENERIC_ONLY_INJECTIVE_FOOTPRINT', curve=list(map(str,E.a_invariants())),
        basis=[point_record(P) for P in basis], inherited_rank=r,
        finite_rows=[list(map(int,row)) for row in M.rows()],
        records=records, exposures=exposures, no_two_torsion_prime=torsion_prime,
        prime_cap=prime_cap, candidate_inputs=0)


class Classifier:
    def __init__(self, frame):
        self.E = E = EllipticCurve(QQ, list(map(QQ,frame['curve'])))
        self.basis = [decode_point(E,p) for p in frame['basis']]
        self.groups = []
        rows = []
        for record in frame['records']:
            p = record['prime']
            require(ZZ(p).is_prime(proof=True) and p > 2 and E.discriminant() % p != 0, 'proof prime')
            A, B = [int(a % p) for a in E.a_invariants()[3:]]
            points, _, labels, _, dim = finite_group(p,0,A,B)
            require(len(points) == record['order'] and dim == record['dimension'], 'finite group metadata')
            codes = [labels[reduce_point(P,p)] for P in self.basis]
            require(codes == record['codes'], 'inherited reduction changed')
            rows.extend([[(v>>j)&1 for v in codes] for j in range(dim)])
            self.groups.append((p,labels,dim))
            if p == frame['no_two_torsion_prime']:
                require(len(points)%2 == 1, 'odd-order reduction required')
        require(frame['no_two_torsion_prime'] in [r['prime'] for r in frame['records']], 'missing torsion proof')
        self.M = matrix(GF(2),rows)
        require(self.M.rank() == len(self.basis), 'inherited mod2 independence')
        require(rows == frame['finite_rows'], 'frame rows')

    def code(self,P):
        out = []
        for p,labels,dim in self.groups:
            c = labels[reduce_point(P,p)]
            out.extend((c>>j)&1 for j in range(dim))
        return vector(GF(2),out)

    def word_point(self,w):
        return sum((ZZ(n)*P for n,P in zip(w,self.basis)), self.E(0))

    def halves(self,P):
        """Complete rational halving by the explicit degree-four duplication equation."""
        E = self.E
        if P.is_zero():
            return [E(0)], dict(kind='ZERO_TARGET_NO_RATIONAL_TWO_TORSION')
        R = PolynomialRing(QQ,'x'); x = R.gen()
        A,B = E.a4(), E.a6(); a = P[0]
        f = x**3+A*x+B
        poly = (3*x*x+A)**2-4*f*(a+2*x)
        require(poly == x**4-4*a*x**3-2*A*x*x-(8*B+4*A*a)*x+A*A-4*B*a, 'duplication polynomial')
        fac = poly.factor()
        require(fac.prod() == poly, 'factorization product')
        candidates, result = [], []
        for g,e in fac:
            if g.degree() != 1:
                continue
            xx = -g[0]/g[1]; yy2 = f(xx)
            candidate = dict(x=str(xx), y_squared=str(yy2), square=bool(yy2.is_square()))
            if yy2.is_square():
                for yy in sorted(set([yy2.sqrt(),-yy2.sqrt()])):
                    Q = E([xx,yy])
                    if 2*Q == P:
                        result.append(Q)
            candidates.append(candidate)
        result = sorted(set(result), key=lambda Q:tuple(Q))
        require(len(result) <= 1, 'frame excludes rational 2-torsion')
        return result, dict(kind='EXACT_DUPLICATION_QUARTIC', polynomial=list(map(str,poly.list())),
            unit=str(fac.unit()), factors=[dict(coefficients=list(map(str,g.list())), exponent=int(e)) for g,e in fac],
            rational_abscissas=candidates, rational_halves=[point_record(Q) for Q in result])

    def classify(self, P, max_steps=8):
        E=self.E; original=P; r=len(self.basis)
        seen={}; history=[]; W=vector(ZZ,[0]*r)
        for n in range(max_steps+1):
            require(original == 2**n*P+self.word_point(W), 'chain identity')
            key=tuple(P)
            if key in seen:
                a,Wa=seen[key]; factor=ZZ(2**(n-a)-1)
                relation=2**(n-a)*Wa-W
                require(factor*original == self.word_point(relation), 'cycle dependence relation')
                return dict(status='INHERITED_RATIONAL_SPAN', reason='EXACT_HALVING_CYCLE', steps=n,
                    original=point_record(original), history=history, cycle_start=a,
                    relation_multiplier=str(factor), relation_word=list(map(str,relation)))
            if P.is_zero():
                return dict(status='INHERITED_RATIONAL_SPAN', reason='ZERO_REACHED', steps=n,
                    original=point_record(original), history=history,
                    relation_multiplier='1', relation_word=list(map(str,W)))
            seen[key]=(n,vector(ZZ,W))
            c=self.code(P)
            try:
                bits=self.M.solve_right(c)
            except ValueError:
                return dict(status='NEW_INDEPENDENT_DIRECTION', reason='FINITE_FOOTPRINT_ESCAPE', steps=n,
                    original=point_record(original), terminal=point_record(P), history=history,
                    terminal_code=list(map(int,c)), chain_word=list(map(str,W)))
            require(self.M*bits == c, 'finite parity solution')
            if n == max_steps:
                return dict(status='UNKNOWN_STEP_CAP', steps=n, original=point_record(original), history=history)
            word=list(map(int,bits)); target=P-self.word_point(word)
            halves,certificate=self.halves(target)
            record=dict(step=n, point=point_record(P), code=list(map(int,c)), parity_word=word,
                        target=point_record(target), halving=certificate)
            history.append(record)
            if not halves:
                return dict(status='NEW_INDEPENDENT_DIRECTION', reason='GLOBAL_NONHALVING_ESCAPE', steps=n,
                    original=point_record(original), terminal=point_record(P), history=history,
                    chain_word=list(map(str,W)))
            record['next']=point_record(halves[0])
            W += 2**n*vector(ZZ,word)
            P=halves[0]
        raise AssertionError('unreachable')
