"""Exact relation algebra and a GRH sufficient test for Cl(K)/2 generation.

The low-level RelationSpan and quadratic_margin functions are mathematical
components: their caller must certify the field, primes and principal rows.
Only verify_document performs that audit and emits a class-rank upper bound.
No elliptic-curve rank or class-group independence is inferred here.
"""
from hashlib import sha256
from fractions import Fraction
import json
from sage.all import QQ, ZZ, PolynomialRing, RealIntervalField, RealBallField, pari, prime_range

SCHEMA = 'number-fields.class-span-input.v1'
ASSUMPTION = ('GRH for nontrivial quadratic characters of the ordinary ideal class '
              'group of the specified field that are trivial on the proposed span.')


def digest(data):
    return sha256(json.dumps(data, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def integer(x):
    if isinstance(x, bool) or not isinstance(x, (int, str, type(ZZ(0)))):
        raise ValueError('exact integer required')
    value = int(x)
    if str(value) != str(x):
        raise ValueError('canonical integer required')
    return value


def rational(x):
    if not isinstance(x, (int, str)) or isinstance(x, bool):
        raise ValueError('exact rational string or integer required')
    if isinstance(x, str) and not all(c in '0123456789-/' for c in x):
        raise ValueError('rational literal required')
    return QQ(x)


def interval_record(x):
    return {'lower':str(x.lower()), 'upper':str(x.upper())}


def sparse_row(factors, width):
    row, used = 0, set()
    for c, e in factors:
        c, e = integer(c), integer(e)
        if not 0 <= c < width or c in used or e == 0:
            raise ValueError('duplicate, zero or out-of-range factor')
        used.add(c)
        if e % 2:
            row ^= 1 << c
    return row


class RelationSpan:
    """Incremental F2 algebra. All coordinates, including outside primes, persist.

    Rows supplied to this class are not automatically arithmetic certificates.
    A membership includes coordinates in the supplied anchor list. The returned
    dimension is an upper bound for their image modulo actual principal ideals;
    independence in an incomplete presentation is not actual class independence.
    """
    def __init__(self, width):
        self.width = integer(width)
        if self.width < 0:
            raise ValueError('negative width')
        self.basis = {}
        self.rows = []

    def validate(self, row):
        row = integer(row)
        if row < 0 or row.bit_length() > self.width:
            raise ValueError('row exceeds declared prime coordinates')
        return row

    @staticmethod
    def reduce(row, basis, tag=0):
        while row:
            pivot = row.bit_length()-1
            if pivot not in basis:
                break
            other, label = basis[pivot]
            row ^= other
            tag ^= label
        return row, tag

    def add_relation(self, row):
        row = self.validate(row)
        self.rows.append(row)
        reduced, _ = self.reduce(row, self.basis)
        if reduced:
            self.basis[reduced.bit_length()-1] = (reduced, 0)
        return bool(reduced)

    def analyze(self, anchors):
        anchors = [self.validate(row) for row in anchors]
        augmented = dict(self.basis)
        for j, row in enumerate(anchors):
            reduced, tag = self.reduce(row, augmented, 1 << j)
            if reduced:
                augmented[reduced.bit_length()-1] = (reduced, tag)
        known = {}
        for c in range(self.width):
            residual, coordinates = self.reduce(1 << c, augmented)
            if residual:
                continue
            represented = 1 << c
            for j, row in enumerate(anchors):
                if (coordinates >> j) & 1:
                    represented ^= row
            if self.reduce(represented, self.basis)[0]:
                raise ArithmeticError('anchor membership witness does not replay')
            known[c] = coordinates
        return {'known_coordinates':known, 'relation_rank':len(self.basis),
                'formal_quotient_dimension':self.width-len(self.basis),
                'anchor_image_dimension_upper_bound':len(augmented)-len(self.basis),
                'anchor_count':len(anchors)}


def quadratic_margin(discriminant, signature, prime_norms, known_columns, cutoff, bits=192):
    """Interval sufficient test, conditional on caller-audited arithmetic inputs.

    Unknown signs are -1 at odd powers and +1 at even powers. The archimedean
    upper bound omits a negative tail in BDyDF2008 equation(13). This supports
    every number-field signature, not just totally real cubics.
    """
    d, T, bits = integer(discriminant), integer(cutoff), integer(bits)
    r1, r2 = map(integer, signature)
    n = r1+2*r2
    if d == 0 or T <= 1 or bits < 64 or r1 < 0 or r2 < 0 or n < 1:
        raise ValueError('invalid field or interval parameters')
    if (d < 0) != bool(r2 % 2):
        raise ValueError('discriminant sign disagrees with signature')
    norms = [integer(q) for q in prime_norms]
    if any(q < 2 for q in norms):
        raise ValueError('prime ideal norms must be at least2')
    known = {integer(c) for c in known_columns}
    if any(c < 0 or c >= len(norms) for c in known):
        raise ValueError('known prime coordinate out of range')
    R = RealIntervalField(bits)
    L, total, signed, penalty = R(T).log(), R(0), R(0), R(0)
    terms, unknown_odd, contributing, unresolved = 0, 0, 0, []
    for c, q in enumerate(norms):
        if q >= T:
            continue
        contributing += 1
        power, exponent, weight_sum = q, 1, R(0)
        logq = R(q).log()
        while power < T:
            weight = logq*(1-exponent*logq/L)/R(power).sqrt()
            total += weight
            if c not in known and exponent % 2:
                signed -= weight
                penalty += 4*weight
                weight_sum += 4*weight
                unknown_odd += 1
            else:
                signed += weight
            power *= q
            exponent += 1
            terms += 1
        if c not in known:
            unresolved.append({'column':c, 'norm':q, 'penalty_interval':interval_record(weight_sum)})
    rhs = (R(abs(d)).log()-n*(R.euler_constant()+(8*R.pi()).log())-r1*R.pi()/2
           +(n*R.pi()**2/2+4*r1*R(RealBallField(bits).catalan_constant()))/L)
    margin = 2*signed-rhs
    if not margin.overlaps(2*total-rhs-penalty):
        raise ArithmeticError('signed-sum and penalty evaluations disagree')
    # This ranking is a scheduling aid; only the full interval decides the gate.
    unresolved.sort(key=lambda x:(-Fraction(x['penalty_interval']['lower']), x['column']))
    return {'cutoff':T, 'precision_bits':bits, 'prime_ideals_contributing':contributing,
            'prime_power_terms':terms, 'unknown_odd_prime_power_terms':unknown_odd,
            'all_positive_sum_interval':interval_record(total),
            'worst_signed_sum_interval':interval_record(signed),
            'archimedean_rhs_upper_bound_interval':interval_record(rhs),
            'penalty_interval':interval_record(penalty), 'margin_interval':interval_record(margin),
            'positive_margin_certified':bool(margin.lower()>0),
            'unresolved_primes_by_penalty':unresolved,
            'claim':'Numerical sufficient test only; field, prime coverage and memberships require arithmetic audit.'}


class FieldAudit:
    """Certified maximal order and exact prime/element conversion; no bnfinit."""
    def __init__(self, coefficients, prime_hints=()):
        coefficients = [integer(c) for c in coefficients]
        if len(coefficients) < 3 or coefficients[-1] != 1:
            raise ValueError('monic integral polynomial of degree>=2 required')
        ring = PolynomialRing(QQ, 'x')
        polynomial = ring(coefficients)
        if not polynomial.is_irreducible():
            raise ValueError('irreducible number-field polynomial required')
        hints = [integer(q) for q in prime_hints]
        if len(set(hints)) != len(hints) or any(q < 2 or not pari.isprime(q) for q in hints):
            raise ValueError('distinct proved prime hints required')
        self.polynomial = pari(polynomial)
        self.nf = pari.nfinit([self.polynomial, hints]) if hints else pari.nfinit(self.polynomial)
        if len(pari.nfcertify(self.nf)):
            raise ArithmeticError('maximal order is not certified')
        self.degree = len(coefficients)-1
        self.theta = pari.Mod('x', self.polynomial)
        self.decompositions = {}
        self.description = {
            'polynomial':list(map(str, coefficients)), 'prime_hints':hints,
            'discriminant':str(self.nf.disc()), 'signature':[int(x) for x in self.nf[1]],
            'integral_basis':[[str(pari.polcoef(b, i)) for i in range(self.degree)] for b in self.nf.nf_get_zk()]}

    def element(self, coefficients):
        if len(coefficients) != self.degree:
            raise ValueError('element must have degree-many power-basis coordinates')
        value = sum(pari(rational(c))*self.theta**i for i, c in enumerate(coefficients))
        if value == 0:
            raise ValueError('principal generator must be nonzero')
        return value

    def column(self, P):
        H = pari.idealhnf(self.nf, P)
        return {'p':int(P[0]), 'e':int(P[2]), 'f':int(P[3]),
                'hnf':[[int(H[i,j]) for j in range(self.degree)] for i in range(self.degree)]}

    def primes_above(self, p):
        p = integer(p)
        if p not in self.decompositions:
            if p < 2 or not pari.isprime(p):
                raise ValueError('rational prime is not proved')
            primes = list(pari.idealprimedec(self.nf, p))
            product = pari.idealhnf(self.nf, 1)
            if sum(int(P[2])*int(P[3]) for P in primes) != self.degree:
                raise ArithmeticError('prime decomposition has wrong degree')
            for P in primes:
                if pari.idealnorm(self.nf, P) != p**int(P[3]):
                    raise ArithmeticError('prime ideal norm differs')
                product = pari.idealmul(self.nf, product, pari.idealpow(self.nf, P, int(P[2])))
            if product != pari.idealhnf(self.nf, p):
                raise ArithmeticError('prime decomposition does not factor(p)')
            self.decompositions[p] = primes
        return self.decompositions[p]


def verify_document(document):
    """Replay a portable exact principal-relation ledger and its generation test.

    Invalid input raises; a valid ledger with a nonpositive/uncertain margin
    returns UNKNOWN and no upper bound. Missing prime coverage is invalid.
    """
    if document['schema'] != SCHEMA:
        raise ValueError('unsupported class-span schema')
    T = integer(document['cutoff'])
    if T <= 1:
        raise ValueError('cutoff must exceed1')
    field = document['field']
    audit = FieldAudit(field['polynomial'], field['prime_hints'])
    if audit.description != field:
        raise ValueError('field discriminant, signature or integral basis differs')
    columns = document['columns']
    lookup, ideals = {}, []
    for c, col in enumerate(columns):
        p = integer(col['p'])
        integer(col['e']); integer(col['f'])
        if len(col['hnf']) != audit.degree or any(len(row) != audit.degree for row in col['hnf']):
            raise ValueError('prime HNF has wrong dimensions')
        for row in col['hnf']:
            for entry in row:
                integer(entry)
        matches = [P for P in audit.primes_above(p) if audit.column(P) == col]
        if len(matches) != 1:
            raise ValueError('column is not the claimed prime ideal in the certified basis')
        key = digest(col)
        if key in lookup:
            raise ValueError('duplicate prime ideal column')
        lookup[key] = c
        ideals.append(matches[0])
    for p in prime_range(T):
        for P in audit.primes_above(int(p)):
            if int(p)**int(P[3]) < T and digest(audit.column(P)) not in lookup:
                raise ValueError('missing prime ideal below the analytic cutoff')
    span = RelationSpan(len(columns))
    canonical_count = 0
    # Add a rational-prime row only when all its prime factors are represented.
    for p in sorted({integer(c['p']) for c in columns}):
        group = audit.primes_above(p)
        if all(digest(audit.column(P)) in lookup for P in group):
            factors = [[lookup[digest(audit.column(P))], int(P[2])] for P in group]
            span.add_relation(sparse_row(factors, len(columns)))
            canonical_count += 1
    for relation in document['relations']:
        factors = relation['factorization']
        row = sparse_row(factors, len(columns))
        beta = audit.element(relation['element'])
        product = pari.idealhnf(audit.nf, 1)
        for c, e in factors:
            product = pari.idealmul(audit.nf, product, pari.idealpow(audit.nf, ideals[integer(c)], integer(e)))
        if product != pari.idealhnf(audit.nf, beta):
            raise ArithmeticError('principal ideal factorization is incomplete or incorrect')
        span.add_relation(row)
    anchors = [sparse_row(f, len(columns)) for f in document['anchors']]
    algebra = span.analyze(anchors)
    norms = [c['p']**c['f'] for c in columns]
    analytical = quadratic_margin(integer(field['discriminant']), field['signature'], norms,
                                  algebra['known_coordinates'], T, document['precision_bits'])
    success = analytical['positive_margin_certified']
    return {'schema':'number-fields.class-span-certificate.v1',
            'status':'CERTIFIED_UNDER_GRH' if success else 'UNKNOWN',
            'input_sha256':digest(document), 'field_sha256':digest(field),
            'maximal_order_certified':True, 'complete_prime_coverage_below':T,
            'principal_relations_replayed':len(document['relations']),
            'canonical_rational_relations':canonical_count,
            'relation_rank':algebra['relation_rank'],
            'formal_quotient_dimension':algebra['formal_quotient_dimension'],
            'anchor_count':algebra['anchor_count'],
            'anchor_image_dimension_upper_bound':algebra['anchor_image_dimension_upper_bound'],
            'known_prime_coordinates':[[c, hex(bits)] for c,bits in sorted(algebra['known_coordinates'].items())],
            'analytic_test':analytical,
            'class_two_rank_upper_bound_under_grh':algebra['anchor_image_dimension_upper_bound'] if success else None,
            'class_two_rank_lower_bound':None, 'exact_class_two_rank':None,
            'elliptic_curve_rank':None, 'assumption':ASSUMPTION,
            'claim_boundary':'A conditional class-2-rank upper bound only. Positive margin certifies generation; nonpositive margin is inconclusive. Formal independence never establishes actual class independence. Elliptic rank needs separately certified points and Selmer bounds.'}
