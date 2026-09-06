"""Small exact fixtures; full class groups are used only as calibration oracles.

The production verifier does not call bnfinit, bnfisprincipal or bnfcertify.
All oracle output is converted into replayable principal-element witnesses.
"""
from collections import defaultdict
from sage.all import pari, prime_range
from class_span_grh import FieldAudit, SCHEMA, digest

CASES = [
    ('imaginary_trivial', [1,0,1]),
    ('imaginary_c2', [5,0,1]),
    ('imaginary_c2_squared', [21,0,1]),
    ('real_quadratic', [-5,0,1]),
    ('mixed_cubic', [-1,-1,0,1]),
    ('real_cubic', [1,-3,0,1]),
    ('complex_quartic', [1,0,0,0,1]),
]


def fixture(coefficients, cutoff=100):
    audit = FieldAudit(coefficients)
    bnf = pari.bnfinit(audit.nf, 1)
    if int(pari.bnfcertify(bnf)) != 1:
        raise ArithmeticError('small-field oracle not certified')
    cyclic = [int(x) for x in bnf.bnf_get_cyc()]
    generators = list(bnf.bnf_get_gen())
    columns, lookup = [], {}

    def insert(P):
        column = audit.column(P)
        key = digest(column)
        if key not in lookup:
            lookup[key] = len(columns)
            columns.append(column)
        return lookup[key]

    def factor(I):
        factors = pari.idealfactor(audit.nf, I)
        return [[insert(factors[j,0]), int(factors[j,1])] for j in range(factors.nrows())]

    primes = []
    for p in prime_range(cutoff):
        for P in audit.primes_above(int(p)):
            insert(P)
            if int(p)**int(P[3]) < cutoff:
                primes.append(P)
    generator_factors = [factor(I) for I in generators]
    relations = []
    for P in primes:
        exponents, coordinates = pari.bnfisprincipal(bnf, P, 1)
        beta = pari.nfbasistoalg(audit.nf, coordinates)
        factors = defaultdict(int)
        factors[insert(P)] = 1
        for exponent, group in zip(exponents, generator_factors):
            for c, e in group:
                factors[c] -= int(exponent)*e
        relations.append({'element':[str(pari.polcoef(pari.lift(beta), i)) for i in range(audit.degree)],
                          'factorization':[[c,e] for c,e in sorted(factors.items()) if e]})
    # Include every cyclic generator, also odd-order ones. The upper bound need
    # not be sharp; no knowledge of their orders enters the production checker.
    document = {'schema':SCHEMA, 'field':audit.description, 'cutoff':cutoff,
                'precision_bits':192, 'columns':columns, 'relations':relations,
                'anchors':generator_factors}
    return document, {'certified_cyclic_invariants':cyclic,
                      'unconditional_class_two_rank':sum(n%2 == 0 for n in cyclic)}
