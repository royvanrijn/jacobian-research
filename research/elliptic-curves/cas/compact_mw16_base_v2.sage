#!/usr/bin/env sage-python
"""Bounded coordinate compactification of five anonymous A1/MW16 families.

Only the sanitized equation template selects one presentation per fibration.
Auxiliary hyperelliptic reduction proposes maps; exact elliptic identities
decide validity. No parameters, exceptional points, or rank labels are read.
"""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import PolynomialRing, QQ, ZZ, pari, Matrix, gcd, lcm, prime_range

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / 'elliptic-curves/cas'
sys.path.insert(0, str(CAS))
import certify_compact_r17_candidates as cert
import audit_r17_constant_scaling as scaling
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture, Limits
helpers = SourceFileLoader('mw16_base_helpers', str(CAS/'reduce_r17_family_base.sage')).load_module()
TEMPLATE = ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json'
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-mw16-base-v2'


def sources():
    paths = (Path(__file__).resolve(), TEMPLATE, Path(scaling.__file__).resolve(),
             CAS/'reduce_r17_family_base.sage', Path(cert.__file__).resolve(),
             CAS/'research_runtime/store.py', CAS/'research_runtime/supervisor.py')
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}


def representatives():
    template = cert.read(TEMPLATE)
    if template['status'] != 'PASS_TARGET_FREE_A1_MW16_FAMILY_PRESENTATIONS':
        raise ArithmeticError('template status changed')
    selected = {}
    for row in template['presentations']:
        selected.setdefault(row['fibration_id'], row)
    if len(selected) != 5:
        raise ArithmeticError('expected five distinct fibration classes')
    return list(selected.values())


def prepare(directory):
    if (directory/'protocol.json').exists():
        raise FileExistsError('preserve frozen protocol')
    checkpoint(directory/'protocol.json', {
        'schema': 'elliptic-curves.compact-mw16-base.v2', 'sources': sources(),
        'presentations': [r['presentation_id'] for r in representatives()],
        'selection': 'First template presentation of each anonymous fibration, in template order.',
        'prime_trial_bound': 10000, 'maximum_power_exponent': 132,
        'composite_factor_bit_gate': 192, 'prime_root_bit_gate': 512,
        'precision_bits': 8192, 'pari_stack_bytes': 256000000,
        'worker_wall_seconds': 120, 'worker_rss_bytes': 1073741824,
        'maximum_workers': 5, 'maximum_concurrent_workers': 1,
        'gate': 'First representative must improve coefficient bits by at least 25 percent before the other four run.',
        'continuation': 'The v1 first-family 673-bit composite failed the 160-bit factor gate. Exact gcds with sanitized template rational coefficients split it into factors of 17, 95 and 188 bits (with multiplicity). V2 permits coefficient-gcd splitting, recursive perfect-power extraction, and proved factorization of pieces through 192 bits within the unchanged worker caps.',
        'method': 'Coefficient-gcd splitting and recursive perfect-power extraction; weighted constant scaling; exact perfect-power extraction from gcd of primitive A/B discriminants and resultant; bounded proved factorization; prime-local auxiliary minimization and reduction; exact weighted coefficient identities.',
        'scope': 'Five bounded coordinate computations, no parameter or point search. Full section transport is a separate prerequisite for usable search export.'})


def template_integers(value):
    if isinstance(value, dict):
        for item in value.values(): yield from template_integers(item)
    elif isinstance(value, list):
        for item in value: yield from template_integers(item)
    else:
        try: q = QQ(value)
        except (TypeError, ValueError): return
        yield q.numerator()
        yield q.denominator()


def bounded_factor(root, protocol, template):
    pieces = [root]
    for value in template_integers(template):
        for part in list(pieces):
            g = gcd(part, value)
            if 1 < g < part:
                pieces.remove(part); pieces.extend([g, part//g])
    records = []; factors = {}
    def add(part, multiplicity=1):
        if part == 1: return
        for k in range(protocol['maximum_power_exponent'], 1, -1):
            candidate, exact = part.nth_root(k, truncate_mode=True)
            if exact:
                add(candidate, multiplicity*k); return
        if part.nbits() <= protocol['prime_root_bit_gate'] and part.is_prime(proof=True):
            ff = [(part, 1)]
        elif part.nbits() <= protocol['composite_factor_bit_gate']:
            ff = list(part.factor(proof=True))
        else:
            raise ArithmeticError('split composite exceeds frozen factor gate')
        records.append({'piece': str(part), 'bits': int(part.nbits()), 'multiplicity': multiplicity,
                        'factors': [[str(p), int(e)] for p,e in ff]})
        for p,e in ff: factors[p] = factors.get(p,0) + e*multiplicity
    for piece in pieces: add(piece)
    return sorted(factors.items()), records


def run(directory, presentation):
    protocol = cert.read(directory/'protocol.json')
    if protocol['sources'] != sources() or presentation not in protocol['presentations']:
        raise ArithmeticError('frozen binding changed')
    output = directory/(presentation+'.json')
    if output.exists():
        raise FileExistsError('preserve previous coordinate attempt')
    row = next(r for r in representatives() if r['presentation_id'] == presentation)
    R = PolynomialRing(QQ, 't')
    original_A, original_B = (R(list(map(QQ, row['pencil'][k]))) for k in
                              ('A_coefficients_low_to_high', 'B_coefficients_low_to_high'))
    values = lambda f: [cert.F(str(x)) for x in f.list()]
    initial_u, _ = scaling.scale_for(values(original_A), values(original_B))
    A, B = original_A/QQ(str(initial_u))**4, original_B/QQ(str(initial_u))**6
    def primitive(f):
        g = f*lcm(q.denominator() for q in f)
        return g/gcd(ZZ(q) for q in g)
    P, Q = primitive(A), primitive(B)
    invariants = [ZZ(P.discriminant()), ZZ(Q.discriminant()), ZZ(P.resultant(Q))]
    if any(x == 0 for x in invariants):
        raise ArithmeticError('degenerate auxiliary invariant')
    common = abs(gcd(invariants)); remainder = common; small = []
    for p in prime_range(protocol['prime_trial_bound']):
        e = remainder.valuation(p)
        if e:
            small.append([int(p), int(e)]); remainder //= p**e
    exponent, root = 1, remainder
    if remainder != 1:
        for k in range(2, protocol['maximum_power_exponent']+1):
            candidate, exact = remainder.nth_root(k, truncate_mode=True)
            if exact:
                exponent, root = k, candidate
    diagnostic = {'common_invariant_gcd': str(common), 'small_prime_factors': small,
                  'remaining_cofactor': str(remainder), 'largest_power_exponent': exponent,
                  'power_root': str(root), 'power_root_bits': int(root.nbits())}
    base = {'presentation_id': presentation, 'fibration_id': row['fibration_id'],
            'protocol_sha256': cert.hashed(directory/'protocol.json'), 'diagnostic': diagnostic}
    checkpoint(output, {**base, 'status': 'RUNNING_BAD_PRIME_ANALYSIS'})
    factors, pieces = bounded_factor(root, protocol, cert.read(TEMPLATE))
    diagnostic['coefficient_gcd_factorization'] = pieces
    restored = ZZ(1)
    for p, e in factors:
        if not p.is_prime(proof=True):
            raise ArithmeticError('unproved prime')
        restored *= p**e
    if restored != root or root**exponent != remainder:
        raise ArithmeticError('power factorization failed')
    restored = root**exponent
    for p, e in small:
        restored *= ZZ(p)**e
    if restored != common:
        raise ArithmeticError('factor reconstruction failed')
    primes = sorted(set([ZZ(p) for p, e in small]+[p for p, e in factors]))
    diagnostic['root_factorization'] = [[str(p), int(e)] for p, e in factors]
    checkpoint(output, {**base, 'status': 'RUNNING_AUXILIARY_REDUCTION'})
    pari.allocatemem(protocol['pari_stack_bytes'], silent=True)
    pari.set_real_precision_bits(protocol['precision_bits'])
    reducer = pari('(P,L)->{my(m,n);my(Q=hyperellminimalmodel(P,&m,L));my(S=hyperellred(Q,&n));[Q,m,S,n]}')
    minimal, m, reduced, n = reducer(pari(P), pari(primes))
    matrices = [Matrix(QQ, 2, 2, [z[1][i,j] for i,j in ((0,0),(0,1),(1,0),(1,1))]) for z in (m,n)]
    M = matrices[0]*matrices[1]; a,b,c,d = M.list()
    if not M.det():
        raise ArithmeticError('singular base map')
    AA, BB = helpers.homogeneous(A,8,a,b,c,d), helpers.homogeneous(B,12,a,b,c,d)
    extra_u, content = scaling.scale_for(values(AA), values(BB))
    new_A, new_B = AA/QQ(str(extra_u))**4, BB/QQ(str(extra_u))**6
    u = QQ(str(initial_u*extra_u))
    if helpers.homogeneous(original_A,8,a,b,c,d) != u**4*new_A or helpers.homogeneous(original_B,12,a,b,c,d) != u**6*new_B:
        raise ArithmeticError('exact elliptic coefficient identity failed')
    if helpers.homogeneous(new_A,8,d,-b,-c,a)*u**4 != M.det()**8*original_A or helpers.homogeneous(new_B,12,d,-b,-c,a)*u**6 != M.det()**12*original_B:
        raise ArithmeticError('inverse coefficient identity failed')
    result = {**base, 'status': 'PASS_EXACT_BASE_CHANGE', 'primes': list(map(str,primes)),
              'before_bits': scaling.bits(values(original_A)+values(original_B)),
              'constant_scaled_bits': scaling.bits(values(A)+values(B)),
              'after_bits': scaling.bits(values(new_A)+values(new_B)),
              'base_matrix_a_b_c_d': list(map(str,M.list())), 'scale_u': str(u),
              'A_coefficients_low_to_high': list(map(str,new_A.list())),
              'B_coefficients_low_to_high': list(map(str,new_B.list())),
              'auxiliary_minimization_transform': str(m), 'auxiliary_reduction_transform': str(n),
              'constant_scaling': content}
    checkpoint(output, result)
    print('MW16 BASE', presentation, result['before_bits'], '->', result['after_bits'], flush=True)


def batch(directory):
    protocol = cert.read(directory/'protocol.json')
    if protocol['sources'] != sources():
        raise ArithmeticError('sources changed')
    if (directory/'ledger.json').exists():
        raise FileExistsError('preserve prior bounded campaign')
    ledger = {'protocol_sha256': cert.hashed(directory/'protocol.json'), 'rows': []}
    for i, presentation in enumerate(protocol['presentations']):
        if i and not ledger['rows'][0].get('gate_passed'):
            ledger['rows'].append({'presentation_id': presentation, 'status': 'UNRUN_FIRST_FAMILY_GATE'})
            checkpoint(directory/'ledger.json',ledger)
            continue
        try:
            result = capture([sys.executable, str(Path(__file__).resolve()), 'run', '--directory', str(directory), '--presentation', presentation],
                             limits=Limits(protocol['worker_wall_seconds'], protocol['worker_rss_bytes']), log_path=directory/(presentation+'.log'))
            row = cert.read(directory/(presentation+'.json'))
            if row['status'] != 'PASS_EXACT_BASE_CHANGE':
                raise ArithmeticError('nonterminal worker output')
            entry = {'presentation_id': presentation, 'status': row['status'], 'before_bits': row['before_bits'],
                     'after_bits': row['after_bits'], 'gate_passed': 4*row['after_bits'] <= 3*row['before_bits'],
                     'result_sha256': cert.hashed(directory/(presentation+'.json')), 'supervision': result.supervision}
        except Exception as error:
            entry = {'presentation_id': presentation, 'status': 'FAILED_OR_CENSORED', 'error': str(error)}
        ledger['rows'].append(entry); checkpoint(directory/'ledger.json', ledger)
        print('MW16 BATCH', entry, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage', choices=['prepare','run','batch'])
    p.add_argument('--directory', type=Path, default=DIRECTORY)
    p.add_argument('--presentation', default='a1-presentation-01')
    args = p.parse_args()
    if args.stage == 'prepare': prepare(args.directory)
    elif args.stage == 'run': run(args.directory,args.presentation)
    else: batch(args.directory)
