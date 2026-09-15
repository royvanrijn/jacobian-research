"""Join retained binary and partial RR orders without discarding integrality.

Exactly two fixed fields, 36 products per case, 60-second per-case alarm.
No factorization, number-field initialization, class group, or Selmer call.
"""
import hashlib
import json
import resource
import signal
import sys
import time
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, matrix, vector, lcm

sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OLD = ART / 'det1092_rr_global_pair_v3'
OUT = ART / 'det1092_rr_order_join_v1'
R = PolynomialRing(QQ, 'x')
x = R.gen()


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def coords(b):
    return vector(QQ, [b[j] for j in range(6)])


def check_order(f, basis):
    B = matrix(QQ, [coords(b) for b in basis])
    inv = B.inverse()
    assert all(v in ZZ for v in coords(R(1)) * inv)
    multiplication = []
    for b in basis:
        M = matrix(QQ, [coords(b*c % f)*inv for c in basis])
        assert all(v in ZZ for v in M.list())
        multiplication.append(M)
    traces = matrix(QQ, 6, 6, lambda i,j: (multiplication[i]*multiplication[j]).trace())
    disc = traces.det()
    assert disc == f.discriminant()*B.det()**2 and disc in ZZ and disc
    return B, ZZ(disc)


def run(i):
    signal.alarm(60)
    start = time.monotonic()
    paths = [OLD / ('case-%02d-%s.json' % (i, s)) for s in
             ['binary-order', 'reduced-generator', 'order']]
    binary, reduced, partial = [json.loads(p.read_text()) for p in paths]
    f = R(reduced['defining_polynomial'])
    alpha = R(reduced['monic_root_in_reduced_field'])
    assert R(binary['monic_polynomial'])(alpha) % f == 0
    assert matrix(QQ, [coords(alpha**j % f) for j in range(6)]).det()
    assert f == R(partial['defining_polynomial'])
    basis1 = [R(b)(alpha) % f for b in binary['basis']]
    basis2 = [R(b) for b in partial['order_basis']]
    B1, d1 = check_order(f, basis1)
    B2, d2 = check_order(f, basis2)
    assert d1 == ZZ(binary['order_discriminant'])
    assert d2 == ZZ(partial['order_discriminant'])
    # For commutative orders O1,O2, sum Z*(b_i*c_j) is already a ring.
    products = matrix(QQ, [coords(b*c % f) for b in basis1 for c in basis2])
    den = lcm([v.denominator() for v in products.list()])
    H = matrix(ZZ, den*products).hermite_form(include_zero_rows=False)
    assert H.nrows() == 6
    joined = matrix(QQ, H)/den
    basis = [sum(row[j]*x**j for j in range(6)) for row in joined]
    B, disc = check_order(f, basis)
    assert B == joined
    inv = B.inverse()
    assert all(v in ZZ for v in (B1*inv).list()+(B2*inv).list())
    indices = [abs(ZZ((old*inv).det())) for old in [B1,B2]]
    assert d1 == disc*indices[0]**2 and d2 == disc*indices[1]**2
    result = {
        'case_index': i, 'status': 'PASS_JOINED_INTEGRAL_ORDER',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in paths},
        'source_sha256': sha(Path(__file__)),
        'limits': {'wall_seconds': 60, 'address_space_bytes': 2*1024**3,
                   'pairwise_products': 36, 'factorization_calls': 0},
        'defining_polynomial': list(map(str, f.list())),
        'basis': [list(map(str, b.list())) for b in basis],
        'order_discriminant': str(disc),
        'order_discriminant_bits': abs(disc).nbits(),
        'input_order_discriminant_bits': [abs(d1).nbits(), abs(d2).nbits()],
        'input_order_indices': list(map(str, indices)),
        'elapsed_seconds': time.monotonic()-start,
        'maximal_order_certified': False, 'Selmer_dimension': None,
    }
    OUT.mkdir(exist_ok=True)
    path = OUT / ('case-%02d.json' % i)
    with path.open('x') as stream:
        stream.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['case_index', 'status',
        'order_discriminant_bits', 'input_order_discriminant_bits', 'elapsed_seconds']}), flush=True)
    signal.alarm(0)


if __name__ == '__main__':
    assert len(sys.argv) == 2 and sys.argv[1] in ['8', '9']
    run(int(sys.argv[1]))
