"""Check gcd-supported RR order outputs without nfinit or nfcertify."""
import hashlib
import json
import runpy
import signal
from pathlib import Path
from sage.all import QQ, ZZ, matrix, prod, gcd

BASE = Path(__file__).resolve().parent
env = runpy.run_path(str(BASE/'verify_det1092_rr_global_pair.sage'))
ART, ROOT, R, coords = [env[k] for k in ['ART','ROOT','R','coords']]
OUT = ART/'det1092_rr_coprime_support_v1'
signal.alarm(60)
bindings = {}


def read(p):
    bindings[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text())


results = []
for i in [8,9]:
    inp, order, cert = [read(OUT/('case-%02d-%s.json' % (i,s)))
                        for s in ['input','order','certification']]
    for d in [inp,order,cert]:
        assert d['source_sha256'] == hashlib.sha256((BASE/'refine_rr_coprime_support.sage').read_bytes()).hexdigest()
        assert d['case_index'] == i and d['Selmer_dimension'] is None
    for p,h in inp['inputs'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest() == h
    f = R(inp['defining_polynomial'])
    assert f == R(order['defining_polynomial'])
    divisors = list(map(ZZ,inp['coprime_divisors']))
    assert all(c>1 for c in divisors)
    assert all(gcd(c,d)==1 for j,c in enumerate(divisors) for d in divisors[:j])
    assert prod(c**e for c,e in zip(divisors,inp['exponents'])) == abs(f.discriminant())
    basis = [R(b) for b in order['basis']]
    bits = env['check_order'](f,basis,order['order_discriminant'])
    assert bits == order['order_discriminant_bits']
    old = read(ART/'det1092_rr_order_join_v1'/('case-%02d.json'%i))
    old_basis = [R(b) for b in old['basis']]
    env['check_order'](f,old_basis,old['order_discriminant'])
    B = matrix(QQ,[coords(b) for b in basis]).transpose()
    inclusion = B.inverse()*matrix(QQ,[coords(b) for b in old_basis]).transpose()
    contained = all(v in ZZ for v in inclusion.list())
    if contained:
        assert ZZ(old['order_discriminant']) == ZZ(order['order_discriminant'])*inclusion.det()**2
    assert cert['status']=='UNRESOLVED_MAXIMAL_ORDER' and cert['maximal_order_certified'] is False
    assert cert['unresolved_cofactors'] and cert['unresolved_bits']==[ZZ(v).nbits() for v in cert['unresolved_cofactors']]
    results.append({'case_index':i,'order_discriminant_bits':bits,
                    'joined_order_contained':contained,'unresolved_bits':cert['unresolved_bits']})
print(json.dumps({'status':'PASS_EXACT_ORDER_AND_SUPPORT_REPLAY','cases':results,
    'inputs':bindings,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'dependency_sha256':hashlib.sha256((BASE/'verify_det1092_rr_global_pair.sage').read_bytes()).hexdigest(),
    'boundary':'Verifies integral orders, exact coprime decomposition and containment; certification failure is an execution receipt, not a proof of infeasibility or a Selmer result.'},indent=2))
