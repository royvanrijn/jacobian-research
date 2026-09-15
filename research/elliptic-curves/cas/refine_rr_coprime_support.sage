"""One bounded composite-support order attempt on each frozen RR field.

Use only gcd splitting and primes through1000 for support preparation.
PARI's composite-list nfinit contract does not assume the entries prime.
"""
import hashlib
import json
import resource
import signal
import sys
import time
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, pari, prime_range, gcd

sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'det1092_rr_coprime_support_v1'
R = PolynomialRing(QQ, 'x')
start = time.monotonic()


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(i, stage, data):
    OUT.mkdir(exist_ok=True)
    data.update(case_index=i, stage=stage, source_sha256=sha(Path(__file__)),
                elapsed_seconds=time.monotonic()-start, Selmer_dimension=None)
    with (OUT/('case-%02d-%s.json' % (i,stage))).open('x') as s:
        s.write(json.dumps(data, indent=2)+'\n')
    print(json.dumps({k: data[k] for k in data if k in
        ['case_index','stage','elapsed_seconds','status','order_discriminant_bits','unresolved_bits']}), flush=True)


def run(i):
    def expired(signum, frame):
        raise TimeoutError('60-second whole-attempt wall limit')
    signal.signal(signal.SIGALRM, expired)
    signal.alarm(60)
    try:
        path = ART/'det1092_rr_order_join_v1'/('case-%02d.json' % i)
        joined = json.loads(path.read_text())
        f = R(joined['defining_polynomial'])
        D = abs(ZZ(f.discriminant()))
        d = abs(ZZ(joined['order_discriminant']))
        index = (D//d).isqrt()
        assert D == d*index**2
        # Every value is a divisor of D; gcd refinement loses no prime support.
        initial = [D, d, index]+[gcd(D,ZZ(v)) for v in joined['input_order_indices']]
        initial += [ZZ(p) for p in prime_range(1001) if D%p == 0]
        chunks = sorted(set(v for v in initial if v > 1))
        while True:
            found = False
            for j,a in enumerate(chunks):
                for k in range(j):
                    b = chunks[k]
                    g = gcd(a,b)
                    if g > 1:
                        chunks = sorted(set(chunks[:k]+chunks[k+1:j]+chunks[j+1:]+[g,a//g,b//g])-{ZZ(1)})
                        found = True
                        break
                if found:
                    break
            if not found:
                break
        rest = D
        powers = []
        for c in chunks:
            e = 0
            while rest%c == 0:
                rest //= c
                e += 1
            powers.append(e)
        assert rest == 1 and all(gcd(a,b)==1 for j,a in enumerate(chunks) for b in chunks[:j])
        write(i,'input',{'inputs':{str(path.relative_to(ROOT)):sha(path)},
            'defining_polynomial':list(map(str,f.list())),
            'coprime_divisors':list(map(str,chunks)), 'exponents':powers,
            'divisor_bits':[c.nbits() for c in chunks],
            'limits':{'wall_seconds':60,'address_space_bytes':8*1024**3,'pari_stack_bytes':512*1024**2,
                      'nfinit_calls':1,'nfcertify_calls':1,'trial_prime_bound':1000},
            'boundary':'Coprime divisors are not assumed prime or squarefree. No supplied-basis certification.'})
        pari.allocatemem(32*1024**2,512*1024**2,silent=True)
        nf = pari.nfinit(pari([pari(f),pari(chunks)]))
        write(i,'order',{'status':'PARTIAL_ORDER_REQUIRES_CERTIFICATION',
            'defining_polynomial':list(map(str,R(nf.nf_get_pol()).list())),
            'basis':[list(map(str,R(b).list())) for b in nf.nf_get_zk()],
            'order_discriminant':str(nf.disc()), 'order_discriminant_bits':abs(ZZ(nf.disc())).nbits()})
        pending = list(pari.nfcertify(nf))
        write(i,'certification',{'status':'CERTIFIED_MAXIMAL_ORDER' if not pending else 'UNRESOLVED_MAXIMAL_ORDER',
            'unresolved_cofactors':list(map(str,pending)), 'unresolved_bits':[ZZ(v).nbits() for v in pending],
            'maximal_order_certified':not pending,'supported_squareclasses_complete':False})
    except (TimeoutError, Exception) as exc:
        write(i,'failure',{'status':'INCOMPLETE', 'exception':type(exc).__name__, 'message':str(exc)})
        raise
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    assert len(sys.argv)==2 and sys.argv[1] in ['8','9']
    run(int(sys.argv[1]))
