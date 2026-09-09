#!/usr/bin/env python3
"""Checkpoint exact local conductors; replay uses no factor discovery."""
import argparse
from pathlib import Path
from sage.all import EllipticCurve, QQ, ZZ, pari, prime_range, proof, NumberField, polygen
from sage.version import version
from v3_warm_support import read, atomic, sha, require
proof.all(True)


def local(E, ep, remaining, p, certificate=None):
    pc = pari(certificate) if certificate is not None else pari(p).primecert()
    root = ZZ(pc) if pc.type() == 't_INT' else ZZ(pc[0][0])
    require(root == p and pari.primecertisvalid(pc), 'prime certificate failed')
    e = remaining.valuation(p)
    require(e > 0, 'listed prime does not divide remaining discriminant')
    sd, pd = E.local_data(E.base_ring().ideal(p), algorithm='generic', proof=True), ep.elllocalred(p)
    f = int(sd.conductor_valuation())
    require(f == int(pd[0]), 'Sage and PARI local conductor disagreement')
    return {'prime': str(p), 'prime_certificate': str(pc), 'discriminant_valuation': int(e),
            'conductor_exponent': f, 'minimal_discriminant_valuation': int(sd.discriminant_valuation()),
            'kodaira': str(sd.kodaira_symbol()), 'pari_local_reduction': str(pd)}


def run(packet, output, check=False):
    source = read(packet); curve = list(map(QQ, source['curve']))
    require(all(a.denominator() == 1 for a in curve), 'integral model required')
    # Sage 10.9's generic Tate path can call QQ.uniformizer on nonminimal
    # models. Q[a]/(a-1) is exactly Q and supplies the number-field interface
    # that path expects, without using PARI to minimize the curve first.
    K = NumberField(polygen(QQ)-1, 'a')
    E, ep = EllipticCurve(K, curve), pari.ellinit(curve)
    delta = abs(ZZ(E.discriminant())); require(delta == abs(ZZ(ep.disc())), 'discriminant differs')
    remaining, known, rows = delta, ZZ(1), []
    def result():
        require(remaining * ZZ.prod(ZZ(r['prime'])**r['discriminant_valuation'] for r in rows) == delta,
                'discriminant reconstruction differs')
        return {'id': source['id'], 'curve': source['curve'], 'rank_lower_bound': source['rank_lower_bound'],
                'input_sha256': sha(packet), 'worker_sha256': sha(Path(__file__)),
                'status': 'EXACT' if remaining == 1 else 'UNKNOWN', 'local_data': rows,
                'discriminant_abs': str(delta), 'remaining_cofactor': str(remaining),
                'conductor_divisor': str(known), 'conductor_upper_bound': str(known * remaining),
                'exact_conductor': str(known) if remaining == 1 else None,
                'sage_version': version, 'pari_version': str(pari.version()),
                'argument': 'Exact primality certificates and agreement of Sage generic Tate over the degree-one presentation Q[a]/(a-1) with PARI elllocalred over Q. The local product divides N; N is at most that product times the unprocessed discriminant cofactor. Exact only when the cofactor is one. No rank, minimal-model metric or record claim.'}
    def consume(p, pc=None):
        nonlocal remaining, known
        r = local(E, ep, remaining, p, pc)
        remaining //= p**r['discriminant_valuation']; known *= p**r['conductor_exponent']; rows.append(r)
        if not check: atomic(output, result())
    if check:
        saved = read(output)
        for row in saved['local_data']:
            consume(ZZ(row['prime']), row['prime_certificate'])
        require(saved == result(), 'independent local replay differs')
        print('CONDUCTOR_REPLAY_PASS', source['id'], saved['status'], flush=True)
        return
    require(not output.exists(), 'preserve existing conductor checkpoint')
    atomic(output, result())
    candidates = set(map(ZZ, source['known_primes'])) | set(prime_range(2, 1001))
    for p in sorted(candidates):
        if remaining % p == 0: consume(p)
    # Retain earlier discovered factors, including composite splits, instead of
    # repeating factorization of the full discriminant from scratch.
    chunks = list(map(ZZ, source['factor_hints']))
    chunks.append(remaining)
    for candidate in sorted(chunks, key=lambda n: (len(str(n)), n)):
        n = remaining.gcd(candidate)
        if n == 1: continue
        factors = n.factor(proof=True)
        require(ZZ.prod(p**e for p,e in factors) == n, 'factor product differs')
        for p,e in factors:
            if remaining % p == 0: consume(p)
    print('CONDUCTOR_BUILD_FINISHED', source['id'], result()['status'], flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--packet',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--check',action='store_true')
    a=p.parse_args();run(a.packet,a.output,a.check)
