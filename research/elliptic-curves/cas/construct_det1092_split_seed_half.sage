#!/usr/bin/env sage-python
"""One exact halving attempt on the unresolved conic seed, with no search box.

A mod2 kernel only proposes a coset to halve. A rational doubling identity
and an independent18-point certificate are required to export a new seed.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, GF, EllipticCurve, matrix, gcd, lcm
import det1092_funnel as f
from det1092_funnel_worker import checked_protocol, prepare, conic_points, first_seed, seal_seed
from v3_warm_support import atomic, bindings, read, require, sha

KERNEL = Path(__file__).with_name('verify_det1092_funnel_small_conic_seed.sage')
finite_group = runpy.run_path(str(KERNEL))['finite_group']


def construct(parent_run, source, output):
    parent = read(parent_run/'protocol.json')
    bindings(f.ROOT, parent['inputs']); bindings(f.ROOT, parent['sources'])
    require(read(source/'result.json')['status'] == 'CONIC_POINT_INDEPENDENCE_UNRESOLVED', 'expected unresolved seed')
    row = read(source/'intake.json')
    paths = [source/p for p in ('intake.json', 'm17.json', 'conic-points.json', 'result.json')]
    p = dict(parent)
    p.update(schema='det1092-split-seed-one-half.v1', profile='one_exact_half',
        sources=f.source_bindings(), inputs={**parent['inputs'], **{str(q.relative_to(f.ROOT)):sha(q) for q in paths}},
        population=0, maximum_draws=0, domain='det1092-split-seed-one-half-v1',
        halving=dict(branch_index=0, maximum_layers=1, maximum_cosets=1, wall_seconds=120,
                     rss_bytes=1024**3, quartic_search_charts=0),
        scope='One exact rational half of one finite-mod2 coset proposal on the existing unresolved conic fibre. Exact doubling and independent M17+Q rank18 required. Separate version; original seed results unchanged.')
    atomic(output/'protocol.json', p, immutable=True)
    checked_protocol(output)
    atomic(output/'selection.json', dict(status='SEALED_SINGLE_HALVING_ATTEMPT',
        parameters=1, seed_inputs=[row], arithmetic_candidates=[row],
        protocol_digest=f.digest(f.packed(p))), immutable=True)
    seed = output/'seeds'/row['id']
    model, base, proof = prepare(seed, row)
    points = conic_points(row)
    require([list(map(str, P)) for P in points] == read(source/'conic-points.json'), 'conic recomputation differs')
    P = points[0]
    E = EllipticCurve(QQ, model)
    rows, primes = [], [s['prime'] for s in proof['signatures']]
    for prime in primes:
        A, B = [int(QQ(v.numerator)/QQ(v.denominator) % prime) for v in model[3:]]
        _, _, labels, _, dim = finite_group(prime, 0, A, B)
        codes = []
        for xx, yy in (*base, P):
            x, y = QQ(str(xx)), QQ(str(yy))
            z = lcm(x.denominator(), y.denominator())
            coords = [ZZ(x*z), ZZ(y*z), ZZ(z)]
            common = gcd(coords)
            X, Y, Z = [int((v//common) % prime) for v in coords]
            if not Z:
                require(X == 0 and bool(Y), 'bad infinity reduction')
                point = None
            else:
                point = (X*pow(Z, -1, prime) % prime, Y*pow(Z, -1, prime) % prime)
            codes.append(labels[point])
        rows.extend([[(code >> bit) & 1 for code in codes] for bit in range(dim)])
    M = matrix(GF(2), rows)
    require(M[:, :17].rank() == 17, 'generic finite proof differs')
    kernel = M.right_kernel()
    require(kernel.dimension() == 1, 'no unique halving coset proposal')
    word = list(map(int, kernel.basis()[0]))
    require(word[17] == 1, 'kernel lost candidate coordinate')
    rational_base = [E(list(map(QQ, q))) for q in base]
    target = E(list(map(QQ, P)))-sum((word[i]*Q for i, Q in enumerate(rational_base)), E(0))
    proposal = dict(finite_primes=primes, word=word, conic_point=list(map(str, P)),
                    target=list(map(str, target.xy())), proposal_only=True)
    atomic(output/'halving-proposal.json', proposal, immutable=True)
    halves = sorted(target.division_points(2), key=lambda Q:Q.xy())
    require(all(2*Q == target for Q in halves), 'exact doubling failed')
    atomic(output/'halves.json', [list(map(str, Q.xy())) for Q in halves], immutable=True)
    found = first_seed(model, base, proof, [tuple(F(str(v)) for v in Q.xy()) for Q in halves])
    if found:
        seal_seed(seed, row, model, base, *found,
            evidence=dict(kind='one_exact_conic_half', proposal_sha256=sha(output/'halving-proposal.json'),
                          halves_sha256=sha(output/'halves.json')))
    else:
        atomic(seed/'result.json', dict(status='HALVING_SEED_UNRESOLVED', parameter=row['parameter'],
            rank_lower_bound=17, rational_halves=len(halves)), immutable=True)
    result = dict(status=read(seed/'result.json')['status'], parameter=row['parameter'],
        rational_halves=len(halves), one_layer_only=True, quartic_search_charts=0,
        boundary='Exact rational doubling and certified seed if present. No failed halving/certificate is a whole-curve rank upper bound.')
    atomic(output/'constructive-result.json', result, immutable=True)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent-run', type=Path, required=True)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    construct(a.parent_run.resolve(), a.source.resolve(), a.output.resolve())
