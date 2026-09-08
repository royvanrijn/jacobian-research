#!/usr/bin/env sage-python
"""Standalone rational-point proof from one completed funnel V3 epoch.

The rest of the cascade may be live. This checks subgroup independence only;
the separate V3 replay checks chart provenance and the frozen search policy.
"""
import argparse
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, gcd, lcm
from v3_warm_support import atomic, bindings, read, require, sha

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
KERNEL = Path(__file__).with_name('verify_det1092_funnel_small_conic_seed.sage')
finite_group = runpy.run_path(str(KERNEL))['finite_group']


def verify(run, case, stage_path, output):
    protocol = read(run/'protocol.json')
    bindings(ROOT, protocol['inputs'])
    bindings(ROOT, protocol['sources'])
    folder = run/'amplifiers'/case
    require(stage_path.is_relative_to(folder/'replay-M17'), 'stage outside selected amplifier')
    stage = read(stage_path)
    audit_path = stage_path.parent/stage['audit']
    require(audit_path.resolve().is_relative_to(folder), 'audit outside amplifier')
    require(sha(audit_path) == stage['audit_sha256'], 'stage/audit binding differs')
    audit = read(audit_path)
    require(audit['status'] == 'COMPLETE_DECLARED_FINITE_AUDIT', 'incomplete audit')
    cloud = ROOT/audit['input_path']
    require(cloud.resolve().is_relative_to(folder) and sha(cloud) == audit['input_sha256'], 'cloud binding differs')
    seed = run/'seeds'/case
    packet = read(seed/'m18.json')
    intake = read(seed/'intake.json')
    require(packet['status'] == 'CERTIFIED_M18' and packet['intake_sha256'] == sha(seed/'intake.json'), 'unbound seed')
    parent_path = ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
    chart_path = ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
    parent, chart = read(parent_path), read(chart_path)
    require(chart['export_sha256'] == sha(parent_path), 'parent/chart binding differs')
    R = PolynomialRing(QQ, 't')
    s = QQ(intake['parameter'])

    def val(value):
        return R(value['numerator'])(s)/R(value['denominator'])(s)

    ai = [val(v)*s.denominator()**k for v, k in zip(parent['a_invariants'], [2, 4, 6, 8, 12])]
    require(ai == list(map(QQ, audit['curve'])) == list(map(QQ, packet['curve'])) == list(map(QQ, intake['model'])), 'specialization differs')
    require(ai[:3] == [0, 0, 0], 'expected short model')
    E = EllipticCurve(QQ, ai)
    require(bool(E.discriminant()), 'singular fibre')
    points = [tuple(map(QQ, p)) for p in audit['independent_points']]
    rank = int(stage['after'])
    require(18 <= rank <= 100 and len(points) == rank, 'invalid retained rank')
    require(all(E(list(P)) for P in points), 'off-curve point')
    inherited = [(val(x)*s.denominator()**4, val(y)*s.denominator()**6)
                 for x, y in parent['basis_weierstrass_coordinates']]
    require(points[:17] == inherited and points[:18] == [tuple(map(QQ, P)) for P in packet['points']], 'initial subgroup prefix differs')
    a, b, c, d = map(ZZ, chart['parameter_matrix'])
    tau = (a*s+b)/(c*s+d)
    require(str(tau) == packet['original_parameter'] == intake['original_parameter'], 'original parameter differs')
    proof = audit['rank_certificate']
    primes = sorted(set([int(v['prime']) for v in proof['signatures']]
                        + [int(proof['no_rational_2_torsion_prime'])]))
    require(all(2 < p <= 10000 for p in primes), 'finite proof prime cap exceeded')
    rows, certificates, torsion = [], [], False
    for p in primes:
        require(ZZ(p).is_prime(proof=True) and bool(E.discriminant() % p), 'not a good prime')
        A, B = [int(v % p) for v in ai[3:]]
        allpoints, doubled, labels, reps, dim = finite_group(p, 0, A, B)
        reductions = []
        for x, y in points:
            z = lcm(x.denominator(), y.denominator())
            coords = [ZZ(x*z), ZZ(y*z), ZZ(z)]
            common = gcd(coords)
            X, Y, Z = [int((v//common) % p) for v in coords]
            require(bool(X or Y or Z), 'zero reduction vector')
            if not Z:
                require(X == 0 and bool(Y), 'bad point at infinity')
                reductions.append(None)
            else:
                reductions.append((X*pow(Z, -1, p) % p, Y*pow(Z, -1, p) % p))
        codes = [labels[P] for P in reductions]
        rows.extend([[(code >> bit) & 1 for code in codes] for bit in range(dim)])
        if p == proof['no_rational_2_torsion_prime']:
            require(len(allpoints) % 2 == 1, 'torsion exclusion failed')
            torsion = True
        certificates.append(dict(prime=p, finite_order=len(allpoints), dimension=dim,
            all_points=allpoints,
            doubled_points=sorted(doubled, key=lambda P: (-1, -1) if P is None else P),
            representatives=reps, point_reductions=reductions, column_codes=codes))
    M = matrix(GF(2), rows)
    require(torsion and M.ncols() == rank and M.rank() == rank, 'independence unresolved')
    require(M[:, :17].rank() == 17 and M[:, :18].rank() == 18, 'inherited prefix not certified')
    paths = [run/'protocol.json', stage_path, audit_path, cloud, seed/'m18.json', seed/'intake.json', parent_path, chart_path, KERNEL, Path(__file__)]
    result = dict(status='PASS_STANDALONE_FUNNEL_EPOCH', case=case, epoch=stage['epoch'],
        parameter=str(s), original_parameter=str(tau), rank_lower_bound=rank,
        curve=list(map(str, ai)), points=[list(map(str, P)) for P in points],
        proof_primes=primes, finite_group_certificates=certificates,
        matrix_rank=int(M.rank()), initial_rank=18, inherited_rank=17,
        inputs={str(p.relative_to(ROOT)): sha(p) for p in paths},
        full_search_policy_replayed=False,
        argument='Exact rational point identities and generic17/seed18 prefixes. Complete finite group quotients by doubling have full column rank. Odd-order good reduction excludes rational2-torsion; infinite descent proves the stated subgroup lower bound. Search-policy replay is separate. No exact rank, conductor or literature-wide novelty claim.')
    atomic(output, result, immutable=True)
    print('PASS_STANDALONE_FUNNEL_EPOCH', case, rank, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--case', required=True)
    p.add_argument('--stage', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    verify(a.run.resolve(), a.case, a.stage.resolve(), a.output.resolve())
