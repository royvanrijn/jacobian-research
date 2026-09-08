#!/usr/bin/env sage-python
"""Standalone equation and finite-group proof for any retained funnel M18.

The complete finite-group kernel comes from the standalone small-conic checker.
No point constructor, Kummer implementation or repository rank backend is used.
This proves the point subgroup, not the seed search's ordering or completeness.
"""
import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, gcd, lcm

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
KERNEL = Path(__file__).with_name('verify_det1092_funnel_small_conic_seed.sage')
finite_group = runpy.run_path(str(KERNEL))['finite_group']


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(seed, output):
    packet = read(seed/'m18.json')
    intake = read(seed/'intake.json')
    assert packet['status'] == 'CERTIFIED_M18'
    assert packet['intake_sha256'] == sha(seed/'intake.json')
    assert packet['parameter'] == intake['parameter']
    paths = [seed/'m18.json', seed/'intake.json', KERNEL, Path(__file__)]
    parent_path = ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
    chart_path = ART/'det1092_reduced_parameter_chart_v1/generic-proof.json'
    parent, chart = read(parent_path), read(chart_path)
    assert chart['export_sha256'] == sha(parent_path)
    paths += [parent_path, chart_path]
    R = PolynomialRing(QQ, 't')

    def val(value, t):
        return R(value['numerator'])(t)/R(value['denominator'])(t)

    s = QQ(packet['parameter'])
    ai = [val(v, s)*s.denominator()**k
          for v, k in zip(parent['a_invariants'], [2, 4, 6, 8, 12])]
    assert ai == list(map(QQ, packet['curve'])) == list(map(QQ, intake['model']))
    assert ai[:3] == [0, 0, 0]
    E = EllipticCurve(QQ, ai)
    assert E.discriminant()
    points = [tuple(map(QQ, p)) for p in packet['points']]
    assert len(points) == 18 and all(E(list(p)) for p in points)
    expected = [(val(x, s)*s.denominator()**4, val(y, s)*s.denominator()**6)
                for x, y in parent['basis_weierstrass_coordinates']]
    assert points[:17] == expected
    a, b, c, d = map(ZZ, chart['parameter_matrix'])
    tau = (a*s+b)/(c*s+d)
    assert str(tau) == packet['original_parameter'] == intake['original_parameter']
    conic_checked = packet['evidence']['kind'] in ('frozen_conic', 'primitive_conic_solve', 'exact_split_followup')
    if conic_checked:
        original_path = ART/'curve302_recovered_mw17_parent_v1.json'
        cover_path = ART/'det1092_orbit8044_rank18_base_change_v2.json'
        original, cover = read(original_path), read(cover_path)
        paths += [original_path, cover_path]
        old = EllipticCurve(QQ, [val(v, tau) for v in original['a_invariants']])
        h = c*s.numerator()+d*s.denominator()
        w = QQ(chart['weierstrass_u'])
        X, Y = points[17]
        x = X*w**2/h**4-old.b2()/12
        y = Y*w**3/h**6-(old.a1()*x+old.a3())/2
        assert old([x, y])
        rc, rb, ra = [val(v, tau) for v in cover['lift']['residual_coefficients']]
        f0, f1, f2 = [R(v)(tau) for v in cover['lift']['line_coefficients']]
        assert ra*x*x+rb*x+rc == 0 and f0+f1*x+f2*y == 0
    proof = packet['proof']
    primes = sorted(set([int(v['prime']) for v in proof['signatures']]
                        + [int(proof['no_rational_2_torsion_prime'])]))
    rows, certificates, torsion = [], [], False
    for p in primes:
        assert ZZ(p).is_prime(proof=True) and p > 2 and E.discriminant() % p
        A, B = [int(v % p) for v in ai[3:]]
        allpoints, doubled, labels, reps, dim = finite_group(p, 0, A, B)
        reductions = []
        for x, y in points:
            z = lcm(x.denominator(), y.denominator())
            coords = [ZZ(x*z), ZZ(y*z), ZZ(z)]
            common = gcd(coords)
            X, Y, Z = [int((v//common) % p) for v in coords]
            assert X or Y or Z
            if not Z:
                assert X == 0 and Y
                reductions.append(None)
            else:
                reductions.append((X*pow(Z, -1, p) % p, Y*pow(Z, -1, p) % p))
        codes = [labels[P] for P in reductions]
        rows.extend([[(code >> bit) & 1 for code in codes] for bit in range(dim)])
        if p == proof['no_rational_2_torsion_prime']:
            assert len(allpoints) % 2 == 1
            torsion = True
        certificates.append(dict(prime=p, finite_order=len(allpoints), dimension=dim,
            all_points=allpoints,
            doubled_points=sorted(doubled, key=lambda P: (-1, -1) if P is None else P),
            representatives=reps, point_reductions=reductions, column_codes=codes))
    M = matrix(GF(2), rows)
    assert torsion and M.ncols() == 18 and M.rank() == 18 and M[:, :17].rank() == 17
    j = E.j_invariant()
    result = dict(status='PASS_STANDALONE_FUNNEL_M18', parameter=str(s),
        original_parameter=str(tau), rank_lower_bound=18, inherited_rank=17,
        j_numerator_bits=int(abs(j.numerator()).nbits()),
        j_denominator_bits=int(j.denominator().nbits()), conic_incidence_checked=conic_checked,
        proof_primes=primes, finite_group_certificates=certificates, matrix_rank=18,
        inputs={str(p.relative_to(ROOT)): sha(p) for p in paths},
        argument='Exact specialized parent and generic17 prefix; all18 rational points satisfy the equation. Complete finite group quotients by doubling give18 independent columns. Odd-order good reduction excludes rational2-torsion; infinite descent proves independence. No full-rank upper bound, novelty or conductor claim.')
    payload = json.dumps(result, indent=2, sort_keys=True, default=int)+'\n'
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        assert output.read_text() == payload
    else:
        output.write_text(payload)
    print('PASS_STANDALONE_FUNNEL_M18', str(s), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    verify(a.seed.resolve(), a.output.resolve())
