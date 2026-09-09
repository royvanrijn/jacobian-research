#!/usr/bin/env sage-python
"""Independent cold-construction and admission replay, <=25 seconds.

Manual reversed rational-function group law verifies the trace; Sage finite
groups verify the producer's manual finite quotient engine independently.
Root-free reduction, not the producer's factorization, certifies nonhalving.
"""
import argparse
import csv
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, gcd, lcm, prime_range

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'det1092_euclidean_seed_admission_v1'
R = PolynomialRing(QQ, 't')
K = R.fraction_field()
PROOF_PRIME_CAP = 257


def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def decode(v): return K(R(v['numerator']))/R(v['denominator'])
def point(E, v): return E(0) if v is None else E(list(map(QQ, v)))
def word(E, basis, v): return sum((ZZ(n)*P for n, P in zip(v, basis)), E(0))
def key(P): return tuple(int(a) for a in P)
def save(name, value):
    p = OUT/name
    if p.exists(): assert read(p) == value, 'immutable replay changed'
    else:
        with p.open('x') as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write('\n')


def replay():
    protocol = read(OUT/'protocol.json')
    for name, digest in protocol['inputs'].items(): assert sha(ROOT/name) == digest
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in
              [OUT/'protocol.json', OUT/'construction.json',
               *[OUT/r['label']/n for r in protocol['roster']
                 for n in ['generic-frame.json', 'result.json']]]}
    save('replay-protocol.json', dict(checker_sha256=sha(Path(__file__)), inputs=inputs,
        seconds=25, nonhalving_proof_prime_cap=PROOF_PRIME_CAP,
        method='Manual reversed Q(t) group law; independent Sage finite quotients; '
               'exact square bounds; exact cycle relation or root-free duplication quartic.'))
    parent = read(ART/'curve302_recovered_mw17_parent_v1.json')
    original = EllipticCurve(K, [decode(a) for a in parent['a_invariants']])
    A, B = R(-original.c4()/48), R(-original.c6()/864)
    basis = []
    for row in parent['basis_weierstrass_coordinates']:
        x, y = map(decode, row)
        X, Y = x+original.b2()/12, y+(original.a1()*x+original.a3())/2
        assert Y*Y == X**3+A*X+B
        basis.append((X, Y))
    c = read(OUT/'construction.json')
    assert c['protocol_sha256'] == sha(OUT/'protocol.json') and c['mask'] == protocol['mask'] == 8044
    with (ART/'curve302_parent_degree2_multisection_orbits_v1.tsv').open() as stream:
        table = next(row for row in csv.DictReader(stream, delimiter='\t')
                     if int(row['orbit_mask']) == 8044)
    assert table['category'] == 'rational' and list(map(int, table['parent_MW17_w'].split())) == c['word']
    w = vector(ZZ, c['word'])
    assert w*matrix(ZZ, parent['generic_height_gram'])*w == 10
    def add(P, Q):
        if P is None: return Q
        if Q is None: return P
        x, y = P; u, v = Q
        if x == u and y == -v: return None
        slope = (v-y)/(u-x) if x != u else (3*x*x+A)/(2*y)
        z = slope*slope-x-u
        return (z, slope*(x-z)-y)
    trace = None
    for n, P in reversed(list(zip(c['word'], basis))):
        for unused in range(abs(n)):
            trace = add(trace, (P[0], -P[1] if n > 0 else P[1]))
    h, N, V, m, g, b, k, q = [R(c[name]) for name in ['h', 'nx', 'ny', 'm', 'g', 'b', 'k', 'q']]
    assert trace == (K(N)/h**2, K(V)/h**3)
    assert h.is_monic() and h.degree() == 3 and gcd(h, N) == 1
    assert [f.degree() for f in [N, V, m, g, b, k]] <= [10, 15, 5, 9, 4, 3]
    assert m*N+V == h*h*g and m*m-N == h*h*b and m*b-2*g == h*h*k
    assert 4*m*k-3*b*b-4*A == h*h*q and 1 <= q.degree() <= 2
    assert gcd(q, q.derivative()) == 1
    x0, x1, y0, y1 = [R(v) for v in c['maps']]
    assert (x0, x1, y0, y1) == (b/2, h/2, -h*k/2, -m/2)
    assert y0*y0+y1*y1*q == x0**3+3*x0*x1*x1*q+A*x0+B
    assert 2*y0*y1 == 3*x0*x0*x1+x1**3*q+A*x1
    chart = read(ART/'det1092_reduced_parameter_chart_v1/generic-proof.json')
    aa, bb, cc, dd = map(QQ, chart['parameter_matrix'])
    rows = []
    for entry in protocol['roster']:
        label = entry['label']; folder = OUT/label
        if 'original_parameter' in entry: tau = QQ(entry['original_parameter'])
        else:
            s = QQ(entry['reduced_parameter']); tau = (aa*s+bb)/(cc*s+dd)
        frame, result = read(folder/'generic-frame.json'), read(folder/'result.json')
        assert result['protocol_sha256'] == sha(OUT/'protocol.json')
        assert result['construction_sha256'] == sha(OUT/'construction.json')
        assert result['frame_sha256'] == sha(folder/'generic-frame.json')
        assert QQ(result['original_parameter']) == tau and QQ(result['q_value']) == q(tau)
        E = EllipticCurve(QQ, [A(tau), B(tau)])
        specialized = [E([v(tau) for v in P]) for P in basis]
        assert list(map(QQ, frame['curve'])) == list(E.a_invariants())
        assert [point(E, v) for v in frame['basis']] == specialized and frame['candidate_inputs'] == 0
        assert frame['prime_cap'] == protocol['limits']['prime_cap']
        groups = []; seen_odd = False; boundaries = []
        for record in frame['records']:
            p = record['prime']
            assert ZZ(p).is_prime(proof=True) and 2 < p <= frame['prime_cap']
            assert all(a.denominator() % p for a in E.a_invariants()) and E.discriminant() % p
            ep = EllipticCurve(GF(p), E.a_invariants()); points = sorted(ep.points(), key=key)
            doubles = {key(2*P): 2*P for P in points}
            labels = {v: 0 for v in doubles}; reps = [ep(0)]; dim = 0
            for P in points:
                if key(P) in labels: continue
                new = []
                for i, Q in enumerate(reps):
                    V0 = P+Q; new.append(V0)
                    for D in doubles.values():
                        v = key(V0+D); assert v not in labels
                        labels[v] = i | (1 << dim)
                reps += new; dim += 1
            assert len(labels) == len(points) == len(doubles)*2**dim == record['order']
            assert dim == record['dimension']
            seen_odd |= bool(len(points) % 2)
            if p == frame['no_two_torsion_prime']: assert len(points) % 2
            groups.append((p, ep, labels, dim))
            boundaries.append((sum(v[3] for v in groups), seen_odd))
        def code(P):
            out = []
            for p, ep, labels, dim in groups:
                den = lcm([v.denominator() for v in P]); coords = [ZZ(v*den) for v in P]
                common = gcd(coords); red = ep([int(v/common % p) for v in coords])
                value = labels[key(red)]
                out.extend((value >> j) & 1 for j in range(dim))
            return vector(GF(2), out)
        M = matrix(GF(2), [code(P) for P in specialized]).transpose()
        assert M.rank() == 17 and seen_odd
        assert frame['no_two_torsion_prime'] in [r['prime'] for r in frame['records']]
        for nrows, odd in boundaries[:-1]: assert not (odd and M[:nrows].rank() == 17)
        last = frame['records'][-1]['prime']
        good = [int(p) for p in prime_range(3, last+1)
                if all(a.denominator() % p for a in E.a_invariants()) and E.discriminant() % p]
        assert good == [r['prime'] for r in frame['records']]
        row = dict(label=label, status=result['status'], footprint_rank=17, last_footprint_prime=last)
        value = q(tau); cert = result['square_certificate']
        if value < 0: assert result['status'] == 'EXACT_NONSPLIT' and cert['negative']
        else:
            n, d = value.numerator(), value.denominator()
            aa0, bb0 = ZZ(cert['numerator_floor_sqrt']), ZZ(cert['denominator_floor_sqrt'])
            assert (ZZ(cert['numerator']), ZZ(cert['denominator'])) == (n, d)
            assert aa0*aa0 <= n < (aa0+1)**2 and bb0*bb0 <= d < (bb0+1)**2
            if aa0*aa0 != n or bb0*bb0 != d: assert result['status'] == 'EXACT_NONSPLIT'
            else:
                W = QQ(result['root']); assert W >= 0 and W*W == value
                P = E([x0(tau)+x1(tau)*W, y0(tau)+y1(tau)*W])
                other = E([x0(tau)-x1(tau)*W, y0(tau)-y1(tau)*W])
                assert point(E, result['point']) == P and point(E, result['companion']) == other
                assert P+other == word(E, specialized, c['word'])
                decision = result['admission']; current = P
                assert decision['status'] == result['status'] and point(E, decision['original']) == P
                for step, hop in enumerate(decision['history']):
                    assert step == hop['step'] and current == point(E, hop['point'])
                    assert len(hop['parity_word']) == 17 and set(hop['parity_word']) <= {0, 1}
                    target = current-word(E, specialized, hop['parity_word'])
                    assert target == point(E, hop['target']) and not any(code(target))
                    if 'next' in hop:
                        current = point(E, hop['next']); assert 2*current == target
                    else: assert step == len(decision['history'])-1
                row.update(reason=decision['reason'], steps=decision['steps'])
                if decision['status'] == 'INHERITED_RATIONAL_SPAN':
                    multiplier = ZZ(decision['relation_multiplier'])
                    assert multiplier > 0 and multiplier % 2
                    assert multiplier*P == word(E, specialized, decision['relation_word'])
                    row.update(relation_multiplier=str(multiplier), relation_word=decision['relation_word'])
                else:
                    assert decision['status'] == 'NEW_INDEPENDENT_DIRECTION'
                    assert current == point(E, decision['terminal'])
                    assert P == 2**decision['steps']*current+word(E, specialized, decision['chain_word'])
                    if decision['reason'] == 'FINITE_FOOTPRINT_ESCAPE':
                        assert M.augment(matrix(GF(2), len(code(current)), 1, list(code(current)))).rank() == 18
                    else:
                        assert decision['reason'] == 'GLOBAL_NONHALVING_ESCAPE' and not target.is_zero()
                        S = PolynomialRing(QQ, 'x'); x = S.gen(); a, a4, a6 = target[0], E.a4(), E.a6()
                        f = x**4-4*a*x**3-2*a4*x*x-(8*a6+4*a4*a)*x+a4*a4-4*a6*a
                        assert f == S(decision['history'][-1]['halving']['polynomial'])
                        den = lcm([v.denominator() for v in f]); integers = [ZZ(v*den) for v in f]
                        content = gcd(integers); integers = [v/content for v in integers]
                        witness = None
                        for p in prime_range(3, PROOF_PRIME_CAP+1):
                            if not integers[-1] % p: continue
                            residues = [int(v % p) for v in integers]
                            values = [sum(v*pow(z, j, int(p)) for j, v in enumerate(residues)) % p for z in range(p)]
                            if all(values):
                                witness = dict(prime=int(p), residues=residues, values=list(map(int, values)))
                                break
                        assert witness is not None, 'UNKNOWN_PROOF_PRIME_CAP'
                        row['nonhalving_certificate'] = witness
        rows.append(row)
    return dict(status='PASS_INDEPENDENT_COLD_CONSTRUCTION_AND_ADMISSION', cases=rows,
        replay_protocol_sha256=sha(OUT/'replay-protocol.json'),
        classification='Verified application, no new rank or prospective302 seed',
        point_searches=0, complete_atlas_runs=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write', action='store_true'); args = parser.parse_args()
    signal.alarm(25)
    result = replay()
    if args.write: save('independent-replay.json', result)
    else: assert result == read(OUT/'independent-replay.json')
    print(result['status'], flush=True)
    for row in result['cases']: print(row, flush=True)
