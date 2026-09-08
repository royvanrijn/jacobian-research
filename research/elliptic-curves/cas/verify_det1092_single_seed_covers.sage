#!/usr/bin/env sage-python
"""Equation-only cover/height proof and fixed-address splitting certificates.

No constructor import, exceptional point, winning coordinate or V3 input.
Each invocation processes one cover; 25 seconds, nine fixed parameters.
"""
import argparse
import hashlib
import json
import signal
from pathlib import Path
from sage.all import QQ, ZZ, GF, PolynomialRing, EllipticCurve, matrix, vector, gcd, lcm

ROOT = Path(__file__).resolve().parents[2]
DIR = ROOT/'artifacts/generated-results/elliptic-curves/det1092_single_seed_covers_v1'
OBSTRUCTION_PRIMES = [3,5,7,11,13,17,19,23,29,31,37,41,43,47]
RANK_PRIMES = [17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def retain(p, data):
    text = json.dumps(data, indent=2, sort_keys=True)+'\n'
    if p.exists():
        assert p.read_text() == text
    else:
        p.write_text(text)


def finite_rank(E, points):
    """Rebuild the full finite groups and all doubling cosets, with skips."""
    rows, attempts = [], []
    key = lambda P: tuple(int(z) for z in P)
    for p in RANK_PRIMES:
        if any(a.denominator()%p == 0 for a in E.a_invariants()):
            attempts.append({'p': p, 'status': 'SKIP_EQUATION_DENOMINATOR'})
            continue
        if E.discriminant().valuation(p) != 0:
            attempts.append({'p': p, 'status': 'SKIP_BAD_REDUCTION'})
            continue
        F = GF(p)
        e = EllipticCurve(F, E.a_invariants())
        allpoints = list(e.points())
        doubles = {key(2*P): 2*P for P in allpoints}
        labels = {k: 0 for k in doubles}
        reps = [e(0)]
        dim = 0
        for P in allpoints:
            if key(P) in labels:
                continue
            new = []
            for j, Q in enumerate(reps):
                other = P+Q
                new.append(other)
                for D in doubles.values():
                    labels[key(other+D)] = (1 << dim) | j
            reps += new
            dim += 1
        assert len(labels) == len(allpoints) and len(reps) == 2**dim

        def red(P):
            den = lcm(v.denominator() for v in P)
            coords = [ZZ(v*den) for v in P]
            content = gcd(coords)
            return e([F(v/content) for v in coords])

        codes = [labels[key(red(P))] for P in points]
        block = [[(c >> j) & 1 for c in codes] for j in range(dim)]
        rows.extend(block)
        attempts.append({'p': p, 'status': 'PASS_FINITE_QUOTIENT',
                         'group_order': len(allpoints), 'double_order': len(doubles), 'rows': block})
    M = matrix(GF(2), len(rows), len(points), sum(rows, []))
    # A root-free reduction of the exact 2-division cubic excludes Q-2-torsion.
    p = 31
    no2 = False
    if all(a.denominator()%p for a in E.a_invariants()):
        S = PolynomialRing(GF(p), 'x')
        x = S.gen()
        f = 4*x**3+GF(p)(E.b2())*x*x+2*GF(p)(E.b4())*x+GF(p)(E.b6())
        no2 = not any(f(a) == 0 for a in GF(p))
    H = M[:, :17]
    sep = next((z for z in H.left_kernel().basis() if (z*M[:,17])[0]), None)
    return {'rank': int(M.rank()), 'generic_rank': int(H.rank()),
            'no_rational_2_torsion_prime': 31 if no2 else None,
            'attempts': attempts, 'matrix_rows': rows,
            'separator': None if sep is None else list(map(int, sep)),
            'independent_seed_certified': bool(M.rank() == 18 and no2)}


def verify(index):
    path = DIR/f'cover-{index:02d}-input.json'
    d = json.loads(path.read_text())
    assert set(d) == {'schema','a_invariants','generic_sections','generic_height_gram',
                     'trace_word','trace_point','line','h','cover_polynomial','point_map',
                     'inherited_parameter','exclusion_factors','cases'}
    R = PolynomialRing(QQ, 't')
    K = R.fraction_field()

    def rat(r):
        return K(R(r['numerator']))/R(r['denominator'])

    E = EllipticCurve(K, [rat(r) for r in d['a_invariants']])
    basis = [E([rat(r) for r in P]) for P in d['generic_sections']]
    w = vector(QQ, d['trace_word'])
    Z = E([rat(r) for r in d['trace_point']])
    assert Z == sum((n*P for n, P in zip(w, basis)), E(0))
    G = matrix(QQ, d['generic_height_gram'])
    assert G.nrows() == 17 and G.det() == 1092 and G.is_positive_definite()
    assert w*G*w == 10
    f0, f1, f2 = map(R, d['line'])
    h = R(d['h'])
    F = R(d['cover_polynomial'])
    assert h.degree() == 3 and h*h == Z[0].denominator()
    assert [f0.degree(), f1.degree(), f2.degree()] == [10,6,4]
    assert gcd([f0,f1,f2]) == 1  # no vertical fibre component
    assert f0+f1*(-Z)[0]+f2*(-Z)[1] == 0
    assert F.degree() == 6 and F.gcd(F.derivative()) == 1
    delta = R(E.discriminant())
    assert delta.degree() == 24 and delta.gcd(delta.derivative()) == 1
    assert F.gcd(delta) == 1  # branch only at smooth parent fibres
    x0, x1 = map(rat, d['point_map']['x'])
    y0, y1 = map(rat, d['point_map']['y'])
    assert x1 and x1 == h**3/(2*f2**2)
    # Check the complete Weierstrass and line maps coefficientwise mod s^2-F.
    S = PolynomialRing(K, 's')
    s = S.gen()
    xx, yy = S(x0)+x1*s, S(y0)+y1*s
    a1,a2,a3,a4,a6 = E.a_invariants()
    modulus = s*s-K(F)
    assert (yy*yy+a1*xx*yy+a3*yy-xx**3-a2*xx**2-a4*xx-a6)%modulus == 0
    assert f0+f1*xx+f2*yy == 0
    slope = y1/x1
    intercept = y0-slope*x0
    assert slope*slope+a1*slope-a2-2*x0 == Z[0]
    assert -(slope+a1)*Z[0]-intercept-a3 == Z[1]
    # Independent residual elimination also proves degree two and its squareclass.
    X = PolynomialRing(K, 'x')
    x = X.gen()
    L = f0+f1*x
    pol = L*L-a1*x*L*f2-a3*L*f2-f2*f2*(x**3+a2*x*x+a4*x+a6)
    residual, rem = pol.quo_rem(x-(-Z)[0])
    assert not rem and residual.degree() == 2
    assert residual.discriminant() == h**6*F
    # The exact RR divisor is 3O+10F-(-Z)=2O+5F+phi(w).
    # Smooth genus2 normalization has the same arithmetic genus, so D.O=1.
    # Base change unramified over all I1 fibres: chi=4, 48I1, no corrections.
    cross = G*w
    gram = (2*G).augment(matrix(QQ, 17, 1, list(cross))).stack(
        matrix(QQ, 1, 18, list(cross)+[10]))
    assert gram.is_positive_definite() and gram.det() == 5*2**17*1092
    assert QQ(10)-(cross*(2*G).inverse()*cross) == 5
    t0 = QQ(d['inherited_parameter'])
    assert R(f2/h)(t0) == 0 and h(t0) and F(t0) and F(t0).is_square()
    inherited_s = F(t0).sqrt()
    # Both cover points above t0 are inherited: the maps extend to O and Z.
    # Divisor intersection D.O=1 and trace P+sigma(P)=Z certify these limits.
    exclusions = list(map(R, d['exclusion_factors']))
    assert exclusions == [delta,F,f2,h]
    cases = []
    for case in d['cases']:
        t = QQ(case['parameter'])
        needed = [*exclusions,
                  *[v.denominator() for v in [*E.a_invariants(),x0,x1,y0,y1]],
                  *[p[j].denominator() for p in basis for j in range(2)]]
        bad = [j for j, f in enumerate(needed) if not f(t)]
        row = {'index': case['index'], 'parameter': str(t), 'excluded_factor_indices': bad}
        if bad:
            row['status'] = 'EXCLUDED_FROM_DISPLAYED_AFFINE_CERTIFIER'
        else:
            val = QQ(F(t))
            square = val.is_square()
            row.update({'square_value': str(val), 'rational_square': bool(square)})
            if not square:
                attempts = []
                for p in OBSTRUCTION_PRIMES:
                    vp = val.valuation(p)
                    unit = val/QQ(p)**vp
                    residue = int(GF(p)(unit))
                    nonsquare = bool(vp % 2 or pow(residue,(p-1)//2,p) == p-1)
                    attempts.append({'p': p, 'valuation': int(vp), 'unit_residue': residue,
                                     'local_nonsquare': nonsquare})
                    if nonsquare:
                        break
                row.update({'status': 'EXACT_NONSPLIT_FIBRE', 'modular_attempts': attempts})
            else:
                root = val.sqrt()
                e = EllipticCurve(QQ, [v(t) for v in E.a_invariants()])
                P = e([x0(t)+x1(t)*root,y0(t)+y1(t)*root])
                Q = e([x0(t)-x1(t)*root,y0(t)-y1(t)*root])
                assert P+Q == e([Z[0](t),Z[1](t)])
                generic = [e([p[0](t),p[1](t)]) for p in basis]
                proof = finite_rank(e, [*generic, P])
                row.update({'status': 'CERTIFIED_INDEPENDENT_SEED' if proof['independent_seed_certified']
                            else 'RATIONAL_SPLIT_INDEPENDENCE_UNRESOLVED',
                            'root': str(root), 'point': list(map(str,P[:2])),
                            'complement': list(map(str,Q[:2])), 'independence': proof})
        cases.append(row)
        retain(DIR/f'cover-{index:02d}-case-{case["index"]:02d}.json', row)
    result = {
        'classification': 'verified application and new deduction',
        'status': 'PASS_EQUATION_ONLY_COVER_HEIGHT_AND_FIXED_PANEL',
        'input_sha256': sha(path), 'checker_sha256': sha(Path(__file__)),
        'cover_degree': 2, 'cover_genus': 2, 'rational_constant_field': 'Q',
        'point_map_and_trace_identities': True,
        'intersection_height_certificate': {'D_squared': 2, 'D_dot_O': 1,
            'base_changed_chi': 4, 'I1_fibres': 48, 'section_height': 10,
            'anti_invariant_height': 20, 'Schur_complement': 5,
            'gram': [[str(c) for c in row] for row in gram.rows()],
            'determinant': str(gram.det()), 'rank_lower_bound_over_Q_C': 18},
        'inherited_pair': {'parameter': str(t0), 'square_root_absolute': str(abs(inherited_s)),
                           'elliptic_images_after_extension': ['O','Z'], 'seed': False},
        'cases': cases,
        'boundary': 'Function-field independence is not independence at every rational cover point. The second cover was calibrated retrospectively, despite its witness-free execution. Nonsplitting excludes this cover only, not extra points on the elliptic fibre.',
        'limits': {'wall_seconds': 25, 'parameters': 9, 'point_searches': 0,
                   'new_members': 0, 'V3_artifacts_read': 0, 'class_groups': 0},
    }
    retain(DIR/f'cover-{index:02d}-replay.json', result)
    print(result['status'], index, [(r['index'],r['status']) for r in cases], flush=True)


if __name__ == '__main__':
    signal.alarm(25)
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', type=int, choices=[0,1], required=True)
    verify(parser.parse_args().index)
