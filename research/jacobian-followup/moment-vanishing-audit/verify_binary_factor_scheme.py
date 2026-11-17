#!/usr/bin/env python3
"""Exact inputs for the constructive binary-cubic/factor scheme classification."""
from pathlib import Path
from itertools import combinations, permutations
import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
import sympy as S
from sympy.matrices.normalforms import smith_normal_form

P = Path(__file__).resolve().parent
RUNS = ['binary_cubic_coefficients', 'binary_cubic_primary_audit',
        'binary_factor_socles', 'binary_factor_weighted_geometry', 'binary_factor_torsion']
COORDINATE_CASES = [
    ('MIN_ZERO', 'b,c,e,f', 'a,d'), ('MIN_PAIR', 'c,d,e,f', 'a,b'),
    ('FACE_FULL_PAIR', 'a,b,d', 'c,e,f'), ('FACE_SAME', 'a,b,c', 'd,e,f'),
    ('FACE_MIXED', 'a,b,f', 'c,d,e'), ('EDGE_FULL_PAIR', 'a,d', 'b,c,e,f'),
    ('EDGE_SAME', 'a,b', 'c,d,e,f'), ('EDGE_MIXED', 'a,e', 'b,c,d,f'),
    ('VERTEX', 'a', 'b,c,d,e,f'), ('ORIGIN', '', 'a,b,c,d,e,f')]
WEIGHTED_CASES = [
    ('FACE_SAME', 'abc', 'def', [1, 1, 1], 13, 84),
    ('EDGE_FULL_PAIR', 'ad', 'bcef', [1, 1, 1, 1], 18, 62),
    ('EDGE_SAME', 'ab', 'cdef', [1, 1, 1, 2], 21, 8),
    ('VERTEX', 'a', 'bcdef', [1, 1, 1, 2, 2], 32, 2),
    ('ORIGIN', '', 'abcdef', [1, 1, 1, 1, 1, 1], 30, 13)]


def sha(name):
    return hashlib.sha256((P/name).read_bytes()).hexdigest()


def read_blocks(text):
    # Singular comments announce benign variable/library redefinitions.
    clean = '\n'.join(line for line in text.splitlines() if not line.startswith('//'))
    markers = list(re.finditer(r'(?m)^([A-Z][A-Z0-9_]*)$', clean))
    blocks = {}
    for i, marker in enumerate(markers):
        body = clean[marker.end():markers[i+1].start() if i+1 < len(markers) else len(clean)].strip()
        blocks.setdefault(marker[1], []).append(body)
    return blocks


def one(blocks, name):
    assert len(blocks[name]) == 1, name
    return blocks[name][0]


def vector(blocks, name):
    out = list(map(int, re.findall(r'-?\d+', one(blocks, name))))
    while out and out[-1] == 0:
        out.pop()
    return out


def parse_polynomial(raw, env):
    answer = 0
    for token in re.findall(r'[+-]?[^+-]+', raw.strip().strip(',').replace(' ', '')):
        sign = -1 if token.startswith('-') else 1
        token = token.lstrip('+-')
        number = re.match(r'\d+(?:/\d+)?', token)
        coefficient = S.Rational(number[0]) if number else S.Integer(1)
        rest = token[len(number[0]):] if number else token
        monomial = 1
        while rest:
            part = re.match(r'([a-z])(\d*)', rest)
            assert part and part[1] in env, (raw, rest)
            monomial *= env[part[1]]**(int(part[2]) if part[2] else 1)
            rest = rest[len(part[0]):]
        answer += sign*coefficient*monomial
    return S.expand(answer)


def symmetry_maps():
    for perm in permutations(range(3)):
        for flip in [False, True]:
            yield {str('abcdef'[i]): 'abcdef'[perm[i % 3] + 3*((i//3) ^ flip)] for i in range(6)}


def main():
    started = time.monotonic()
    ap = argparse.ArgumentParser()
    ap.add_argument('--audit-existing', action='store_true')
    args = ap.parse_args()
    if not args.audit_existing:
        for stem in RUNS:
            run = subprocess.run([sys.executable, str(P/('run_'+stem+'.py'))], text=True, capture_output=True)
            assert run.returncode == 0, (stem, run.stdout, run.stderr)
    logs, sources, provenance = {}, {}, {}
    for stem in RUNS:
        rec = json.loads((P/(stem+'_run.json')).read_text())
        assert rec['terminal'] and rec['completed'] and rec['exit_code'] == 0 and not rec['timed_out']
        assert rec['wall_limit_seconds'] == 50 and rec['memory_limit_MiB'] == 768
        assert rec['input_sha256'] == sha(stem+'.sing') and rec['log_sha256'] == sha(stem+'_stdout.txt')
        raw = (P/(stem+'_stdout.txt')).read_text()
        assert '?' not in raw and raw.rstrip().endswith('END')
        logs[stem] = read_blocks(raw)
        sources[stem] = (P/(stem+'.sing')).read_text()
        for name in [stem+'.sing', stem+'_stdout.txt', 'run_'+stem+'.py']:
            provenance[name] = sha(name)
    x, y, z, w, T = S.symbols('x y z w T')
    a, b, c, d, e, f = S.symbols('a b c d e f')
    env = dict(zip('abcdefxyzw', [a, b, c, d, e, f, x, y, z, w]))
    data = json.loads((P/'binary_factor_moment_inputs.json').read_text())
    binary = [x*w+y*z, x*z**3+y**3*w-3*y*y*z*z, y**3*z**3]
    raw_moments = [S.Poly(S.expand((x+y*T+z*T*T+w*T**3)**(2*m)), T).coeff_monomial(T**(3*m)) for m in range(1, 4)]
    assert list(map(str, raw_moments)) == data['binary_moments_without_factorials']
    G = S.groebner(binary, x, y, z, w, domain=S.QQ)
    GM = S.groebner(raw_moments, x, y, z, w, domain=S.QQ)
    assert all(G.reduce(v)[1] == 0 for v in raw_moments)
    assert all(GM.reduce(v)[1] == 0 for v in binary)
    expected_binary_source = 'ideal I='+','.join(str(v).replace('**', '^') for v in binary)+';'
    assert expected_binary_source in sources['binary_cubic_coefficients']
    assert expected_binary_source in sources['binary_cubic_primary_audit']
    mapping = {x: d*e*f, y: a*e*f+d*b*f+d*e*c,
               z: a*b*f+a*e*c+d*b*c, w: a*b*c}
    factors = [S.expand(v.subs(mapping, simultaneous=True)) for v in binary]
    assert factors == [S.sympify(v, locals=env) for v in data['factor_generators']]
    # Balanced monomials have generators s=xw,t=yz,r=xz^3,q=y^3w.
    s, t, u, r, q = S.symbols('s t u r q')
    reducer = S.groebner([r*r-u*r+s*t**3], r, domain=S.QQ.poly_ring(s, t, u))
    invariant_moments = []
    for moment in raw_moments:
        expr = 0
        for (aa, bb, cc, dd), coefficient in S.Poly(moment, x, y, z, w).terms():
            common = min(aa, dd)
            aa -= common
            dd -= common
            if aa:
                assert cc == 3*aa+bb
                monomial = s**common*r**aa*t**bb
            elif dd:
                assert bb == cc+3*dd
                monomial = s**common*q**dd*t**cc
            else:
                assert bb == cc
                monomial = s**common*t**bb
            expr += coefficient*monomial
        rem = S.expand(reducer.reduce(S.expand(expr.subs(q, u-r)))[1])
        assert r not in rem.free_symbols
        invariant_moments.append(rem)
    assert [S.expand(v.subs(s, -t).subs(u, 3*t*t)) for v in invariant_moments] == [0, 0, -90*t**3]
    # The five-component auxiliary decomposition is checked, including primaryness.
    four = logs['binary_cubic_coefficients']
    assert one(four, 'PRIMARY_COUNT') == '5'
    assert one(four, 'INTERSECTION_FORWARD') == one(four, 'INTERSECTION_REVERSE') == '0'
    four_components = []
    for rawq, rawp in zip(four['PRIMARY_IDEAL'], four['ASSOCIATED_PRIME']):
        qgens = [parse_polynomial(v, env) for v in rawq.split(',')]
        pgens = [parse_polynomial(v, env) for v in rawp.split(',')]
        assert all(S.Poly(v, x, y, z, w).is_homogeneous for v in qgens)
        four_components.append({'generators': list(map(str, qgens)), 'prime': list(map(str, pgens))})
    assert [set(v['prime']) for v in four_components] == [
        {'w', 'z'}, {'y', 'x'}, {'w', 'z', 'y'}, {'z', 'y', 'x'}, {'w', 'z', 'y', 'x'}]
    for index, rawq in enumerate(four['PRIMARY_IDEAL'], 1):
        assert 'ideal Q'+str(index)+'='+rawq+';' in sources['binary_cubic_primary_audit']
    audit = logs['binary_cubic_primary_audit']
    for label in ['AXIS_X_RADICAL_FORWARD', 'AXIS_X_RADICAL_REVERSE', 'AXIS_X_REGULAR',
                  'AXIS_W_RADICAL_FORWARD', 'AXIS_W_RADICAL_REVERSE', 'AXIS_W_REGULAR', 'ORIGIN_RADICAL_REVERSE']:
        assert one(audit, label) == '0'
    assert [list(map(int, v.split())) for v in audit['IRREDUNDANT_COUNT']] == [[1, 1], [2, 1], [3, 3], [4, 2], [5, 1]]
    def h(coefficients):
        return sum(v*T**i for i, v in enumerate(coefficients))
    def hf(poly):
        assert S.denom(poly) == 1
        values = [int(S.expand(poly).coeff(T, i)) for i in range(S.degree(poly, T)+1)]
        assert all(v >= 0 for v in values)
        return values
    four_torsion = S.cancel((h(vector(audit, 'FULL_HILBERT'))-h(vector(audit, 'ORIGIN_SAT_HILBERT')))/(1-T)**4)
    assert four_torsion == T**8
    assert [list(map(int, v.split())) for v in audit['RADICAL_POWER']][-2:] == [[5, 36, 2], [6, 49, 0]]
    # Complete Smith matrices of the normalized zero-factor family.
    zz = S.symbols('zz')
    norm = [S.expand(v.subs({b: 1, c: 1, e: 1, f: zz})) for v in factors]
    D = zz**6+4*zz**5+14*zz**4+18*zz**3+14*zz**2+4*zz+1
    assert S.gcd(D, S.diff(D, zz)) == 1 and D.subs(zz, 0) == 1
    assert S.expand(zz**6*D.subs(zz, 1/zz)-D) == 0
    assert S.factor(S.resultant(norm[0].subs(d, 1), norm[1].subs(d, 1), a)-zz**2*D**2) == 0
    smith = []
    for degree in range(8):
        columns = []
        for generator in norm:
            deg = S.Poly(generator, a, d).total_degree()
            for i in range(degree-deg+1):
                poly = S.Poly(S.expand(generator*a**i*d**(degree-deg-i)), a, d)
                columns.append([poly.coeff_monomial(a**j*d**(degree-j)) for j in range(degree+1)])
        matrix = S.Matrix.hstack(*(S.Matrix(col) for col in columns)) if columns else S.zeros(degree+1, 0)
        snf = smith_normal_form(matrix, domain=S.QQ.poly_ring(zz)) if columns else matrix
        diagonal = [S.factor(snf[i, i]) for i in range(min(snf.shape)) if snf[i, i] != 0]
        stripped = []
        for value in diagonal:
            poly = S.Poly(value, zz)
            while poly.eval(0) == 0:
                poly = S.quo(poly, S.Poly(zz, zz))
            stripped.append(poly.monic().as_expr())
        assert all(v in [S.Integer(1), D] for v in stripped)
        assert stripped.count(D) == {4: 1, 5: 2, 6: 1}.get(degree, 0)
        assert degree+1-len(diagonal) == [1, 2, 2, 2, 1, 0, 0, 0][degree]
        smith.append({'A_degree': degree, 'matrix_shape': list(matrix.shape),
                      'free_rank': degree+1-len(diagonal), 'diagonal': list(map(str, diagonal))})
    # The nine minimal coordinate primes form the edge set of a triangular prism.
    variables = 'abcdef'
    edges = [set(v) for v in ['ad', 'be', 'cf', 'ab', 'ac', 'bc', 'de', 'df', 'ef']]
    covers = []
    for k in range(7):
        for subset in combinations(variables, k):
            support = set(subset)
            if all(support & edge for edge in edges) and not any(set(v) <= support for v in covers):
                covers.append(''.join(subset))
    assert set(covers) == {'bcdf', 'bcde', 'acef', 'acde', 'abef', 'abdf'}
    radical_generators = [S.prod(env[v] for v in cover) for cover in covers]
    for edge in edges:
        assert all(S.expand(v.subs({env[name]: 0 for name in edge})) == 0 for v in factors)
    facets = {''.join(v for v in variables if v not in edge) for edge in edges}
    faces = {''.join(v) for F in facets for k in range(5) for v in combinations(F, k)}
    torus_weights = S.Matrix([[1, 0, 0, 1, 0, 0], [0, 1, 0, 0, 1, 0],
                             [0, 0, 1, 0, 0, 1], [1, 1, 1, -1, -1, -1]])
    for F in faces:
        if len(F) <= 3:
            assert torus_weights[:, [variables.index(v) for v in F]].rank() == len(F)
    symmetries = list(symmetry_maps())
    def orbit(F):
        return {''.join(sorted(mapping[v] for v in F)) for mapping in symmetries}
    candidate_orbits = {frozenset(orbit(F)) for F in faces}
    assert len(candidate_orbits) == 10
    assert {frozenset(orbit(params.replace(',', ''))) for _, params, _ in COORDINATE_CASES} == candidate_orbits
    soc = logs['binary_factor_socles']
    expected_socles = [1, 1, 0, 1, 0, 5, 2, 0, 2, 4]
    coordinate_associated, coordinate_tests = set(), []
    factor_source = 'ideal I='+','.join(str(v).replace('**', '^') for v in factors)+';'
    assert sources['binary_factor_socles'].count(factor_source) == 10
    for (label, params, remaining), expected in zip(COORDINATE_CASES, expected_socles):
        field = '(0,'+params+')' if params else '0'
        assert 'ring r'+label+'='+field+',('+remaining+'),dp;' in sources['binary_factor_socles']
        assert one(soc, label+'_FORWARD') == one(soc, label+'_REVERSE') == '0'
        polynomial = S.cancel((h(vector(soc, label+'_I_HILBERT'))-h(vector(soc, label+'_K_HILBERT')))/(1-T)**len(remaining.split(',')))
        assert S.denom(polynomial) == 1 and polynomial.subs(T, 1) == expected
        if expected:
            coordinate_associated |= orbit(params.replace(',', ''))
        coordinate_tests.append({'free_face': params, 'socle_dimension': expected,
                                 'filtered_Hilbert': str(S.expand(polynomial))})
    assert len(coordinate_associated) == 27
    assert one(soc, 'MIN_ZERO_LENGTH') == '8' and one(soc, 'MIN_PAIR_LENGTH') == '2'
    # The five embedded coordinate orbits admit positive gradings after normalization.
    weighted = logs['binary_factor_weighted_geometry']
    torsion = logs['binary_factor_torsion']
    component_orbits, embedded_components = [], {}
    for label, free, remaining, weights, bound, length in WEIGHTED_CASES:
        normalized = [S.expand(v.subs({env[name]: 1 for name in free})) for v in factors]
        expected_ring = 'ring r'+label+'=0,('+','.join(remaining)+'),wp('+','.join(map(str, weights))+');'
        expected_ideal = 'ideal I='+','.join(str(v).replace('**', '^') for v in normalized)+';'
        for stem in ['binary_factor_weighted_geometry', 'binary_factor_torsion']:
            source = sources[stem]
            assert expected_ring in source and expected_ideal in source
            assert 'intvec weights='+','.join(map(str, weights))+';' in source
        case_block = sources['binary_factor_torsion'].split(expected_ring)[1].split('ring r')[0]
        assert 'ideal Sat=1;' in case_block
        for name in remaining:
            assert 'sat(G,ideal('+name+'))' in case_block
        for generator in normalized:
            poly = S.Poly(generator, *(env[name] for name in remaining))
            assert len({sum(n*w for n, w in zip(exp, weights)) for exp, _ in poly.terms()}) == 1
        denominator = S.prod(1-T**weight for weight in weights)
        hs = S.cancel((h(vector(weighted, label+'_I_HILBERT'))-h(vector(weighted, label+'_K_HILBERT')))/denominator)
        ht = S.cancel((h(vector(torsion, label+'_I_HILBERT'))-h(vector(torsion, label+'_SAT_HILBERT')))/denominator)
        assert S.degree(hs, T) == bound and S.degree(ht, T) == bound
        assert ht.subs(T, 1) == length
        assert one(weighted, label+'_FORWARD') == one(weighted, label+'_REVERSE') == '0'
        powers = {name: bound//weight+1 for name, weight in zip(remaining, weights)}
        assert all(powers[name]*weight > bound for name, weight in zip(remaining, weights))
        for mapping in symmetries:
            F = ''.join(sorted(mapping[name] for name in free))
            exponent_map = dict(sorted((mapping[name], exponent) for name, exponent in powers.items()))
            if F in embedded_components:
                assert embedded_components[F] == exponent_map
            embedded_components[F] = exponent_map
        component_orbits.append({'free_face': free, 'remaining_variables': remaining,
                                 'weights': weights, 'torsion_Hilbert_function': hf(ht),
                                 'intrinsic_torsion_length': length, 'top_weight': bound,
                                 'socle_Hilbert': str(S.expand(hs)), 'powers': powers,
                                 'orbit_size': len(orbit(free))})
    assert len(embedded_components) == 18
    assert coordinate_associated == facets | set(embedded_components)
    # Complete radical-power upper test and an explicit surviving eighth power.
    powers = [list(map(int, value.split())) for value in weighted['RADICAL_POWER']]
    exponent_set = {tuple(int(S.Poly(v, a, b, c, d, e, f).monoms()[0][i]) for i in range(6)) for v in radical_generators}
    current = {(0,)*6}
    for index, count, remainder_count in powers:
        current = {tuple(x+y for x, y in zip(left, right)) for left in current for right in exponent_set}
        assert count == len(current)
        if index < 9:
            assert remainder_count > 0
        else:
            assert index == 9 and remainder_count == 0
    witness = [parse_polynomial(v, env) for v in weighted['LOWER_WITNESS'][-1].splitlines()]
    assert witness[0] == (a*c*e*f)**8 and witness[1] != 0
    full_hilbert = vector(soc, 'ORIGIN_I_HILBERT')
    assert S.cancel(h(full_hilbert)/(1-T)**2).subs(T, 1) == 36
    assert 3*8+6*2 == 36
    failed = json.loads((P/'binary_three_factor_run.json').read_text())
    assert failed['terminal'] and failed['timed_out'] and not failed['completed']
    assert failed['input_sha256'] == sha('binary_three_factor.sing')
    assert failed['log_sha256'] == sha('binary_three_factor_stdout.txt')
    provenance['binary_factor_moment_inputs.json'] = sha('binary_factor_moment_inputs.json')
    cert = {'claim': 'G2CUBICSCHEME20261114', 'field': 'C, characteristic zero',
            'binary_cubic': {'all_order_generators': list(map(str, binary)),
                             'balanced_ring': 'Q[s,t,r,q]/(r*q-s*t^3); reversal invariants Q[s,t,u=r+q]',
                             'first_three_invariant_moments': list(map(str, invariant_moments)),
                             'primary_components': four_components, 'nilradical_index': 6,
                             'degree': 2, 'intrinsic_origin_torsion_Hilbert': 'T^8'},
            'factor_generators': list(map(str, factors)), 'factor_coefficient_map': {'x': str(d*e*f), 'y': str(a*e*f+d*b*f+d*e*c), 'z': str(a*b*f+a*e*c+d*b*c), 'w': str(a*b*c)},
            'radical': {'minimal_coordinate_primes': ['('+','.join(sorted(edge))+')' for edge in edges],
                        'squarefree_generators': list(map(str, radical_generators))},
            'normalized_zero_factor': {'parameter': 'zz=b*f/(c*e)', 'D': str(D), 'Smith': smith,
                                       'free_rank': 8, 'torsion_ranks_by_A_degree': {'4': 1, '5': 2, '6': 1},
                                       'intrinsic_torsion_per_root': 4, 'chosen_root_component_length': 12,
                                       'decomposition': 'I:D^infinity intersect all(I+(zz-lambda)), for the six distinct roots lambda of D; zz is inverted'},
            'associated_primes': {'coordinate_free_faces': sorted(coordinate_associated),
                                  'coordinate_count': 27, 'noncoordinate_count': 18, 'total_over_C': 45,
                                  'noncoordinate_family': '(a,d,b*f-lambda*c*e), lambda a root of D, and the two factor-label permutations'},
            'coordinate_candidate_tests': coordinate_tests, 'embedded_coordinate_orbits': component_orbits,
            'primary_construction': {'minimal_zero_factor': 'I:(H*Delta)^infinity, H=b*c*e*f, Delta=homogenization of D(b*f/(c*e))',
                                     'minimal_pair': 'I:(c*d*e*f)^infinity for (a,b), and factor/color permutations',
                                     'noncoordinate': '(I+(b*f-lambda*c*e)):(b*c*e*f)^infinity and permutations',
                                     'embedded_coordinate': [{'free_face': F, 'powers': powers} for F, powers in sorted(embedded_components.items())],
                                     'equality': 'The intersection of all45 explicitly described primary ideals is I, by the written graded-separation and minimal-associated-prime argument; no expanded45-way CAS replay is claimed.'},
            'factor_scheme': {'minimal_primes': 9, 'degree': 36, 'reduced_degree': 9,
                              'generic_minimal_lengths': {'zero_factor': 8, 'isotropic_pair': 2},
                              'nilradical_index': 9, 'power_tests': powers,
                              'eighth_power_lower_witness': list(map(str, witness)),
                              'origin_torsion_length': 13, 'origin_torsion_degrees': [28, 29, 30],
                              'origin_torsion_Hilbert': '3*T^28+6*T^29+4*T^30'},
            'provenance': provenance,
            'retained_failure': 'The generic six-variable primdecGTZ run is a terminal50-second timeout. Its standard-basis/Hilbert checkpoint is retained; no full decomposition is attributed to it.',
            'assurance': 'Internal constructive written proof and exact symbolic/CAS inputs. The all-order tail and45-way equality use explicit proofs, not finite extrapolation. No independent replay, external review or formal verification. N=3,4 coefficient schemes, affine factors and arbitrary cubic GMC remain separate.'}
    path = P/'binary_factor_scheme_certificate.json'
    path.write_text(json.dumps(cert, indent=2)+'\n')
    receipt = {'claim': cert['claim'], 'completed': True, 'primary_components_over_C': 45,
               'minimal_primes': 9, 'noncoordinate_embedded_primes': 18,
               'nilradical_index': 9, 'degree': 36, 'intrinsic_origin_torsion_length': 13,
               'checker_sha256': sha(Path(__file__).name), 'certificate_sha256': sha(path.name),
               'audit_existing_only': args.audit_existing, 'elapsed_seconds': time.monotonic()-started,
               'assurance': cert['assurance']}
    (P/'binary_factor_scheme_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
