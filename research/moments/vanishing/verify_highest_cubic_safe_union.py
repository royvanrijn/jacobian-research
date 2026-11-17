"""Verify the complete-chart 31-generator safe-union target ideal.

This proves a presentation of known safe families, not their exhaustion
of any moment zero locus. All31 moment radical memberships remain open.
Elimination/intersection have an explicit exact QQ Singular dependency;
the generator-change identities and family substitutions are multiplied
independently in SymPy.
"""
from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess
import time

import sympy as S

import highest_cubic_chart as H
from highest_cubic_fibre import exact_moment

HERE = Path(__file__).resolve().parent
VARIABLES = H.VARIABLES
LOCALS = {str(v): v for v in VARIABLES}


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def polynomial(text):
    assert re.fullmatch(r'[A-Za-z0-9_+*/^() .\-]+', text), text
    value = S.sympify(text.replace('^', '**'), locals=LOCALS)
    assert value.free_symbols <= set(VARIABLES)
    return S.Poly(value, *VARIABLES, domain=S.QQ)


def primitive_singular(f, variables):
    value = S.Poly(f, *variables, domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr()
    return str(value).replace('**', '^').replace(' ', '')


def extract(text, begin, end, prefix):
    block = text.split(begin+'\n')[1].split(end)[0]
    return [polynomial(line.split('=', 1)[1]) for line in block.splitlines()
            if line.startswith(prefix+'[')]


def source_text():
    s = S.Symbol('s')
    image = {
        H.LW: -H.LZ*s*s/2, H.LT: -H.LZ*s,
        H.BM2: -s**3*H.B1/2-3*s**4*H.B2/4,
        H.BM1: -3*s*s*H.B1/2-2*s**3*H.B2, H.B0: -s*H.B1-s*s*H.B2,
        H.CM2: -S.Rational(5, 2)*H.C1*s**3-S.Rational(9, 2)*s**5,
        H.CM1: -3*H.C1*s*s-S.Rational(9, 2)*s**4,
    }
    equations = [s**3+H.C1*s+H.C0]+[x-image[x] for x in H.DEPENDENT]
    rod = [H.CM2, H.CM1, H.C0, H.C1, H.BM2, H.BM1, H.B0,
           H.LW+H.B1*H.B1/2, H.LT-H.B1*H.B2/2, H.LZ-H.B2*H.B2/4]
    names = ','.join(map(str, VARIABLES))
    text = ('ring R=0,(s,'+names+'),dp;\noption(redSB);\nideal I='
            +','.join(primitive_singular(f, (s, *VARIABLES)) for f in equations)
            +';\nideal E=eliminate(I,s);\nring T=0,('+names+'),dp;\n'
            +'option(redSB);map phi=R,0,'+names+';ideal J=phi(E);\n'
            +'print("NULLCONE_BEGIN");J;print("NULLCONE_END");ideal ROD='
            +','.join(primitive_singular(f, (s, *VARIABLES)) for f in rod)
            +';\nideal K=intersect(J,ROD);print("SAFE_UNION_BEGIN");K;'
            +'print("SAFE_UNION_END");print("NULLCONE_DIM");print(dim(std(J)));'
            +'print("SAFE_UNION_DIM");print(dim(std(K)));quit;\n')
    return text, image, rod


def replay(source_name, log_name):
    result = subprocess.run(['Singular', '-q', str(HERE/source_name)],
                            capture_output=True, text=True, timeout=35)
    assert result.returncode == 0 and not result.stderr
    assert '?' not in result.stdout and 'error' not in result.stdout.lower()
    assert result.stdout == (HERE/log_name).read_text()


def check_vectors(polynomials, data):
    assert data['variables'] == list(map(str, VARIABLES))
    assert data['weights'] == list(H.WEIGHTS)
    assert len(polynomials) == len(data['generators'])
    for f, row in zip(polynomials, data['generators']):
        assert f == polynomial(row['polynomial'])
        assert {sum(e*w for e, w in zip(pp, H.WEIGHTS)) for pp, c in f.terms()} == {row['weight']}
        assert f.total_degree() == row['ordinary_degree']
        assert len(f.terms()) == row['terms']


def power_obstruction():
    """An exact power lower bound, not a new support-classification search."""
    a, b = S.symbols('A B')
    moments = []
    for n in (5, 10):
        raw = sum(S.Rational(c.numerator, c.denominator)*a**pp[2]*b**pp[5]
                  for pp, c in exact_moment(n).items()
                  if not any(pp[j] for j in (0, 1, 3, 4, 6)))
        moments.append(S.Poly(raw, a, b).clear_denoms()[1].primitive()[1].as_expr())
    assert moments[0] == 5*a**3+a*b*b
    assert moments[1] == 8091*a**6+4495*a**4*b*b+465*a*a*b**4+5*b**6
    g = [5*a**3+a*b*b, 2759*a*a*b**4-125*b**6, a*b**6, b**8]
    quotient, remainder = S.div(moments[1], g[0], a, b)
    scale = S.Poly(remainder, a, b).LC()/2759
    assert S.expand(remainder-scale*g[1]) == 0
    rows = [[S.Integer(1), S.Integer(0)], [-quotient/scale, 1/scale]]
    rows.append([S.expand((b**4*rows[0][j]-5*a*rows[1][j]/2759)*2759/3384)
                 for j in range(2)])
    rows.append([S.expand((b*b*rows[1][j]-2759*a*rows[2][j])/(-125))
                 for j in range(2)])
    for f, row in zip(g, rows):
        assert S.expand(sum(c*m for c, m in zip(row, moments))-f) == 0
    assert S.expand(quotient*g[0]+scale*g[1]-moments[1]) == 0
    # Verify Buchberger's criterion on the supplied four polynomials by
    # ordinary exact division, rather than trusting a Groebner finder.
    for i, left in enumerate(g):
        for right in g[:i]:
            lp, rp = S.Poly(left, a, b), S.Poly(right, a, b)
            le, re = lp.monoms()[0], rp.monoms()[0]
            lcm = tuple(max(x, y) for x, y in zip(le, re))
            lm = a**(lcm[0]-le[0])*b**(lcm[1]-le[1])/lp.LC()
            rm = a**(lcm[0]-re[0])*b**(lcm[1]-re[1])/rp.LC()
            assert S.reduced(S.expand(lm*left-rm*right), g, a, b)[1] == 0
    assert S.reduced(a**6, g, a, b)[1] == 5*b**6/2759
    assert S.reduced(a**7, g, a, b)[1] == 0
    leading = [S.Poly(f, a, b).monoms()[0] for f in g]
    footprint = [(i, j) for i in range(3) for j in range(8)
                 if not any(i >= x and j >= y for x, y in leading)]
    assert len(footprint) == 18 and max(i+j for i, j in footprint) == 7
    # On F=z^3+(A+B*t)*w^2, circular balance forces n=5s and
    # coefficient degree3s. Thus moments beyond10 have degree>=9 and
    # already belong to this Artinian homogeneous ideal, whose top is7.
    data = {
        'status': 'PASS, EXACT_POWER_OBSTRUCTION_ONLY',
        'phase': 'z^3+(A+B*t)*w^2', 'variables': ['A', 'B'],
        'nonzero_moment_indices': 'Multiples of5 only; coefficient degree3s at n=5s',
        'full_selected_moment_ideal_generators': [str(x) for x in moments],
        'supplied_verified_standard_basis': [str(x) for x in g],
        'standard_basis_in_moment_generators': [[str(c) for c in row] for row in rows],
        'length': 18, 'top_ordinary_degree': 7,
        'A_nilpotence_index': 7, 'A6_normal_form': '5*B^6/2759',
        'global_target1_restriction': '-6*A',
        'consequence': 'The first global safe-union target has no power e<=6 in the FULL ordinary moment ideal. Its radical membership and any power e>=7 remain open.',
    }
    (HERE/'highest_cubic_safe_union_power_obstruction.json').write_text(json.dumps(data, indent=2)+'\n')
    return data


def run(replay_cas=True):
    started = time.monotonic()
    expected, image, rod = source_text()
    assert (HERE/'highest_cubic_safe_union_input.sing').read_text() == expected
    if replay_cas:
        replay('highest_cubic_safe_union_input.sing', 'highest_cubic_safe_union_output.log')
        replay('highest_cubic_safe_union_reduce.sing', 'highest_cubic_safe_union_reduce.log')
    text = (HERE/'highest_cubic_safe_union_output.log').read_text()
    assert text.rstrip().endswith('NULLCONE_DIM\n5\nSAFE_UNION_DIM\n5')
    nullcone = extract(text, 'NULLCONE_BEGIN', 'NULLCONE_END', 'J')
    raw_union = extract(text, 'SAFE_UNION_BEGIN', 'SAFE_UNION_END', 'K')
    reduced = extract((HERE/'highest_cubic_safe_union_reduce.log').read_text(),
                      'MINIMAL_BEGIN', 'MINIMAL_END', 'M')
    assert (len(nullcone), len(raw_union), len(reduced)) == (67, 77, 31)
    for kind, values in [('nullcone', nullcone), ('union', raw_union)]:
        check_vectors(values, json.loads((HERE/f'highest_cubic_{kind}_targets.json').read_text()))
    gate = json.loads((HERE/'highest_cubic_safe_union_gate.json').read_text())
    check_vectors(reduced, gate)
    assert gate['verified_global_moment_memberships'] == 0

    # Independent direct evaluation on every parameter of both families.
    s = S.Symbol('s')
    null_subs = dict(image, **{})
    null_subs[H.C0] = -s*H.C1-s**3
    rod_subs = {H.LW: -H.B1*H.B1/2, H.LT: H.B1*H.B2/2, H.LZ: H.B2*H.B2/4,
                H.BM2: 0, H.BM1: 0, H.B0: 0,
                H.CM2: 0, H.CM1: 0, H.C0: 0, H.C1: 0}
    for f in nullcone:
        assert S.expand(f.as_expr().subs(null_subs, simultaneous=True)) == 0
    for f in raw_union+reduced:
        assert S.expand(f.as_expr().subs(null_subs, simultaneous=True)) == 0
        assert S.expand(f.as_expr().subs(rod_subs, simultaneous=True)) == 0

    # Complete explicit generator-change identities in BOTH directions.
    witness = (HERE/'highest_cubic_safe_union_generator_witness.log').read_text()
    witness_entries = {}
    for label, source, target in [('A', raw_union, reduced), ('B', reduced, raw_union)]:
        sums = [S.Poly(0, *VARIABLES, domain=S.QQ) for _ in target]
        entries = []
        for line in witness.split(label+'_BEGIN\n')[1].split(label+'_END')[0].splitlines():
            i, j, coefficient = line.split('|')
            i, j = int(i)-1, int(j)-1
            assert 0 <= i < len(source) and 0 <= j < len(target)
            multiplier = polynomial(coefficient)
            sums[j] += source[i]*multiplier
            entries.append([i+1, j+1, str(multiplier.as_expr())])
        assert sums == target
        witness_entries[label] = entries
    (HERE/'highest_cubic_safe_union_generator_change.json').write_text(
        json.dumps({'status': 'PASS, EXACT_POLYNOMIAL_IDENTITIES',
                    'A_convention': '77-generator row times A equals31-generator row',
                    'B_convention': '31-generator row times B equals77-generator row',
                    'entries': witness_entries}, indent=2)+'\n')
    histogram = Counter(row['weight'] for row in gate['generators'])
    assert min(histogram) == 5 and max(histogram) == 10
    assert max(row['ordinary_degree'] for row in gate['generators']) == 4
    obstruction = power_obstruction()
    assert S.expand(reduced[0].as_expr().subs({v: (S.Symbol('A') if v == H.BM2
                                                  else S.Symbol('B') if v == H.CM2 else 0)
                                             for v in VARIABLES})+6*S.Symbol('A')) == 0
    receipt = {
        'claim': 'SIC2B33GLOBALSAFEGATE20261117', 'status': 'PASS',
        'checker_sha256': digest(Path(__file__).name),
        'input_sha256': digest('highest_cubic_safe_union_input.sing'),
        'elimination_and_intersection_output_sha256': digest('highest_cubic_safe_union_output.log'),
        'gate_sha256': digest('highest_cubic_safe_union_gate.json'),
        'generator_change_sha256': digest('highest_cubic_safe_union_generator_change.json'),
        'complete_nullcone_generators': 67, 'raw_union_generators': 77,
        'preferred_union_generators': 31, 'maximum_ordinary_degree': 4,
        'positive_weight_histogram': [[k, v] for k, v in sorted(histogram.items())],
        'safe_union_dimension': 5,
        'complete_chart_boundary': 'All H3-nonzero phases after seven fixed rotations; no D or base localization; lower H3=0 handled separately.',
        'direct_generator_change_identity_counts': {key: len(value) for key, value in witness_entries.items()},
        'all_parameter_family_substitutions': 'PASS',
        'first_target_power_lower_bound': obstruction['A_nilpotence_index'],
        'power_obstruction_sha256': digest('highest_cubic_safe_union_power_obstruction.json'),
        'Singular_dependency': 'Exact QQ elimination and ideal intersection; generator-change identities independently multiplied',
        'CAS_replayed_this_run': replay_cas,
        'verified_global_moment_memberships': 0, 'global_safe_exhaustion': 'OPEN',
        'radical_exponents': 'UNKNOWN', 'sufficient_prefix': 'UNKNOWN',
        'full_Image33': 'OPEN',
        'proof_boundary': 'This presents the ideal of a known safe union and a sufficient global radical gate. It neither proves those memberships nor excludes additional safe pure components. No independent specialist/external/formal verification.',
        'elapsed_seconds': round(time.monotonic()-started, 3),
    }
    (HERE/'highest_cubic_safe_union_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--retained-cas-output', action='store_true')
    args = parser.parse_args()
    run(not args.retained_cas_output)
