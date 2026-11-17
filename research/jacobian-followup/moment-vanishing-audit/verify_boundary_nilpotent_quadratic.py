"""Replay the complete zero-linear/rank-two nilpotent quadratic boundary.

The largest chart retains an explicit exact-QQ Singular dependency. Its
direct unit-multiplier export is not claimed complete. Other chart unit
witnesses are checked by multiplication, with no unit-ideal oracle.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import time

import sympy as S

from boundary_nilpotent_quadratic import (
    CHARGES, COEFFICIENTS, primitive_expression, raw_expression,
)
from conformal_euler_cubic import harmonic_basis, laplace, q, sphere_average, z, w, t

HERE = Path(__file__).resolve().parent


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def run():
    started = time.monotonic()
    cert = json.loads((HERE/'boundary_nilpotent_quadratic_certificate.json').read_text())
    m3, m2, m1, y, x, high2, high3 = COEFFICIENTS
    basis = harmonic_basis()[3]
    assert tuple(basis) == CHARGES
    assert all(laplace(f) == 0 for f in basis.values())
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    B = S.Matrix([[0, 0, 0], [S.Rational(1, 2), 0, 0], [0, S.Rational(1, 2), 0]])
    spatial = S.Matrix([z, t, w])
    assert (spatial.T*metric*B*spatial)[0] == z*t
    assert metric*B == B.T*metric and B.rank() == 2 and B**3 == S.zeros(3)
    assert B**2 != S.zeros(3)

    # A full Jordan-chain Gram normalization, rather than a selected slice.
    c = S.Symbol('c', nonzero=True)
    d, e = S.symbols('d e')
    gram = S.Matrix([[0, 0, c], [0, c, d], [c, d, e]])
    alpha = -d/(2*c)
    beta = -(e+2*alpha*d+alpha**2*c)/(2*c)
    adapter = S.Matrix([[1, alpha, beta], [0, 1, alpha], [0, 0, 1]])
    assert (adapter.T*gram*adapter-c*metric).applyfunc(S.simplify) == S.zeros(3)
    jordan = S.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    assert jordan*adapter == adapter*jordan
    # c!=0 can always be chosen since B^2 is nonzero and self-adjoint.
    # Scale the chain, swap its endpoints and apply a diagonal isometry:
    # 2*t*w becomes z*t. Orientation can be corrected by -Id.
    assert S.expand((2*t*w).subs({w: z/2, z: 2*w}, simultaneous=True)) == z*t
    assert S.expand(q.subs({w: z/2, z: 2*w}, simultaneous=True)-q) == 0

    rho = S.Symbol('rho', nonzero=True)
    transform = {z: rho*z, w: w/rho, t: t}
    assert S.expand((z*t).subs(transform, simultaneous=True)/rho-z*t) == 0
    for j, f in basis.items():
        assert S.expand(f.subs(transform, simultaneous=True)/rho-rho**(j-1)*f) == 0
    assert all(j-1 != 0 for j in [-3, -2, -1, 0])

    numeric = dict(zip(COEFFICIENTS, range(1, 14, 2)))
    F = z*t+sum(numeric[coefficient]*basis[j]
                for j, coefficient in zip(CHARGES, COEFFICIENTS))
    for n in range(1, 7):
        assert S.expand(raw_expression(n).subs(numeric)-sphere_average(F**n)) == 0
    assert primitive_expression(2) == 15*y**2+5*x*m1+4*high2*m2+12*high3*m3
    assert primitive_expression(3) == y*m1+x*m2+2*high2*m3

    expected = {
        'lowest_minus3': {m3: 1, high2: -(y*m1+x*m2)/2,
                          high3: -(15*y*y+5*x*m1-2*(y*m1+x*m2)*m2)/12},
        'lowest_minus2': {m3: 0, m2: 1, x: -y*m1,
                          high2: -(15*y*y-5*y*m1*m1)/4},
        'lowest_minus1': {m3: 0, m2: 0, m1: 1, y: 0, x: 0},
    }
    checks = {}
    for name, substitution in expected.items():
        row = cert['charts'][name]
        variables = sorted(set(COEFFICIENTS)-set(substitution), key=str)
        local = {str(v): v for v in COEFFICIENTS}
        assert row['variables'] == list(map(str, variables))
        assert {local[k]: S.sympify(v, locals=local) for k, v in row['substitution'].items()} == substitution
        for n in [2, 3]:
            assert S.expand(primitive_expression(n).subs(substitution, simultaneous=True)) == 0
        equations = []
        for r in row['moments']:
            f = S.Poly(S.expand(primitive_expression(r['n']).subs(substitution, simultaneous=True)),
                       *variables).clear_denoms()[1].primitive()[1]
            assert f.set_domain(S.QQ) == S.Poly(S.sympify(r['polynomial'], locals=local),
                                               *variables, domain=S.QQ)
            equations.append(f.set_domain(S.QQ))
        if row['unit_multipliers'] is not None:
            multipliers = [S.Poly(S.sympify(s.replace('^', '**'), locals=local),
                                  *variables, domain=S.QQ) for s in row['unit_multipliers']]
            assert len(multipliers) == len(equations)
            value = sum((a*b for a, b in zip(multipliers, equations)),
                        S.Poly(0, *variables, domain=S.QQ))
            assert value.as_expr() == 1
            checks[name] = {'prefix': row['prefix'], 'method': 'direct polynomial unit identity',
                            'multiplier_term_counts': [len(f.terms()) for f in multipliers]}
        else:
            assert name == 'lowest_minus3' and row['prefix'] == 8
            # Keep the actual exact input and output as audit artifacts.
            script = 'ring R=0,('+','.join(map(str, variables))+'),dp;\noption(redSB);\n'
            script += 'ideal I='+','.join(str(f.as_expr()).replace('**', '^') for f in equations)+';\n'
            script += 'ideal G=std(I);\nprint("UNIT_CHECK_BEGIN");\nG;\nprint(dim(G));\nprint(vdim(G));\nprint("UNIT_CHECK_END");\nquit;\n'
            source = HERE/'boundary_nilpotent_largest_chart_verified_input.sing'
            source.write_text(script)
            cas_started = time.monotonic()
            result = subprocess.run(['Singular', '-q', str(source)], capture_output=True,
                                    text=True, timeout=30)
            log = result.stdout+result.stderr
            (HERE/'boundary_nilpotent_largest_chart_verified_output.log').write_text(log)
            assert result.returncode == 0 and '?' not in log, log[-3000:]
            body = log.split('UNIT_CHECK_BEGIN\n')[1].split('UNIT_CHECK_END')[0].strip().splitlines()
            assert body == ['G[1]=1', '-1', '0'], body
            checks[name] = {'prefix': 8, 'method': 'exact QQ Singular standard-basis unit calculation',
                            'elapsed_seconds': round(time.monotonic()-cas_started, 3),
                            'input_sha256': digest(source.name),
                            'output_sha256': digest('boundary_nilpotent_largest_chart_verified_output.log'),
                            'direct_unit_multiplier_export': 'INCOMPLETE; three45-second timeouts retained'}

    # All negative coefficients zero: mu2 alone forces the zero charge zero.
    assert primitive_expression(2).subs({m3: 0, m2: 0, m1: 0}) == 15*y*y
    # With no nonpositive cubic charge, every phase monomial has charge>=1.
    for n in range(1, 9):
        assert raw_expression(n).subs({m3: 0, m2: 0, m1: 0, y: 0}) == 0
    false = z*t+basis[-1]+S.Rational(39, 2)*z**3
    values = [sphere_average(false**n) for n in range(1, 9)]
    assert values[:7] == [0]*7 and values[7] == S.Rational(111856, 9561123)
    false_sub = dict(zip(COEFFICIENTS, [0, 0, 1, 0, 0, 0, S.Rational(39, 2)]))
    assert all(S.expand(raw_expression(n).subs(false_sub)-values[n-1]) == 0 for n in range(1, 9))
    assert sphere_average((z*t+z**3)**2*w**2) == S.Rational(2, 105)
    # Exact equations of the whole closed safe container, with no chart units.
    quadratics = S.symbols('h2_w2 h2_wt h2_t2minuszw h2_zt h2_z2')
    a, b, c, d, e = quadratics
    Bg = S.Matrix([[-c/2, b/2, a], [d/2, c, b/2], [e, d/2, -c/2]])
    H2 = a*w*w+b*w*t+c*(t*t-z*w)+d*z*t+e*z*z
    assert S.expand((spatial.T*metric*Bg*spatial)[0]-H2) == 0
    p = S.expand(S.trace(Bg**2)/2)
    r = S.expand(Bg.det())
    assert (Bg**3-p*Bg-r*S.eye(3)).applyfunc(S.expand) == S.zeros(3)
    container = {
        'coefficient_field': 'QQ', 'linear_harmonic': 'h1_z*z+h1_t*t+h1_w*w',
        'quadratic_harmonic': str(H2), 'cubic_harmonic': 'unrestricted',
        'defining_equations': ['h1_z', 'h1_t', 'h1_w', str(p), str(r)],
        'interpretation': 'H1=0 and B nilpotent; rank2 new theorem, rank1 previous exact theorem, B=0 odd theorem.',
        'scope': 'Every all-order pure characteristic-zero FIELD point of this closed container is full-Image safe. The container need not be pure; it is not asserted to exhaust the pure locus.',
        'full_Image_cutoff': 'm>coordinate_degree Q',
        'nonzero_B_minimal_consecutive_set_theoretic_prefix': 10,
        'B_zero_prefix': 'No new prefix assertion; existing odd theorem supplies safety.',
        'ordinary_ideal_stabilization': 'NOT ASSERTED',
        'assurance': 'Rank2 largest exclusion retains exact QQ Singular dependency; B=0 imports odd theorem assurance. No msolve.',
    }
    (HERE/'boundary_nilpotent_safe_container.json').write_text(json.dumps(container, indent=2)+'\n')
    version = subprocess.run(['Singular', '--version'], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    receipt = {
        'claim': cert['claim'], 'status': 'PASS', 'checker_sha256': digest(Path(__file__).name),
        'model_sha256': digest('boundary_nilpotent_quadratic.py'),
        'certificate_sha256': digest('boundary_nilpotent_quadratic_certificate.json'),
        'safe_container_sha256': digest('boundary_nilpotent_safe_container.json'),
        'complete_stratum': 'H1=0, nonzero rank2 nilpotent self-adjoint quadratic operator; unrestricted H3',
        'classification': 'After H2=z*t normalization, exactly cubic charges1,2,3 are pure.',
        'least_consecutive_set_theoretic_prefix': 8,
        'full_original_Image_cutoff': 'm>coordinate_degree Q',
        'checks': checks, 'Singular_version': version, 'msolve_dependency': False,
        'first7_false_positive_mu8': str(values[7]),
        'ordinary_ideal_stabilization': 'NOT ASSERTED', 'full_Image33': 'OPEN',
        'boundary': 'Largest chart retains the displayed exact QQ Singular dependency. Two other exclusions have directly checked rational unit witnesses. No specialist, external or formal verification.',
        'elapsed_seconds': round(time.monotonic()-started, 3),
    }
    (HERE/'boundary_nilpotent_quadratic_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    run()
