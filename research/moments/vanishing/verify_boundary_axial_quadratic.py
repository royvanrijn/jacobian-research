"""Exact replay of the complete zero-linear semisimple axial exclusion.

The five-variable chart retains an exact QQ Singular dependency. Its
direct polynomial unit witness is unexported. No finite-field unit is
used to infer a characteristic-zero exclusion.
"""
from pathlib import Path
import hashlib
import itertools
import json
import subprocess
import time

import sympy as S

from boundary_axial_quadratic import CHARGES, COEFFICIENTS, primitive_expression, raw_expression
from conformal_euler_cubic import harmonic_basis, laplace, q, sphere_average, z, w, t

HERE = Path(__file__).resolve().parent


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def run():
    started = time.monotonic()
    cert = json.loads((HERE/'boundary_axial_quadratic_certificate.json').read_text())
    a, b, c, y, x, d, e = COEFFICIENTS
    basis = harmonic_basis()[3]
    assert tuple(basis) == CHARGES and all(laplace(f) == 0 for f in basis.values())
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    B = S.diag(-S.Rational(1, 2), 1, -S.Rational(1, 2))
    spatial = S.Matrix([z, t, w])
    assert (spatial.T*metric*B*spatial)[0] == t*t-z*w
    assert metric*B == B.T*metric and B.trace() == 0 and B.rank() == 3
    assert B**2-B/2-S.eye(3)/2 == S.zeros(3)
    rho = S.Symbol('rho', nonzero=True)
    torus = {z: rho*z, t: t, w: w/rho}
    reversal = {z: w, w: z, t: -t}
    reversal_matrix = S.Matrix([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
    assert reversal_matrix.det() == 1
    for transform in [torus, reversal]:
        assert S.expand(q.subs(transform, simultaneous=True)-q) == 0
        assert S.expand((t*t-z*w).subs(transform, simultaneous=True)-(t*t-z*w)) == 0
    for j, f in basis.items():
        assert S.expand(f.subs(torus, simultaneous=True)-rho**j*f) == 0
        assert S.expand(f.subs(reversal, simultaneous=True)-S.Integer(-1)**(j+1)*basis[-j]) == 0

    numeric = dict(zip(COEFFICIENTS, range(1, 14, 2)))
    F = t*t-z*w+sum(numeric[v]*basis[j] for j, v in zip(CHARGES, COEFFICIENTS))
    for n in range(1, 7):
        assert S.expand(raw_expression(n).subs(numeric)-sphere_average(F**n)) == 0
    M2 = 15*y*y+5*x*c+4*d*b+12*e*a+21
    M3 = 4*y*y+x*c-4*e*a+2
    assert primitive_expression(2) == M2 and primitive_expression(3) == M3
    # Highest absolute charge zero, including the zero cubic.
    central = {a: 0, b: 0, c: 0, x: 0, d: 0, e: 0}
    assert S.expand((4*M2-15*M3).subs(central)/54) == 1

    expected = {
        'highest3': ({e: 1, a: -S.Rational(1, 2)-S.Rational(5, 108)*x*c-S.Rational(4, 27)*d*b},
                     27*y*y+8*x*c+4*d*b+27),
        'highest2': ({e: 0, a: 0, d: 1, b: -(5*x*c+54)/16}, 4*y*y+x*c+2),
        'highest1': ({e: 0, a: 0, d: 0, b: 0, x: 1, c: -S.Rational(54, 5)}, 5*y*y-11),
    }
    checks = {}
    for name, (substitution, relation) in expected.items():
        row = cert['charts'][name]
        variables = sorted(set(COEFFICIENTS)-set(substitution), key=str)
        local = {str(v): v for v in COEFFICIENTS}
        assert row['variables'] == list(map(str, variables))
        assert {local[k]: S.sympify(v, locals=local) for k, v in row['substitution'].items()} == substitution
        assert S.sympify(row['quadratic_relation'], locals=local) == relation
        # The retained relation plus substitution is equivalent to M2,M3:
        # all pivots are displayed nonzero rational constants.
        if name == 'highest3':
            raw3 = S.expand(M3.subs({e: 1}))
            solved_a = y*y+x*c/4+S.Rational(1, 2)
            assert S.expand(raw3.subs(a, solved_a)) == 0
            assert S.expand(M2.subs({e: 1, a: solved_a}, simultaneous=True)) == relation
            assert S.expand((solved_a-substitution[a])-relation/27) == 0
        elif name == 'highest2':
            assert S.expand(M3.subs(substitution, simultaneous=True)) == relation
            assert S.expand(M2.subs(substitution, simultaneous=True)-S.Rational(15, 4)*relation) == 0
        else:
            assert S.expand(M2.subs(substitution, simultaneous=True)-3*relation) == 0
            assert S.expand(M3.subs(substitution, simultaneous=True)-S.Rational(4, 5)*relation) == 0
        equations = [S.Poly(relation, *variables, domain=S.QQ)]
        for r in row['moments']:
            substituted = S.expand(primitive_expression(r['n']).subs(substitution, simultaneous=True))
            quotient, remainder = S.div(substituted, relation, y)
            assert S.expand(substituted-quotient*relation-remainder) == 0
            f = S.Poly(remainder, *variables).clear_denoms()[1].primitive()[1].set_domain(S.QQ)
            assert f == S.Poly(S.sympify(r['polynomial'], locals=local), *variables, domain=S.QQ)
            equations.append(f)
        if name == 'highest1':
            assert equations[1].as_expr() == 1
            checks[name] = {'method': 'fourth residual is a nonzero rational constant', 'prefix': 4}
        elif row['unit_multipliers'] is not None:
            ms = [S.Poly(S.sympify(s.replace('^', '**'), locals=local), *variables, domain=S.QQ)
                  for s in row['unit_multipliers']]
            assert len(ms) == len(equations)
            assert sum((m*f for m, f in zip(ms, equations)), S.Poly(0, *variables, domain=S.QQ)).as_expr() == 1
            used = [index+3 for index, m in enumerate(ms) if index and not m.is_zero]
            checks[name] = {'method': 'direct rational polynomial unit identity',
                            'prefix': max(used), 'multiplier_term_counts': [len(m.terms()) for m in ms]}
        else:
            assert name == 'highest3'
            script = 'ring R=0,('+','.join(map(str, variables))+'),dp;\noption(redSB);\n'
            script += 'ideal I='+','.join(str(f.as_expr()).replace('**', '^') for f in equations)+';\n'
            script += 'ideal G=std(I);print("EXACT_UNIT_BEGIN");G;print(dim(G));print(vdim(G));print("EXACT_UNIT_END");quit;\n'
            source = HERE/'boundary_axial_largest_chart_verified_input.sing'
            source.write_text(script)
            cas_started = time.monotonic()
            result = subprocess.run(['Singular', '-q', str(source)], capture_output=True, text=True, timeout=60)
            output = result.stdout+result.stderr
            (HERE/'boundary_axial_largest_chart_verified_output.log').write_text(output)
            assert result.returncode == 0 and '?' not in output, output[-2000:]
            body = output.split('EXACT_UNIT_BEGIN\n')[1].split('EXACT_UNIT_END')[0].strip().splitlines()
            assert body == ['G[1]=1', '-1', '0'], body
            checks[name] = {'method': 'exact QQ Singular standard-basis unit calculation',
                            'prefix': 10, 'elapsed_seconds': round(time.monotonic()-cas_started, 3),
                            'input_sha256': digest(source.name),
                            'output_sha256': digest('boundary_axial_largest_chart_verified_output.log'),
                            'direct_unit_multiplier_export': 'INCOMPLETE;45-second tracked export timeout retained'}

    # The closed harmonic-square quadratic cone, including its null boundary.
    container = json.loads((HERE/'boundary_axial_safe_container.json').read_text())
    aa, bb, cc, dd, ee, lam = S.symbols('h2_w2 h2_wt h2_t2minuszw h2_zt h2_z2 lam')
    Bg = S.Matrix([[-cc/2, bb/2, aa], [dd/2, cc, bb/2], [ee, dd/2, -cc/2]])
    A = Bg+lam*S.eye(3)
    minors = [S.expand(A.extract(i, j).det()) for i in itertools.combinations(range(3), 2)
              for j in itertools.combinations(range(3), 2)]
    assert container['lifted_rank_equations'] == list(map(str, minors))
    principal = sum(A.extract(i, i).det() for i in itertools.combinations(range(3), 2))
    assert S.expand(S.trace(Bg**2)-6*lam*lam+2*principal) == 0
    gb = S.groebner(minors, lam, aa, bb, cc, dd, ee, order='lex')
    elimination = [f.as_expr() for f in gb.polys if not f.as_expr().has(lam)]
    assert container['base_generators'] == ['h1_z', 'h1_t', 'h1_w']+list(map(str, elimination))
    assert len(elimination) == 7 and all(S.Poly(f, aa, bb, cc, dd, ee).total_degree() == 3 for f in elimination)
    # The non-null projection of a square is axial, the null square is rank1.
    nz, nt, nw = S.symbols('nz nt nw')
    normal = S.Matrix([nz, nt, nw])
    norm = (normal.T*metric*normal)[0]
    square_B = normal*normal.T*metric-norm*S.eye(3)/3
    square_coeffs = {aa: nz*nz, bb: 2*nz*nt, cc: S.Rational(2, 3)*(nt*nt-nz*nw),
                     dd: 2*nw*nt, ee: nw*nw}
    assert (Bg.subs(square_coeffs)-square_B).applyfunc(S.expand) == S.zeros(3)
    assert all(S.expand(f.subs(square_coeffs)) == 0 for f in elimination)
    version = subprocess.run(['Singular', '--version'], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    receipt = {
        'claim': cert['claim'], 'status': 'PASS', 'checker_sha256': digest(Path(__file__).name),
        'model_sha256': digest('boundary_axial_quadratic.py'),
        'certificate_sha256': digest('boundary_axial_quadratic_certificate.json'),
        'safe_container_sha256': digest('boundary_axial_safe_container.json'),
        'result': 'No pure H1=0 phase with nonzero semisimple repeated-eigenvalue B, H3 unrestricted.',
        'sufficient_consecutive_prefix': 10, 'minimal_prefix': 'NOT CLAIMED', 'checks': checks,
        'safe_container': 'H1=0, seven cubic equations for the closed harmonic-square quadratic cone; pure points force rank B<=1.',
        'full_original_Image_cutoff_on_pure_container': 'm>coordinate_degree Q',
        'Singular_version': version, 'msolve_dependency': False, 'ordinary_ideal_stabilization': 'NOT ASSERTED',
        'unrestricted_Image33': 'OPEN',
        'boundary': 'Largest exclusion retains exact rational Singular dependency. No inference from finite-field units. Distinct-eigenvalue and nonsemisimple repeated-eigenvalue B are not covered. No independent specialist, external or formal verification.',
        'elapsed_seconds': round(time.monotonic()-started, 3),
    }
    (HERE/'boundary_axial_quadratic_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    run()
