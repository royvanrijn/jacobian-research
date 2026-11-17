"""Replay nilpotent-quadratic contraction-zero safety and three global targets.

This classifies the full pure locus in a CLOSED CONDITIONAL CONTAINER.
It does not prove that every pure phase lies in that container. The
remaining four-variable exclusion has a retained exact QQ Singular
dependency; the smaller chart identities are multiplied directly.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import subprocess
import time

import sympy as S

from conformal_euler_cubic import harmonic_basis, laplace, q, sphere_average, z, w, t
import highest_cubic_chart as H
from nilpotent_contraction_zero import (
    A, B, C, D, E, F, G, VARIABLES, exact_moment, raw_expression,
    primitive, d_nonzero, d_zero_c_nonzero,
)

HERE = Path(__file__).resolve().parent


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def lower(f, parameter):
    return S.expand(f.subs({z: z-parameter*t-parameter**2*w/2,
                           t: t+parameter*w, w: w}, simultaneous=True))


def coefficient(f, monomial):
    return S.Poly(f, z, w, t).coeff_monomial(monomial)


def multiplication(value, data):
    ma = S.Matrix([[S.Rational(c) for c in row] for row in data['multiplication_a']])
    mb = S.Matrix([[S.Rational(c) for c in row] for row in data['multiplication_b']])
    result = S.zeros(17)
    for c, (i, j) in zip(value, data['quotient_basis']):
        if c:
            result += S.Rational(c)*(ma**i)*(mb**j)
    return result


def run(replay_cas=True):
    started = time.monotonic()
    hb = harmonic_basis()
    cubic_coefficients = S.symbols('gm3 gm2 gm1 g0 g1 g2 g3')
    full3 = S.expand(sum(c*hb[3][j] for c, j in zip(cubic_coefficients, range(-3, 4))))
    assert S.expand(laplace(laplace(w*w*full3))
                    +8*cubic_coefficients[4]*w-16*cubic_coefficients[5]*t
                    -48*cubic_coefficients[6]*z) == 0
    assert S.expand(laplace(laplace(z*t*full3))
                    +24*cubic_coefficients[3]*z-16*cubic_coefficients[2]*t
                    -16*cubic_coefficients[1]*w) == 0

    # Rank1: after the contraction equations only cubic charges<=0
    # remain. The third moment forces the positive linear coefficient0.
    lz, lt, lw, gm3, gm2, gm1, g0 = S.symbols('lz lt lw gm3 gm2 gm1 g0')
    rank1 = w*w+lz*z+lt*t+lw*w+sum(v*hb[3][j]
                    for v, j in zip((gm3, gm2, gm1, g0), range(-3, 1)))
    assert S.expand(sphere_average(rank1**3)-2*lz*lz/5) == 0
    central = json.loads((HERE/'nilpotent_contraction_zero_rank1_central_unit.json').read_text())
    r = S.Symbol('r')
    central_phase = r*t+hb[3][0]
    equations = [primitive(sphere_average(central_phase**n), (r,)) for n in (2, 4)]
    assert equations == [S.sympify(f, locals={'r': r}) for f in central['moments']]
    assert S.expand(sum(S.sympify(a, locals={'r': r})*f
                        for a, f in zip(central['multipliers'], equations))) == 1

    # Regular nilpotent: NO restriction on the three linear coefficients.
    phase = A*z+B*t+C*w+z*t+D*w**3+E*hb[3][1]+F*z*z*t+G*z**3
    assert S.expand(raw_expression(2)-(2*A*C+B*B)/3-4*D*G/35) == 0
    assert S.expand(raw_expression(3)-2*B*C/5-4*D*F/105) == 0
    point = dict(zip(VARIABLES, range(1, 8)))
    for n in range(1, 7):
        computed = sum(S.Rational(c.numerator, c.denominator)
                       *math.prod(point[v]**e for v, e in zip(VARIABLES, pp))
                       for pp, c in exact_moment(n).items())
        assert computed == sphere_average(phase.subs(point)**n)
    rho = S.Symbol('rho', nonzero=True)
    rescaled = S.expand(phase.subs({z: rho*z, w: w/rho}, simultaneous=True)/rho)
    assert S.expand(rescaled-phase.subs({B: B/rho, C: C/rho**2,
                    D: D/rho**4, F: rho*F, G: rho**2*G}, simultaneous=True)) == 0

    data = json.loads((HERE/'nilpotent_contraction_zero_charts.json').read_text())
    largest, rod = [], []
    for row in data['d_nonzero']['moments']:
        value = d_nonzero(row['n'])
        assert value == S.sympify(row['primitive'], locals={str(v): v for v in VARIABLES})
        largest.append(value)
    for row in data['d_zero_c_nonzero']['moments']:
        value = d_zero_c_nonzero(row['n'])
        assert value == S.sympify(row['primitive'], locals={str(v): v for v in VARIABLES})
        rod.append(value)
    assert [r['n'] for r in data['d_nonzero']['moments']] == list(range(4, 13))
    assert [r['n'] for r in data['d_zero_c_nonzero']['moments']] == list(range(4, 13))

    # Independently multiply all three Rodrigues-family ideal identities.
    log = (HERE/'nilpotent_contraction_zero_rodrigues_witness.log').read_text()
    sums = [S.Poly(0, E, F, G, domain=S.QQ) for _ in range(3)]
    for line in log.split('WITNESS_BEGIN\n')[1].split('WITNESS_END')[0].splitlines():
        i, j, text = line.split('|')
        multiplier = S.Poly(S.sympify(text.replace('^', '**'),
                            locals={'e': E, 'f': F, 'g': G}), E, F, G, domain=S.QQ)
        sums[int(j)-1] += multiplier*S.Poly(rod[int(i)-1], E, F, G, domain=S.QQ)
    assert [v.as_expr() for v in sums] == [E, F, 2*G+1]
    assert all(S.expand(value.subs({E: 0, F: 0, G: -S.Rational(1, 2)})) == 0 for value in rod)

    # Complete four-variable exclusion, retaining the exact CAS dependency.
    source = ('ring R=0,(a,b,c,e),dp;\noption(redSB);\nideal I='
              +','.join(str(f).replace('**', '^').replace(' ', '') for f in largest)
              +';\nideal K=std(I);print("BASIS_BEGIN");K;print("BASIS_END");'
              +'print("DIM");print(dim(K));print("VDIM");print(vdim(K));quit;\n')
    source_name = 'nilpotent_contraction_zero_lowest_minus3_prefix12.sing'
    log_name = 'nilpotent_contraction_zero_lowest_minus3_prefix12.log'
    assert (HERE/source_name).read_text() == source
    expected = (HERE/log_name).read_text()
    assert expected.strip() == 'BASIS_BEGIN\nK[1]=1\nBASIS_END\nDIM\n-1\nVDIM\n0'
    if replay_cas:
        result = subprocess.run(['Singular', '-q', str(HERE/source_name)],
                                capture_output=True, text=True, timeout=45)
        assert result.returncode == 0 and not result.stderr and result.stdout == expected

    # The three global targets in the COMPLETE highest-cubic chart.
    bvars = (H.BM2, H.BM1, H.B0, H.B1, H.B2)
    quadratic = S.expand(sum(v*hb[2][j] for v, j in zip(bvars, range(-2, 3))))
    operator = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])*S.hessian(quadratic, (z, t, w))/2
    P = S.expand(2*S.trace(operator*operator))     # 4p
    R = S.expand(4*operator.det())                # 4r
    unnormalized_T = -coefficient(laplace(laplace(quadratic*full3)), z)/8
    target_T = S.expand(unnormalized_T.subs({cubic_coefficients[0]: H.eliminated_cminus3(),
                        cubic_coefficients[1]: H.CM2, cubic_coefficients[2]: H.CM1,
                        cubic_coefficients[3]: H.C0, cubic_coefficients[4]: H.C1,
                        cubic_coefficients[5]: 0, cubic_coefficients[6]: 1}, simultaneous=True))
    gate31 = json.loads((HERE/'highest_cubic_safe_union_gate.json').read_text())
    assert target_T == S.sympify(gate31['generators'][0]['polynomial'], locals={str(v): v for v in H.VARIABLES})
    assert P == S.sympify(gate31['generators'][3]['polynomial'], locals={str(v): v for v in H.VARIABLES})
    for target, weight, degree in [(P, 6, 2), (R, 9, 3), (target_T, 5, 2)]:
        polynomial = S.Poly(target, *H.VARIABLES)
        assert polynomial.total_degree() == degree
        assert {sum(e*d for e, d in zip(pp, H.WEIGHTS)) for pp, c in polynomial.terms()} == {weight}
    assert S.expand(operator.charpoly().as_expr().subs({operator.charpoly().gen: S.Symbol('lambda')})
                    -(S.Symbol('lambda')**3-P*S.Symbol('lambda')/4-R/4)) == 0

    # Covariance and a FINITE nine-rotation inference, without silently
    # assuming three good rotations among the earlier seven-chart cover.
    s = S.Symbol('s')
    covariant = laplace(laplace(quadratic*full3))
    assert laplace(laplace(lower(quadratic, s)*lower(full3, s))) == lower(covariant, s)
    assert coefficient(lower(covariant, s), z) == coefficient(covariant, z)
    low3 = [coefficient(lower(hb[3][j], s), w**3) for j in range(-3, 4)]
    assert S.Matrix([[f.subs(s, k) for f in low3] for k in range(7)]).det() != 0
    vector = S.symbols('v_z v_t v_w')
    lowest_vector = coefficient(lower(vector[0]*z+vector[1]*t+vector[2]*w, s), w)
    assert S.expand(lowest_vector-(-vector[0]*s*s/2+vector[1]*s+vector[2])) == 0
    # A degree<=6 nonzero cubic evaluation polynomial has at least THREE
    # nonzero values on0..8. A degree<=2 vector polynomial vanishing at
    # those three values is identically zero.
    assert len(range(9))-6 == 3

    # Local lower bounds use ACTUAL multiplication, not weighted initial
    # multiplication. Amplitude normalization multiplies p,r by units.
    local = json.loads((HERE/'rodrigues_local_length17_certificate/certificate_data.json').read_text())
    targets = json.loads((HERE/'rodrigues_global_seven_target_certificate/global_target_nilindices.json').read_text())
    for name, exponent in [('p', 7), ('r', 5)]:
        matrix = multiplication(targets['normalized_target_vectors'][name], local)
        assert matrix**exponent == S.zeros(17) and matrix**(exponent-1) != S.zeros(17)
    obstruction = json.loads((HERE/'highest_cubic_safe_union_power_obstruction.json').read_text())
    assert obstruction['A_nilpotence_index'] == 7

    program = {
        'claim': 'SIC2B33THREEGATE20261117', 'status': 'CONDITIONAL_CLASSIFICATION_PROVED, ALL_THREE_GLOBAL_MEMBERSHIPS_OPEN',
        'variables': list(map(str, H.VARIABLES)), 'positive_weights': list(H.WEIGHTS),
        'targets': [{'name': name, 'polynomial': str(target), 'degree': degree, 'weight': weight,
                     'known_exponent_lower_bound': exponent}
                    for name, target, degree, weight, exponent in
                    [('4p', P, 2, 6, 7), ('4r', R, 3, 9, 5), ('spin1_highest', target_T, 2, 5, 7)]],
        'unit_test_convention': 'Three SEPARATE1 in I_N+(1-v*f), any finite N and powers; all are unproved.',
        'fixed_power_ansatz': 'Target of weight d at power e uses n<=floor(e*d/3) and cofactor weight e*d-3*n only.',
        'finite_rotation_parameters_for_covariant_inference': list(range(9)),
        'conditional_nonzero_B_prefix': 12,
        'B_zero': 'All-order odd theorem, not a claim that twelve odd moments suffice.',
        'verified_global_memberships': 0,
        'full_Image33': 'OPEN',
    }
    (HERE/'highest_cubic_three_target_gate.json').write_text(json.dumps(program, indent=2)+'\n')
    local_cyclic = json.loads((HERE/'nilpotent_cyclic_slice_verification.json').read_text())
    assert local_cyclic['status'] == 'PASS'
    assert local_cyclic['actual17_algebra_quotient_length'] == 1
    assert local_cyclic['minimal_consecutive_ordinary_local_prefix_on_slice'] == 10
    assert S.Rational(local_cyclic['Rodrigues_raw_Jacobian_determinant']) != 0
    cyclic = {
        'claim': program['claim'], 'source': 'nilpotent_cyclic_slice_prefix14.json.gz',
        'source_sha256': digest('nilpotent_cyclic_slice_prefix14.json.gz'),
        'initial_prefix': 14, 'allow_later_prefixes': True, 'valid_equivalence_for_N_at_least': 12,
        'augmented_ordinary_identity': 'J12+(gminus2,gminus1,g0)=Rodrigues maximal ideal',
        'reduced_local_input': 'nilpotent_cyclic_slice_verification.json',
        'reduced_local_input_sha256': digest('nilpotent_cyclic_slice_verification.json'),
        'targets': [{'polynomial': name, 'degree': 1, 'character_mod2': character}
                    for name, character in [('gminus2', 1), ('gminus1', 0), ('g0', 1)]],
        'equivalence': 'Three ordinary memberships OR three radical memberships iff full normalized cyclic nilpotent branch length1, by the certified reduced local point.',
        'ordinary_witness_degree_bound': 'UNKNOWN', 'verified_unrestricted_memberships': 0,
        'global_branch_length1': 'OPEN', 'normalized_three_branch_length3': 'OPEN',
    }
    (HERE/'nilpotent_cyclic_three_target_program.json').write_text(json.dumps(cyclic, indent=2)+'\n')
    receipt = {
        'claim': program['claim'], 'status': 'PASS', 'checker_sha256': digest(Path(__file__).name),
        'model_sha256': digest('nilpotent_contraction_zero.py'),
        'largest_input_sha256': digest(source_name), 'largest_output_sha256': digest(log_name),
        'chart_input_sha256': digest('nilpotent_contraction_zero_charts.json'),
        'global_gate_sha256': digest('highest_cubic_three_target_gate.json'),
        'cyclic_gate_sha256': digest('nilpotent_cyclic_three_target_program.json'),
        'unrestricted_linear_coefficients': 'ALL THREE RETAINED',
        'ranks_of_nilpotent_B': [0, 1, 2],
        'independent_spatial_moments_checked': list(range(1, 7)),
        'rank1_central_unit': 'PASS, direct exact Bezout multiplication',
        'regular_nilpotent_d0_Rodrigues_ideal': 'PASS, three directly multiplied ideal identities through8',
        'regular_nilpotent_d_nonzero_exclusion': 'PASS, QQ standard basis1 through12',
        'Singular_dependency': 'Largest four-variable exact QQ exclusion; optional independent unit export timed out at40s',
        'msolve_dependency': False, 'CAS_replayed_this_run': replay_cas,
        'conditional_marked_cutoff': 'm>3 coordinate-degree Q; sharp by Rodrigues',
        'nine_rotation_global_three_target_reduction': 'PASS, exact covariant/weight adapters plus polynomial root counts',
        'global_three_memberships': 'ALL OPEN', 'cyclic_three_memberships': 'ALL OPEN',
        'full_Image33': 'OPEN',
        'boundary': 'A complete pure-safe CLOSED CONDITIONAL CONTAINER and a reduced global gate. It does not prove p,r or the contraction vanish at every pure phase. Twelve moments suffice only on the NONZERO nilpotent quadratic part under the contraction condition. B0 imports all-order odd classification. No independent specialist/external/formal verification.',
        'elapsed_seconds': round(time.monotonic()-started, 3),
    }
    (HERE/'nilpotent_contraction_zero_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--retained-cas-output', action='store_true')
    args = parser.parse_args()
    run(not args.retained_cas_output)
