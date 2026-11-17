"""Verify the complete cyclic nilpotent slice and its five-target gate.

This checks exact input, the polynomial moving-frame adapter, the reduced
Rodrigues local ring, and an augmented GLOBAL ideal identity. It does not
prove any of the five unrestricted global memberships or global length1.
"""
from fractions import Fraction
from pathlib import Path
import gzip
import hashlib
import itertools
import json
import math
import time

import sympy as S

from conformal_euler_cubic import harmonic_basis, laplace, q, sphere_average, z, w, t
from nilpotent_cyclic_slice import (
    CHARGES, CUBIC, LINEAR, RODRIGUES, VARIABLES,
    exact_moment, raw_expression, rodrigues_linearization,
)

HERE = Path(__file__).resolve().parent


def digest(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def run():
    started = time.monotonic()
    with gzip.open(HERE/'nilpotent_cyclic_slice_prefix14.json.gz', 'rt') as stream:
        data = json.load(stream)
    assert data['variables'] == list(map(str, VARIABLES))
    assert [r['n'] for r in data['moments']] == list(range(2, 15))
    counts = {}
    parities = [0, 1, 0, 1, 0, 1, 0, 1, 0]
    for row in data['moments']:
        scale = S.Rational(row['raw_equals_primitive_times'])
        fscale = Fraction(int(scale.p), int(scale.q))
        decoded = {tuple(k): fscale*c for k, c in row['terms']}
        assert len(decoded) == len(row['terms']) and decoded == exact_moment(row['n'])
        assert all(sum(k*v for k, v in zip(powers, parities)) % 2 == row['n'] % 2
                   for powers in decoded)
        # The scout primes must not kill any raw/primitive scaling.
        assert all(int(scale.p) % prime and int(scale.q) % prime for prime in [65521, 1073741827])
        value = sum((c*Fraction(-1, 2)**powers[-1] for powers, c in decoded.items()
                     if not any(powers[:-1])), Fraction(0))
        assert value == 0
        counts[str(row['n'])] = len(decoded)

    basis = harmonic_basis()[3]
    a, b = LINEAR
    F0 = w+z*t-z**3/2
    directions = [z, t, *(basis[j] for j in CHARGES)]
    numeric = dict(zip(VARIABLES, range(1, 18, 2)))
    F = w+z*t+numeric[a]*z+numeric[b]*t+sum(numeric[v]*basis[j] for j, v in zip(CHARGES, CUBIC))
    for n in range(1, 7):
        assert S.expand(raw_expression(n).subs(numeric)-sphere_average(F**n)) == 0
    jacobian = S.Matrix([rodrigues_linearization(n) for n in range(2, 11)])
    independent = S.Matrix([[n*sphere_average(F0**(n-1)*direction) for direction in directions]
                            for n in range(2, 11)])
    assert jacobian == independent and jacobian.rank() == 9
    assert jacobian[:8, :].rank() == 8
    determinant = jacobian.det()
    assert determinant == -S.Rational(134217728, 94732977184601460705)
    local = json.loads((HERE/'nilpotent_cyclic_slice_local_jacobian.json').read_text())
    assert S.Matrix(local['raw_jacobian']) == jacobian and S.Rational(local['raw_determinant']) == determinant
    assert all(rodrigues_linearization(n) == [0]*9 for n in range(11, 15))

    # Complete global polynomial adapter for p=r=0 and s2=1/4.
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    B = S.Matrix([[0, 0, 0], [S.Rational(1, 2), 0, 0], [0, S.Rational(1, 2), 0]])
    v = S.Matrix([1, b, a])
    frame = S.Matrix.hstack(v, B*v, B**2*v)
    assert frame.det() == S.Rational(1, 8)
    gram = frame.T*metric*frame
    expected_gram = S.Matrix([[2*a+b*b, b, S.Rational(1, 4)],
                              [b, S.Rational(1, 4), 0], [S.Rational(1, 4), 0, 0]])
    assert gram == expected_gram and gram.det() == -S.Rational(1, 64)
    companion = S.Matrix([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    assert B*frame == frame*companion
    yy = S.symbols('y0 y1 y2')
    forward_y = [(frame.T*metric*S.Matrix([z, t, w]))[i] for i in range(3)]
    inverse_y = {z: 4*yy[2], t: 2*yy[1]-4*b*yy[2],
                 w: yy[0]-2*b*yy[1]+4*(b*b-a)*yy[2]}
    assert all(S.expand(f.subs(inverse_y, simultaneous=True)-yi) == 0 for f, yi in zip(forward_y, yy))
    triples = list(itertools.combinations_with_replacement(range(3), 3))
    tensor = {ijk: S.Symbol('U'+''.join(map(str, ijk))) for ijk in triples}
    trace = [sum(gram[i, j]*tensor[tuple(sorted((i, j, k)))] for i in range(3) for j in range(3))
             for k in range(3)]
    solved = S.solve(trace, [tensor[0, 1, 1], tensor[1, 1, 1], tensor[1, 1, 2]], dict=True)[0]
    assert not any(S.denom(S.cancel(value)).free_symbols for value in solved.values())
    generic_cubic = 0
    for ijk in triples:
        multiplicity = math.factorial(3)//math.prod(math.factorial(ijk.count(i)) for i in range(3))
        generic_cubic += multiplicity*tensor[ijk]*math.prod(forward_y[i] for i in ijk)
    generic_cubic = S.expand(generic_cubic.subs(solved, simultaneous=True))
    assert laplace(generic_cubic) == 0
    polynomial = S.Poly(generic_cubic, z, w, t)
    monomials = [w**3, w*w*t, w*t*t, t**3, z*t*t, z*z*t, z**3]
    coefficients = [polynomial.coeff_monomial(mon) for mon in monomials]
    assert S.expand(generic_cubic-sum(c*basis[j] for j, c in zip(CHARGES, coefficients))) == 0

    # Check the restriction of the actual17-dimensional algebra, without jets.
    vendor = HERE/'rodrigues_global_seven_target_certificate'
    base_data = HERE/'rodrigues_local_length17_certificate/certificate_data.json'
    assert base_data.read_bytes() == (vendor/'rodrigues_local_length17_certificate/certificate_data.json').read_bytes()
    old = json.loads(base_data.read_text())
    adapter = json.loads((vendor/'global_target_nilindices.json').read_text())
    Ma = S.Matrix(old['multiplication_a'])
    Mb = S.Matrix(old['multiplication_b'])
    operators = [Ma**i*Mb**j for i, j in old['quotient_basis']]
    def operator(vector):
        return sum((S.Rational(c)*M for c, M in zip(vector, operators)), S.zeros(17))
    Mp = operator(adapter['normalized_target_vectors']['p'])
    Mr = operator(adapter['normalized_target_vectors']['r'])
    columns = Mp.row_join(Mr)
    assert columns[0, :] == S.zeros(1, 34)
    _, pivots = columns.rref()
    assert len(pivots) == 16
    rank_minor = columns.extract(range(1, 17), pivots).det()
    assert rank_minor != 0

    # Exact GLOBAL augmented ideal certificate after five linear equations.
    cert = json.loads((HERE/'nilpotent_cyclic_slice_five_target_certificate.json').read_text())
    A, beta, C, D = S.symbols('A B C D')
    selected_variables = (A, beta, C, D)
    loc = dict(zip(cert['variables'], selected_variables))
    assert cert['five_targets'] == ['ell_z', 'ell_t', 'gminus3', 'gminus1', 'g1']
    equations = []
    for row in cert['moments']:
        raw = S.Integer(0)
        for pp, coefficient in exact_moment(row['n']).items():
            if any(pp[i] for i in [0, 1, 2, 4, 6]):
                continue
            raw += S.Rational(coefficient.numerator, coefficient.denominator)*A**pp[3]*beta**pp[5]*C**pp[7]*D**pp[8]
        polynomial = S.sympify(row['primitive'], locals=loc)
        assert S.expand(raw-S.Rational(row['raw_equals_primitive_times'])*polynomial) == 0
        equations.append(S.Poly(polynomial, *selected_variables, domain=S.QQ))
    targets = {'A': A, 'B': beta, 'C': C, 'D+1/2': D+S.Rational(1, 2)}
    for name, values in cert['witnesses'].items():
        multipliers = [S.Poly(S.sympify(c, locals=loc), *selected_variables, domain=S.QQ) for c in values]
        assert len(multipliers) == len(equations) == 8
        identity = sum((m*f for m, f in zip(multipliers, equations)),
                       S.Poly(0, *selected_variables, domain=S.QQ))
        assert identity.as_expr() == targets[name]
    assert all(f.as_expr().subs({A: 0, beta: 0, C: 0, D: -S.Rational(1, 2)}) == 0 for f in equations)

    program = {
        'claim': 'SIC2B33NILPCYCLIC20261117', 'input_prefix': 14, 'allow_later_prefixes': True,
        'source': 'nilpotent_cyclic_slice_prefix14.json.gz', 'source_sha256': digest('nilpotent_cyclic_slice_prefix14.json.gz'),
        'source_variables': list(map(str, VARIABLES)), 'source_moment_indices': list(range(2, 15)),
        'known_local_length_for_N_at_least10': 1,
        'known_augmented_global_identity': 'I9+(ell_z,ell_t,gminus3,gminus1,g1)=Rodrigues maximal ideal',
        'ordinary_target': 'dim B_N=1; equivalent to five exact degree1 memberships for N>=10',
        'radical_target': 'Five separate unit tests1 in I_N+(1-u*f); equivalent to length1 for fixed N>=10 by the local Jacobian',
        'cyclic_fixed_D_comparison': 'A_N/(p,r) over C is a product of3 isomorphic branch algebras B_N, so the optional normalized nilpotent target is length3.',
        'tests': [{'target': f, 'degree': 1, 'character_mod2': ch,
                   'auxiliary_character_mod2': ch, 'unit_target': '1'}
                  for f, ch in zip(cert['five_targets'], [0, 1, 0, 0, 0])],
        'unit_certificates_verified': 0, 'unrestricted_global_memberships_verified': 0,
        'witness_degree_bound': 'UNKNOWN', 'global_length': 'UNCOMPUTED',
        'full_Image33': 'OPEN',
    }
    (HERE/'nilpotent_cyclic_slice_five_target_program.json').write_text(json.dumps(program, indent=2)+'\n')
    receipt = {
        'claim': program['claim'], 'status': 'PASS', 'checker_sha256': digest(Path(__file__).name),
        'model_sha256': digest('nilpotent_cyclic_slice.py'),
        'input_sha256': program['source_sha256'],
        'selected_certificate_sha256': digest('nilpotent_cyclic_slice_five_target_certificate.json'),
        'moments_regenerated': list(range(2, 15)), 'term_counts': counts,
        'independent_spatial_checks': list(range(1, 7)),
        'independent_Rodrigues_derivative_checks': list(range(2, 11)),
        'Rodrigues_raw_Jacobian_determinant': str(determinant),
        'minimal_consecutive_ordinary_local_prefix_on_slice': 10,
        'actual17_algebra_p_r_image_rank': 16,
        'actual17_algebra_quotient_length': 1, 'rank_minor_columns': list(pivots),
        'rank_minor_determinant': str(rank_minor),
        'polynomial_cyclic_adapter': 'PASS; all cubic coefficients and constant trace pivots retained',
        'augmented_I9_identity': 'PASS, via four exact polynomial multiplier rows after the five linear targets',
        'unrestricted_global_five_memberships': 'OPEN', 'global_length1': 'OPEN',
        'normalized_nilpotent_length3': 'OPEN', 'full_Image33': 'OPEN',
        'msolve_dependency_for_proved_results': False, 'Singular_dependency_for_proved_results': False,
        'boundary': 'Global scouts have no certified output. Complete normalized slice and local/augmented identities are proved, not all pure-point exhaustion or generic marked safety. No independent specialist or formal verification.',
        'elapsed_seconds': round(time.monotonic()-started, 3),
    }
    (HERE/'nilpotent_cyclic_slice_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    run()
