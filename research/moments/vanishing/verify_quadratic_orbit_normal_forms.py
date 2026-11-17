"""Exact adapters for a complete field-valued quadratic orbit normalization.

Exhaustiveness is proved by self-adjoint Jordan decomposition in the linked
note. These matrix checks certify the displayed adapters and D factors;
they do not classify purity or identify ordinary quotient scheme lengths.
"""
from pathlib import Path
import hashlib
import json

import sympy as S

from conformal_euler_cubic import harmonic_basis, laplace, q, z, w, t

HERE = Path(__file__).resolve().parent


def run():
    eta = S.Symbol('eta')
    aa, bb, cc = S.symbols('ell_z ell_t ell_w')
    metric = S.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    spatial = S.Matrix([z, t, w])
    v = S.Matrix([cc, bb, aa])
    cubic = harmonic_basis()[3]
    forms = {
        'zero': S.zeros(3),
        'nilpotent_rank1': S.Matrix([[0, 0, 1], [0, 0, 0], [0, 0, 0]]),
        'nilpotent_rank2': S.Matrix([[0, 0, 0], [S.Rational(1, 2), 0, 0], [0, S.Rational(1, 2), 0]]),
        'nonnilpotent_Jordan': S.Matrix([[-S.Rational(1, 2), 0, 1], [0, 1, 0], [0, 0, -S.Rational(1, 2)]]),
        'semisimple': S.Matrix([[-eta/2, 0, 1], [0, eta, 0], [1, 0, -eta/2]]),
    }
    expected_H2 = {
        'zero': S.Integer(0), 'nilpotent_rank1': w*w, 'nilpotent_rank2': z*t,
        'nonnilpotent_Jordan': t*t-z*w+w*w,
        'semisimple': z*z+w*w+eta*(t*t-z*w),
    }
    expected_D = {
        'zero': 0, 'nilpotent_rank1': 0, 'nilpotent_rank2': cc**3/8,
        'nonnilpotent_Jordan': S.Rational(9, 4)*aa*aa*bb,
        'semisimple': bb*(aa-cc)*(aa+cc)*(9*eta*eta-4)/4,
    }
    assert forms['nilpotent_rank1']**2 == S.zeros(3)
    assert forms['nilpotent_rank1'].rank() == 1
    assert forms['nilpotent_rank2']**3 == S.zeros(3)
    assert forms['nilpotent_rank2'].rank() == 2
    J = forms['nonnilpotent_Jordan']
    assert J.det() == S.Rational(1, 4)
    assert (J+S.eye(3)/2)**2*(J-S.eye(3)) == S.zeros(3)
    assert (J+S.eye(3)/2)*(J-S.eye(3)) != S.zeros(3)
    semi = forms['semisimple']
    eigenvectors = S.Matrix([[1, 0, 1], [0, 1, 0], [1, 0, -1]])
    eigenvalues = S.diag(1-eta/2, eta, -1-eta/2)
    assert semi*eigenvectors == eigenvectors*eigenvalues
    assert eigenvectors.det() != 0

    rows = []
    for name, B in forms.items():
        assert B.trace() == 0 and metric*B == B.T*metric
        H2 = S.expand((spatial.T*metric*B*spatial)[0])
        assert S.expand(H2-expected_H2[name]) == 0 and laplace(H2) == 0
        D = S.factor(S.Matrix.hstack(v, B*v, B**2*v).det())
        assert S.expand(D-expected_D[name]) == 0
        p = S.expand(S.trace(B**2)/2)
        r = S.factor(B.det())
        disc = S.factor(4*p**3-27*r*r)
        assert (B**3-p*B-r*S.eye(3)).applyfunc(S.expand) == S.zeros(3)
        rows.append({'name': name, 'quadratic': str(H2),
                     'B': [[str(x) for x in B.row(i)] for i in range(3)],
                     'p_trace_B2_over2': str(p), 'r_det_B': str(r),
                     'characteristic_discriminant': str(disc), 'D': str(D),
                     'parameters': ['eta'] if name == 'semisimple' else [],
                     'linear_harmonic': str(aa*z+bb*t+cc*w),
                     'linear_coefficients': ['ell_z', 'ell_t', 'ell_w'],
                     'cubic_coefficients': ['gamma_'+str(j) for j in range(-3, 4)],
                     'cubic_basis': {str(j): str(f) for j, f in cubic.items()},
                     'coefficient_restrictions': 'None on H1 or H3.',
                     'full_purity_implies_Image_status': 'PROVED by odd theorem for pure points' if name == 'zero' else 'OPEN',
                     'known_subcases': ({
                         'nilpotent_rank1': 'H1=0 completely classified and safe.',
                         'nilpotent_rank2': 'H1=0 completely classified and safe; Rodrigues component safe when H1 is nonzero.',
                         'nonnilpotent_Jordan': 'Zero-linear case and full family open.',
                         'semisimple': 'H1=0 and eta=+/-2/3 excluded by the axial theorem. Other zero-linear eta and full mixed-linear family open.',
                     }.get(name, 'All pure points safe by the odd theorem.'))})
    output = {
        'claim': 'SIC2B33QUADFORMS20261117', 'status': 'EXACT_FIELD_COVERAGE_INPUT',
        'equivalence': 'Oriented linear isometry and nonzero phase amplitude, after algebraic field extension.',
        'coverage': 'Every trace-free self-adjoint three-dimensional B; every H1 and H3 retained.',
        'normal_forms': rows,
        'zero_constant': 'The first moment removes H0.',
        'generic_semisimple_boundaries': 'eta=+/-2/3 repeated eigenvalue; eta=0,+/-2 rank2. No eta value is omitted.',
        'maximal_phase_parameter_count': 11,
        'D0_coverage': 'All D=0 points retained. Factorizations give complete D0 equations in every normal form.',
        'scheme_boundary': 'Field-valued coverage only. No isomorphism or length comparison with the fixed-D cyclic ordinary algebra; no division of51 or17 by a residual symmetry.',
        'unrestricted_Image33': 'OPEN',
    }
    path = HERE/'quadratic_orbit_normal_form_input.json'
    path.write_text(json.dumps(output, indent=2)+'\n')
    receipt = {
        'claim': output['claim'], 'status': 'PASS',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'input_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'proof': 'QUADRATIC_ORBIT_NORMAL_FORMS.md, self-adjoint Jordan decomposition and gap normalization.',
        'forms': 5, 'largest_family_parameters': 11, 'D_factor_identities': 'PASS',
        'msolve_dependency': False, 'Singular_dependency_for_normalization': False,
        'purity_classification_or_global_exhaustion': 'NOT CLAIMED',
    }
    (HERE/'quadratic_orbit_normal_forms_verification.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    run()
