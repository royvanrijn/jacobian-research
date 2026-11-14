"""Load the certified fixed1020 Weierstrass model and known section points.

Use with sage -python; no search, specialization, or rank computation.
"""
from pathlib import Path
import json
import sys
from sage.all import QQ, PolynomialRing, FunctionField, EllipticCurve

PACKET = Path(__file__).resolve().parents[2] / 'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1'
sys.set_int_max_str_digits(0)


def load_model():
    data = json.loads((PACKET / 'rr-weierstrass-invariants-QQ.json').read_text())
    assert data['status'] == 'PASS_EXACT_QQ_POINTED_CUBIC_JACOBIAN_INVARIANTS'
    base = FunctionField(QQ, 'u')
    polynomials = PolynomialRing(QQ, 'u')
    curve = EllipticCurve(base, [base(polynomials(data['A'])), base(polynomials(data['B']))])
    return curve


def load_known_sections(curve=None):
    curve = load_model() if curve is None else curve
    data = json.loads((PACKET / 'rr-component-sections-QQ.json').read_text())
    assert data['status'] == 'PASS_EXACT_QQ_FOUR_COMPONENT_SECTION_COORDINATES'
    ring = PolynomialRing(QQ, 'u')
    return {row['component']: curve([curve.base_ring()(ring(row['x'])), curve.base_ring()(ring(row['y']))]) for row in data['records']}


def load_marked_basis(curve=None):
    """Return all17 points in the order of the certified physical marking."""
    curve = load_model() if curve is None else curve
    data = json.loads((PACKET / 'rr-seventeen-sections-QQ.json').read_text())
    assert data['status'] == 'PASS_EXACT_QQ_SEVENTEEN_SECTION_IDENTITIES_HEIGHTS_AND_SATURATION'
    marking = json.loads((PACKET / 'rr-seventeen-marking-identification.json').read_text())
    assert marking['status'] == 'PASS_EXACT_GOOD_REDUCTION_IDENTIFICATION_OF_ALL17_MARKED_SECTIONS'
    ring = PolynomialRing(QQ, 'u')
    return [curve([curve.base_ring()(ring(row['x'])), curve.base_ring()(ring(row['y']))]) for row in data['records']]


def source_to_weierstrass(t, x, y, field=QQ):
    """Evaluate the original-source forward map on its affine domain."""
    sys.path.insert(0, str(PACKET))
    from determinantal_source_inverse import forward
    obj, w, v = forward(field, field(t), field(x), field(y))
    return obj.u, w, v


def weierstrass_to_source(u, w, v, field=QQ):
    """Evaluate the fixed determinant inverse; poles raise arithmetic errors.

    These are exact rational formulas. QQ evaluation may be expensive because
    the coefficients are large; finite-field evaluation is also supported.
    """
    sys.path.insert(0, str(PACKET))
    from determinantal_source_inverse import Inverse
    certificate = json.loads((PACKET / 'rr-source-inverse-certificate.json').read_text())
    assert certificate['status'] == 'PASS_EXACT_INVERSE_NONZERO_MINORS_AND_MODULAR_COMPOSITIONS'
    return Inverse(field, field(u)).recover(field(w), field(v), certificate['t_matrix_rows'], certificate['x_matrix_rows'])


if __name__ == '__main__':
    E = load_model()
    points = load_marked_basis(E)
    known = load_known_sections(E)
    assert len(points) == 17 and points[0] == -known['I10_theta2']
    print('Loaded determinant1020 model over QQ(u), coefficient degrees8,12.')
    print('Loaded all17 exact sections in the retained saturated physical basis order.')
    print('Source forward and determinantal inverse evaluation APIs are available.')
    print('No rank computation performed; rank and saturation use the retained certificates.')
