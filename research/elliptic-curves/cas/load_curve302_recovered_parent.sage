"""Load the explicit MW17 parent; usable with sage -python or Sage load().

    E, basis, t0 = load_curve302_recovered_parent()

Coefficients are exact rational functions, stored low degree first in the
generated JSON. The fibre at t0=0 is literally the public curve302 model.
"""
import json
from pathlib import Path
from sage.all import QQ, PolynomialRing, EllipticCurve


def load_curve302_recovered_parent(data_path=None):
    if data_path is None:
        # Sage load() may inherit the caller's __file__. Support a repository
        # working directory as well as direct Python/runpy execution.
        here = Path(globals().get('__file__', '.')).resolve()
        relative = Path('artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json')
        candidates = [Path.cwd(), *Path.cwd().parents, here.parent, *here.parents]
        data_path = next((root/relative for root in candidates if (root/relative).is_file()), None)
        if data_path is None:
            raise FileNotFoundError('Pass the explicit path to curve302_recovered_mw17_parent_v1.json')
    data = json.loads(Path(data_path).read_text())
    R = PolynomialRing(QQ, 't')
    K = R.fraction_field()
    def value(record):
        return K(R(record['numerator']) / R(record['denominator']))
    E = EllipticCurve(K, [value(a) for a in data['a_invariants']])
    basis = [E([value(c) for c in point]) for point in data['basis_weierstrass_coordinates']]
    return E, basis, QQ(data['specialization_parameter'])


if __name__ == '__main__':
    E, basis, t0 = load_curve302_recovered_parent()
    print('Loaded 17 exact sections; specialization parameter:', t0)
    print('Specialized a-invariants:', [a(t0) for a in E.a_invariants()])
