"""Small exact norm matrices for the two canonical immersed r7 image charts."""
from pathlib import Path
import importlib.util
import sympy as s

HERE = Path(__file__).resolve().parent
FAMILY_PATH = HERE.parent/'jc2-r7-image-family/family.py'


def families():
    spec = importlib.util.spec_from_file_location('r7_canonical_family', FAMILY_PATH)
    family = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(family)
    return [(name, parameter, s.expand(family.p.subs(seed)), s.expand(family.q.subs(seed)), family.t)
            for name, parameter, seed in [('L', family.d, family.line),
                                          ('K', family.x, family.laurent)]]


def norm_matrix(p, q, t, U, V):
    """Matrix of q(t)-V on the rank-six algebra C[t]/(p(t)-U)."""
    pp = s.Poly(p, t)
    assert pp.degree() == 6 and pp.LC() == 1
    relation = [s.expand((U if i == 0 else 0)-pp.nth(i)) for i in range(6)]

    def times_t(vec):
        top = vec[-1]
        return [s.expand((vec[i-1] if i else 0)+top*relation[i]) for i in range(6)]

    reduced = [s.Integer(0)]*6
    for coefficient in s.Poly(q, t).all_coeffs():
        reduced = times_t(reduced)
        reduced[0] = s.expand(reduced[0]+coefficient)
    columns = []
    for j in range(6):
        columns.append([s.expand(value-(V if i == j else 0)) for i, value in enumerate(reduced)])
        reduced = times_t(reduced)
    return s.Matrix.hstack(*[s.Matrix(col) for col in columns])


def matching_rows(c, grouping):
    """All necessary branch-count completions, not claimed global realizations."""
    rows = []
    if grouping.startswith('A'):
        for b in range(4):
            for h in range(7):
                if grouping == 'A_H=C' and h != c:
                    continue
                if grouping == 'A_H=B' and h != b:
                    continue
                retained = 6-2*c-2*b-h
                if retained >= 0:
                    rows.append(dict(b=b, h=h, cusp=False, retained_length=retained))
    elif grouping == 'B':
        for b in range(4):
            for cusp in [False, True]:
                if cusp and b == 0:
                    continue
                retained = 6-2*c-2*b-int(cusp)
                if retained >= 0:
                    rows.append(dict(b=b, h=0, cusp=cusp, retained_length=retained))
    elif grouping == 'C':
        for b in range(3):
            retained = 6-2*c-3*b
            if retained >= 0:
                rows.append(dict(b=b, h=0, cusp=False, retained_length=retained))
    else:
        raise ValueError(grouping)
    return rows
