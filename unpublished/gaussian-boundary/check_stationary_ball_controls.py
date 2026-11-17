"""Exact controls for the written theorem; not a completeness proof.

New v003 check. Preserve its result; do not replay sealed checks merely
to regenerate a status receipt.
"""
from pathlib import Path
import hashlib
import json
import sympy as S

r, x, y, t, b, f = S.symbols('r x y t b f', nonzero=True)
def valuation(a):
    a = S.cancel(a)
    if a == 0:
        return S.oo
    n, d = S.fraction(a)
    return min(k[0] for k in S.Poly(n, r).monoms()) - min(k[0] for k in S.Poly(d, r).monoms())

controls = []
def balls(name, phase, roots):
    records = []
    for center in roots:
        shifted = S.Poly(S.expand(phase.subs(x, center+x)-phase.subs(x, center)), x)
        rho = max(-S.Rational(valuation(a), j[0]) for j, a in shifted.terms() if j[0] >= 2)
        scaled = S.Poly(S.expand(shifted.as_expr().subs(x, r**rho*y)), y)
        assert all(valuation(a) >= 0 for a in scaled.all_coeffs())
        reduction = S.expand(sum(S.limit(a, r, 0)*y**j[0] for j, a in scaled.terms()))
        records.append(dict(center=str(center), rho=str(rho), reduction=str(reduction), rank=S.degree(reduction, y)-1))
    # Group identical balls by exact containment; all roots in a ball have the same scale.
    groups = []
    for i, root in enumerate(roots):
        if any(i in g for g in groups):
            continue
        rho = S.Rational(records[i]['rho'])
        group = [j for j, other in enumerate(roots) if valuation(root-other) >= rho]
        assert all(S.Rational(records[j]['rho']) == rho for j in group)
        assert len(group) == records[i]['rank']
        groups.append(group)
    assert sum(records[g[0]]['rank'] for g in groups) == S.degree(phase, x)-1
    controls.append(dict(name=name, phase=str(phase), roots=records, disjoint_balls=groups))

balls('colliding non-Morse block', x**4/4-r**2*x**2/2, [0, r, -r])
balls('finite and escaping blocks', r*x**3/3-x**2/2, [0, 1/r])
phase = x**4/4-(2/r**2+1)*x**3/3+(1/r**4+1/r**2)*x**2/2
balls('unequal radii with nested center clusters', phase, [0, 1/r**2, 1/r**2+1])

u = x-1/r**2
filters = [1-r**2*x, (1-u)*r**2*x, u*r**2*x]
centers = [0, 1/r**2, 1/r**2+1]
scales = [r**2, r, r]
matrix = []
for filt in filters:
    row = []
    for center, delta in zip(centers, scales):
        value = S.Poly(S.expand(filt.subs(x, center+delta*y)), y)
        assert all(valuation(a) >= 0 for a in value.all_coeffs())
        row.append(S.expand(sum(S.limit(a, r, 0)*y**j[0] for j, a in value.terms())))
    matrix.append(row)
assert S.Matrix(matrix) == S.eye(3)

gaussian = []
for k in (1, 2, 3):
    bound = 1-S.Rational(k, 2)
    gaussian.append(dict(phase=f'-r^{k}*x^2/2', h='r', lower_bound=str(bound), normalized_seed_power=str(bound)))
assert [a['lower_bound'] for a in gaussian] == ['1/2', '0', '-1/2']

q = 1-b*t
flattened = S.cancel((-b/(2*f))*t*f/q-S.Rational(1, 2))
assert S.cancel(flattened+1/(2*q)) == 0

record = dict(version='v003', date='2026-11-17', status='passed',
    scope='Exact critical-ball, degree/multiplicity, uneven-scale filter and strict-threshold controls only.',
    not_proved_by_this_check=['General polynomial separation lemma', 'Formal functional horizontality',
        'Full-frame theorem', 'Semisimplicity', 'All-coefficient rank-one classification', 'Logarithmic cut estimates'],
    critical_ball_controls=controls, uneven_scale_reduction_matrix=[[str(a) for a in row] for row in matrix],
    strict_gaussian_threshold_controls=gaussian, flattened_quadratic_coefficient=str(flattened),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    independent_replay=False, formal_verification=False, external_review=False)
out = Path(__file__).with_name('stationary_ball_controls.json')
assert not out.exists(), 'Preserve the successful receipt instead of overwriting it.'
out.write_text(json.dumps(record, indent=2, default=str)+'\n')
print(json.dumps(dict(status='passed', critical_ball_families=3, uneven_scale_matrix='identity3', strict_threshold_cases=3,
    mathematical_completeness_proved_by_checks=False)))
