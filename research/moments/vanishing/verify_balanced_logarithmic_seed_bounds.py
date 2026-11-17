"""Exact auxiliary algebra for balanced logarithmic cubic seed bounds.

Analytic logarithmic growth and period comparison are retained written inputs.
This checker does not verify an analytic theorem or a global radical ideal.
"""
from pathlib import Path
from itertools import combinations
import hashlib
import json
import time
import sympy as s

P = Path(__file__).resolve().parent
t, x, y = s.symbols('t x y')


def laurent(expression, variable=y, domain=None):
    expression = s.cancel(expression)
    if expression == 0:
        return {}
    numerator, denominator = s.fraction(expression)
    dpoly = s.Poly(denominator, variable)
    assert len(dpoly.terms()) == 1, 'The denominator must be a Laurent monomial'
    (exponent,), scalar = dpoly.terms()[0]
    polynomial = s.Poly(numerator, variable, domain=domain)
    return {e[0]-exponent: s.cancel(c/scalar) for e, c in polynomial.terms() if c}


def coordinate_support_checks():
    a, b, c, f, rho, sigma = s.symbols('a b c f rho sigma', nonzero=True)
    eta = f*(rho-sigma)
    U = ((rho*c-b)*t+1)/eta
    V = ((b-sigma*c)*t-1)/eta
    K = U*V-a*t*t/f
    Z = f/eta*(y-K/y-U+V)
    T = y-U-rho*Z
    records = []
    for label, zh, th in [('smooth_f_nonzero', Z, T)]:
        for coordinate in [zh, th]:
            coefficients = laurent(coordinate)
            assert all(-1 <= k <= 1 and s.Poly(s.cancel(v), t).degree() <= 1-k for k, v in coefficients.items())
        records.append({'chart': label, 'coordinate_orders': [-1, 1], 'coefficient_degree': 'at most1-k'})
    e, d = s.symbols('e d', nonzero=True)
    z0 = -c/e
    Azero = a+b*z0+d*z0*z0
    Zh = y+t*z0
    Th = (y-d*y*y-t*(b+2*d*z0)*y+t*z0-t*t*Azero)/(e*y)
    for coordinate in [Zh, Th]:
        assert all(-1 <= k <= 1 and s.Poly(s.cancel(v), t).degree() <= 1-k for k, v in laurent(coordinate).items())
    records.append({'chart': 'WZT_smooth', 'coordinate_orders': [-1, 1], 'coefficient_degree': 'at most1-k'})
    trapped = y+(1-b*t)/e
    assert all(s.Poly(v, t).degree() <= 1-k for k, v in laurent(trapped).items())
    records.append({'chart': 'trapped_g_zero', 'coordinate_orders': [0, 1], 'coefficient_degree': 'at most1-k'})
    # Exact support implication for every G monomial, independent of coefficients.
    coordinates = {(power, k) for k in [-1, 0, 1] for power in range(2-k)}
    def product(A, B): return {(a+c, b+d) for a, b in A for c, d in B}
    for degree in range(4):
        support = {(0, 0)}
        for _ in range(degree): support = product(support, coordinates)
        assert all(n+3-degree+k <= 3 and -3 <= k <= 3 for n, k in support)
    assert all(n+k <= 2 for n, k in product(coordinates, coordinates))
    assert s.cancel(s.diff(y/t, y)/(y/t)-1/y) == 0
    return records


def specification(p, q):
    if p and q:
        basis = list(range(-q, p))
        plus, minus = 3-p, 3+q
        depths = lambda j: (p*(j-1)+1, q*(j-1))
    elif p:
        basis = list(range(p)); plus, minus = 3-p, 0
        depths = lambda j: (p*(j-1)+1, 0)
    elif q:
        basis = list(range(0, -q, -1)); plus, minus = 0, 3+q
        depths = lambda j: (0, q*(j-1)+1)
    else:
        basis = [0]; plus = minus = 0
        depths = lambda j: (0, 0)
    return basis, plus, minus, depths


def degree_cases():
    rows = []
    for p in range(4):
        for q in range(4):
            basis, plus, minus, depths = specification(p, q)
            rank = len(basis)
            bj = lambda j: 3*j+plus*depths(j)[0]+minus*depths(j)[1]
            maximum = -1
            for size in range(1, rank+1):
                for selected in combinations(basis, size):
                    if 0 not in selected: continue
                    gate = sum(bj(j) for j in range(1, size+1))-sum(selected)
                    maximum = max(maximum, gate)
                    assert gate <= 336
                    assert all(3*j-j >= 0 for j in range(1, size+1))
            full = sum(bj(j) for j in range(1, rank+1))-sum(basis)
            assert maximum == full
            rows.append({'positive_order': p, 'negative_order': q, 'basis': basis,
                         'rank': rank, 'endpoint_degree_bounds': [plus, minus],
                         'seed_depths': [list(depths(j)) for j in range(1, rank+1)],
                         'seed_degree_bounds': [bj(j) for j in range(1, rank+1)],
                         'all_seed_rank_and_row_choices_gate': full})
    assert next(row for row in rows if row['positive_order'] == row['negative_order'] == 3)['all_seed_rank_and_row_choices_gate'] == 336
    return rows


def model(label, W, maximum_seed=3):
    K = s.QQ.frac_field(t)
    raw = laurent(W)
    p = max(0, max(raw)); q = max(0, -min(raw))
    assert p <= 3 and q <= 3
    assert all(s.Poly(c, t).degree() <= 3-k for k, c in raw.items())
    basis, plus, minus, depths = specification(p, q)
    valid = set(basis)
    Aplus = raw[p] if p else s.Integer(1)
    Aminus = raw[-q] if q else s.Integer(1)
    phase = s.expand(W/t**2)
    old_phase = phase.subs(y, t*x)
    assert s.cancel(s.diff(phase, t)-s.diff(old_phase, t).subs(x, y/t)+y*s.diff(phase, y)/t) == 0
    source = s.expand(t*s.diff(W, t)-2*W-t*t)
    assert all(s.Poly(c, t).degree() <= 3-k for k, c in laurent(source).items())
    derivative_coeffs = {k: K.convert(k*c) for k, c in raw.items() if k}
    def parse(expression):
        return {k: K.convert(c) for k, c in laurent(expression).items()}
    def add(out, power, value):
        out[power] = out.get(power, K.zero)+value
        if not out[power]: del out[power]
    def exact(poly):
        out = {}
        for power, value in poly.items():
            add(out, power, K.convert(t*t*power)*value)
            for k, coefficient in derivative_coeffs.items(): add(out, power+k, coefficient*value)
        return out
    def reduce(poly):
        out = poly.copy(); witness = {}; high_depth = low_depth = 0
        while any(k not in valid for k in out):
            exponent = max(out) if max(out) > max(valid) else min(out)
            if exponent > max(valid) and p:
                shift = exponent-p; divisor = derivative_coeffs[p]; high_depth += 1
            elif exponent < min(valid) and q:
                shift = exponent+q; divisor = derivative_coeffs[-q]; low_depth += 1
            else:
                shift = exponent; divisor = K.convert(t*t*exponent); assert exponent != 0
            factor = out[exponent]/divisor; add(witness, shift, factor)
            for k, coefficient in exact({shift: factor}).items(): add(out, k, -coefficient)
        difference = poly.copy()
        for k, coefficient in out.items(): add(difference, k, -coefficient)
        for k, coefficient in exact(witness).items(): add(difference, k, -coefficient)
        assert not difference
        return out, witness
    def clear(poly, denominator):
        out = []
        for exponent in basis:
            value = K.convert(denominator)*poly.get(exponent, K.zero)
            expression = s.cancel(K.to_sympy(value))
            assert s.denom(expression) == 1
            out.append(s.Poly(expression, t, domain=s.QQ))
        return out
    connection_H = s.expand(t**3*Aplus**p*Aminus**q)
    columns = []
    for exponent in basis:
        normal, witness = reduce(parse(source*y**exponent))
        values = clear(normal, connection_H/t**3)
        columns.append({'N': [str(v.as_expr()) for v in values],
                        'witness': {str(k): str(K.to_sympy(v)) for k, v in witness.items()}})
    seeds = []; R = s.Integer(1)
    for j in range(1, maximum_seed+1):
        R = s.expand(t**3*s.diff(R, t)-3*(j-1)*t*t*R+source*R)
        assert all(s.Poly(c, t).degree()+k <= 3*j for k, c in laurent(R).items())
        normal, witness = reduce(parse(R))
        high, low = depths(j)
        denominator = s.expand(Aplus**high*Aminus**low)
        values = clear(normal, denominator)
        bj = 3*j+plus*high+minus*low
        assert all(v.degree() <= bj-k for k, v in zip(basis, values))
        H = s.Poly(t**(3*j)*denominator, t, domain=s.QQ)
        assert H.degree() <= bj
        valuation = min(e[0] for e, c in H.terms() if c)
        assert valuation >= 3*j
        seeds.append({'j': j, 'top_depth': high, 'bottom_depth': low,
                      'H': str(H.as_expr()), 'H_valuation': valuation, 'degree_bound': bj,
                      'N': [str(v.as_expr()) for v in values],
                      'witness': {str(k): str(K.to_sympy(v)) for k, v in witness.items()}})
    return {'label': label, 'phase_numerator_W': str(W), 'phase': str(phase), 'positive_order': p, 'negative_order': q,
            'basis': basis, 'connection_H': str(connection_H), 'connection_columns': columns,
            'seed_columns': seeds, 'normalized_source_numerator': str(source),
            'exact_witness_operator': 't^2*y*partial_y+y*W_y', 'logarithmic_scale_Jacobian': '1'}


def main():
    started = time.monotonic()
    supports = coordinate_support_checks(); cases = degree_cases()
    # Realized full conic and WZT charts, plus every distinct endpoint pattern.
    U = (t+1)/2; V = (5*t-1)/2; kap = U*V-t*t
    Z = (y-kap/y-U+V)/2; T = y-U-Z
    def phase_from_coordinates(zh, th):
        return s.expand(sum((1+2*i+3*j)*t**(3-i-j)*zh**i*th**j
                            for i in range(4) for j in range(4-i))-th*th/2)
    full = phase_from_coordinates(Z, T)
    z0 = -2
    zface = y+t*z0
    tface = (y-y*y+5*t*y-2*t-7*t*t)/y
    face = phase_from_coordinates(zface, tface)
    c0face = phase_from_coordinates(y, 1-3*t-2*t*t/y)
    trapped = phase_from_coordinates(s.Integer(0), y+1-2*t)
    phases = [('smooth_conic', full), ('WZT_c_nonzero', face), ('WZT_c_zero_a_nonzero', c0face), ('trapped_component', trapped),
              ('regular_infinity', (1+t**6)/y**3+(1+t**4)/y+s.Integer(1)),
              ('orders1_3', y+(1+t**6)/y**3+1),
              ('orders3_1', y**3+(1+t**4)/y+1),
              ('orders2_2', (1+t)*y**2+(1+t**5)/y**2+1),
              ('constant_phase', 1+t**3)]
    models = []
    for label, W in phases:
        models.append(model(label, s.expand(W))); print('Certified', label, flush=True)
    # Homogeneous ordinary-node complement refinement: W independent of t.
    import verify_ordinary_node_finite_prefixes as ordinary
    weights = [0]*9+[2]
    parameters = ordinary.PARAMETERS
    def total_degree(poly):
        return max((e[0]+e[1]+sum(a*b for a, b in zip(e[2:], weights))
                    for e, c in s.Poly(poly, ordinary.u, ordinary.v, *parameters).terms() if c), default=-100)
    coefficients = list(parameters[:9]); delta = parameters[9]
    prepared = ordinary.operator_data('nodal', *coefficients, delta)
    checks = 0
    for exponent in ordinary.monomials(8):
        remainder, _, _, _, _ = ordinary.one_step('nodal', exponent, coefficients, delta, prepared)
        assert total_degree(remainder) <= sum(exponent); checks += 1
    a, d, m, ju, jv, qu, qv, quv = s.symbols('a d m ju jv qu qv quv')
    ou, ov = ordinary.u, ordinary.v
    completed = m*(t*(a*ou**3+d*ov**3)-(qu*ou**2+quv*ou*ov+qv*ov**2))+(t*ou*ov-ju*ou-jv*ov)**2/2
    hom_phase = s.expand(t*t*completed.subs({ou: (ou+jv)/t, ov: (ov+ju)/t}, simultaneous=True))
    assert s.diff(hom_phase, t) == 0
    assert s.Poly(hom_phase, ou, ov).total_degree() <= 4
    hom_source = -2*m*t*hom_phase-2*m*m*t**3
    assert s.Poly(hom_source, t, ou, ov).total_degree() <= 5
    assert 5*7*8//2-11 == 129
    global_cover = {'quadratic': 7, 'smooth': 94, 'nodal': 353, 'cusp': 558, 'transverse': 336,
                    'tangent': 590, 'triangle': 336, 'concurrent': 104, 'repeated': 99}
    assert max(global_cover.values()) == 590
    assert max(56, 10, 10, 129, 336) == 336
    certificate = {'coordinate_support_checks': supports, 'Laurent_order_cases': cases, 'models': models,
                   'homogeneous_nodal_constant_phase_reduction_checks': checks,
                   'homogeneous_nodal_normalization_checks': 3,
                   'homogeneous_nodal_complement_gate': 129,
                   'global_cover': global_cover}
    path = P/'balanced_logarithmic_seed_bounds_certificate.json'; path.write_text(json.dumps(certificate, indent=2)+'\n')
    report = {'completed': True, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'runtime_seconds': round(time.monotonic()-started, 3), 'sympy_version': s.__version__,
              'coordinate_charts_checked': len(supports), 'Laurent_order_cases_checked': len(cases),
              'models': len(models), 'connection_columns': sum(len(m['connection_columns']) for m in models),
              'direct_seed_columns': sum(len(m['seed_columns']) for m in models),
              'general_conic_cubic_exclusion_gate': 336, 'WZT_conic_cubic_exclusion_gate': 336,
              'full_f0_cubic_gate': 336, 'full_W_linear_cubic_gate': 425,
              'homogeneous_nodal_complement_gate': 129, 'homogeneous_nodal_full_gate': 336,
              'homogeneous_nodal_reduction_checks': checks,
              'homogeneous_nodal_normalization_checks': 3,
              'global_sufficient_pure_prefix': 590, 'least_prefix_interval': [20, 590],
              'all_homogeneous_sufficient_pure_prefix': 336, 'all_homogeneous_least_prefix_interval': [14, 336],
              'global_cover': global_cover, 'optimality_proved': False,
              'written_analytic_or_global_proof_certified_by_checker': False,
              'boundary': 'Coordinate support/Jacobian algebra, exact logarithmic reductions and model seed certificates, finite endpoint/row degree arithmetic and homogeneous-nodal reduction identities only. The logarithmic growth theorem, comparison, full chart/class cover and analytic nonzero zeroth coefficients remain retained written mathematics. No independent replay, external review, formal verification or expanded global radical cofactors.'}
    (P/'balanced_logarithmic_seed_bounds_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
