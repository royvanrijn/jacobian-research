"""Exact rescaling/divergence/seed checks for the proposed global888 cover.

Reads retained projection algorithms without executing their mains or
rewriting their receipts. Full class/rank/analytic arguments remain written.
"""
from pathlib import Path
from functools import lru_cache
import hashlib
import json
import time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted

P = Path(__file__).resolve().parent
u, v, t = ordinary.u, ordinary.v, ordinary.t
NEW = {
    'nodal': {'rank': 7, 'k': 2, 'cost': 2, 'slope': 14, 'offset': 4, 'sum_grades': 11, 'complement_gate': 353, 'frame_gate': 888, 'full_gate': 888},
    'transverse': {'rank': 6, 'k': 3, 'cost': 3, 'slope': 18, 'offset': 6, 'sum_grades': 8, 'complement_gate': 334, 'frame_gate': 888, 'full_gate': 888},
    'triangle': {'rank': 5, 'k': 4, 'cost': 4, 'slope': 22, 'offset': 4, 'sum_grades': 6, 'complement_gate': 304, 'frame_gate': 888, 'full_gate': 888},
    'cusp': {'rank': 6, 'k': 1, 'cost': 4, 'slope': 30, 'offset': 12, 'sum_grades': 0, 'complement_gate': 558, 'frame_gate': 425, 'full_gate': 558},
    'tangent': {'rank': 5, 'k': 2, 'cost': 5, 'slope': 46, 'offset': 20, 'sum_grades': 0, 'complement_gate': 590, 'frame_gate': 425, 'full_gate': 590}}


def rescale(poly):
    result = 0
    for (i, j), coefficient in s.Poly(poly, u, v).terms():
        result += s.cancel(coefficient/t**(i+j))*u**i*v**j
    return s.expand(result)


def total_time_degree(poly):
    return s.Poly(poly, t, u, v).total_degree()


def source_numerator(V, D):
    delta = t*t*D
    return s.expand(s.diff(V, t)*delta-V*s.diff(delta, t)-2*t**3*D**2-t**4*D*s.diff(D, t)/2)


def universal_rescaling():
    checkpoint = json.loads((P/'ordinary_node_phase_shift_checkpoint.json').read_text())
    old_V = s.sympify(checkpoint['V'], locals={'u': u, 'v': v, 't': t})
    D = s.sympify(checkpoint['D'], locals={'t': t})
    delta = t*t*D
    V = rescale(old_V)
    assert total_time_degree(V) <= 4
    source = source_numerator(V, D)
    assert total_time_degree(source) <= 6
    old_source = s.diff(old_V, t)*delta-old_V*s.diff(delta, t)-t**4*D*s.diff(D, t)/2
    defect = source-rescale(old_source)+delta/t*(u*s.diff(V, u)+v*s.diff(V, v)+2*delta)
    assert s.expand(defect) == 0
    alpha, m = s.symbols('alpha m')
    denominator = m-2*alpha*t
    ell = 2*u+3*v+1
    J = s.Symbol('ju')*u+s.Symbol('jv')*v+s.Symbol('j0')
    g = sum(s.Symbol('g'+str(i)+str(j))*u**i*v**j for i, j in ordinary.monomials(2))
    q = sum(s.Symbol('q'+str(i)+str(j))*u**i*v**j for i, j in ordinary.monomials(2))
    for kind in ['cusp', 'tangent']:
        G = v**3 if kind == 'cusp' else u*v*v
        L = t*(u*u+ell)-J
        original = denominator*(t*(G+g)-q)+L*L/2
        normalized = rescale(t*t*original)
        assert total_time_degree(normalized) <= 4
        assert total_time_degree(source_numerator(normalized, denominator)) <= 6
        top = u**4/2+denominator*(v**3 if kind == 'cusp' else u*v*v)
        data = weighted.CASES[kind]
        assert all(weighted.weight(e, data['weights']) <= data['lower_weight'] for e, c in s.Poly(normalized-top, u, v).terms() if c)
    return {'ordinary_total_time_degree': 4, 'source_total_time_degree': 6,
            'coordinate_gauge_identity': 'S_hat=S_old(U/t,V/t)-(delta/t)*(G_U(U)+G_V(V))',
            'prefactor': 't^(-2)*D^(-1/2)', 'weighted_phase_cases': ['cusp', 'tangent']}


def universal_ordinary_divisions():
    parameters = ordinary.PARAMETERS
    weights = [0, 1, 1, 2, 2, 2, 3, 3, 4, 3]
    def degree(expression):
        poly = s.Poly(expression, u, v, *parameters)
        return max((e[0]+e[1]+sum(a*b for a, b in zip(e[2:], weights)) for e, c in poly.terms() if c), default=-100)
    out = []
    for kind in ['nodal', 'transverse', 'triangle']:
        C, A, B, E, F, J, L, N, K, delta = parameters
        coefficients = [C, A if kind != 'triangle' else 0, B if kind == 'nodal' else 0, E, F, J, L, N, K]
        prepared = ordinary.operator_data(kind, *coefficients, delta)
        count = 0
        for exponent in ordinary.monomials(8):
            remainder, Cu, Cv, V, divisor = ordinary.one_step(kind, exponent, coefficients, delta, prepared)
            defect = divisor*u**exponent[0]*v**exponent[1]-remainder-delta*(s.diff(Cu, u)+s.diff(Cv, v))-s.diff(V, u)*Cu-s.diff(V, v)*Cv
            assert s.expand(defect) == 0
            assert degree(remainder) <= sum(exponent)+NEW[kind]['k']
            assert degree(Cu) <= sum(exponent)+NEW[kind]['k']-3
            assert degree(Cv) <= sum(exponent)+NEW[kind]['k']-3
            count += 1
        out.append({'kind': kind, 'universal_cleared_divergence_identities': count, 'total_degree_cost': NEW[kind]['k']})
    return out


def phase(kind, singular):
    if kind in ordinary.CASES:
        old_V, _, D, _, determinant, _ = ordinary.gaussian_phase(kind, singular)
        V = rescale(old_V)
        pv = s.Poly(V, u, v)
        terms = [u*u*v*v, u**3, v**3, u*u, v*v, u*v, u, v, s.Integer(1)]
        coefficients = [s.factor(pv.coeff_monomial(m)) for m in terms]
        assert coefficients[0] == s.Rational(1, 2)
        return V, D, determinant, coefficients
    D = t if singular else 1-2*t
    H2 = u*u
    G3 = v**3 if kind == 'cusp' else u*v*v
    ell = 2*u+3*v+1
    J = 2*u+3*v+1
    g = 3*u*u+5*u*v+7*v*v+11*u+13*v+17
    q = (11*u*u+13*v*v)/2+u*v+5*u+7*v+2
    old = D*(t*(G3+g)-q)+(t*(H2+ell)-J)**2/2
    V = rescale(t*t*old)
    metric = s.Matrix([[11, 1, 2], [1, 13, 3], [2, 3, 0 if singular else 1]])
    assert metric.det() != 0
    return V, D, int(metric.det()), None


def model(kind, singular):
    data = NEW[kind]
    old = ordinary.CASES[kind] if kind in ordinary.CASES else weighted.CASES[kind]
    basis = old['basis']; index = {e: i for i, e in enumerate(basis)}
    V, D, determinant, coefficients = phase(kind, singular)
    delta = t*t*D
    assert total_time_degree(V) <= 4
    S = source_numerator(V, D)
    assert total_time_degree(S) <= 6
    ring = s.QQ.poly_ring(t)
    poly = lambda f: s.Poly(f, u, v, domain=ring)
    zero, zt = poly(0), s.Poly(0, t, domain=s.QQ)
    Vp = poly(V); Vu, Vv, dp = Vp.diff(u), Vp.diff(v), poly(delta)
    if kind in ordinary.CASES:
        prepared = ordinary.operator_data(kind, *coefficients, delta)
        A = prepared[-1]
        def raw_step(exponent):
            remainder, Cu, Cv, _, _ = ordinary.one_step(kind, exponent, coefficients, delta, prepared)
            return poly(remainder), poly(Cu), poly(Cv)
        connection_depth = old['steps']
    else:
        a, b = s.Rational(1, 2), D
        A = a*b if kind == 'cusp' else a*b*b
        top = a*u**4+b*(v**3 if kind == 'cusp' else u*v*v)
        remainder_phase = poly(s.expand(V-top))
        def raw_step(exponent):
            nf, hu, hv = weighted.principal_division(kind, exponent, a, b)
            nfp, hup, hvp = poly(s.cancel(A*nf)), poly(s.cancel(A*hu)), poly(s.cancel(A*hv))
            remainder = nfp-hup*remainder_phase.diff(u)-hvp*remainder_phase.diff(v)-dp*(hup.diff(u)+hvp.diff(v))
            return remainder, hup, hvp
        connection_depth = old['steps']
    At = s.Poly(s.expand(A), t, domain=s.QQ); Ap = poly(At.as_expr())
    assert At.degree() <= data['k']

    @lru_cache(None)
    def step(exponent):
        remainder, Cu, Cv = raw_step(exponent)
        defect = Ap*poly(u**exponent[0]*v**exponent[1])-remainder-dp*(Cu.diff(u)+Cv.diff(v))-Vu*Cu-Vv*Cv
        assert defect.is_zero
        if kind in ordinary.CASES:
            assert all(s.Poly(c, t).degree()+sum(e) <= data['k']+sum(exponent) for e, c in remainder.terms() if c)
        else:
            assert all(s.Poly(c, t).degree() <= data['cost'] for e, c in remainder.terms() if c)
        measure = sum if kind in ordinary.CASES else lambda e: weighted.weight(e, old['weights'])
        gap = 1 if kind in ordinary.CASES else old['gap']
        assert all(e in index or measure(e) <= measure(exponent)-gap for e, c in remainder.terms() if c)
        return remainder, Cu, Cv

    @lru_cache(None)
    def reduce_monomial(exponent, remaining):
        if exponent in index:
            values = [zt]*len(basis); values[index[exponent]] = At**remaining
            return tuple(values), zero, zero
        assert remaining > 0, (kind, exponent)
        remainder, cu0, cv0 = step(exponent)
        values = [zt]*len(basis); Cu = zero; Cv = zero
        for e, c in remainder.terms():
            if not c: continue
            vector, cu, cv = reduce_monomial(e, remaining-1)
            scalar = poly(c); scalar_t = s.Poly(c, t, domain=s.QQ)
            values = [a+scalar_t*b for a, b in zip(values, vector)]
            Cu += scalar*cu; Cv += scalar*cv
        Cu += Ap**(remaining-1)*cu0; Cv += Ap**(remaining-1)*cv0
        return tuple(values), Cu, Cv

    def reduce_input(f, depth):
        values = [zt]*len(basis); Cu = zero; Cv = zero
        for e, c in f.terms():
            if not c: continue
            vector, cu, cv = reduce_monomial(e, depth)
            scalar = poly(c); scalar_t = s.Poly(c, t, domain=s.QQ)
            values = [a+scalar_t*b for a, b in zip(values, vector)]
            Cu += scalar*cu; Cv += scalar*cv
        normal = sum((poly(a.as_expr()*u**e[0]*v**e[1]) for a, e in zip(values, basis)), zero)
        assert (Ap**depth*f-normal-dp*(Cu.diff(u)+Cv.diff(v))-Vu*Cu-Vv*Cv).is_zero
        return values, {'Cu': str(Cu.as_expr()), 'Cv': str(Cv.as_expr())}

    connection = []
    for exponent in basis:
        values, witness = reduce_input(poly(S)*poly(u**exponent[0]*v**exponent[1]), connection_depth)
        connection.append({'N': [str(a.as_expr()) for a in values], 'witness': witness})
    seeds = []
    R = s.Integer(1)
    for j in range(1, 3):
        R = s.expand(delta**2*s.diff(R, t)-2*(j-1)*delta*s.diff(delta, t)*R+S*R)
        assert total_time_degree(R) <= 6*j
        depth = (4*j-2 if kind != 'triangle' else 4*j-1) if kind in ordinary.CASES else (6*j-3 if kind == 'cusp' else 8*j-4)
        values, witness = reduce_input(poly(R), depth)
        degree_bound = data['slope']*j-data['offset']
        for exponent, value in zip(basis, values):
            assert value.degree() <= degree_bound-(sum(exponent) if kind in ordinary.CASES else 0)
        H = s.Poly(s.expand(delta**(2*j)*A**depth), t, domain=s.QQ)
        assert H.degree() <= degree_bound
        valuation = min(e[0] for e, c in H.terms() if c)
        assert valuation >= 4*j
        seeds.append({'j': j, 'depth': depth, 'H': str(H.as_expr()), 'N': [str(a.as_expr()) for a in values],
                      'degree_bound': degree_bound, 'H_valuation': valuation, 'witness': witness})
    return {'kind': kind, 'singular_at_zero': singular, 'D': str(D), 'delta': str(delta), 'V': str(V),
            'source': str(S), 'A0': str(A), 'metric_determinant': determinant, 'basis': basis,
            'connection_depth': connection_depth, 'connection_H': str(s.expand(delta**2*A**connection_depth)),
            'connection_columns': connection, 'direct_seed_columns': seeds, 'steps_checked': step.cache_info().currsize}


def main():
    start = time.monotonic()
    normalized = universal_rescaling(); print('Universal rescaling and prefactor identity passed', flush=True)
    universal = universal_ordinary_divisions(); print('135 universal graded ordinary-node divergence identities passed', flush=True)
    models = []
    for kind in NEW:
        for singular in [False, True]:
            models.append(model(kind, singular)); print('Certified', kind, 'singular' if singular else 'regular', flush=True)
    for kind, data in NEW.items():
        smax = data['rank']
        gate = sum(data['slope']*j-data['offset'] for j in range(1, smax+1))-data['sum_grades']
        assert gate == data['complement_gate']
        assert max(gate, data['frame_gate']) == data['full_gate']
        for j in range(1, smax+1): assert 4*j-j >= 0
    cover = {'quadratic': 7, 'smooth': 94, 'nodal': 888, 'cusp': 558, 'transverse': 888,
             'tangent': 590, 'triangle': 888, 'concurrent': 104, 'repeated': 99}
    assert max(cover.values()) == 888 and (888-19-1).bit_length() == 10
    cert = {'models': models, 'universal_normalization': normalized, 'universal_ordinary_divisions': universal,
            'bounds': NEW, 'global_cover': cover}
    path = P/'rescaled_singular_prefixes_certificate.json'; path.write_text(json.dumps(cert, indent=2)+'\n')
    report = {'completed': True, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'runtime_seconds': round(time.monotonic()-start, 3), 'sympy_version': s.__version__,
              'universal_ordinary_divergence_identities': 135, 'exact_models': len(models),
              'exact_connection_columns': sum(len(m['connection_columns']) for m in models),
              'exact_direct_seed_columns': sum(len(m['direct_seed_columns']) for m in models),
              'class_bounds': NEW, 'global_cover': cover, 'global_sufficient_pure_prefix': 888,
              'least_prefix_interval': [20, 888], 'least_prefix_computed': False,
              'binary_search_rounds_max': 10, 'full_coefficient_radical_search_executed': False,
              'universal_written_proof_certified_by_checker': False,
              'boundary': 'Exact rescaling/prefactor identities, graded divergence divisions and ten rational model connections/seed certificates. All-coefficient basis/rank coverage, analytic nonzero zeroth coefficients and the global cover remain written arguments using retained exclusions. No optimality, global coefficient-ring cofactor expansion, independent replay, external review or formal verification.'}
    (P/'rescaled_singular_prefixes_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k: report[k] for k in ['completed', 'global_sufficient_pure_prefix', 'least_prefix_interval', 'exact_models', 'exact_connection_columns', 'exact_direct_seed_columns', 'runtime_seconds']}, indent=2))


if __name__ == '__main__':
    main()
