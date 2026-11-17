"""Declared-weight certificates and small exact controls of the written transfer."""
from pathlib import Path
from fractions import Fraction
import hashlib
import json
import platform
import sympy as s


def weight_part(poly, variables, weights, degree):
    return s.expand(sum(c*s.prod(v**a for v, a in zip(variables, mon))
                        for mon, c in s.Poly(poly, *variables).terms()
                        if sum(w*a for w, a in zip(weights, mon)) == degree))


def split_certificate(symbol, phase, variables, derivatives, weights, degree):
    assert len(variables) == len(derivatives) == len(weights)
    assert all(isinstance(w, int) and w > 0 for w in weights) and degree > 0
    assert set(variables).isdisjoint(s.sympify(symbol).free_symbols)
    ps, ls = s.Poly(phase, *variables), s.Poly(symbol, *derivatives)
    assert all(sum(w*a for w, a in zip(weights, mon)) <= degree for mon, c in ps.terms() if c)
    assert all(sum(w*a for w, a in zip(weights, mon)) >= degree for mon, c in ls.terms() if c)
    f = weight_part(phase, variables, weights, degree)
    b = weight_part(symbol, derivatives, weights, degree)
    g, rr = s.expand(phase-f), s.expand(symbol-b)
    gp = [degree-sum(w*a for w, a in zip(weights, mon)) for mon, c in s.Poly(g, *variables).terms() if c]
    gl = [sum(w*a for w, a in zip(weights, mon))-degree for mon, c in s.Poly(rr, *derivatives).terms() if c]
    return b, f, rr, g, min(gp) if gp else 1, min(gl) if gl else 1


def apply_symbol(symbol, poly, variables, derivatives):
    value = 0
    for mon, c in s.Poly(symbol, *derivatives).terms():
        term = poly
        for v, count in zip(variables, mon):
            if count:
                term = s.diff(term, v, count)
        value += c*term
    return s.expand(value)


def run():
    x, y, z = s.symbols('x y z')
    dx, dy, dz = s.symbols('dx dy dz')
    alpha, beta, u, v, ell, mm, cc, gamma, eta = s.symbols('alpha beta u v ell mm cc gamma eta')
    def iterate(symbol, poly, xs, ds, count):
        for _ in range(count):
            poly = apply_symbol(symbol, poly, xs, ds)
        return poly
    examples = [
        ('isotropic_principal', dx**2-dy**2+dz**3, (x+y)**2+x-z+1,
         (x, y, z), (dx, dy, dz), (1, 1, 1), 2, x+z**2),
        ('weighted_drift', dy+dx**3, x**2+alpha*x+beta,
         (x, y), (dx, dy), (1, 2), 2, y**2+x),
        ('nonpure_principal', dx**2+dz**3, x**2+y+1,
         (x, y, z), (dx, dy, dz), (1, 1, 1), 2, z+x**2),
    ]
    rows = []
    budget_controls = 0
    for label, symbol, phase, xs, ds, ws, degree, q in examples:
        b, f, rr, g, dp, dl = split_certificate(symbol, phase, xs, ds, ws, degree)
        e = max(sum(w*a for w, a in zip(ws, mon)) for mon, c in s.Poly(q, *xs).terms() if c)
        qe = weight_part(q, xs, ws, e)
        for m in range(1, 4):
            assert s.expand(iterate(symbol, phase**m, xs, ds, m)-iterate(b, f**m, xs, ds, m)) == 0
            full = iterate(symbol, q*phase**m, xs, ds, m)
            leading = iterate(b, qe*f**m, xs, ds, m)
            assert s.expand(weight_part(full, xs, ws, e)-leading) == 0
            for k in range(m+1):
                for l in range(m+1):
                    if k*dp+l*dl <= e:
                        continue
                    value = iterate(b, q*g**k*f**(m-k), xs, ds, m-l)
                    value = iterate(rr, value, xs, ds, l)
                    assert value == 0
                    budget_controls += 1
        rows.append({'case': label, 'weights': ws, 'r': degree,
                     'phase_gap': dp, 'operator_gap': dl,
                     'pure_equality_controls': 3, 'marked_top_controls': 3})

    # Universal ordinary pure ideal independence in a small parameter family.
    symbol = alpha*dx**2+beta*dy**2+gamma*dx**3+eta*dy**4
    phase = (u*x+v*y)**2+ell*x+mm*y+cc
    b, f, rr, g, dp, dl = split_certificate(symbol, phase, (x, y), (dx, dy), (1, 1), 2)
    ideal_rows = []
    for m in range(1, 4):
        pure = iterate(symbol, phase**m, (x, y), (dx, dy), m)
        expected = s.factorial(2*m)*(alpha*u**2+beta*v**2)**m
        assert s.expand(pure-expected) == 0
        assert not (set((ell, mm, cc, gamma, eta)) & pure.free_symbols)
        ideal_rows.append({'m': m, 'lower_phase_and_higher_symbol_parameters_absent': True})

    # Equality is literal over coefficients,including a nontrivial nilpotent specialization.
    eps = s.symbols('eps')
    for m in range(1, 4):
        value = iterate(eps*dx**2+dy**3, (eps*x**2+y+1)**m, (x, y), (dx, dy), m)
        expected = iterate(eps*dx**2, (eps*x**2)**m, (x, y), (dx, dy), m)
        assert s.expand(value-expected) == 0
        assert s.Poly(s.expand(value), eps).nth(0) == s.Poly(s.expand(value), eps).nth(1) == 0

    # A bounded DVG1 perturbation regression;the all-order input remains imported.
    aa, bb, c, d, tt = s.symbols('a b c d t')
    da, db, dc, dd, dt = s.symbols('da db dc dd dt')
    xs, ds = (aa, bb, c, d, tt), (da, db, dc, dd, dt)
    seed = (tt+c)*(aa*d+bb*tt)
    base = dt*(da*dd-db*dc)
    marker = -c
    # Deterministically choose one derivative monomial extracting a degree-four term.
    one_step = iterate(base, marker*seed**2, xs, ds, 1)
    mon = next(mon for mon, coeff in s.Poly(one_step, *xs).terms() if coeff and sum(mon) == 4)
    higher = s.prod(q**power for q, power in zip(ds, mon))
    newphase = seed+aa**2+c*tt+1
    newsymbol = base+higher
    split_certificate(newsymbol, newphase, xs, ds, (1,)*5, 3)
    bad_rows = []
    for m in (2, 3):
        assert iterate(newsymbol, newphase**m, xs, ds, m) == 0
        value = iterate(newsymbol, marker*newphase**m, xs, ds, m)
        top = weight_part(value, xs, (1,)*5, 1)
        expected = (-1)**(m+1)*s.factorial(m+1)*s.factorial(m)**2*tt
        assert s.expand(top-expected) == 0
        bad_rows.append({'m': m, 'degree_one_marker': str(top), 'lower_marker_part': str(s.expand(value-top))})

    # Quantitative inequality checked with exact fractions,without a cutoff campaign.
    bound_controls = 0
    for leading_c in (1, 2, 3):
        for dp, dl in [(1, 1), (1, 2), (2, 1), (2, 3)]:
            for dg, df in [(0, 2), (1, 3), (2, 5)]:
                for wmax in (1, 2, 3):
                    cprime = leading_c+wmax*(Fraction(1, min(dp, dl))+leading_c*max(Fraction(dg, dp), Fraction(df, dl)))
                    for degq in (1, 2):
                        for e in range(wmax*degq+1):
                            for k in range(e//dp+1):
                                for l in range((e-k*dp)//dl+1):
                                    threshold = max(k, l)+leading_c*(degq+k*dg+max(l-k, 0)*df)
                                    assert threshold <= cprime*degq
                                    bound_controls += 1
    # Negative/zero weights are not accepted as certificates.
    rejected = False
    try:
        split_certificate(dx**2, x**2, (x, y), (dx, dy), (1, 0), 2)
    except AssertionError:
        rejected = True
    assert rejected
    return {'completed': True, 'python_version': platform.python_version(),
            'sympy_version': s.__version__, 'filtered_identity_cases': rows,
            'outside_defect_budget_zero_controls': budget_controls,
            'ordinary_pure_ideal_parameter_controls': ideal_rows,
            'nilpotent_coefficient_identity_controls': 3,
            'bounded_DVG1_perturbation_controls': bad_rows,
            'chosen_DVG1_order_four_perturbation': str(higher),
            'quantitative_transfer_inequality_controls': bound_controls,
            'nonpositive_weight_mutation_rejected': True,
            'declared_weight_certificate_checker_implemented': True,
            'all_order_transfer_proof_inferred_from_samples': False,
            'DVG1_all_order_proof_reestablished': False,
            'ternary_Delta6_producer_executed': False,
            'arbitrary_principal_pair_classification_proved': False,
            'general_Gaussian_uniformity_proved': False,
            'ordinary_mixed_cutoff_ideals_claimed_unchanged': False,
            'zero_or_negative_weights_covered': False,
            'protected_ternary_frontier_changed': False,
            'independent_review': False, 'formal_verification': False,
            'novelty_claim': False}


if __name__ == '__main__':
    path = Path(__file__).resolve()
    result = run()
    result['checker_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('filtered_principal_transfer_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
