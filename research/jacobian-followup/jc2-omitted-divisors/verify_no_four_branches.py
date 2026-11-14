"""Reconstruct quartic-divisor equations and replay y^6 membership with Fractions."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import platform
import re
import time
import sympy as s
from image_norm import families

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHECKS = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode(text):
    """Only signed rational monomials in L,M,y; no symbolic evaluator."""
    if not re.fullmatch(r'[0-9LMy+*/^\-]+', text):
        raise ValueError('unrecognized polynomial grammar')
    out = {}
    for token in re.findall(r'[+-]?[^+-]+', text):
        coefficient = Fraction(-1 if token.startswith('-') else 1)
        exponents = [0, 0, 0]
        for factor in token.lstrip('+-').split('*'):
            match = re.fullmatch(r'([LMy])(?:\^(\d+))?', factor)
            if match:
                exponents['LMy'.index(match[1])] += int(match[2] or 1)
            else:
                coefficient *= Fraction(factor)
        key = tuple(exponents)
        out[key] = out.get(key, Fraction(0))+coefficient
    return {key: value for key, value in out.items() if value}


def main():
    started = time.monotonic()
    certificate = json.loads((HERE/'four_branch_certificate.json').read_text())
    _, x, p, q, t = next(row for row in families() if row[0] == 'K')
    y, L, M = s.symbols('y L M')

    def scaled(poly, n):
        value = s.Poly(s.expand(x**n*poly.subs(t, t/x)), t, x)
        check(f'weight{n} scaling has only nonnegative powers of x^3',
              all(j >= 0 and j % 3 == 0 for (i, j), c in value.terms()))
        return s.expand(sum(c*t**i*y**(j//3) for (i, j), c in value.terms()))

    P, Q = scaled(p, 6), scaled(q, 10)
    check('scaled family P reconstructed', str(P) == certificate['P'])
    check('scaled family Q reconstructed', str(Q) == certificate['Q'])
    pp = s.Poly(P, t)
    check('scaled P monic centered degree6', pp.degree() == 6 and pp.LC() == 1 and pp.nth(5) == 0)
    check('scaled Q monic degree10', s.Poly(Q, t).degree() == 10 and s.Poly(Q, t).LC() == 1)
    g0 = pp.nth(4)-M+L*L
    h1 = s.expand(pp.nth(3)+L*M-L*g0)
    h0 = s.expand(pp.nth(2)+L*h1-M*g0)
    h = t**4+L*t**3+M*t*t+h1*t+h0
    raw_f0 = s.expand(-L*h0+h1*g0-pp.nth(1))
    check('whole degree6 cofactor identity',
          s.expand(P-h*(t*t-L*t+g0)+raw_f0*t+h0*g0) == 0)
    quotient, rem = s.div(Q, h, t)
    check('whole monic degree10 division identity', s.expand(Q-quotient*h-rem) == 0)
    raw = [raw_f0]+[s.expand(rem).coeff(t, i) for i in [3, 2, 1]]
    generators = [s.Poly(f, L, M, y).clear_denoms()[1].primitive()[1].as_expr() for f in raw]
    check('all four original necessary generators reconstructed', [str(f) for f in generators] == certificate['generators'])
    equations = [{tuple(e): Fraction(int(c.p), int(c.q)) for e, c in s.Poly(f, L, M, y).terms()}
                 for f in generators]
    check('generator complete supports', list(map(len, equations)) == [17, 42, 67, 85])
    multipliers = [decode(text) for text in certificate['multipliers']]
    check('four complete multiplier supports', list(map(len, multipliers)) == [386, 254, 162, 147])
    total = {}
    products = 0
    for generator, multiplier in zip(equations, multipliers):
        for a, ca in generator.items():
            for b, cb in multiplier.items():
                key = tuple(i+j for i, j in zip(a, b))
                total[key] = total.get(key, Fraction(0))+ca*cb
                products += 1
    total = {k: v for k, v in total.items() if v}
    check('literal exact rational identity y^6=sum C_i F_i', total == {(0, 0, 6): Fraction(1)})
    check('complete uncropped product support replay', products == 40579)
    check('certificate target power6', certificate['target'] == 'y^6')
    strata_path = ROOT/'jacobian-followup/jc2-r7-collisions/strata.json'
    strata = json.loads(strata_path.read_text())
    check('L nonimmersed type imported exactly',
          [(r['cusps'], r['A5'], r['nodes']) for r in strata['singularities'] if r['cusps']] == [(1, 1, 8)])
    check('L all immersed singularity types at most triple',
          all(r['cusps'] == 0 and r['A5'] == 0 for r in strata['singularities'] if not r['cusps']))
    image_ledger = json.loads((HERE/'image_norm_ledger.json').read_text())
    check('combined bound attained by saved L triple-rank certificate',
          image_ledger['triple_values_at_d0']['R_omitted_local_lengths'] == [2, 2, 2])
    ledger = dict(as_of='2026-11-14', claim='PF2OMIT12',
        hypotheses='Complete normalized degree(6,10) gap21 coefficient families PF2R7FAM1; no immersion required for parametrization fibre length bound',
        Laurent_scaling='y=x^3, tau=x*t, target scales x^6,x^10; x nonzero',
        Laurent_parameter_fibre_scheme_length_at_most=3,
        polynomial_parameter_fibre_scheme_length_at_most=3,
        repeated_quartic_common_divisors_excluded=True,
        certificate_target='y^6', exact_rational_product_count=products,
        first_image_norm_matrix_rank_at_least=3,
        previously_open_Laurent_rank_at_most2_locus='empty for every nonzero parameter; removes no parameter',
        immersed_first_omitted_prime=dict(e=2, f=1, omitted=True, ramifies=True,
            possible_normalization_point_counts=[0, 1, 2, 3],
            possible_total_ambient_lengths=[0, 2, 4, 6],
            ambient_length_per_point=2),
        nonimmersed_corank_is_branch_count=False,
        combined_bound_sharp=True, separate_Laurent_sharpness_claimed=False,
        Laurent_triple_and_tangency_strata_classified=False,
        actual_second_image_equation=None, actual_unramified_image_equation=None,
        actual_global_normalization_constructed=False, full_goal_complete=False)
    target = HERE/'four_branch_ledger.json'
    target.write_text(json.dumps(ledger, indent=2)+'\n')
    inputs = [HERE/n for n in ['NO_FOUR_BRANCHES.md', 'four_branch_certificate.json', 'image_norm.py', 'image_norm_ledger.json']]
    inputs += [ROOT/'jacobian-followup/jc2-r7-image-family'/n for n in ['family.py', 'families.json', 'PROOF.md']]
    inputs += [ROOT/'jacobian-followup/jc2-r7-collisions'/n for n in ['PROOF.md', 'strata.json']]
    inputs += [ROOT/'archive/non-elliptic/scripts/verify_f2_degree_6_10_gap.py']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), checks=CHECKS,
        checker_sha256=digest(Path(__file__)), input_sha256={str(p.relative_to(ROOT)): digest(p) for p in inputs},
        ledger_sha256=digest(target), runtime_seconds=round(time.monotonic()-started, 6),
        python=platform.python_version(), sympy=s.__version__,
        proof_boundary='Original quartic-divisor equations reconstructed; all40579 rational term products replayed independently of Singular. Cofactor necessity and imported family/collision completeness remain written arguments. No global cover or additional image is constructed.',
        formal_verification=False, external_review=False, full_goal_complete=False)
    (HERE/'four_branch_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
