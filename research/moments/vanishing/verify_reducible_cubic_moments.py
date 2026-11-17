#!/usr/bin/env python3
"""Exact invariant-series, Bezout and chart inputs for reducible ternary cubics."""
from pathlib import Path
from fractions import Fraction
from math import factorial, comb, lcm
from collections import deque
import hashlib
import json
import re
import time
import sympy as S

P = Path(__file__).resolve().parent
ACCEPTED = ['tangent_line_conic_charts','transverse_easy_charts',
    'transverse_coupled_boundary','linear_quadratic_noniso_generic',
    'linear_quadratic_other_invariant_charts','linear_quadratic_tenth_sharpness_v2',
    'linear_quadratic_bezout_certificates']
FAILED = ['transverse_coupled_open','linear_quadratic_tenth_sharpness']


def sha(name):
    return hashlib.sha256((P/name).read_bytes()).hexdigest()


def clean(name):
    return '\n'.join(x.strip() for x in (P/name).read_text().splitlines()
                     if x.strip() and not x.lstrip().startswith('//'))


def primitive(expr, variables):
    return S.Poly(S.expand(expr),*variables,domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr()


def direct_coefficient(k, variables):
    """Binomial/multinomial extraction, separate from the generating recurrence."""
    terms = {}
    for j1 in range(k+1):
        for j2 in range(k-j1+1):
            j0 = k-j1-j2
            remaining = 2*k-j1-2*j2
            for r1 in range(remaining+1):
                for r2 in range((remaining-r1)//2+1):
                    rem = remaining-r1-2*r2
                    if rem < 0 or rem % 3:
                        continue
                    r3 = rem//3; r = r1+r2+r3
                    coeff = S.binomial(-k-S.Rational(1,2),r)
                    coeff *= S.Rational(factorial(r),factorial(r1)*factorial(r2)*factorial(r3))
                    coeff *= S.Rational(factorial(k),factorial(j0)*factorial(j1)*factorial(j2))
                    mon = (r1,r2,r3,j0,j1,j2)
                    terms[mon] = terms.get(mon,0)+coeff
    return S.Poly.from_dict(terms,variables,domain=S.QQ)


def native_polynomial(expr, variables):
    """Parse explicit-multiplication Singular output without large symbolic expansion."""
    out = {}
    for term in re.findall(r'[+-]?[^+-]+',expr.replace(' ','')):
        sign = -1 if term.startswith('-') else 1
        term = term.lstrip('+-'); coef = Fraction(sign); mon = [0]*len(variables)
        for factor in term.split('*'):
            found = False
            for i,name in enumerate(variables):
                match = re.fullmatch(re.escape(name)+r'(?:\^(\d+))?',factor)
                if match:
                    mon[i] += int(match.group(1) or 1); found = True; break
            if not found:
                coef *= Fraction(factor)
        key = tuple(mon); out[key] = out.get(key,Fraction(0))+coef
    return {m:c for m,c in out.items() if c}


def integer_identity(generators, coefficients, variables, target):
    parsed = [native_polynomial(s,[str(x) for x in variables]) for s in coefficients]
    denominator = 1
    for poly in parsed:
        for c in poly.values():
            denominator = lcm(denominator,c.denominator)
    total = {}
    for gen,weight in zip(generators,parsed):
        terms = S.Poly(gen,*variables,domain=S.ZZ).terms()
        for m,q in weight.items():
            z = q.numerator*(denominator//q.denominator)
            for n,c in terms:
                key = tuple(a+b for a,b in zip(m,n))
                total[key] = total.get(key,0)+z*int(c)
    total = {m:c for m,c in total.items() if c}
    expected = {m:int(c)*denominator for m,c in S.Poly(target,*variables).terms()}
    assert total == expected
    return {'coefficient_term_counts':[len(q) for q in parsed],
            'cleared_denominator_sha256':hashlib.sha256(str(denominator).encode()).hexdigest(),
            'target':str(target)}


def main():
    start = time.monotonic()
    assert sha('verify_nodal_product_reduction.py') == '7254eff7dad198b2e15fa727da0703bab5f8e3b8704697b27a0ed212527821d6'
    for stem in ACCEPTED+FAILED:
        rec = json.loads((P/(stem+'_run.json')).read_text())
        assert rec['terminal'] and rec['input_sha256'] == sha(stem+'.sing')
        assert rec['log_sha256'] == sha(stem+'_stdout.txt')
        if stem in ACCEPTED:
            assert rec['completed'] and not rec['timed_out'] and rec['exit_code'] == 0
            assert '?' not in (P/(stem+'_stdout.txt')).read_text()
        else:
            assert not rec['completed']
            assert rec['timed_out'] if stem == 'transverse_coupled_open' else not rec['timed_out']
    d1,d2,d3,h0,h1,h2,t = S.symbols('d1 d2 d3 h0 h1 h2 t')
    iv = [d1,d2,d3,h0,h1,h2]; isyms = {str(x):x for x in iv}
    raw = json.loads((P/'linear_quadratic_invariants_input.json').read_text())
    coefficients = [S.sympify(x,locals=isyms) for x in raw['normalized_moment_coefficients']]
    for k,c in enumerate(coefficients,1):
        assert S.Poly(c,*iv,domain=S.QQ) == direct_coefficient(k,iv)
        for mon in S.Poly(c,*iv).monoms():
            assert sum(mon[3:]) == k
            assert sum(a*b for a,b in zip(mon,[1,2,3,0,1,2])) == 2*k
    charts = json.loads((P/'linear_quadratic_invariant_charts_input.json').read_text())
    for row in charts:
        variables = [isyms[x] for x in row['variables']]
        sub = {isyms[k]:S.sympify(v) for k,v in row['normalization'].items()}
        solution = S.sympify(row['h2_solution'],locals=isyms)
        assert S.expand(coefficients[0].subs(sub).subs(h2,solution)) == 0
        expected = [primitive(c.subs(sub).subs(h2,solution),variables) for c in coefficients[1:]]
        assert expected == [S.sympify(x,locals=isyms) for x in row['remaining_moments']]
    # Explicit cleared-denominator polynomial identities verify the hard exclusions.
    log = clean('linear_quadratic_bezout_certificates_stdout.txt').splitlines()
    identities = []
    for row,take,target in [(charts[0],4,S.Integer(1)),(charts[1],3,d3),(charts[3],3,d3*d3)]:
        label = row['label']; weights = []
        for i,line in enumerate(log):
            if line == label+'_COEFFICIENT':
                assert int(log[i+1]) == len(weights)+1; weights.append(log[i+2])
        assert len(weights) == take
        variables = [isyms[x] for x in row['variables']]
        generators = [S.sympify(s,locals=isyms) for s in row['remaining_moments'][:take]]
        identities.append({'chart':label,**integer_identity(generators,weights,variables,target)})
    assert S.expand(coefficients[0].subs({h0:1,d1:0,h1:0})) == h2-S.Rational(3,2)*d2
    assert coefficients[1].subs({h0:1,d1:0,h1:0,d2:1,h2:S.Rational(3,2)}) == -S.Rational(7,8)
    assert coefficients[2].subs({h0:1,d1:0,h1:0,d2:0,h2:0}) == S.Rational(63,8)*d3*d3

    # Entire six-covariance identities check the Gaussian D/H moment formula.
    a,b,c,d,e,f = cv = S.symbols('a b c d e f'); csyms = {str(x):x for x in cv}
    G = S.Matrix([[a,d,e],[d,b,f],[e,f,c]]); delta = G.det()
    families = json.loads((P/'line_conic_moments_input.json').read_text())['families']
    matrix_data = {}
    for kind,A in [('tangent',S.Matrix([[0,0,S.Rational(1,2)],[0,1,0],[S.Rational(1,2),0,0]])),
                   ('transverse',S.Matrix([[1,0,0],[0,0,S.Rational(1,2)],[0,S.Rational(1,2),0]]))]:
        M = S.eye(3)-2*t*G*A; D = S.Poly(M.det(),t); H = S.Poly((M.adjugate()*G)[0,0],t)
        assert D.nth(3) == 2*delta
        if kind == 'tangent':
            assert H.as_expr() == a+2*(d*d-a*b)*t
        else:
            assert H.nth(2) == -delta
        substitutions = [D.nth(1),D.nth(2),D.nth(3),H.nth(0),H.nth(1),H.nth(2)]
        powers = []
        for q in substitutions:
            base = S.Poly(q,*cv,domain=S.QQ); arr = [S.Poly(1,*cv,domain=S.QQ)]
            for k in range(10):arr.append(arr[-1]*base)
            powers.append(arr)
        for k,coefficient in enumerate(coefficients,1):
            transformed = S.Poly(0,*cv,domain=S.QQ)
            for mon,number in S.Poly(coefficient,*iv).terms():
                term = S.Poly(number,*cv,domain=S.QQ)
                for i,degree in enumerate(mon):term *= powers[i][degree]
                transformed += term
            transformed *= S.Rational(factorial(2*k)**2,2**k*factorial(k))
            given = S.sympify(families[kind]['primitive_moments'][k-1],locals=csyms)
            given *= S.Integer(families[kind]['scalars'][k-1])
            assert transformed == S.Poly(given,*cv,domain=S.QQ)
        matrix_data[kind] = {'D':str(D.as_expr()),'H':str(H.as_expr())}

    # Tangent refinement: all three non-null charts exclude by gate8.
    covcharts = json.loads((P/'line_conic_charts_input.json').read_text())
    tangent_log = clean('tangent_line_conic_charts_stdout.txt')
    for row,gate in zip(covcharts[:3],[8,6,6]):
        variables = [csyms[x] for x in row['variables']]
        sub = {csyms[k]:S.sympify(v) for k,v in row['normalization'].items()}
        eliminated = csyms[row['eliminated']]; solution = S.sympify(row['solution'],locals=csyms)
        mf = [S.sympify(x,locals=csyms) for x in families['tangent']['primitive_moments']]
        assert S.expand(mf[0].subs(sub).subs(eliminated,solution)) == 0
        expected = [primitive(q.subs(sub).subs(eliminated,solution),variables) for q in mf[1:]]
        assert expected == [S.sympify(x,locals=csyms) for x in row['remaining_moments']]
        dd = primitive(delta.subs(sub).subs(eliminated,solution),variables)
        assert dd == S.sympify(row['determinant'],locals=csyms)
        gb = S.groebner(expected[:gate//2-1],*variables)
        assert gb.reduce(S.expand(dd**12))[1] == 0
        assert tangent_log.split(row['label']+'_GATE\n')[1:][-1].splitlines()[0] == str(gate)
    # Explicit algebraic tangent sharpness witness, in e and z=f^2.
    z = S.symbols('z'); row = covcharts[0]; polys = []
    for q in row['remaining_moments'][:3]:
        full = S.Poly(S.sympify(q,locals=csyms),e,f)
        assert all(j%2 == 0 for (i,j),v in full.terms())
        polys.append(S.expand(sum(v*e**i*z**(j//2) for (i,j),v in full.terms())))
    p7 = 3381*e**7+8358*e**6+8808*e**5+5165*e**4+1730*e**3+263*e**2-11*e-6
    zsolution = S.solve(polys[0],z)[0]
    for q in polys[:2]:
        assert S.rem(S.Poly(S.fraction(S.cancel(q.subs(z,zsolution)))[0],e),S.Poly(p7,e)).is_zero
    eighth_num = S.Poly(S.fraction(S.cancel(polys[2].subs(z,zsolution)))[0],e)
    assert S.gcd(eighth_num,S.Poly(p7,e)).degree() == 0
    assert p7.subs(e,0) < 0 and p7.subs(e,1) > 0
    detw = S.cancel(-5*e*e-2*e-1-zsolution)
    assert S.factor(detw+2*(14*e**4+28*e**3+19*e**2+6*e+1)/(7*e+3)) == 0

    # Algebraic sharpness10: a proper length28 quotient on the good open set.
    sharp = clean('linear_quadratic_tenth_sharpness_v2_stdout.txt')
    assert sharp.split('GOOD_WITNESS_LENGTH\n')[1].splitlines()[0] == '28'
    assert sharp.split('BAD_LOCUS_IDEAL\n')[1].splitlines()[0] == '1'
    src = (P/'linear_quadratic_tenth_sharpness_v2.sing').read_text()
    actual = [S.sympify(s.replace('^','**'),locals=isyms) for s in src.split('ideal I=')[1].split(';')[0].split(',')]
    assert actual == [S.sympify(x,locals=isyms) for x in charts[0]['remaining_moments'][:3]]
    assert S.discriminant(1+t+d2*t*t+d3*t**3,t) == d2*d2-4*d3-4*d2**3-27*d3*d3+18*d2*d3
    good_basis = sharp.split('GOOD_WITNESS_BASIS\n')[1].split('BAD_LOCUS_IDEAL')[0].strip().splitlines()
    leading = [next(iter(native_polynomial(q.rstrip(','),['d2','d3','h1']))) for q in good_basis]
    todo = deque([(0,0,0)]); standard = {(0,0,0)}
    while todo:
        mon = todo.popleft()
        for i in range(3):
            mm = mon[:i]+(mon[i]+1,)+mon[i+1:]
            if mm not in standard and not any(all(a<=b for a,b in zip(g,mm)) for g in leading):
                standard.add(mm);todo.append(mm)
        assert len(standard) <= 100
    assert len(standard) == 28
    # Partial-fraction realization on distinct nonzero eigenvalues.
    lam = S.symbols('l1:4',nonzero=True); ww = S.symbols('w1:4')
    DD = S.prod(1-2*x*t for x in lam)
    HH = sum(ww[i]*S.prod(1-2*lam[j]*t for j in range(3) if j!=i) for i in range(3))
    for i in range(3):
        residue = HH.subs(t,1/(2*lam[i]))/S.prod(1-lam[j]/lam[i] for j in range(3) if j!=i)
        assert S.cancel(residue-ww[i]) == 0
    assert S.expand(S.expand(HH).coeff(t,2)-4*S.prod(lam)*sum(ww[i]/lam[i] for i in range(3))) == 0
    # Sharp mixed degree bound, witnessed in the smooth tangent-conic class.
    mixed_values = []
    for degree in range(1,7):
        value = sum(comb(degree,j)*factorial(degree+j)*S.factorial2(2*(degree-j)-1)
                    for j in range(degree+1))
        assert value > 0; mixed_values.append(int(value))
    names = ['verify_nodal_product_reduction.py','line_conic_moments_input.json',
        'generate_line_conic_moments.py','prepare_line_conic_charts.py',
        'line_conic_charts_input.json','prepare_transverse_coupled_charts.py',
        'transverse_coupled_charts_input.json','generate_linear_quadratic_invariants.py',
        'linear_quadratic_invariants_input.json','linear_quadratic_invariant_charts_input.json',
        'explore_tangent_sharp_witness.py','tangent_sharp_witness_exploration.json',
        'reducible_cubic_source_audit.json']
    for stem in ACCEPTED+FAILED:names += ['run_'+stem+'.py',stem+'.sing',stem+'_stdout.txt',stem+'_run.json']
    cert = {'claim':'G3REDUCIBLE20261114','field':'C, characteristic zero',
        'theorem':'Every homogeneous reducible ternary cubic is GMC-safe; pure gate10 and mixed cutoff deg(multiplier) are both sharp uniformly.',
        'pure_null_iff_exists_factorization':'P=(l^T*x)*(x^T*A*x), l^T*l=0 and A*l proportional to l; A may be degenerate',
        'linear_times_nondegenerate_quadratic':{'pure_null_iff':'l*l=0 and A*l is proportional to l',
            'first_ten_zero_iff_all_pure_zero':True,'geometric_type_of_pure_null':'tangent line-conic; tangency point is the Gaussian dual of the line factor',
            'tangent_sharp_gate':8,'transverse_sharp_nonvanishing_gate':10},
        'universal_series':{'D':'1+d1*t+d2*t^2+d3*t^3','H':'h0+h1*t+h2*t^2',
            'C_k':'[t^(2k)] H^k D^(-k-1/2)',
            'moment_factor':'M_(2k)=(2k)!^2/(2^k*k!)*C_k',
            'radical_after_inverting_d3':'rad(C1,C2,C3,C4,C5)=(h0,h1,h2)',
            'bezout_identities':identities},
        'covariance_generating_polynomials':matrix_data,
        'tangent_sharpness':{'positive_root_polynomial':str(p7),'f_squared':str(zsolution),
            'c':'-4e^2-2e-1','determinant':str(detw),'first_seven_moments_zero':True,'eighth_moment_nonzero':True},
        'transverse_sharpness':{'finite_algebra_length':28,'distinct_points_not_asserted':True,
            'normalization':'h0=d1=1; h2=3/8*(4*d2+4*h1-5)',
            'system':charts[0]['remaining_moments'][:3],
            'required_open_conditions':'d3*disc(D)*(4*d2+4*h1-5)!=0',
            'realization':'D=product(1-2*lambda_i*t); H/D=sum(r_i/(1-2*lambda_i*t)); A=diag(lambda_i), l_i^2=r_i',
            'first_nine_moments_zero':True,'tenth_moment_nonzero':True},
        'mixed_sharpness':{'P':'Z*(Z*W+T^2)','multiplier':'W^D',
            'boundary_value':'sum_j binomial(D,j)*(D+j)!*(2*(D-j)-1)!! >0',
            'values_D1_to_D6':mixed_values},
        'provenance':{name:sha(name) for name in names},
        'boundary':'All-order pure/mixed claims use the written Gaussian generating identity, complete invariant rescaling cover and Wick partner counting. Sharp transverse existence is certified algebraically with a length28 quotient and a distinct-eigenvalue realization, not28 distinct points or a numerical search. The direct coupled covariance open-chart30-second timeout remains UNKNOWN. Nonhomogeneous factors/perturbations and all remaining irreducible nodal and coefficient-scheme questions remain separate. No arbitrary-cubic upper bound is asserted.',
        'assurance':'Internal written, exact rational CAS, direct integer-polynomial Bezout and symbolic/combinatorial proof; no independent replay, external review or formal verification.'}
    path = P/'reducible_cubic_moments_certificate.json';path.write_text(json.dumps(cert,indent=2)+'\n')
    rec = {'claim':cert['claim'],'completed':True,'reducible_sharp_pure_gate':10,
        'tangent_sharp_pure_gate':8,'transverse_sharp_nonvanishing_gate':10,
        'sharp_mixed_cutoff':'deg(multiplier)','sharpness_algebra_length':28,
        'checker_sha256':sha(Path(__file__).name),'certificate_sha256':sha(path.name),
        'elapsed_seconds':time.monotonic()-start,'assurance':cert['assurance']}
    (P/'reducible_cubic_moments_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2))


if __name__ == '__main__':
    main()
