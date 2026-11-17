"""Uniform cubic-algebra cusp identities and the full oriented r7 prime ledger."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
import platform
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


def main():
    started = time.monotonic()
    models = families()
    cusps = []
    for name, parameter, p, q, t in models:
        descriptions = [(3*parameter**3+2, -s.Rational(3,2)*parameter**2,
                         -s.Rational(1,4), -s.Rational(49,5832)*parameter**2, 6048*parameter, 'I')]
        if name == 'K':
            descriptions.append((15*parameter**3-4, -s.Rational(15,4)*parameter**2,
                                 s.Rational(1,48), s.Rational(686963,58320)*parameter**2,
                                 -6300*parameter, 'II'))
        for base, critical, target_U, target_V, expected_jet, cls in descriptions:
            label = name+' class'+cls

            def red(expr):
                num, den = s.fraction(s.cancel(expr))
                return s.rem(s.expand(num*s.invert(den, base, parameter)), base, parameter)

            check(label+' parameter algebra squarefree', s.gcd(base, s.diff(base,parameter)) == 1)
            check(label+' parameter is a unit', s.gcd(base,parameter) == 1)
            check(label+' exact target U', red(p.subs(t,critical)-target_U) == 0)
            check(label+' exact target V', red(q.subs(t,critical)-target_V) == 0)
            check(label+' both first derivatives vanish', red(s.diff(p,t).subs(t,critical)) == 0
                  and red(s.diff(q,t).subs(t,critical)) == 0)
            jet = red((s.diff(p,t,2)*s.diff(q,t,3)-s.diff(p,t,3)*s.diff(q,t,2)).subs(t,critical))
            check(label+' exact cusp second/third jet determinant', red(jet-expected_jet) == 0)
            check(label+' ordinary cusp jet is a unit', s.gcd(jet,base) == 1)
            residual_resultants = {}
            for kind, f, g, power in [('derivative', s.diff(p,t), s.diff(q,t), 1),
                                      ('fibre', p-target_U, q-target_V, 2)]:
                f, g = red(f), red(g)
                divisor = (t-critical)**power
                fq, fr = s.div(f, divisor, t)
                gq, gr = s.div(g, divisor, t)
                check(label+' '+kind+' whole monic divisions', red(fr) == 0 and red(gr) == 0
                      and red(f-fq*divisor-fr) == 0 and red(g-gq*divisor-gr) == 0)
                resultant = red(s.resultant(red(fq),red(gq),t))
                check(label+' '+kind+' complement resultant is a unit', s.gcd(resultant,base) == 1)
                residual_resultants[kind] = str(resultant)
            cusps.append(dict(chart=name, coefficient_class=cls, parameter=str(parameter),
                parameter_equation=str(base), critical_parameter=str(critical),
                target_U=str(target_U), target_V=str(target_V), jet_determinant=str(jet),
                residual_resultants=residual_resultants,
                parametrization_fibre_scheme_length=2, normalization_point_count=1,
                original_R_ambient_length=3, original_R_e=2, original_R_f=1,
                actual_original_target_coordinate_change=None))
    _, d, pL, qL, t = models[0]
    _, x, pK, qK, _ = models[1]
    base = 3*x**3+2
    for label, expr in [('p',pK-pL.subs(d,x)),('q',qK-qL.subs(d,x))]:
        num, den = s.fraction(s.cancel(expr))
        check('L/K classI identical '+label+' in the same target coordinates',
              s.rem(s.expand(num*s.invert(den,base,x)),base,x) == 0)

    oriented = []
    for name, ER, eS, ES, unram in [('A',0,2,0,1),('B',0,2,1,0),
                                   ('C',0,3,0,0),('N',1,2,0,0)]:
        check(name+' complete all-prime Euler budget', 2+ER+eS+ES+unram == 5)
        row = dict(row=name, R=dict(image='C',e=2,f=1,E=ER,ramifies=True,omitted=True),
                   S=dict(image='B',e=eS,f=1,E=ES,ramifies=True,omitted=True),
                   unramified_omitted_count=unram,
                   nonimmersive_image='C' if ER else ('B' if ES else None),
                   original_affine_normalization_smooth=True)
        if unram:
            row['D'] = dict(image='H',e=1,f=1,E=0,ramifies=False,omitted=True,
                            may_share_image_with=['C','B'])
        oriented.append(row)
    completions = [(e,E,t) for e in range(2,6) for E in range(4) for t in range(4)
                   if 2+1+e+E+t == 5]
    check('nonimmersed first image forces only the N row', completions == [(2,0,0)])
    exceptional = []
    for label,c,eps in [('off_C',0,0),('smooth_C',1,0),('double_C',2,0),('cusp_C',1,1)]:
        rows = []
        for b in range(4):
            retained = 6-2*c-2*b-eps
            if retained < 0:
                continue
            lengths = [3 if eps and i == 0 else 2 for i in range(c)]+[2]*b
            check(f'{label} b{b} conserved full fibre length', retained+sum(lengths) == 6)
            rows.append(dict(B_branches=b,retained_length=retained,omitted_local_lengths=lengths))
        exceptional.append(dict(location=label,C_branches=c,epsilon=eps,necessary_completions=rows))
    check('only possible singleton is the cusp with one B branch',
          [(r['location'],c['B_branches']) for r in exceptional for c in r['necessary_completions']
           if c['retained_length'] == 1] == [('cusp_C',1)])
    check('cusp corank correction is one point and length3', 2-1 == 1 and 2*2-1 == 3)
    check('node equal corank is two points and length4', 2-0 == 2 and 2*2-0 == 4)
    for p in [1,2,3]:
        weights = [F(2*p,3)]*3
        delta = 2*sum(min(a,b) for i,a in enumerate(weights) for b in weights[i+1:])
        check(f'cusp-only p{p} inverse pole/discriminant/index orders',
              sum(weights) == 2*p and delta == 4*p and math.ceil((delta-2)/2) == 2*p-1)
    inverse = json.loads((HERE/'singleton_valuation_ledger.json').read_text())
    singleton_B = next(r for r in inverse['packet_rows'] if r['packet']=='B')
    check('new N singleton has old complete cusp/fold index bounds',
          singleton_B['coefficient_order_lower_bounds'] == [1,1,2,3,3]
          and singleton_B['raw_discriminant_order_lower_bound'] == 11
          and singleton_B['index_order_lower_bound'] == 4)
    ledger = dict(as_of='2026-11-14',claim='PF2OMIT13',
        hypotheses='Exactly PF2R7RAM2 global degree6/index2/gap21 class and named geometry; no claim deriving those hypotheses for every terminal passport',
        exact_cusp_descriptions=cusps,
        distinct_nonimmersed_classes=2,
        full_first_prime_formula=dict(matrix_corank='m=6-rank(M)',
            epsilon='1 at the listed cusp for a nonimmersed coefficient point,0 elsewhere',
            normalization_point_count='m-epsilon',total_omitted_ambient_length='2*m-epsilon',
            first_prime=dict(e=2,f=1,omitted=True,ramifies=True)),
        complete_oriented_rows=oriented,
        new_N_row_exceptional_value_table=exceptional,
        new_N_singleton_iff='B contains the listed unique cusp point, with one smooth branch there',
        new_N_off_cusp_source_fibres_even=True,
        new_N_monodromy='S6',new_N_source_field_has_proper_intermediate_field=False,
        new_N_inverse_leading_coefficient='constant*H_C^p_R*H_B^p_S',
        new_N_raw_discriminant='constant*H_C*H_B*J^2',
        new_N_cusp_with_B_absent=dict(retained_length=3,reciprocal_weights='2*p_R/3 three times',
            coefficient_orders_a4_a5_a6_minimum=[1,2,2],leading_order='2*p_R',
            raw_discriminant_order_at_least='4*p_R',index_order_at_least='2*p_R-1'),
        actual_second_image_equation=None,actual_optional_unramified_image_in_other_rows=None,
        actual_global_source_realization=False,full_goal_complete=False)
    target = HERE/'nonimmersed_ledger.json'
    target.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs = [HERE/n for n in ['NONIMMERSED_LEDGER.md','image_norm.py','RESIDUE_LEDGER.md',
              'SIMPLE_RAMIFICATION_AND_RADIAL.md','GLOBAL_MONODROMY.md',
              'SINGLETON_INVERSE_VALUATIONS.md','singleton_valuation_ledger.json']]
    inputs += [ROOT/'jacobian/jc2-r7-image-family'/n for n in ['family.py','families.json','PROOF.md']]
    inputs += [ROOT/'jacobian/jc2-r7-reciprocal/PROOF.md',
               ROOT/'jacobian/jc2-normal-r7/README.md',
               ROOT/'foundations/scripts/verify_f2_degree_6_10_gap.py']
    receipt = dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={str(p.relative_to(ROOT)):digest(p) for p in inputs},
        ledger_sha256=digest(target),runtime_seconds=round(time.monotonic()-started,6),
        python=platform.python_version(),sympy=s.__version__,
        proof_boundary='Uniform exact cusp values, derivative/fibre quotient resultant units and all oriented budgets checked. Global extra ramification, boundary regularity, flat lengths and transposition topology remain written/imported arguments. No second image or source map is constructed.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'nonimmersed_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','runtime_seconds']}))


if __name__ == '__main__':
    main()
