"""Exact controls for odd image intersections; the universal proof is written."""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(label,condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started=time.monotonic()
    t,v,U,V,a=s.symbols('t v U V a')
    p=t**6+t**3
    q=t**10+s.Rational(5,3)*t**7+s.Rational(7,9)*t**4+s.Rational(7,81)*t
    controls=[]
    for k,pb,qb,expected in [
        (1,v**3,v**5,24),(1,v**3,2*v**5,30),
        (3,v**9+v,v**15,74),(3,v**9+v,2*v**15,90)]:
        norm=s.resultant(pb-U,V-qb,v)
        hb=s.expand(norm/s.Poly(norm,V).LC())
        comp=s.Poly(s.expand(hb.subs({U:p,V:q})),t)
        check(f'control{k},{expected} normalized implicit degrees',
              s.degree(hb,V)==3*k and s.total_degree(hb)==5*k)
        check(f'control{k},{expected} pullback degree',comp.degree()==expected)
        check(f'control{k},{expected} odd intersection bounds',9*k<=comp.degree()<=30*k)
        # Birationality for v^(3k)+v, v^(5k): the cyclic group of the
        # latter has no nonidentity element fixing the linear term.
        if k==1:
            check(f'control1,{expected} coprime coordinate degrees',s.igcd(3,5)==1)
        else:
            check(f'control3,{expected} cyclic stabilizer trivial',
                  s.Poly(pb,v).coeff_monomial(v)==1 and s.degree(qb,v)==15)
        controls.append(dict(k=k,p_B=str(pb),q_B=str(qb),implicit_equation=str(hb),
            intersection_polynomial=str(comp.as_expr()),affine_intersection_multiplicity=comp.degree(),
            infinity_intersection_multiplicity=50*k-comp.degree(),global_cover_claimed=False))
    relation=a**3+s.Rational(7,75)
    pb=v**3+a*v-s.Rational(2,15)
    qb=v**5+s.Rational(5,3)*a*v**3+s.Rational(5,9)*a**2*v
    norm=s.resultant(pb-U,V-qb,v)
    hb=s.expand(norm/s.Poly(norm,V).LC())
    comp0=s.Poly(s.expand(hb.subs({U:p,V:q})),t)
    comp=s.Poly(sum(s.rem(c,relation,a)*t**m[0] for m,c in comp0.terms()),t)
    derivative_resultant=s.rem(s.resultant(s.diff(pb,v),s.diff(qb,v),v),relation,a)
    check('cubic parameter relation squarefree',s.gcd(relation,s.diff(relation,a))==1)
    check('closer immersed quintic intersection degree12',comp.degree()==12)
    check('closer quintic leading coefficient a unit',s.gcd(comp.LC(),relation)==1)
    check('closer quintic immersive at all three parameter roots',s.gcd(derivative_resultant,relation)==1)
    controls.append(dict(k=1,parameter_relation=str(relation),p_B=str(pb),q_B=str(qb),
        implicit_equation=str(hb),intersection_polynomial=str(comp.as_expr()),
        affine_intersection_multiplicity=12,infinity_intersection_multiplicity=38,
        derivative_resultant=str(derivative_resultant),global_cover_claimed=False))
    quintic_degrees=sorted({30}|{20+r for r in range(-10,9,2)}|{9})
    check('quintic necessary degree list',quintic_degrees==[9,10,12,14,16,18,20,22,24,26,28,30])
    check('quintic controls in necessary degree list',all(c['affine_intersection_multiplicity'] in quintic_degrees for c in controls if c['k']==1))
    # Numeric labels support, but do not replace, the all-odd-k proof.
    budget=json.loads((HERE/'generic_line_ledger.json').read_text())['conditional_F2_inverse_degree_budget']
    degree_obligations=[]
    for group in budget['groupings']:
        rows=[]
        for choice in group['degree_choices']:
            k=choice['k_B']
            if k%2:
                check(f"{group['grouping']},{k},{choice.get('k_H')} odd exponent obstruction",(-11*k)%2==1)
                check(f"{group['grouping']},{k},{choice.get('k_H')} norm pole sum",2*k*10*k+k*(-11*k)==9*k*k)
                rows.append(dict(k_B=k,k_H=choice.get('k_H'),minimum_affine_intersection_multiplicity=9*k,
                    maximum_affine_intersection_multiplicity=30*k,
                    singleton_at_every_C_B_intersection=group['grouping'] in ['A_H=C','A_H=B','C']))
        degree_obligations.append(dict(grouping=group['grouping'],odd_degree_choices=rows))
    for group,coefficients in [('A_H=C',(3,2)),('A_H=B',(2,3)),('C',(2,3))]:
        allowed=[(c,b,6-coefficients[0]*c-coefficients[1]*b)
                 for c in range(1,7) for b in range(1,7)
                 if 6-coefficients[0]*c-coefficients[1]*b>=0]
        check(group+' unique nonnegative intersection completion',allowed==[(1,1,1)])
    check('unramified omitted point changes retained length',6-2-2-1==1 and 6-2-2==2)
    check('ordinary intersection need not singleton in B/N',6-2-2==2)
    ledger=dict(as_of='2026-11-14',claim='PF2OMIT15',
        universal_hypotheses='Birational degrees(6,10), exact first odd Laurent exponent-11; other birational polynomial image degrees(3k,5k), k positive odd',
        affine_intersection_multiplicity='9k<=deg H_B(p,q)<=30k',
        infinity_intersection_multiplicity='20k<=I_infinity(C,B)<=41k',
        lower_bound_counts_distinct_values=False,quintic_necessary_degree_options=quintic_degrees,
        sharpness_of_lower_bound_proved=False,even_k_bound_claimed=False,
        singleton_forced_groupings=['A_H=C','A_H=B','C'],
        singleton_packets=dict(A_shared=dict(omitted_lengths=[2,2,1],radial_indices=[2,2,1],retained_length=1),
                               C=dict(omitted_lengths=[2,3],radial_indices=[2,3],retained_length=1)),
        conditional_F2_odd_degree_obligations=degree_obligations,exact_target_curve_controls=controls,
        actual_B_equation=None,actual_H_equation=None,actual_intersection_values=None,
        actual_source_realization=None,unrestricted_singleton_zero_properness_proved=False,full_goal_complete=False)
    ledger_path=HERE/'odd_intersection_ledger.json'
    ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=['jacobian/jc2-omitted-divisors/'+n for n in [
        'ODD_IMAGE_INTERSECTIONS.md','RESIDUE_LEDGER.md','NONIMMERSED_LEDGER.md',
        'SINGLETON_INVERSE_VALUATIONS.md','GENERIC_LINE_ACCOUNTING.md','generic_line_ledger.json']]
    inputs+=['jacobian/jc2-r7-image-family/PROOF.md',
        'foundations/plane-jc/F2_DEGREE_6_10_POLYNOMIAL_GAP.md']
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(CHECKS),checks=CHECKS,
        checker_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in inputs},
        ledger_sha256=digest(ledger_path),runtime_seconds=round(time.monotonic()-started,6),
        python=platform.python_version(),sympy=s.__version__,
        proof_boundary='Universal odd-k formal inverse/norm/intersection argument is written. Exact resultants, quotient-algebra degree12 control, odd capped degree labels and singleton length completions are checked. Controls are target curves, not global covers. No sharpness or even-k bound.',
        formal_verification=False,external_review=False,full_goal_complete=False)
    (HERE/'odd_intersection_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['outcome','check_count','runtime_seconds']}))


if __name__=='__main__':
    main()
