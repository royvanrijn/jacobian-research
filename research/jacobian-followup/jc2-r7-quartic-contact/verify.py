"""Exact quartic contact minors and independent Fraction unit replay.

No Groebner-basis discovery or numerical parameter samples are needed.
The inherited geometry and polynomial-base reduction are written proof.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib
import json
import time

import sympy as S
from sympy.polys.rings import ring

from contact import R,u,v,w,z,t,contact_matrix,selected_minor,z_coefficient

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(name,value):
    if not value: raise AssertionError(name)
    CHECKS.append(name)


def zero(name,value):
    check(name,S.expand(value)==0)


def trim(p):
    p=list(p)
    while p and p[-1]==0: p.pop()
    return p


def add(a,b):
    out=[Fraction(0)]*max(len(a),len(b))
    for i,c in enumerate(a):out[i]+=c
    for i,c in enumerate(b):out[i]+=c
    return trim(out)


def multiply(a,b):
    if not a or not b:return []
    out=[Fraction(0)]*(len(a)+len(b)-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b):out[i+j]+=c*d
    return trim(out)


def fractions(a):return list(map(Fraction,a))


def unit_value(f,g,a,b):return add(multiply(fractions(f),fractions(a)),multiply(fractions(g),fractions(b)))


def polynomial(coeffs,var):return sum(S.Rational(c)*var**i for i,c in enumerate(coeffs))


def direct_remainder_check(polynomials,raw):
    # A separate multivariate long division implementation, rather than
    # the scalar-coefficient power recurrence used by contact.py.
    Q,tt,uu,vv,ww,zz=ring('t,u,v,w,z',S.QQ)
    divisor=(tt**4+uu*tt**3+vv*tt**2+ww*tt-zz)**2
    for i,(f,expected) in enumerate(zip(polynomials,raw)):
        _,remainder=Q.from_expr(f).div(divisor)
        assembled=sum((S.Symbol('t')**j*c.as_expr() for j,c in enumerate(expected)),S.Integer(0))
        zero('independent long division remainder '+str(i),remainder.as_expr()-assembled)


def main():
    start=time.monotonic()
    certificates=json.loads((HERE/'certificates.json').read_text())
    classes=json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases']
    check('exactly two retained coefficient classes',len(classes)==2)
    check('exactly four frozen branch certificates',len(certificates['branches'])==4)
    k_patterns=[(k,h0,l0,m0) for k,h0,l0,m0 in product(range(1,4),range(6),range(7),range(7))
                if k*(h0+1-l0)==14 and k*(h0+1-m0)==10]
    check('no smaller local degree for quartic base',k_patterns==[])
    weakened=[(k,h0,l0,m0) for k,h0,l0,m0 in product(range(1,4),range(7),range(7),range(7))
              if k*(h0+1-l0)==14 and k*(h0+1-m0)==10]
    check('degree bound mutation restores an escape',(2,6,0,2) in weakened)
    check('three leading coefficient ratios pairwise distinct',len({Fraction(10,3),Fraction(1),Fraction(5,3)})==3)
    output=[]
    for case,row in enumerate(classes,1):
        p,q=S.sympify(row['p']),S.sympify(row['q'])
        matrix,raw=contact_matrix(p,q)
        direct_remainder_check((q,p*p,p*q,q*q),raw)
        spec=certificates['cubic_coefficient_minor']
        leading=z_coefficient(selected_minor(matrix,spec['rows']),spec['z_power']).as_expr()
        zero(f'class{case} cubic coefficient forces u0',leading-S.sympify(spec['identity']))
        check(f'class{case} u5 coefficient is nonzero',leading!=0)
        atzero=[[f.compose(u,R.zero) for f in row] for row in matrix]
        principal=selected_minor(atzero,[0,1,2,3])
        if case==1: expected=-4*(90*S.Symbol('v')+31)*(9*S.Symbol('w')+4)**2/243
        else: expected=-4*(90*S.Symbol('v')+157)*(45*S.Symbol('w')+8)**2/6075
        zero(f'class{case} two complete parameter branches',z_coefficient(principal,10).as_expr()-expected)
        branch_records=[]
        items=[c for c in certificates['branches'] if c['case']==case]
        expected_values=({'fixed_v':('v','w',S.Rational(-31,90)),
                          'fixed_w':('w','v',S.Rational(-4,9))} if case==1 else
                         {'fixed_v':('v','w',S.Rational(-157,90)),
                          'fixed_w':('w','v',S.Rational(-8,45))})
        check(f'class{case} both complete factor branches covered',
              len(items)==2 and {item['branch'] for item in items}==set(expected_values)
              and all((item['fixed_parameter'],item['free_parameter'],S.Rational(item['fixed_value']))
                      ==expected_values[item['branch']] for item in items))
        for item in items:
            variable=v if item['fixed_parameter']=='v' else w
            free=S.Symbol(item['free_parameter'])
            value=R.ground_new(S.Rational(item['fixed_value']))
            specialized=[[f.compose(variable,value) for f in row] for row in atzero]
            det=selected_minor(specialized,item['rows'])
            for i,power in enumerate(item['z_powers']):
                actual=z_coefficient(det,power).as_expr()
                zero(f'class{case} {item["branch"]} exact coefficient Z{power}',
                     actual-polynomial(item['coefficient_arrays'][i],free))
            f,g=item['coefficient_arrays'];aa,bb=item['bezout_arrays']
            check(f'class{case} {item["branch"]} separate Fraction unit',unit_value(f,g,aa,bb)==[Fraction(1)])
            tampered=list(aa);tampered[0]=str(Fraction(tampered[0])+1)
            check(f'class{case} {item["branch"]} tampered certificate rejected',unit_value(f,g,tampered,bb)!=[Fraction(1)])
            branch_records.append(dict(branch=item['branch'],fixed_parameter=item['fixed_parameter'],fixed_value=item['fixed_value'],unit_identity=True))
        output.append(dict(case=case,u_forced_zero=True,branches=branch_records))

    # The rank criterion must accept actual conic identities. These low
    # degree controls are independent of every cusp source coefficient.
    good_p=t;good_q=t*t
    bad_q=t*t+t**3
    columns=lambda pp,qq:S.Matrix([[S.Poly(f,t).nth(j) for f in [1,pp,qq,pp*pp,pp*qq,qq*qq]] for j in range(8)])
    check('positive conic relation has contact rank below six',columns(good_p,good_q).rank()<6)
    check('nonconic low degree mutation has full contact rank',columns(good_p,bad_q).rank()==6)
    zero('positive conic has identically double contact',good_q-good_p**2)

    inputs=[HERE/'PROOF.md',HERE/'README.md',HERE/'REVIEW.md',HERE/'certificates.json',HERE/'contact.py',Path(__file__).resolve(),
            HERE.parent/'jc2-normal-r7/classification.json',
            HERE.parent/'jc2-r7-universal-contact/PROOF.md',
            HERE.parent/'jc2-r7-cubic-cover/PROOF.md']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    receipt=dict(schema='r7-quartic-contact-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),input_sha256=hashes,cases=output,
                 all_complex_quartic_coefficients_excluded_for_both_classes=True,
                 geometry='Written proof of quartic polynomial base and inherited conic/field identification; not an output of this checker.',
                 assurance='Exact symbolic minors plus separate Fraction polynomial-unit replay; no fresh Groebner calculation, numerical roots, modular promotion, source finality, formal verification or independent external review.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} QUARTIC CONTACT CHECKS ({receipt["elapsed_seconds"]}s)')
    print('Both exact cusp pairs: every quartic polynomial parameter excluded.')


if __name__=='__main__':main()
