"""Exact controls for leading-valuation factor degrees and critical-curve gates."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(S.expand(item))==0 for item in values))

    T,V,Z,U,Y,z,x,y,v,s,t,w,r,tau = S.symbols('T V Z U Y z x y v s t w r tau')
    m = 6
    leading = T*T*(T-1)**3
    polynomial = leading*(z**6+T*z**4+z+T*T)+T**4+V

    def order(value,level):
        if value==0:
            return math.inf
        return min(monomial[0] for monomial,coefficient in S.Poly(value.subs(T,level+r),r).terms())

    coefficients = S.Poly(polynomial,z)
    for level,q in [(0,2),(1,3)]:
        check(f'leading order at level {level} is exactly {q}',order(coefficients.LC(),level)==q)
        check(f'all positive coefficients satisfy the single-slope inequalities at level {level}',
            all(order(coefficients.nth(i),level)>=S.Rational(i*q,m) for i in range(1,m)))
        check(f'the generic constant term has valuation zero at level {level}',order(coefficients.nth(0),level)==0)
        zero(f'the binomial root realizes the exact slope minus {q}/6',
            (T**q*Z**m-1).subs({T:tau**m,Z:tau**(-q)},simultaneous=True))
    check('the first valuation alone permits only degrees three below six',
        [e for e in range(1,6) if e*2%6==0]==[3])
    check('the second valuation alone permits only degrees two and four below six',
        [e for e in range(1,6) if e*3%6==0]==[2,4])
    check('the two valuation constraints together permit no proper factor degree',
        not [e for e in range(1,6) if e*2%6==0 and e*3%6==0])
    zero('a noncoprime single slope really can split into the permitted degree-three factors',
        (T*Z**3-1)*(T*Z**3+1)-(T*T*Z**6-1))
    check('the uniform arithmetic bound uses the gcd of all leading orders',
        math.lcm(6//math.gcd(6,2),6//math.gcd(6,3))==6//math.gcd(6,2,3)==6)

    a,b,c,d = S.symbols('a b c d')
    cubic = a*z**3+b*z*z+c*z+d+V
    Q = S.cancel(a*a*cubic.subs(z,U/a))
    zero('the normalized cubic root model is monic with no denominators',
        Q-(U**3+b*U*U+a*c*U+a*a*(d+V)))
    check('its U derivative is independent of the second target parameter',not S.diff(Q,U).has(V))
    critical_D = b*b-3*a*c
    zero('the normalized critical curve is exactly Y squared equals the critical discriminant',
        3*S.diff(Q,U)-((3*U+b)**2-critical_D))
    root_discriminant = S.discriminant(cubic,z)
    zero('the V discriminant of the cubic root discriminant is sixteen times the critical discriminant cubed',
        S.discriminant(root_discriminant,V)-16*critical_D**3)
    check('the two discriminants are different objects',root_discriminant.has(V) and not critical_D.has(V))

    P_partial = T*z**3+z
    zero('the partial-drop S3 cubic has its exact nonsquare root discriminant',
        S.discriminant(P_partial+V,z)-(-4*T-27*T*T*V*V))
    check('that root discriminant has two simple V zeros over the parameter field',
        S.discriminant(S.discriminant(P_partial+V,z),V)!=0)
    check('the partial-drop cubic meets the independent square-free-leading top-coefficient gate',
        S.gcd(T,S.diff(T,T))==1 and S.rem(S.Poly(P_partial,z).nth(2),T,T)==0)
    check('that example does not meet the full-degree-drop valuation inequalities',
        order(S.Poly(P_partial,z).nth(1),0)<S.Rational(1,3))

    P = T*z**3+z*z-T**3*z/3
    D = 1+T**4
    zero('the complementary cubic has the positive-genus critical discriminant',
        S.Poly(P,z).nth(2)**2-3*S.Poly(P,z).LC()*S.Poly(P,z).nth(1)-D)
    check('its critical quartic has four distinct roots',S.degree(D,T)==4 and S.gcd(D,S.diff(D,T))==1)
    check('its quadratic coefficient fails the independent divisible-top-coefficient condition',S.rem(1,T,T)==1)
    check('its leading zero is only a drop to quadratic degree',S.degree(P.subs(T,0),z)==2)
    zero('the new S3 cubic has a nonzero exact V discriminant of its root discriminant',
        S.discriminant(S.discriminant(P+V,z),V)-16*D**3)
    zero('the degree-two critical curve has a nonconstant unit obstruction to A1 parametrization',
        (Y-T)*(Y+T)-1-(Y*Y-T*T-1))
    zero('the degree-one critical curve instead has an actual polynomial parametrization',
        (Y*Y-T).subs({T:r*r,Y:r},simultaneous=True))
    zero('the singular square factor does not change the normalized critical function field',
        (Y*Y-T*T*(T*T+1)).subs(Y,T*U)-T*T*(U*U-T*T-1))

    f_bad = -s+s**3*t
    g_bad = -t
    h_bad = 1/s
    P_bad = T*z**3+z*z
    normalized_label = S.cancel(f_bad*h_bad)
    zero('a genuine rational pole can occur when the top-coefficient gate fails outside Keller',
        P_bad.subs({T:f_bad,z:h_bad})+g_bad)
    check('the pole is not silently called a polynomial root',S.denom(h_bad)==s)
    check('the label becomes polynomial after multiplying by the actual leading coefficient',
        normalized_label==-1+s*s*t)
    J_bad = S.Matrix([f_bad,g_bad]).jacobian([s,t]).det()
    check('this rational-pole control explicitly fails the unit-Jacobian hypothesis',J_bad==1-3*s*s*t)
    Q_bad = U**3+U*U+T*T*V
    zero('the global vertical-Jacobian identity holds with its complete non-Keller right side',
        S.Matrix([f_bad,normalized_label]).jacobian([s,t]).det()
        *S.diff(Q_bad,U).subs({T:f_bad,U:normalized_label,V:g_bad})+f_bad*f_bad*J_bad)
    zero('the square critical discriminant has two distinct rational normalized critical factors',
        S.diff(Q_bad,U)-3*U*(U+S.Rational(2,3)))
    check('the second required critical unit really fails away from the leading divisor in the non-Keller control',
        (normalized_label+S.Rational(2,3)).subs({s:1,t:S.Rational(1,3)})==0
        and f_bad.subs({s:1,t:S.Rational(1,3)})!=0)
    check('that second-unit failure is exactly a zero of the non-Keller Jacobian',
        J_bad.subs({s:1,t:S.Rational(1,3)})==0)
    f_nonreduced = -s*s-s**3*t
    zero('dropping square-free source pullback also permits a rational pole even with zero quadratic coefficient',
        P_partial.subs({T:f_nonreduced,z:h_bad})+t)
    check('that second control has a multiple leading-zero source component and is not Keller',
        S.factor(f_nonreduced)==-s*s*(s*t+1)
        and S.Matrix([f_nonreduced,t]).jacobian([s,t]).det().subs(s,0)==0)

    relation = x*y-P-V
    gb = S.groebner([relation],x,y,z,v,domain=S.QQ.poly_ring(T,V),order='lex')
    check('the newly closed cubic retains the global monic xy basis',
        gb.polys[0].LM().exponents==(1,1,0,0) and gb.polys[0].LC()==1)
    zero('its leading-zero height-one fibre is nonzero even after the cubic drops degree',
        relation.subs(T,0)-(x*y-z*z-V))
    image = {T:s,x:1,y:t+P.subs({T:s,z:w}),z:w,v:v}
    zero('the new S3 family has a genuine complete polynomial source inverse',
        S.Matrix([T,y-P,z,v]).subs(image,simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the carried output pair is exactly the identity Keller base',
        S.Matrix([T,x*y-P]).subs(image,simultaneous=True)-S.Matrix([s,t]))
    actual = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('the independent nondegenerate ambient quadratic has its literal pullback',
        (T*v+y*z+x*x).subs(image,simultaneous=True)-actual)
    check('the actual new S3 control has unit Hessian',S.hessian(actual,[s,t,w,v]).det()==1)
    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    evaluate = {T:p4,z:p2}
    inverse = S.Matrix([p4,p3-P.subs(evaluate)-p2*S.diff(P,z).subs(evaluate),
        p2,p1-p2*S.diff(P,T).subs(evaluate)])
    gradient = S.Matrix([S.diff(actual,q) for q in [s,t,w,v]])
    zero('the displayed actual gradient inverse is a left inverse',
        inverse.subs(dict(zip([p1,p2,p3,p4],gradient)),simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the displayed actual gradient inverse is a right inverse',
        gradient.subs(dict(zip([s,t,w,v],inverse)),simultaneous=True)-S.Matrix([p1,p2,p3,p4]))

    retained = {
        'jacobian-followup/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md':'3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
        'jacobian-followup/hc4-triangular-graphs/LEADING_FIBRE_UNIT_GATE.md':'1d054ba8762ea657be879b6e3568aa58c9ccb65e14ab642d45e0865e8f6382f0',
        'jacobian-followup/hc4-triangular-graphs/NORMAL_MONOGENIC_COVERS.md':'223b81738c823c3658079344f57e74a1405b3fcb52eb8adc76c9ee170226226f',
        'jacobian-followup/hc4-triangular-graphs/sources/stacks-leading-riemann-hurwitz-0C1B.html':'616e526d775d71ae3dd706bccfdc38115a846c7eafac1e7080d53691bc2249a3'
    }
    context = json.loads((HERE/'sources/valuation-critical-import.json').read_text())
    for item in context['sources']:
        retained['jacobian-followup/hc4-triangular-graphs/sources/'+item['retained_file']] = item['sha256']
    for name,digest in retained.items():
        check('retained proof or primary input byte pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_valuation_critical.py','LEADING_VALUATIONS_AND_CRITICAL_CURVES.md','sources/valuation-critical-import.json']},
        dependencies=retained,
        written_not_finite_exhaustion=[
            'leading valuation orders constrain every source factor degree by m/gcd(m,q1,...,qk)',
            'coprime combined leading orders force full root transitivity and close every carried output embedding',
            'a rational root has a polynomial normalized label with an exact vertical-Jacobian identity',
            'every normalized critical component must be polynomially parametrizable by the explicit nonproperness input',
            'distinct rational normalized critical roots require independent localized-unit defect classes',
            'square-free leading coefficient dividing the next coefficient excludes rational roots in all degrees',
            'cubic outputs close if the critical discriminant has square-free part degree at least two or the independent no-pole gate holds'],
        nonproperness_parametric_curve_theorem_fully_reproved=False,
        primary_nonproperness_degree_bound_used=False,
        critical_discriminant_confused_with_cubic_root_discriminant=False,
        rational_root_exclusion_asserted_full_irreducibility_in_arbitrary_degree=False,
        all_nonmonic_cubic_outputs_asserted_closed=False,
        all_Keller_component_fibres_asserted_irreducible=False,
        full_splitting_group_preservation_asserted=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-valuation-critical.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
