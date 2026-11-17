"""Exact controls for regular-parameter class monodromy and all polynomial outputs.

The closed-polynomial theorem is an explicit primary import. Neither that
theorem nor arbitrary-degree Galois disjointness is proved by finite controls.
"""
import argparse
import hashlib
import itertools
import json
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

    s,t,w,v,T,V,x,y,z,Y,c = S.symbols('s t w v T V x y z Y c')
    f = s+s*s*t
    zero('a nonsingular parameter has an explicit global gradient Bezout identity',
        (1-2*s*t)*S.diff(f,s)+4*t*t*S.diff(f,t)-1)
    zero('its special zero fibre is the stated product',f-s*(1+s*t))
    zero('the two special components are disjoint by a polynomial ideal certificate',
        (1+s*t)-t*s-1)
    zero('its generic fibre has literal Laurent recovery',f.subs(t,(c-s)/(s*s))-c)
    zero('the generic inverse of s is explicit on the fibre',
        (s*(1+s*t)/c-1).subs(t,(c-s)/(s*s)))
    h = s+t*t
    composite = h**3-h
    zero('a nonlinear polynomial composition creates an actual critical point',
        S.Matrix([S.diff(composite,u) for u in [s,t]]).subs({s:1/S.sqrt(3),t:0}))

    cubic = z**3-3*z+V
    zero('the noncyclic cubic has the exact discriminant with two simple branch values',
        S.discriminant(cubic,z)-27*(4-V*V))
    factors = S.factor_list(S.discriminant(cubic,z),V)[1]
    check('both discriminant zeros have odd valuation over the complex rational field',
        len(factors)==2 and all(S.degree(factor,V)==1 and exponent==1 for factor,exponent in factors))
    check('the root-cover equation is primitive linear in V',S.Poly(cubic,V).LC()==1 and S.degree(cubic,V)==1)
    zero('that triangular equation has a literal polynomial inverse',
        (V-(z**3-3*z)).subs(V,cubic)-V)
    action = []
    for p in itertools.permutations(range(3)):
        columns = []
        for j in range(2):
            columns.append(S.eye(2)[:,p[j]] if p[j]<2 else S.Matrix([-1,-1]))
        action.append(S.Matrix.hstack(*columns))
    check('the actual quotient-lattice root action is nonabelian',
        any(g*h0!=h0*g for g in action for h0 in action))
    check('the full S3 root action has no nonzero invariant class over the rationals',
        S.Matrix.vstack(*(g-S.eye(2) for g in action)).rank()==2)
    alpha = S.Matrix([1,0])
    orbit = {tuple(g*alpha) for g in action}
    check('a proper root stabilizer still fixes a nonzero component class',
        len(orbit)==3 and sum(g*alpha==alpha for g in action)==2)
    zero('the nonzero conjugate class orbit sums to zero',
        sum((S.Matrix(value) for value in orbit),S.zeros(2,1)))

    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    polynomials = [z+2,z**3-3*z,z**5+2*z*z+7]
    controls = []
    for P in polynomials:
        m = S.degree(P,z)
        relation = x*y-P-V
        gb = S.groebner([relation],x,y,z,v,domain=S.QQ.poly_ring(T,V),order='lex')
        check(f'polynomial degree {m}: the full-target relation has a monic xy leading term',
            gb.polys[0].LM().exponents==(1,1,0,0) and gb.polys[0].LC()==1)
        standard = x*x*z**(m+1)+y**3*z+v*z
        zero(f'polynomial degree {m}: standard monomials survive in the common free basis',
            gb.reduce(standard)[1]-standard)
        singular = S.groebner([relation,x,y,S.diff(P,z)],x,y,z,domain=S.QQ.frac_field(V),order='lex')
        check(f'polynomial degree {m}: the geometric generic Danielewski equation is smooth',
            singular.contains(S.Integer(1)))
        Pw = P.subs(z,w)
        image = {T:s,x:1,y:t+Pw,z:w,v:v}
        zero(f'polynomial degree {m}: the embedding has literal complete source recovery',
            S.Matrix([T,y-P,z,v]).subs(image,simultaneous=True)-S.Matrix([s,t,w,v]))
        actual = 1+s*w+t*v+Pw*v
        check(f'polynomial degree {m}: an independent ambient quadratic gives actual unit HC4 Hessian',
            S.hessian(actual,[s,t,w,v]).det()==1)
        inverse = S.Matrix([p3-S.diff(P,z).subs(z,p1)*p2,p4-P.subs(z,p1),p1,p2])
        gradient = S.Matrix([S.diff(actual,u) for u in [s,t,w,v]])
        zero(f'polynomial degree {m}: the written gradient inverse is a left inverse',
            inverse.subs(dict(zip([p1,p2,p3,p4],gradient)),simultaneous=True)-S.Matrix([s,t,w,v]))
        zero(f'polynomial degree {m}: the same formula is a right inverse',
            gradient.subs(dict(zip([s,t,w,v],inverse)),simultaneous=True)-S.Matrix([p1,p2,p3,p4]))
        controls.append(dict(polynomial=str(P),degree=int(m),generic_class_group_rank=int(m-1)))

    B = T**7*z**4+(T+1)*z**3+T**5*z+2*T*T+1
    B0 = B.subs(T,0)
    H = S.cancel((B-B0)/T)
    check('the T-dependent lower-term correction is globally polynomial',S.denom(H)==1)
    zero('the relative Y shear removes every T-dependent coefficient literally',
        (T*(Y+H)+x*y-B)-(T*Y+x*y-B0))
    zero('the relative shear has the opposite translation as polynomial inverse',Y+H-H-Y)
    relative_map = S.Matrix([T,Y-H,x,y,z])
    check('its full ambient Jacobian is one',relative_map.jacobian([T,Y,x,y,z]).det()==1)
    image_B = {T:s,Y:w,x:1,y:t-s*w+B.subs({T:s,z:v}),z:v}
    zero('the arbitrary two-variable suspension carries the desired output pair',
        S.Matrix([T,T*Y+x*y-B]).subs(image_B,simultaneous=True)-S.Matrix([s,t]))
    zero('that suspension embedding has literal polynomial source recovery',
        S.Matrix([T,T*Y+x*y-B,Y,z]).subs(image_B,simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('its generic ambient fibre eliminates Y exactly',
        (T*Y+x*y-B-V).subs(Y,(V-x*y+B)/T))
    check('the special T=0 relation has a nonzero right-hand polynomial over C(V)',
        S.Poly(B0+V,z).LC()==1 and S.degree(B0+V,z)==3)
    constant_B = S.Integer(7)
    zero('even a constant B has a literal integral-special-fibre Laurent parameterization',
        (x*y-constant_B-V).subs(y,(constant_B+V)/x))

    a,b = S.symbols('a b')
    non_keller = S.Matrix([a,b*b-a])
    check('the two-parameter countercontrol is explicitly outside Keller',
        non_keller.jacobian([a,b]).det()==2*b)
    zero('its two-parameter cover nevertheless acquires a source square root',
        (T+V).subs({T:a,V:b*b-a})-b*b)
    zero('the excluded constant first-type polynomial leaves nonconstant Laurent units',
        (x*y/(V+7)-1).subs(y,(V+7)/x))

    retained = {
        'jacobian-followup/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md':'3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
        'jacobian-followup/hc4-triangular-graphs/SMOOTH_INFINITY_OUTPUT_OBSTRUCTION.md':'7ee0a9b293ca95a8f037cfdef0b5c0630f0fb169ec728df5b18c4dc2fb5638db',
        'jacobian-followup/hc4-relative-row-orbit/PROOF.md':'d3f528180c732f2a9c65af0cba4a405a37e34437d48e350e78748eaff61ecb36',
        'jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md':'26a6a76eb579de05b8e71245b8b804bdd5f036602f2e3e6ad5d036fd9f5cd572',
        'jacobian-followup/hc4-triangular-graphs/sources/arzhantsev-petravchuk-closed-math-0608157v2.html':'8042645504f3535bae9ce6d656a68015b893dbf259888ca4708e69722d658cc9'
    }
    for name,digest in retained.items():
        check('retained mathematical or primary import pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_regular_parameter.py','REGULAR_PARAMETER_MONODROMY.md','sources/regular-parameter-import.json']},
        dependencies=retained,polynomial_output_controls=controls,
        written_not_finite_exhaustion=[
            'nowhere-critical polynomial parameter regularity in every dimension, with an explicit classical closed-polynomial import',
            'arbitrary finite one-parameter monodromy is preserved by the Keller source field',
            'the torsion-free class-group output criterion does not require abelian action when it descends from that parameter',
            'all nonconstant P(z) and every two-variable suspension T*Y+xy-B(T,z) are covered in all degrees'],
        closed_polynomial_import_fully_reproved=False,
        every_special_parameter_fibre_asserted_irreducible=False,
        entire_Keller_field_extension_asserted_trivial=False,
        arbitrary_two_parameter_class_actions_asserted_covered=False,
        source_or_target_submersion_asserted_coordinate=False,
        generic_Danielewski_rings_asserted_factorial=False,
        arbitrary_low_rank_or_constant_first_type_asserted_covered=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-regular-parameter.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
