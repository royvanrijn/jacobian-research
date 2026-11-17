"""Finite exact controls for the written nilpotent-cover and output theorems."""
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


def finite_group(elements,multiply,identity):
    elements = frozenset(elements)

    def closure(generators):
        result = {identity}
        pending = [identity]
        while pending:
            current = pending.pop()
            for generator in generators:
                value = multiply(current,generator)
                if value not in result:
                    result.add(value)
                    pending.append(value)
        return frozenset(result)

    subgroups = {frozenset({identity})}
    pending = list(subgroups)
    while pending:
        subgroup = pending.pop()
        for element in elements-subgroup:
            larger = closure(subgroup|{element})
            if larger not in subgroups:
                subgroups.add(larger)
                pending.append(larger)
    inverses = {a:next(b for b in elements if multiply(a,b)==identity) for a in elements}

    def normal(subgroup):
        return all(multiply(multiply(g,h),inverses[g]) in subgroup
            for g in elements for h in subgroup)

    prime_index_normal = [h for h in subgroups if h!=elements
        and S.isprime(len(elements)//len(h)) and normal(h)]
    return elements,subgroups,prime_index_normal,closure


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

    compose = lambda a,b:tuple(a[b[i]] for i in range(len(a)))
    identity = tuple(range(4))
    rotation = (1,2,3,0)
    reflection = (0,3,2,1)
    rotations = [identity]
    for _ in range(3):
        rotations.append(compose(rotation,rotations[-1]))
    dihedral = set(rotations)|{compose(g,reflection) for g in rotations}
    group,subgroups,normal_maximal,closure = finite_group(dihedral,compose,identity)
    check('the quartic root permutation group has order eight',len(group)==8)
    check('the quartic group is nonabelian',compose(rotation,reflection)!=compose(reflection,rotation))
    check('all ten D8 subgroups were independently enumerated',len(subgroups)==10)
    check('each proper D8 source image lies below a normal prime-index quotient',
        all(any(h<=n for n in normal_maximal) for h in subgroups if h!=group))
    check('the reflection image is proper and intransitive on the four roots',
        len(closure({reflection}))==2 and {g[0] for g in closure({reflection})}=={0})

    product = {(g,c) for g in group for c in range(3)}
    multiply_product = lambda a,b:(compose(a[0],b[0]),(a[1]+b[1])%3)
    mixed,subgroups_mixed,normal_mixed,_ = finite_group(product,multiply_product,(identity,0))
    check('the nilpotent control need not be a prime-power group',len(mixed)==24)
    check('all twenty subgroups of D8 times C3 were enumerated',len(subgroups_mixed)==20)
    check('every proper mixed-prime nilpotent subgroup still has a cyclic prime quotient above it',
        all(any(h<=n for n in normal_mixed) for h in subgroups_mixed if h!=mixed))

    permutations = set(itertools.permutations(range(3)))
    symmetric,subgroups_s3,normal_s3,closure_s3 = finite_group(permutations,compose,(0,1,2))
    transposition = closure_s3({(1,0,2)})
    check('the solvable boundary has the six S3 subgroups',len(subgroups_s3)==6)
    check('the S3 transposition subgroup escapes every normal prime-index container',
        len(transposition)==2 and not any(transposition<=n for n in normal_s3))
    check('the S3 boundary nevertheless surjects onto its only nontrivial abelian quotient',
        len(normal_s3)==1 and len({compose(h,a) for h in transposition for a in normal_s3[0]})==6)

    def matrix(permutation):
        result = S.zeros(4)
        for i,target in enumerate(permutation):
            result[target,i] = 1
        return result

    equations = (matrix(rotation)-S.eye(4)).col_join(matrix(reflection)-S.eye(4))
    check('the full nonabelian root action has only the diagonal fixed vector',
        equations.rank()==3 and equations*S.ones(4,1)==S.zeros(8,1))
    check('an excluded reflection source image would leave nonzero quotient-lattice invariants',
        4-(matrix(reflection)-S.eye(4)).rank()-1==2)

    T,V,z,a,b,c,x,y,v,s,t,w = S.symbols('T V z a b c x y v s t w')
    P = T*(z*z-T)**2
    polynomial = P+V
    check('the new quartic has nonconstant leading coefficient',S.Poly(P,z).LC()==T)
    check('the quartic is not a translated fourth power by its cubic and quadratic coefficients',
        S.Poly(P,z).nth(3)==0 and S.Poly(P,z).nth(2)==-2*T*T)
    zero('the first square root gives the rational recovery V=-T*b^2',
        (b*b+V/T).subs(V,-T*b*b))
    check('the first quadratic radicand has odd V valuation',
        S.degree(-V/T,V)==1 and not S.denom(-V/T).has(V))
    check('the next two radicands have distinct simple prime divisors',
        S.Poly(T+b,T).degree()==1 and S.Poly(T-b,T).degree()==1
        and S.gcd(T+b,T-b)==1)
    zero('one prime sees the first square class with multiplicity one',S.diff(T+b,T)-1)
    check('the other square class is nonzero at that prime',(T-b).subs(T,-b)==-2*b)
    check('the roles of the square-class valuation witnesses are symmetric',(T+b).subs(T,b)==2*b)
    field = S.QQ.frac_field(T,V)
    gb = S.groebner([a*a-T-b,c*c-T+b,b*b+V/T],a,c,b,domain=field)
    normal = lambda expression:gb.reduce(S.expand(expression))[1]
    check('the radical quotient has the eight expected standard monomials',
        {item.LM().exponents for item in gb.polys}=={(2,0,0),(0,2,0),(0,0,2)})
    zero('the first signed pair is a root of the actual nonmonic quartic',normal(polynomial.subs(z,a)))
    zero('the second signed pair is also a root',normal(polynomial.subs(z,c)))
    tau = {a:c,c:-a,b:-b}
    sigma = {a:-a,c:c,b:b}
    relations = [a*a-T-b,c*c-T+b,b*b+V/T]
    check('the order-four rotation preserves all radical field relations',
        all(normal(q.subs(tau,simultaneous=True))==0 for q in relations))
    check('the reflection also preserves all radical field relations',
        all(normal(q.subs(sigma,simultaneous=True))==0 for q in relations))
    twice = S.Matrix([a,c,b]).subs(tau,simultaneous=True).subs(tau,simultaneous=True)
    zero('the quartic field rotation has square equal to simultaneous root negation',
        twice-S.Matrix([-a,-c,b]))
    zero('the field rotation has fourth power one',twice.subs(tau,simultaneous=True).subs(tau,simultaneous=True)-S.Matrix([a,c,b]))

    q1 = T*z*z+z+T*T
    q2 = (T-1)*z*z+T*z+1
    q3 = (T+1)*z*z+z+T**3
    nested = q3.subs(z,q2.subs(z,q1))
    check('arbitrary quadratic stages need not have polynomial completed-square shifts',
        S.denom(S.cancel(1/(2*T)))==2*T)
    check('three arbitrary quadratic stages produce degree eight with the exact leading coefficient',
        S.degree(nested,z)==8 and S.factor(S.Poly(nested,z).LC())==(T-1)**2*(T+1)*T**4)
    relation = x*y-polynomial
    ambient = S.groebner([relation],x,y,z,v,domain=S.QQ.poly_ring(T,V),order='lex')
    check('the nonabelian nonmonic output retains a globally monic xy basis',
        ambient.polys[0].LM().exponents==(1,1,0,0) and ambient.polys[0].LC()==1)
    zero('the leading-zero affine height-one fibre stays integral through xy=V',relation.subs(T,0)-(x*y-V))
    check('the special closed origin fibre is allowed to be reducible',S.factor(relation.subs({T:0,V:0}))==x*y)
    image = {T:s,x:1,y:t+P.subs({T:s,z:w}),z:w,v:v}
    zero('the explicit HC4-size embedding has its complete polynomial inverse',
        S.Matrix([T,y-P,z,v]).subs(image,simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the two output observables carry exactly the identity Keller pair',
        S.Matrix([T,x*y-P]).subs(image,simultaneous=True)-S.Matrix([s,t]))
    actual = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('the independent nondegenerate ambient quadratic pulls back literally',
        (T*v+y*z+x*x).subs(image,simultaneous=True)-actual)
    check('the actual nonabelian nonmonic HC4 control has unit Hessian',S.hessian(actual,[s,t,w,v]).det()==1)
    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    evaluate = {T:p4,z:p2}
    inverse = S.Matrix([p4,p3-P.subs(evaluate)-p2*S.diff(P,z).subs(evaluate),
        p2,p1-p2*S.diff(P,T).subs(evaluate)])
    gradient = S.Matrix([S.diff(actual,q) for q in [s,t,w,v]])
    zero('the actual gradient inverse is a left inverse',
        inverse.subs(dict(zip([p1,p2,p3,p4],gradient)),simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the actual gradient inverse is a right inverse',
        gradient.subs(dict(zip([s,t,w,v],inverse)),simultaneous=True)-S.Matrix([p1,p2,p3,p4]))
    cubic = z**3-3*z+V
    zero('the mixed-prime solvable boundary has nonsquare cubic discriminant',
        S.discriminant(cubic,z)-27*(4-V*V))
    u = S.symbols('u')
    zero('the mixed-prime boundary has the written Cardano radical recovery',
        cubic.subs(z,u+1/u)-(u**3+u**-3+V))

    retained = {
        'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md':'3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
        'hessian/hc4-triangular-graphs/LEADING_FIBRE_UNIT_GATE.md':'1d054ba8762ea657be879b6e3568aa58c9ccb65e14ab642d45e0865e8f6382f0',
        'hessian/hc4-triangular-graphs/NORMAL_MONOGENIC_COVERS.md':'223b81738c823c3658079344f57e74a1405b3fcb52eb8adc76c9ee170226226f',
        'hessian/hc4-relative-row-orbit/PROOF.md':'d3f528180c732f2a9c65af0cba4a405a37e34437d48e350e78748eaff61ecb36'
    }
    context = json.loads((HERE/'sources/nilpotent-cover-context.json').read_text())
    for item in context['sources']:
        retained['hessian/hc4-triangular-graphs/sources/'+item['retained_file']] = item['sha256']
    for name,digest in retained.items():
        check('retained prior proof or primary context byte pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_nilpotent_covers.py','NILPOTENT_COVER_MONODROMY.md','sources/nilpotent-cover-context.json']},
        dependencies=retained,
        written_not_finite_exhaustion=[
            'every proper finite nilpotent subgroup lies below a normal cyclic prime quotient',
            'Keller base change preserves the full group of every finite nilpotent Galois target cover in every dimension',
            'the filtered class-group output criterion closes for finite nilpotent class monodromy',
            'every positive-z-degree polynomial output with nilpotent generic root monodromy closes without a special-fibre assumption',
            'arbitrary common-prime radical depth, including every composition of quadratics, has p-group root monodromy'],
        full_nilpotent_Galois_group_preservation_proved=True,
        arbitrary_finite_Galois_group_preservation_asserted=False,
        arbitrary_finite_solvable_group_preservation_asserted=False,
        mixed_prime_radical_towers_asserted_nilpotent=False,
        all_possible_wreath_groups_asserted_attained=False,
        all_Keller_component_fibres_asserted_irreducible=False,
        all_closed_target_fibres_asserted_geometrically_integral=False,
        generic_regular_parameter_asserted_localized_unit_descent=False,
        general_primary_wreath_embedding_theorem_imported=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-nilpotent-covers.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
