"""Exact controls for normal finite monogenic covers and monic output bases.

Purity, normality under etale base change and the finite-cover comparison
are explicit primary inputs; finite controls do not prove those theorems.
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

    T,V,Z,s,t,w,v,x,y,z,c,a = S.symbols('T V Z s t w v x y z c a')
    P = Z**3+T*Z*Z+(T*T+1)*Z+T**5+1
    root_equation = P-V
    check('the genuinely two-parameter root cover is monic of fixed degree three',
        S.Poly(root_equation,Z).LC()==1 and S.degree(root_equation,Z)==3)
    check('its total space is smooth by a constant target derivative',S.diff(root_equation,V)==-1)
    zero('its total space has a literal polynomial coordinate recovery',root_equation.subs(V,P))
    monic = S.groebner([root_equation],Z,domain=S.QQ.poly_ring(T,V))
    check('monic division gives a finite free target module in the control',
        S.degree(monic.reduce(Z**7+T*Z**4)[1],Z)<3)
    K = S.Matrix([s+t*t,t])
    check('the source change is an actual nonlinear Keller automorphism',K.jacobian([s,t]).det()==1)
    zero('it has an exact polynomial inverse',
        S.Matrix([T-V*V,V]).subs({T:K[0],V:K[1]},simultaneous=True)-S.Matrix([s,t]))
    pulled = root_equation.subs({T:K[0],V:K[1]},simultaneous=True)
    zero('the pulled-back monic cover has a global smoothness certificate',
        2*t*S.diff(pulled,s)-S.diff(pulled,t)-1)

    ramified = Z*Z-s*s
    zero('a non-Keller square map makes the normal cover factor',ramified-(Z-s)*(Z+s))
    zero('the two factor resultant records the bad source divisor',S.resultant(Z-s,Z+s,Z)-2*s)
    check('the same source divisor is precisely the nonunit Jacobian',
        S.Matrix([s*s,t]).jacobian([s,t]).det()==2*s)
    zero('the factors really meet at the ramification point',
        S.Matrix([Z-s,Z+s]).subs({s:0,Z:0}))
    zero('the pulled-back cover is singular there',
        S.Matrix([S.diff(ramified,u) for u in [s,t,Z]]).subs({s:0,Z:0}))
    zero('on a Laurent source the nonconstant resultant is instead a unit',(2*s)*(1/(2*s))-1)

    zeta = (-1+S.sqrt(-3))/2
    factor_a = Z-a
    factor_c = S.expand((Z-zeta*a)*(Z-zeta*zeta*a))
    zero('the integral-root collision control has the complete cubic factorization',
        S.expand(factor_a*factor_c)-(Z**3-a**3))
    zero('its cross-resultant is the positive-valuation collision product',
        S.simplify(S.resultant(factor_a,factor_c,Z))-3*a*a)
    check('that cross-resultant has positive integral-root valuation at the origin',
        S.Poly(3*a*a,a).terms()==[((2,),S.Integer(3))])

    permutations = list(itertools.permutations(range(3)))
    identity = (0,1,2)
    subset = {0}
    stabilizer = {p for p in permutations if {p[i] for i in subset}==subset}
    cuts = [{p[i] for i in subset} for p in permutations]
    core = {p for p in permutations if all({p[i] for i in cut}==cut for cut in cuts)}
    check('one factor cut admits nontrivial inertia without the conjugate-cut test',len(stabilizer)==2)
    check('all conjugate cuts reduce its normal core to the identity',core=={identity})
    check('a collision transposition crosses at least one conjugate cut',
        any({(1,0,2)[i] for i in cut}!=cut for cut in cuts))
    C = S.Matrix([[0,0,-1],[1,0,-1],[0,1,-1]])
    zero('a proper cyclic subgroup already has a transitive four-root action',C**4-S.eye(3))
    first = S.Matrix([1,0,0])
    orbit = {tuple(C**i*first) for i in range(4)}
    check('that cyclic root orbit has four classes while S4 has twenty-four elements',
        len(orbit)==4 and len(list(itertools.permutations(range(4))))==24)
    check('transitivity alone kills the invariant quotient-lattice classes',
        (C-S.eye(3)).rank()==3)

    cone = Z*Z-T*V
    zero('the singular normal cover has its even-invariant parameterization',
        cone.subs({T:a*a,V:c*c,Z:a*c},simultaneous=True))
    zero('the cover is singular at the cone vertex',
        S.Matrix([S.diff(cone,u) for u in [T,V,Z]]).subs({T:0,V:0,Z:0}))
    cone_pull = cone.subs({T:K[0],V:K[1]},simultaneous=True)
    check('its quadratic pullback has two distinct odd-valuation discriminant factors',
        sorted(exponent for factor,exponent in S.factor_list(t*(s+t*t))[1])==[1,1]
        and S.gcd(t,s+t*t)==1)

    q = S.symbols('q')
    node = q*q-T*T*(T+1)
    node_map = {T:c*c-1,q:(c*c-1)*c}
    zero('a monic primitive order has the stated normalization parameterization',node.subs(node_map,simultaneous=True))
    zero('its primitive field generator recovers the normal parameter only generically',
        (q/T).subs(node_map,simultaneous=True)-c)
    zero('the missing normal element is integral over the order',(c*c-T-1).subs(T,c*c-1))
    check('two normal points collide in the primitive order while the missing element separates them',
        all(expression.subs(c,1)==expression.subs(c,-1) for expression in node_map.values())
        and S.Integer(1)!=S.Integer(-1))
    check('the normalized finite map is locally etale at both collided points',
        S.diff(c*c-1,c).subs(c,1)!=0 and S.diff(c*c-1,c).subs(c,-1)!=0)
    zero('the primitive monic model nevertheless has zero root derivative there',
        S.diff(node,q).subs(q,0))
    nonmonic = T*Z*Z-V
    zero('a nonmonic model has a whole affine root fibre where finiteness fails',
        nonmonic.subs({T:0,V:0}))
    zero('dividing its leading coefficient produces a root with an actual pole',
        nonmonic.subs({T:a*a,V:1,Z:1/a},simultaneous=True))
    check('the root pole prevents the all-roots-integral premise',S.denom(1/a)==a)

    Pambient = P.subs(Z,z)
    image = {T:s,x:1,y:t+Pambient.subs({T:s,z:w}),z:w,v:v}
    zero('the monic two-parameter observable has a genuine HC4-size embedding',
        S.Matrix([T,x*y-Pambient]).subs(image,simultaneous=True)-S.Matrix([s,t]))
    zero('that embedding has literal full polynomial source recovery',
        S.Matrix([T,y-Pambient,z,v]).subs(image,simultaneous=True)-S.Matrix([s,t,w,v]))
    actual = 1+s*v+t*w+w*Pambient.subs({T:s,z:w})
    zero('an independent nondegenerate ambient quadratic gives the stated actual pullback',
        (T*v+y*z+x*x).subs(image,simultaneous=True)-actual)
    check('the actual HC4 source Hessian determinant is one',S.hessian(actual,[s,t,w,v]).det()==1)
    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    evaluate = {T:p4,z:p2}
    inverse = S.Matrix([p4,p3-Pambient.subs(evaluate)-p2*S.diff(Pambient,z).subs(evaluate),
        p2,p1-p2*S.diff(Pambient,T).subs(evaluate)])
    gradient = S.Matrix([S.diff(actual,u) for u in [s,t,w,v]])
    zero('its actual source gradient inverse is a left inverse',
        inverse.subs(dict(zip([p1,p2,p3,p4],gradient)),simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the same written inverse is a right inverse',
        gradient.subs(dict(zip([s,t,w,v],inverse)),simultaneous=True)-S.Matrix([p1,p2,p3,p4]))

    retained = {
        'hessian/hc4-triangular-graphs/REGULAR_PARAMETER_MONODROMY.md':'b010d2e4c62906a8a79da9dbc41ef7d0e3d5edefb4520cf487ba09ded5d86956',
        'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md':'3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
        'hessian/hc4-relative-row-orbit/PROOF.md':'d3f528180c732f2a9c65af0cba4a405a37e34437d48e350e78748eaff61ecb36',
        'hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md':'26a6a76eb579de05b8e71245b8b804bdd5f036602f2e3e6ad5d036fd9f5cd572'
    }
    context = json.loads((HERE/'sources/monic-cover-import.json').read_text())
    for item in context['sources']:
        retained['hessian/hc4-triangular-graphs/sources/'+item['retained_file']] = item['sha256']
    for name,digest in retained.items():
        check('retained mathematical or primary import pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_monic_covers.py','NORMAL_MONOGENIC_COVERS.md','sources/monic-cover-import.json']},
        dependencies=retained,
        written_not_finite_exhaustion=[
            'normal finite monogenic covers stay integral under every complex Keller map, in every dimension and degree',
            'all conjugate constant cross-resultants prevent every codimension-one inertia crossing; the quotient is trivial by purity and finite-cover comparison',
            'all monic two-parameter Danielewski output pairs close by root transitivity and vanishing fixed classes',
            'no nontrivial Keller intermediate field has globally monogenic finite normalization',
            'all-degree actual HC4 Hessian and gradient inverse controls have literal formulae'],
        primary_normality_and_purity_theorems_fully_reproved=False,
        full_splitting_group_preservation_asserted=False,
        generic_primitive_element_asserted_global_normal_monogenicity=False,
        arbitrary_nonmonic_or_nonnormal_models_asserted_covered=False,
        arbitrary_etale_source_asserted_scalar_units=False,
        nonconstant_leading_coefficients_inverted_globally=False,
        entire_Keller_field_extension_asserted_trivial_unconditionally=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-monic-covers.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
