"""Exact controls for the leading-fibre localized-unit gate.

The written all-degree proofs use explicit prior purity inputs and a new
Riemann-Hurwitz input. The audited 2024 article is not a proof input.
"""
import argparse
import hashlib
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

    s,t,w,v,T,V,x,y,z,u,c = S.symbols('s t w v T V x y z u c')
    f = s+s*s*t
    zero('a reducible special-fibre control is nevertheless a polynomial submersion',
        (1-2*s*t)*S.diff(f,s)+4*t*t*S.diff(f,t)-1)
    zero('the special fibre has its two exact reduced factors',f-s*(1+s*t))
    zero('those factors have disjoint zero sets by an ideal certificate',(1+s*t)-t*s-1)
    zero('the first component factor becomes a localized source unit',s*(1+s*t)/f-1)
    zero('the other component factor is also a localized unit',(1+s*t)*s/f-1)
    check('the adjacent fibre is primitive linear in t and irreducible',
        S.degree(f-1,t)==1 and S.gcd(s*s,s-1)==1)
    lattice = S.Matrix([[1,1,0],[0,0,1]])
    check('the two fibre relations leave a rank-one unit defect',lattice.rank()==2)
    check('the fibre relations form a saturated lattice',abs(lattice[:,[0,2]].det())==1)
    check('one component unit represents a nonzero quotient class',
        S.Matrix.vstack(lattice,S.Matrix([[1,0,0]])).rank()==3)
    zero('the generic regular parameter still has literal rational recovery',
        f.subs(t,(c-s)/(s*s))-c)
    check('that localized component unit is not algebraic over C(f)',
        S.Matrix([f,s]).jacobian([s,t]).det()==-s*s)

    K = S.Matrix([s+t*t,t])
    check('the good-fibre source control is genuinely Keller',K.jacobian([s,t]).det()==1)
    for level in [0,1]:
        check(f'the leading-zero fibre {level} is irreducible by its linear source variable',
            S.degree(K[0]-level,s)==1 and S.Poly(K[0]-level,s).LC()==1)
    resultant_unit = (K[0]-1)**2/K[0]**3
    zero('a localized unit in the good-fibre case descends literally to T',
        resultant_unit-((T-1)**2/T**3).subs(T,K[0]))
    check('the descended unit has zeros and poles only on the removed vertical levels',
        S.factor(S.denom((T-1)**2/T**3))==T**3)

    check('the affine punctured-line identity is etale',S.diff(u,u)==1)
    zero('the omitted finite point maps to zero rather than infinity',u.subs(u,0))
    check('a finite power cover instead records its forbidden finite branch point',
        S.diff(u**3,u).subs(u,0)==0)
    beta = S.symbols('beta')
    constant_cover = beta*beta-T
    zero('a nontrivial constant-field cover has no V-dependent branch polynomial',
        S.diff(S.discriminant(constant_cover,beta),V))
    check('it is etale over the generic V-line with nonzero constant discriminant',
        S.discriminant(constant_cover,beta)==4*T)
    check('absorbing that constant square root in a polynomial source fails Keller',
        S.Matrix([s*s,t]).jacobian([s,t]).det()==2*s)

    leading = T*(T-1)
    P = leading*(z**3+(T*T+1)*z*z+T**3*z)+(T*T+1)
    relation = x*y-P-V
    gb = S.groebner([relation],x,y,z,v,domain=S.QQ.poly_ring(T,V),order='lex')
    check('the nonmonic output still has a globally monic xy reduction',
        gb.polys[0].LM().exponents==(1,1,0,0) and gb.polys[0].LC()==1)
    weighted_top = sum(coefficient*x**monomial[0]*y**monomial[1]*z**monomial[2]
        for monomial,coefficient in S.Poly(relation,x,y,z).terms()
        if 3*monomial[0]+3*monomial[1]+2*monomial[2]==6)
    zero('its generic weighted leading relation has the actual nonconstant leading coefficient',
        weighted_top-(x*y-leading*z**3))
    zero('at a leading zero the affine height-one fibre remains a nonzero Laurent equation',
        relation.subs(T,0)-(x*y-1-V))
    zero('that affine height-one fibre has an explicit nonconstant unit',
        (x*y/(1+V)-1).subs(y,(1+V)/x))
    closed = S.groebner([relation.subs({T:0,V:-1})],x,y,z,order='lex')
    check('a closed codimension-two target fibre can still be reducible',
        closed.reduce(x*y)[1]==0 and closed.reduce(x)[1]!=0 and closed.reduce(y)[1]!=0)
    image = {T:s,x:1,y:t+P.subs({T:s,z:w}),z:w,v:v}
    zero('the degree-drop family has a genuine complete HC4-size source inverse',
        S.Matrix([T,y-P,z,v]).subs(image,simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('its two stipulated output observables are exactly the Keller base',
        S.Matrix([T,x*y-P]).subs(image,simultaneous=True)-S.Matrix([s,t]))
    actual = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('an independent ambient nondegenerate quadratic has the stated pullback',
        (T*v+y*z+x*x).subs(image,simultaneous=True)-actual)
    check('the actual HC4 source Hessian is a unit despite leading-degree drops',
        S.hessian(actual,[s,t,w,v]).det()==1)
    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    evaluate = {T:p4,z:p2}
    inverse = S.Matrix([p4,p3-P.subs(evaluate)-p2*S.diff(P,z).subs(evaluate),
        p2,p1-p2*S.diff(P,T).subs(evaluate)])
    gradient = S.Matrix([S.diff(actual,q) for q in [s,t,w,v]])
    zero('the displayed actual gradient inverse is a left inverse',
        inverse.subs(dict(zip([p1,p2,p3,p4],gradient)),simultaneous=True)-S.Matrix([s,t,w,v]))
    zero('the displayed actual gradient inverse is a right inverse',
        gradient.subs(dict(zip([s,t,w,v],inverse)),simultaneous=True)-S.Matrix([p1,p2,p3,p4]))

    quadratic = T*z*z+(T+1)*z+T**3+V
    discriminant = S.discriminant(quadratic,z)
    zero('arbitrary-coefficient quadratic root monodromy has its exact discriminant',
        discriminant-((T+1)**2-4*T*(T**3+V)))
    check('the quadratic discriminant has a single simple V-zero over C(T)',
        S.degree(discriminant,V)==1 and S.diff(discriminant,V)==-4*T)
    power = 5
    shift = T*T
    constant = T**3+1
    cyclic = leading*(z-shift)**power+constant+V
    zero('the arbitrary-leading translated power normalizes to one cyclic equation',
        (cyclic/leading).subs(z,shift+beta)-(beta**power+(V+constant)/leading))
    check('its Kummer parameter has target valuation one at V plus the constant term',
        S.degree(V+constant,V)==1 and S.diff(V+constant,V)==1
        and not S.denom(1/leading).has(V))

    retained = {
        'hessian/hc4-triangular-graphs/NORMAL_MONOGENIC_COVERS.md':'223b81738c823c3658079344f57e74a1405b3fcb52eb8adc76c9ee170226226f',
        'hessian/hc4-triangular-graphs/REGULAR_PARAMETER_MONODROMY.md':'b010d2e4c62906a8a79da9dbc41ef7d0e3d5edefb4520cf487ba09ded5d86956',
        'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md':'3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
        'hessian/hc4-relative-row-orbit/PROOF.md':'d3f528180c732f2a9c65af0cba4a405a37e34437d48e350e78748eaff61ecb36'
    }
    context = json.loads((HERE/'sources/leading-fibre-import.json').read_text())
    for item in context['sources']:
        retained['hessian/hc4-triangular-graphs/sources/'+item['retained_file']] = item['sha256']
    for name,digest in retained.items():
        check('retained proof, primary input or explicitly excluded audit-source pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_leading_fibre.py','LEADING_FIBRE_UNIT_GATE.md','sources/leading-fibre-import.json']},
        dependencies=retained,
        written_not_finite_exhaustion=[
            'the localized-unit defect is exactly the direct sum of reduced special-fibre component lattices',
            'normal monogenic covers over a vertical open target remain irreducible when the inverted source fibres are irreducible',
            'the remaining unramified generic affine-line quotient is a finite constant-field extension by Riemann-Hurwitz',
            'all nonmonic polynomial output equations close under the finite fibre condition or abelian root monodromy',
            'every quadratic-z case and every arbitrary-leading translated power closes without a special-fibre hypothesis'],
        Riemann_Hurwitz_theorem_fully_reproved=False,
        Tian_Sun_2024_claims_imported=False,
        all_Keller_component_fibres_asserted_irreducible=False,
        generic_regular_parameter_asserted_localized_unit_descent=False,
        all_closed_target_fibres_asserted_geometrically_integral=False,
        arbitrary_etale_affine_curve_maps_asserted_finite=False,
        generic_affine_line_geometric_triviality_asserted_no_constant_field_extension=False,
        full_splitting_group_preservation_asserted=False,
        arbitrary_nonmonic_outputs_without_fibre_or_monodromy_gate_asserted_closed=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-leading-fibre.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
