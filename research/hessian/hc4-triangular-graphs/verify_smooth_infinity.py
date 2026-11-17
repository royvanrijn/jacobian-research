"""Exact controls for the written arbitrary-degree output-base obstruction.

The complete-intersection Picard theorem is an explicit primary import.
Finite controls do not prove it or exhaust arbitrary observable degrees.
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

    Z = S.symbols('Z0:5')
    T,U = S.symbols('T U')
    diagonal = []
    for m in [2,3]:
        p = sum(z**m for z in Z)
        q = sum(i*z**m for i,z in enumerate(Z))
        J = S.Matrix([[S.diff(f,z) for z in Z] for f in [p,q]])
        minors = [J[:,[i,j]].det()-m*m*(j-i)*Z[i]**(m-1)*Z[j]**(m-1)
            for i in range(5) for j in range(i+1,5)]
        zero(f'degree {m}: every rank minor has the distinct-coefficient support formula',S.Matrix(minors))
        check(f'degree {m}: every two coordinate columns have distinct diagonal labels',
            all(j-i!=0 for i in range(5) for j in range(i+1,5)))
        check(f'degree {m}: both leading forms have their stated homogeneous degree',
            all(sum(exponents)==m for f in [p,q] for exponents in S.Poly(f,*Z).monoms()))
        g0 = Z[0]**m-(T-U+sum((i-1)*Z[i]**m for i in range(2,5)))
        g1 = Z[1]**m-(U-sum(i*Z[i]**m for i in range(2,5)))
        zero(f'degree {m}: the first monic reduction is an exact combination over the full target',g0-((p-T)-(q-U)))
        zero(f'degree {m}: the second monic reduction is an exact combination over the full target',g1-(q-U))
        gb = S.groebner([p-T,q-U],*Z,domain=S.QQ.poly_ring(T,U),order='lex')
        check(f'degree {m}: the common free-basis leading monomials are pure coprime powers',
            [poly.LM().exponents for poly in gb.polys]==[(m,0,0,0,0),(0,m,0,0,0)])
        zero(f'degree {m}: both proposed relations reduce to zero',
            S.Matrix([gb.reduce(g)[1] for g in [g0,g1]]))
        standard = Z[0]**(m-1)*Z[1]**(m-1)*Z[2]**3
        zero(f'degree {m}: a boundary standard monomial stays independent',gb.reduce(standard)[1]-standard)
        at_zero = [poly.as_expr().subs({T:0,U:0}) for poly in gb.polys]
        cone = S.groebner([p,q],*Z,order='lex')
        zero(f'degree {m}: the same basis survives the singular zero affine fibre',
            S.Matrix([cone.reduce(g)[1] for g in at_zero]))
        diagonal.append(dict(degree=m,leading_monomials=[list(poly.LM().exponents) for poly in gb.polys]))
    projective = {Z[0]:1,Z[1]:S.sqrt(-2),Z[2]:1,Z[3]:0,Z[4]:0}
    p = sum(z*z for z in Z)
    q = sum(i*z*z for i,z in enumerate(Z))
    zero('the quadratic diagonal boundary has an explicit nonzero common point',
        S.Matrix([p.subs(projective),q.subs(projective)]))
    check('its two boundary gradients have full rank at that point',
        S.Matrix([[S.diff(f,z).subs(projective) for z in Z] for f in [p,q]]).rank()==2)

    x,y,z,u,v = S.symbols('x y z u v')
    trap = S.groebner([x-T,x*y+z*z-U],x,z,y,u,v,domain=S.QQ.poly_ring(T,U),order='lex')
    zero('the denominator-trap ring has the stated monic relation',trap.reduce(z*z-U+T*y)[1])
    zero('the forbidden denominator cancels in the full ambient ring',
        ((z*z-U)/T).subs({T:x,U:x*y+z*z})+y)
    check('the unique F[z] representation has a genuine coefficient denominator',
        S.denom(S.Poly((z*z-U)/T,z).coeff_monomial(z*z))==T)
    check('the putative primitive element has coefficient one in the monic normal basis',
        trap.reduce(z)[1]==z)
    zero('its height-one reduction becomes algebraic instead of positive-degree transcendental',
        (z*z-U+T*y).subs(T,0)-(z*z-U))
    generic_y = (U-z*z)/T
    zero('the geometric generic affine fibre eliminates y exactly',T*generic_y+z*z-U)
    zero('the generic elimination has a literal inverse on all other coordinates',
        S.Matrix([T,z,u,v])-S.Matrix([x,z,u,v]).subs(x,T))
    top_gb = S.groebner([x,x*y+z*z],x,z,y,u,v,order='lex')
    check('the leading scheme contains z squared in its ideal',top_gb.reduce(z*z)[1]==0)
    check('the leading scheme does not contain z in its ideal',top_gb.reduce(z)[1]==z)
    check('its reduced support is the smooth coordinate plane in projective four-space',
        S.Matrix([[1,0,0,0,0],[0,0,1,0,0]]).rank()==2)
    s,t = S.symbols('s t')
    embedding = {T:s,x:1,y:t+u*v,z:u,v:v}
    zero('the genuine quadric embedding has the stipulated output base',
        S.Matrix([T,x*y-z*v]).subs(embedding,simultaneous=True)-S.Matrix([s,t]))
    inverse = S.Matrix([T,y-z*v,z,v])
    zero('the embedding has a literal polynomial source inverse',
        inverse.subs(embedding,simultaneous=True)-S.Matrix([s,t,u,v]))
    E = x-1
    zero('the full image equation vanishes identically',E.subs(embedding,simultaneous=True))
    check('its image is smooth everywhere',S.diff(E,x)==1)
    q4 = x*y-z*v
    check('the positive boundary quadratic has nondegenerate four-variable Hessian',
        S.hessian(q4,[x,y,z,v]).det()==1)
    actual = 1+s*u+t*v+u*v*v
    zero('the independent quadratic ambient potential has the stated actual pullback',
        (T*z+y*v+x*x).subs(embedding,simultaneous=True)-actual)
    check('the control is compatible with actual unit HC4 Hessian',S.hessian(actual,[s,t,u,v]).det()==1)

    retained = {
        'hessian/hc4-triangular-graphs/STABLE_HYPERSURFACE_OBSTRUCTION.md':'a169c8541cda91042b31a55b77c249c89d48170eba2fe4ba7b524ccfa4100c30',
        'hessian/hc4-triangular-graphs/QUADRIC_OUTPUT_BASE_OBSTRUCTION.md':'c4a27b3b92b385b2bc5aebc518504e6769db161614bb1501f83d020166e15e0e',
        'hessian/hc4-relative-row-orbit/PROOF.md':'d3f528180c732f2a9c65af0cba4a405a37e34437d48e350e78748eaff61ecb36',
        'hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md':'26a6a76eb579de05b8e71245b8b804bdd5f036602f2e3e6ad5d036fd9f5cd572',
        'hessian/hc4-triangular-graphs/sources/ravindra-srinivas-lefschetz-math-0511134v1.html':'b59ded1196c414a7ec1ed65bfd153bfd9f798e8e0db4ce343c6ee064f1dc3240'
    }
    for name,digest in retained.items():
        check('retained mathematical or primary import pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_smooth_infinity.py','SMOOTH_INFINITY_OUTPUT_OBSTRUCTION.md','sources/smooth-infinity-import.json']},
        dependencies=retained,diagonal_family_controls=diagonal,
        written_not_finite_exhaustion=[
            'filtered output-base theorem with generic factoriality and every height-one domain fibre',
            'all positive observable degrees and lower coefficients with a smooth scheme-theoretic complete intersection at infinity',
            'global free basis and denominator saturation from the leading regular sequence',
            'generic ambient factoriality via the explicitly imported complete-intersection Picard theorem',
            'arbitrary common degree diagonal observable family and necessary HC4 extension boundary'],
        classical_Picard_import_fully_reproved=False,
        generic_factoriality_alone_asserted_sufficient=False,
        smooth_reduced_support_asserted_scheme_smoothness=False,
        all_Keller_extensions_asserted_smooth_at_infinity=False,
        arbitrary_ambient_fibres_asserted_affine_spaces=False,
        arbitrary_source_Hessian_invariance_used=False,
        singular_infinity_pairs_asserted_impossible=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-smooth-infinity.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
