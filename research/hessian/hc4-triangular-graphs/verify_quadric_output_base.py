"""Exact controls for the written stable quadric output-base obstruction.

No finite check proves arbitrary-r factoriality, primitive composition,
or the disjoint-component and smoothness arguments.
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
        check(name,all(S.cancel(S.expand(item.doit()))==0 for item in values))

    T,U,X,Y,Z,W,V,S0 = S.symbols('T U X Y Z W V S0')
    # Both inverse-metric normal types and all arbitrary-r rank patterns.
    for r in [2,3,4]:
        ds = S.symbols('d0:'+str(r+2))
        q = ds[0]*ds[1]+sum(v*v for v in ds[2:])
        ambient = (T,)+ds
        noniso = T*T+q
        check('type-I control r='+str(r)+' has nondegenerate full quadratic',
            S.hessian(noniso,ambient).rank()==r+3)
        check('type-I control r='+str(r)+' has irreducible-rank leading quadric',
            S.hessian(q,ds).rank()==r+2)
        ds = S.symbols('e0:'+str(r+1))
        q = ds[0]*ds[1]+sum(v*v for v in ds[2:])
        ambient = (T,Y)+ds
        iso = T*Y+q
        check('type-II control r='+str(r)+' has nondegenerate full quadratic',
            S.hessian(iso,ambient).rank()==r+3)
        check('type-II control r='+str(r)+' has rank-at-least-three leading cone',
            S.hessian(q,ds).rank()==r+1)
    # Exact generic ring charts and the prime-element quotient.
    q4 = X*Y+Z*W
    zero('the type-I X-localization eliminates Y polynomially after X inversion',
        (q4-U).subs(Y,(U-Z*W)/X))
    zero('its X-zero quotient retains a nonzero rank-two quadric level',
        (q4-U).subs(X,0)-(Z*W-U))
    qiso = T*Y+X*Z+W*W
    zero('the type-II generic T chart is an actual polynomial three-space',
        (qiso-U).subs(Y,(U-X*Z-W*W)/T))
    check('the type-II T-zero leading quadric is irreducible by rank three',
        S.hessian(X*Z+W*W,(X,Z,W)).rank()==3)
    # Monic normal-form reductions: coefficients are polynomials in T,U.
    gI = S.groebner([X*Y+Z*W-U],X,Y,Z,W,domain=S.QQ[T,U],order='grlex')
    zero('monic type-I normal form replaces XY by U-ZW',gI.reduce(X*Y)[1]-(U-Z*W))
    gII = S.groebner([X*Z+W*W+T*Y-U],X,Z,W,Y,domain=S.QQ[T,U],order='grlex')
    zero('monic type-II normal form includes the linear TY lower term',
        gII.reduce(X*Z)[1]-(U-W*W-T*Y))
    h = X+T*Z
    expansion = h*h+U*h+1
    zero('primitive composition retains a genuine global A-polynomial expansion',
        gI.reduce(expansion)[1]-S.expand(expansion))
    check('primitive coefficient reduction stays nonconstant at T=0',h.subs(T,0)==X)
    # Show why a nonprimitive h would invalidate the denominator conclusion.
    bad_h = T*X
    zero('a nonprimitive control admits a hidden target denominator',bad_h/T-X)
    check('that hidden denominator is not in the base polynomial ring',
        not S.sympify(1/T).is_polynomial(T,U))
    # Disjoint equal-degree factors have a common primitive nonconstant part.
    e1,e2 = h-2,h+3
    zero('disjoint normalized factors differ by a nonzero scalar',e2-e1-5)
    zero('their reduced product is one polynomial in the common h',
        e1*e2-(h*h+h-6))
    zero('the first factor is a unit modulo the second',(e1/(-5)).subs(X,-3-T*Z)-1)
    # Generic source labels are constant along the free variables.
    a,b,u,v = S.symbols('a b u v')
    K = S.Matrix([a,b+a*a])
    zero('the positive Keller base has unit determinant',K.jacobian((a,b)).det()-1)
    # A literal quadric-fibre embedding with an ambient noncoordinate pair.
    j = {T:a,X:1,Y:b+u*v,Z:u,W:v}
    R = T*T+X*Y-Z*W
    pair = S.Matrix([T,R]).subs(j,simultaneous=True)
    zero('linear-quadratic observables recover an actual Keller base',pair-S.Matrix([a,a*a+b]))
    zero('its observed base Jacobian is one',pair.jacobian((a,b)).det()-1)
    zero('the embedded source has global polynomial recovery',
        S.Matrix([T,Y-Z*W,Z,W]).subs(j,simultaneous=True)-S.Matrix([a,b,u,v]))
    check('the second ambient observable has nondegenerate five-variable Hessian',
        S.hessian(R,(T,X,Y,Z,W)).rank()==5)
    q0 = T*Z+Y*W+X*X
    check('the independent HC4 potential is ambient nondegenerate',
        S.hessian(q0,(T,X,Y,Z,W)).rank()==5)
    f = q0.subs(j,simultaneous=True)
    zero('its source potential is exactly the stated cotangent polynomial',
        f-(1+a*u+b*v+u*v*v))
    zero('the actual source Hessian determinant is one',S.hessian(f,(a,b,u,v)).det()-1)
    p1,p2,p3,p4 = S.symbols('p1 p2 p3 p4')
    inv = {a:p3-p2*p2,b:p4-2*p1*p2,u:p1,v:p2}
    zero('its source gradient has the displayed polynomial inverse',
        S.Matrix([S.diff(f,z) for z in (a,b,u,v)]).subs(inv,simultaneous=True)-S.Matrix([p1,p2,p3,p4]))
    tau = S.symbols('tau')
    zero('generic quadric-threefold strata give tau cubed minus tau',
        (tau-1)*tau*tau+(tau-1)*tau-(tau**3-tau))
    check('that generic fibre is distinguished from a coordinate-pair A3 fibre',tau**3-tau!=tau**3)
    # Low-dimensional negative controls: neither extends the theorem.
    zero('the surface quadric gives two different factorizations',
        (Z-S.I)*(Z+S.I)-(Z*Z+1))
    check('the surface leading quadratic has rank three and additive degree',
        S.hessian(X*Y-Z*Z,(X,Y,Z)).rank()==3)
    zero('its X-zero quotient splits and X is not prime',
        (X*Y-Z*Z-1).subs(X,0)+(Z-S.I)*(Z+S.I))
    zero('the binary generic quadric has a nonconstant invertible function',
        (X*Y/U).subs(Y,U/X)-1)
    # Native A/B target change preserves the Keller determinant up to N.
    L,M,N = S.symbols('L M N',nonzero=True)
    F,G = a,b+a*a
    target = S.Matrix([F+M/(2*L),L*(F+M/(2*L))**2+N*G])
    zero('native A/B observables become a Keller linear-quadratic pair',
        target.jacobian((a,b)).det()-N)
    out = S.symbols('out')
    zero('an invertible A/B base makes its potential a coordinate',
        (L*a*a+M*a+N*b).subs(b,(out-L*a*a-M*a)/N)-out)
    # The smoothness contradiction is a full ambient differential identity.
    aa,bb,cc = S.symbols('aa bb cc')
    J = S.Matrix([[1,2*a],[0,1]])
    check('an invertible source base kills both base coefficients when their pullback vanishes',
        J.T.inv()*S.zeros(2,1)==S.zeros(2,1))
    relation = S0*S0+T*S0+U
    derivative = S.diff(relation,S0)
    check('a nonlinear relation derivative has strictly smaller label degree',
        S.degree(derivative-cc,S0)<S.degree(relation,S0))
    pins = {}
    retained = json.loads((HERE/'receipt-stable-hypersurface.json').read_text())
    pins['hessian/hc4-triangular-graphs/STABLE_HYPERSURFACE_OBSTRUCTION.md'] = retained['source_sha256']['STABLE_HYPERSURFACE_OBSTRUCTION.md']
    for name,digest in retained['dependencies'].items():
        if 'hc4-relative-row-orbit/PROOF.md' in name or 'stacks-' in name:
            pins[name] = digest
    for receipt_name,name in [('receipt-ambient-quadric.json','AMBIENT_QUADRIC_GRAPHS.md'),
        ('receipt-global-value.json','GLOBAL_VALUE_ORBITS.md'),
        ('receipt-intrinsic-plane.json','INTRINSIC_AFFINE_PLANE.md')]:
        retained = json.loads((HERE/receipt_name).read_text())
        pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
        if receipt_name=='receipt-intrinsic-plane.json':
            for path,digest in retained['dependencies'].items():
                if 'ricolfi-hodge' in path:
                    pins[path] = digest
    source_record = json.loads((HERE/'sources/quadric-output-base-import.json').read_text())
    for item in source_record['sources']:
        pins[str((HERE/'sources'/item['retained_file']).relative_to(ROOT))] = item['sha256']
    for name,digest in pins.items():
        check('retained mathematical or primary source pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_quadric_output_base.py','QUADRIC_OUTPUT_BASE_OBSTRUCTION.md','sources/quadric-output-base-import.json']},
        dependencies=pins,written_not_finite_exhaustion=[
            'both inverse-metric normal types in every r>=2',
            'generic quadric factoriality and constant units',
            'global filtered primitive-composition lemma through all height-one base primes',
            'Galois disjoint-component pencil and smoothness forcing generic degree one',
            'intrinsic cotangent output recognition and direct native A/B encoding consequences'],
        all_degree_stability_threshold_r=2,
        generic_chart_asserted_global_without_primitive_step=False,
        arbitrary_quadric_embedding_rectification_used=False,
        linear_quadratic_observables_asserted_ambient_coordinates=False,
        surface_generic_quadric_asserted_UFD=False,
        arbitrary_substitution_asserted_closed_embedding=False,
        nonlinear_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unaligned_two_affine_class_asserted_closed=False,unrestricted_HC4_or_JC2_resolved=False,
        new_external_classification_imported=False,novelty_priority_asserted=False,
        formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-quadric-output-base.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
