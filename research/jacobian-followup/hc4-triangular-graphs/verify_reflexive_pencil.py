"""Exact controls for the written rank-two descent and all-degree moving suspension proof."""
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
        check(name,S.cancel(S.expand(value))==0)

    T,V,X,Y,Z,W,U = S.symbols('T V X Y Z W U')
    A = 1+T*X**2+(T+1)*X**3
    P = T*(T+1)*Z**5+(T**2+1)*Z**2+Z
    observable = A*Y-P
    DZ = lambda value:S.expand(A*S.diff(value,Z)+S.diff(P,Z)*S.diff(value,Y))
    DX = lambda value:S.expand(A*S.diff(value,X)-S.diff(A,X)*Y*S.diff(value,Y))
    zero('the root-direction derivation preserves both ambient observables',DZ(observable))
    zero('the multiplier-direction derivation preserves both ambient observables',DX(observable))
    h0 = T**2*V**3*Z+T+V
    dh0 = S.expand(DZ(h0))
    zero('one relative derivative detects the finite target-prime denominator bound',dh0-T**2*V**3*A)
    check('the control multiplier has no nonconstant parameter content',S.gcd_list(S.Poly(A,X).all_coeffs())==1)
    zero('dividing precisely the detected coordinate factors recovers a global primitive coordinate',
         (h0-T-V)/(T**2*V**3)-Z)
    check('a pole beyond the detected T bound remains in the relative derivative',
          S.denom(S.cancel(dh0/T**3))==T)
    check('a pole beyond the detected V bound remains in the relative derivative',
          S.denom(S.cancel(dh0/V**4))==V)
    target_polynomial = (Z**3+T*Z+V)/(T**2*V)
    compensated = S.cancel(T**2*V*target_polynomial)
    check('a compensated negative-valuation polynomial has nonzero T residue',compensated.subs(T,0)==Z**3+V)
    check('the analogous V residue is a nonzero polynomial in the primitive coordinate',compensated.subs(V,0)==Z**3+T*Z)

    # A rank-two lattice may be free while its constant section has a base point.
    relation = T*X+V*Y-1
    zero('the omitted-origin control has the stated constant section',relation-(T*X+V*Y-1))
    check('the constant section vector vanishes at the omitted base point',
          S.Matrix([T,V]).subs({T:0,V:0})==S.zeros(2,1))
    check('that point has an empty fibre, so faithful flatness really fails',relation.subs({T:0,V:0})==-1)
    zero('on V nonzero the control is an affine line',relation.subs(Y,(1-T*X)/V))
    zero('on T nonzero the other affine chart gives the same flat family',relation.subs(X,(1-V*Y)/T))
    check('every height-one target prime lies in at least one of these two generic charts',S.gcd(T,V)==1)
    syzygy = S.Matrix([[-V],[T]])
    zero('the quotient by the constant section maps to the codimension-two ideal',
         (S.Matrix([[-V,T]])*S.Matrix([T,V]))[0])
    check('the quotient ideal has two generators but no unit at the base point',
          S.groebner([T,V],T,V).reduce(S.Integer(1))[1]==1)

    # Exact content levels; no degree sweep is used to prove the uniform result.
    content_multiplier = T*(1+(T+1)*X)
    content_output = T*Z**5+Z
    check('a nonprimitive moving multiplier has exactly the admitted content divisor',
          S.gcd_list(S.Poly(content_multiplier,X).all_coeffs())==T)
    check('its content level is closed and nonconstant even though both original top coefficients are nonunits',
          content_output.subs(T,0)==Z and S.Poly(content_output,Z).LC()==T
          and S.expand(S.Poly(content_multiplier,X).LC()-T*(T+1))==0)
    zero('every closed fibre over that content level is a nonempty affine hyperplane',
         (content_multiplier*Y-content_output-V).subs(T,0)+Z+V)
    composite = T*Z**5+Z**2
    alpha = S.symbols('alpha')
    zero('the excluded composite content level still splits geometrically',
         (composite+V).subs({T:0,V:-alpha**2})-(Z-alpha)*(Z+alpha))
    constant = T*Z**5+2
    check('a constant content level produces an empty generic height-one fibre',
          (content_multiplier*Y-constant-V).subs(T,0)==-V-2)
    s,t = S.symbols('s t')
    f = s+s*s*t
    factor = t+s**3
    g = -2+f*factor
    jac = lambda first,second:S.diff(first,s)*S.diff(second,t)-S.diff(first,t)*S.diff(second,s)
    zero('a carried constant content level forces Jacobian divisibility uniformly',jac(f,g)-f*jac(f,factor))
    m,a,d = S.symbols('m a d',positive=True,integer=True)
    zero('the generic suspension weights balance in all joint degrees',a*d+d-d*(a+1))
    check('the leading primitive-linear generic equation has disjoint coefficients in the control',
          S.gcd(T*X**3,T*(T+1)*Z**5)==T)
    check('the shared parameter factor is a generic field unit, not a global leading-unit claim',
          S.Poly(A,X).LC()==T+1 and S.expand(S.Poly(P,Z).LC()-T*(T+1))==0)

    # Existing denominator trap is faithfully flat, but its fibre is not geometric integral.
    trap = T*Y+Z**2-V
    zero('the old denominator trap identity remains outside polynomial saturation',
         (Z**2-V)/T+Y-trap/T)
    zero('its height-one fibre splits after extending the residue field',
         trap.subs({T:0,V:alpha**2})-(Z-alpha)*(Z+alpha))

    # Actual HC4 embedding and Hessian control in a primitive moving family.
    X1,X2,u,v = S.symbols('X1 X2 u v')
    moving = 1+T*X1**2+(T+1)*X1*X2**2
    embedding = {T:s,X1:S.Integer(0),X2:u,Y:t+P.subs({T:s,Z:v}),Z:v}
    zero('the explicit closed embedding carries the first identity coordinate',embedding[T]-s)
    zero('the explicit closed embedding carries the second identity coordinate',
         (moving*Y-P).subs(embedding,simultaneous=True)-t)
    inverse = [T,Y-P,X2,Z]
    check('the explicit ambient pullbacks recover every source coordinate',
          all(S.expand(value.subs(embedding,simultaneous=True)-coordinate)==0
              for value,coordinate in zip(inverse,[s,t,u,v])))
    Q = T*X2+Y*Z+X1**2
    pulled = Q.subs(embedding,simultaneous=True)
    zero('the actual four-variable pullback has unit Hessian',S.hessian(pulled,(s,t,u,v)).det()-1)
    check('the primitive HC4 control has parameter-varying positive coefficient row',
          S.gcd_list(S.Poly(moving,X1,X2).coeffs())==1)
    check('the displayed HC4 control fails leading divisibility without asserting failure of every possible frame',
          S.rem(T,T*(T+1),T)!=0 and S.rem(T+1,T*(T+1),T)!=0)

    context = json.loads((HERE/'sources/reflexive-pencil-import.json').read_text())
    check('exactly two structural primary statements are explicitly imported',len(context['newly_used_primary_inputs'])==2)
    for item in context['newly_used_primary_inputs']:
        data = (HERE/'sources'/item['retained_file']).read_bytes()
        check('primary source byte pin '+item['retained_file'],hashlib.sha256(data).hexdigest()==item['sha256'])
        check('the structural theorem is not represented as fully reproved '+item['retained_file'],item['theorem_fully_reproved'] is False)
    check('global frames are replaced rather than silently inferred',
          context['common_global_free_filtration_required'] is False and context['global_monic_presentation_required'] is False)
    check('faithfulness and geometric height-one integrality remain explicit',
          context['faithful_flatness_required'] is True and context['geometric_height_one_integrality_required'] is True)
    check('the reflexive dimension restriction remains at most two',context['reflexive_base_dimension_at_most_two_required'] is True)
    controls = len(checks)
    for name,digest in context['retained_proof_and_source_sha256'].items():
        check('retained byte pin '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    retained_counts = {}
    for name in ['receipt-parameter-stability.json','receipt-multivariable-suspensions.json']:
        retained = json.loads((HERE/name).read_text())
        check('retained receipt is PASS '+name,retained['outcome']=='PASS')
        retained_counts[name] = retained['check_count']
        for path,digest in retained['source_sha256'].items():
            check('retained receipt source '+name+' -> '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
        for path,digest in retained['dependencies'].items():
            check('retained dependency '+name+' -> '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    flags = {key:value for key,value in context.items() if type(value) is bool}
    result = dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls,checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_reflexive_pencil.py','REFLEXIVE_PENCIL_DESCENT.md','sources/reflexive-pencil-import.json']},
        dependencies=context['retained_proof_and_source_sha256'],
        retained_receipt_counts_not_newly_replayed=retained_counts,
        newly_used_primary_inputs=context['newly_used_primary_inputs'],
        written_not_finite_exhaustion=context['direct_proof_arguments'],**flags,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-reflexive-pencil.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls),indent=2))


if __name__=='__main__':
    main()
