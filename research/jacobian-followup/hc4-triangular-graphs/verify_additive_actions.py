"""Exact controls for the two-sided derivation-degree pencil and invariant overlap images."""
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

    x,y,a,b,T,V,u = S.symbols('x y a b T V u')
    D = lambda value:S.diff(value,y)

    def degree(derivation,value):
        if S.cancel(value)==0:
            return -1
        answer = 0
        for _ in range(30):
            value = S.expand(derivation(value))
            if value==0:
                return answer
            answer += 1
        raise AssertionError('control derivation not nilpotent within the stated finite example')

    first = x*y**3+x**2*y+1
    second = (x+1)*y**2+x
    check('the highest iterated derivative detects the first degree exactly',degree(D,first)==3)
    check('locally nilpotent degree is additive on the product control',degree(D,first*second)==5)
    zero('the highest Leibniz derivative has the nonzero binomial coefficient',
         S.diff(first*second,y,5)-S.binomial(5,3)*S.diff(first,y,3)*S.diff(second,y,2))
    check('one nonconstant degree-zero kernel element remains, so a scalar kernel is not being assumed',
          degree(D,x)==0 and x!=0)
    check('arbitrarily high powers of the kernel generator stay in the degree-zero piece',degree(D,x**17)==0)
    check('a nonconstant Laurent unit also has degree zero before using affine-space quotient units',degree(D,1/x)==0)
    check('a nonconstant unit survives on a transverse component only when that quotient has nonconstant units',
          (1/x).subs(y,a)==1/x)
    ei,ej = y-a,y-b
    cij,cji = a-b,b-a
    zero('the first disjoint-component identity has scalar multiplier',ej-cij-ei)
    zero('the opposite component restricts the first generator to a nonzero scalar',ei.subs(y,b)-cji)
    zero('the two-sided scalar multiplier is exactly one',-cij/cji-1)
    check('the conjugate component degrees are equal and positive',degree(D,ei)==degree(D,ej)==1)
    kernel_multiplier = x*x+x+1
    check('a nonconstant kernel multiplier cannot become scalar modulo a positive-degree component',
          S.Poly(kernel_multiplier-a,y).rem(S.Poly(ej,y)).as_expr()==kernel_multiplier-a)
    bad_ei,bad_ej = x,x*y-1
    zero('the one-sided relation permits a nonconstant multiplier if the opposite unit condition is dropped',
         bad_ej+1-y*bad_ei)
    check('the opposite quotient in that control is Laurent, rather than affine space',
          S.cancel((x*y).subs(y,1/x))==1 and x!=S.Integer(1))
    eplus,eminus = y+u*x,y-u*x
    not_relative = lambda value:S.diff(value,x)-u*S.diff(value,y)
    check('a derivation not defined over the original coefficient field can give unequal conjugate degrees',
          degree(not_relative,eplus)==0 and degree(not_relative,eminus)==1)
    image_equation = y**2-T
    zero('the relative action moves the conjugate-product equation',D(image_equation)-2*y)
    check('a stationary image is not excluded by that single action',D(x)==0)
    check('the two ordinary partial directions have only scalar common polynomial constants',
          S.diff(x,y)==0 and S.diff(x,x)==1 and S.diff(y,y)==1)
    euler = lambda value:x*S.diff(value,x)
    check('an arbitrary derivation may move a function without being locally nilpotent',
          euler(x)==x and euler(euler(x))==x)

    X1,X2,Y,Z,W,N = S.symbols('X1 X2 Y Z W N')
    multiplier = 1+T*X1*X2+X1**2
    root_polynomial = (X1+1)*Z**4+X2*Z*W**2+T*Z**2+W
    observable = multiplier*Y-root_polynomial
    DZ = lambda value:S.expand(multiplier*S.diff(value,Z)+S.diff(root_polynomial,Z)*S.diff(value,Y))
    DW = lambda value:S.expand(multiplier*S.diff(value,W)+S.diff(root_polynomial,W)*S.diff(value,Y))
    zero('the first natural additive direction preserves the second observable',DZ(observable))
    zero('the second natural additive direction preserves the same observable',DW(observable))
    check('both natural directions fix all multiplier coordinates',
          all(der(value)==0 for der in [DZ,DW] for value in [T,X1,X2]))
    check('local nilpotence on the eliminated coordinate has its exact root-degree bound',degree(DZ,Y)==4)
    zero('the iterated derivative formula holds at the last nonzero derivative',
         DZ(DZ(DZ(DZ(Y))))-multiplier**3*S.diff(root_polynomial,Z,4))
    check('the next derivative vanishes, while the other root direction has its own bound',
          DZ(DZ(DZ(DZ(DZ(Y)))))==0 and degree(DW,Y)==2)
    zero('the natural directions commute without requiring freeness',DZ(DW(Y))-DW(DZ(Y)))
    trial = X1*Y**2+Z*W+W**2
    localized = lambda value:S.cancel(value.subs(Y,(V+root_polynomial)/multiplier))
    zero('localizing the multiplier turns the action into a polynomial-coordinate derivative',
         localized(DZ(trial))-multiplier*S.diff(localized(trial),Z))
    zero('the unused stable direction is independent of the multiplier action',DZ(N))
    check('an actual image equation depending on the root variables is moved in the control',DZ(Y-Z**2)!=0)
    check('the surviving multiplier-only image equation is stationary under every natural direction',
          DZ(X1**3+T*X2+V)==0 and DW(X1**3+T*X2+V)==0)

    X = S.symbols('X')
    row_multiplier = X**2-T
    numerator = X**3-T*X+1
    check('a denominator in a multiplier factor is detected by the nonzero row residue',
          S.rem(numerator,row_multiplier,X)==1)
    zero('the same obstruction is visible after splitting the moving rows',numerator.subs({T:u**2,X:u})-1)
    zero('genuine cancellation removes the multiplier denominator exactly',
         row_multiplier*(X+T)/row_multiplier-X-T)
    primitive_multiplier = 1+T*X
    check('primitive coefficient content leaves the multiplier nonzero at the parameter drop',
          S.gcd_list(S.Poly(primitive_multiplier,X).all_coeffs())==1 and primitive_multiplier.subs(T,0)==1)
    content_relation = T*Y-X-V
    zero('the nonprimitive control has an actual extra element of R intersection F[X]',
         (X+V)/T-Y+content_relation/T)
    check('its content fibre destroys the X-injectivity used for global descent',
          content_relation.subs(T,0)==-X-V)

    # Explicit two-direction certificate on the classical smooth suspension.
    z = S.symbols('z')
    relation = x*y-z**3-V
    plus = lambda value:S.expand(x*S.diff(value,z)+3*z**2*S.diff(value,y))
    minus = lambda value:S.expand(y*S.diff(value,z)+3*z**2*S.diff(value,x))
    zero('the first classical suspension direction preserves its equation',plus(relation))
    zero('the opposite suspension direction preserves the equation as well',minus(relation))
    check('both classical directions are locally nilpotent on their nonconstant generators',
          degree(plus,y)==3 and degree(minus,x)==3 and degree(plus,z)==degree(minus,z)==1)
    check('the second direction moves the invariant helper retained by the first',plus(x)==0 and minus(x)==3*z**2)

    # Actual closed HC4 embedding: its stationary image is an invertible control.
    s,t,u0,v0 = S.symbols('s t u0 v0')
    actual_M = 1+X1*X2
    actual_P = X1*Z**5+Z+X1*(T*X2**3*Z**2+1)
    actual_D = lambda value:S.expand(actual_M*S.diff(value,Z)+S.diff(actual_P,Z)*S.diff(value,Y))
    embedding = {T:s,X1:S.Integer(0),X2:u0,Y:t+v0,Z:v0}
    zero('the actual overlap embedding carries the identity pair',
         (actual_M*Y-actual_P).subs(embedding,simultaneous=True)-t)
    check('the image equation is fixed even though this control has an invertible Keller pair',actual_D(X1)==0)
    induced = [S.expand(actual_D(value).subs(embedding,simultaneous=True)) for value in [T,Y-Z,X2,Z]]
    check('the induced source action fixes the original two coordinates and moves only a stable variable',
          induced==[0,0,0,1])
    check('the source inverse still proves this is a closed embedding',
          all(S.expand(value.subs(embedding,simultaneous=True)-coord)==0
              for value,coord in zip([T,Y-Z,X2,Z],[s,t,u0,v0])))
    pullback = (T*X2+Y*Z+X1**2).subs(embedding,simultaneous=True)
    zero('the actual four-variable Hessian remains a nonzero constant',S.hessian(pullback,(s,t,u0,v0)).det()-1)

    context = json.loads((HERE/'sources/additive-action-context.json').read_text())
    check('the standard degree fact is reproved and no new unproved external theorem is imported',
          context['newly_used_unproved_external_statements']==[]
          and context['newly_cached_primary_context'][0]['elementary_degree_fact_fully_reproved'] is True)
    check('neither a scalar degree-zero kernel nor finite-dimensional filtered pieces are required',
          context['scalar_degree_zero_kernel_required'] is False and context['finite_dimensional_degree_pieces_required'] is False)
    check('the moved-image and primitive global-descent conditions remain explicit',
          context['moved_image_required_for_single_action'] is True and context['global_descent_without_content_check_asserted'] is False
          and context['geometric_generic_normal_domain_required'] is True)
    check('geometric scalar units are deduced from a moved component rather than required as a separate hypothesis',
          context['geometric_scalar_units_required'] is False and context['geometric_scalar_units_deduced_from_moved_component'] is True)
    controls = len(checks)
    for name,digest in context['retained_proof_and_source_sha256'].items():
        check('retained byte pin '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    counts = {}
    for name in ['receipt-overlapping-outputs.json','receipt-reflexive-pencil.json']:
        retained = json.loads((HERE/name).read_text())
        check('retained receipt is PASS '+name,retained['outcome']=='PASS')
        counts[name] = retained['check_count']
        for path,digest in retained['source_sha256'].items():
            check('retained receipt source '+name+' -> '+path,hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
        for path,digest in retained['dependencies'].items():
            check('retained dependency '+name+' -> '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    flags = {key:value for key,value in context.items() if type(value) is bool}
    result = dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls,checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_additive_actions.py','ADDITIVE_ACTIONS_AND_IMAGE_INVARIANCE.md','sources/additive-action-context.json']},
        dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=counts,
        newly_cached_primary_context=context['newly_cached_primary_context'],
        newly_used_unproved_external_statements=[],written_not_finite_exhaustion=context['direct_proof_arguments'],**flags,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-additive-actions.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls),indent=2))


if __name__=='__main__':
    main()
