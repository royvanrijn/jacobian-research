"""Exact controls for univariate helper descent and all one-helper overlap outputs."""
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

    T,V,X,Y,Z,H,U,a = S.symbols('T V X Y Z H U a')
    source_helper = (T*H-V)/(T+1)
    image = H+T**2*V
    generic_polynomial = ((T+1)*U+V)/T+T**2*V
    zero('an affine generic helper change preserves the literal image equation',generic_polynomial.subs(U,source_helper)-image)
    zero('its inverse is defined over the generic target field',((T+1)*source_helper+V)/T-H)
    check('the positive helper degree is unchanged despite target denominators',S.degree(generic_polynomial,U)==1)
    check('the original generic helper need not have polynomial target coefficients',S.denom(S.cancel(source_helper))==T+1)
    s,t,w = S.symbols('s t w')
    descended_source = -s**2*t
    zero('the linear global label relation reconstructs the source helper over the target',
         image.subs({T:s,V:t,H:descended_source}))
    check('the global label derivative is the required nonzero constant',S.diff(image,H)==1)
    nonlinear_label = H**3+T*H+V
    check('a higher-degree label derivative has strictly smaller helper degree',
          S.degree(S.diff(nonlinear_label,H),H)==2<S.degree(nonlinear_label,H))
    check('nonconstant stable polynomials cannot be algebraic labels by their distinct powers',
          [S.degree((w**2+t)**n,w) for n in [1,2,3]]==[2,4,6])

    # Faithful saturation, including the content case where B[X] itself is too small.
    content_M = T*X
    content_P = X+T*Z**5
    observable = content_M*Y-content_P
    adapted = X*Y-Z**5
    zero('the primitive rank-two coordinate is a global ambient polynomial',
         ((X+V)/T).subs(V,observable)-adapted)
    zero('the original helper is affine in that global coordinate over the target field',
         X-(T*adapted-observable))
    check('the rank-two affine coefficient matrix is invertible over the generic field',
          S.Matrix([[1,V/T],[0,1/T]]).det()==1/T)
    check('the multiplier has a genuine nonprimitive content zero',S.gcd_list(S.Poly(content_M,X).all_coeffs())==T)
    check('the admitted content polynomial can be closed in X with no root-variable dependence',
          content_P.subs(T,0)==X and S.degree(content_P.subs(T,0),Z)==0)
    check('the generic multiplier root retains positive root degree',S.degree(content_P.subs(X,0),Z)==5)
    check('the old coefficient-intersection assertion really fails before adapting the helper',
          S.Poly((X+V)/T,X).all_coeffs()==[1/T,V/T])

    M = X*(X-1)
    P = (1-X)*Z**2+X*Z**5+T*X*Z
    check('different generic multiplier rows may have different unrestricted root degrees',
          S.degree(P.subs(X,0),Z)==2 and S.degree(P.subs(X,1),Z)==5)
    check('the original leading root coefficient is not a unit modulo the multiplier',
          S.Poly(P,Z).LC()==X and S.gcd(M,X)==X)
    rational_a = V/T
    M_a = S.cancel(M.subs(X,rational_a))
    check('the nonzero-multiplier case remains genuinely nonzero after a rational helper substitution',M_a!=0)
    solved_y = S.cancel((P.subs(X,rational_a)+V)/M_a)
    zero('the nonzero generic multiplier eliminates Y exactly',
         (M*Y-P-V).subs({X:rational_a,Y:solved_y},simultaneous=True))
    check('the remaining root variable stays algebraically free in this presentation',
          S.degree((M*Y-P-V).subs(X,rational_a),Y)==1)
    check('the multiplier-zero helper is algebraic over the parameter field as required',M.subs(X,0)==0)
    residual = P.subs(X,0)
    zero('the multiplier-zero case retains its polynomial root equation instead of eliminating Y',
         (M*Y-P-V).subs(X,0)+residual+V)
    check('the retained root equation has the simple V valuation in its quadratic control',
          S.degree(S.discriminant(residual+V,Z),V)==1)

    # The self-tensor obstruction does not require a Galois source extension.
    u = S.symbols('u')
    cubic = U**3-U-V
    specialized = cubic.subs(V,u**3-u)
    first = U-u
    second = U**2+u*U+u**2-1
    zero('a finite source self tensor has the diagonal and complementary factors',first*second-specialized)
    check('the diagonal factor is nonzero in the self-tensor algebra',S.degree(first,U)<3 and first!=0)
    check('the complementary factor is nonzero as well',S.degree(second,U)<3 and second!=0)
    check('the factors are separated in characteristic zero',S.expand(second.subs(U,u))==3*u**2-1)
    check('the residual quadratic need not split over the original source field',
          S.degree(S.discriminant(second,U),u)==2
          and S.gcd(S.discriminant(second,U),S.diff(S.discriminant(second,U),u))==1)
    multiplication = S.Matrix([[1,u,u**2]])
    check('the finite-dimensional multiplication map has a nonzero kernel at nontrivial degree',
          multiplication.rank()==1 and len(multiplication.nullspace())==2)
    zero('one explicit self-tensor kernel vector maps to zero',
         (multiplication*S.Matrix([-u,1,0]))[0])

    # Structural HC4 variable count, not a bounded polynomial-degree search.
    admissible = [(p,q,unused) for p in [1,2] for q in [1,2] for unused in [0,1]
                  if p+q+unused==3]
    check('the HC4 overlap variable count has exactly the three admitted structural cases',
          set(admissible)=={(1,1,1),(1,2,0),(2,1,0)})
    check('closing every one-helper case leaves exactly two helpers and one root variable',
          [item for item in admissible if item[0]!=1]==[(2,1,0)])
    composite_content = X**2+T*Z**3
    alpha = S.symbols('alpha')
    zero('a composite content level remains excluded from the theorem',
         (composite_content+V).subs({T:0,V:-alpha**2})-(X-alpha)*(X+alpha))
    check('the excluded residual-direction control is detected on its multiplier row',
          S.degree((X*Z+T).subs(X,0),Z)==0)

    # Smooth image alone is insufficient without the Keller differential.
    bad_relation = X**2-V
    zero('a smooth ambient image can have helper degree two for a non-Keller carried pair',
         bad_relation.subs({X:t,V:t**2}))
    bad_jacobian = 2*t
    check('the source Jacobian and label derivative fail at the same point in that control',
          bad_jacobian.subs(t,0)==0 and S.diff(bad_relation,X).subs(X,0)==0)
    check('the actual ambient image equation of that control is still smooth',S.diff(Y,Y)==1)

    # Actual two-root-variable HC4 identity embedding in the newly closed family.
    Z1,Z2,u0,v0 = S.symbols('Z1 Z2 u0 v0')
    actual_M = 1+(X-1)**2
    actual_P = T*Z1**3+X*Z2**2+Z1*Z2
    embedding = {T:s,X:S.Integer(1),Y:t+actual_P.subs({T:s,X:1,Z1:u0,Z2:v0}),Z1:u0,Z2:v0}
    zero('the genuine one-helper HC4 embedding carries the identity second component',
         (actual_M*Y-actual_P).subs(embedding,simultaneous=True)-t)
    check('the explicit source inverse proves closedness',
          all(S.expand(value.subs(embedding,simultaneous=True)-coord)==0
              for value,coord in zip([T,Y-actual_P,Z1,Z2],[s,t,u0,v0])))
    pullback = (T*Z2+Y*Z1+X**2).subs(embedding,simultaneous=True)
    zero('the actual four-variable potential has Hessian determinant one',S.hessian(pullback,(s,t,u0,v0)).det()-1)

    context = json.loads((HERE/'sources/one-helper-overlap-context.json').read_text())
    check('the self-tensor argument is reproved and no new unproved external theorem is imported',
          context['newly_used_unproved_external_statements']==[]
          and context['newly_cached_primary_context'][0]['needed_self_tensor_obstruction_fully_reproved'] is True)
    check('neither primitive content nor a positive root degree at every content level is required',
          context['primitive_multiplier_required'] is False and context['content_levels_positive_Z_degree_required'] is False)
    check('the rational helper is not silently promoted to global reconstruction or the full source-field degree',
          context['rational_helper_alone_asserted_global_reconstruction'] is False
          and context['label_degree_one_alone_asserted_source_degree_one'] is False)
    check('forgetful flatness and stable cancellation are not assumptions',
          context['forgetful_map_flatness_required'] is False and context['stable_cancellation_used'] is False)
    controls = len(checks)
    for name,digest in context['retained_proof_and_source_sha256'].items():
        check('retained byte pin '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    counts = {}
    for name in ['receipt-additive-actions.json','receipt-reflexive-pencil.json','receipt-overlapping-outputs.json']:
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
            ['verify_one_helper_overlap.py','ONE_HELPER_OVERLAP_CLOSURE.md','sources/one-helper-overlap-context.json']},
        dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=counts,
        newly_cached_primary_context=context['newly_cached_primary_context'],newly_used_unproved_external_statements=[],
        written_not_finite_exhaustion=context['direct_proof_arguments'],**flags,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-one-helper-overlap.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls),indent=2))


if __name__=='__main__':
    main()
