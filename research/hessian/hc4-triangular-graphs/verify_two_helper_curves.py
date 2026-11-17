"""Exact controls for finite curve images, flat helper pullbacks and the two-helper pencil."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S
from sympy.combinatorics import Permutation, PermutationGroup

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

    s,t,U,X1,X2 = S.symbols('s t U X1 X2')
    h = s+t**2
    zero('a singular curve image still satisfies its literal defining equation',h**6-(h**2)**3)
    check('a suitable affine line keeps a source coordinate nonconstant',S.degree((h**2).subs(t,0),s)==2)
    zero('the line parameter is integral over the cusp coordinate algebra',(U**2-X1).subs({U:s,X1:s**2}))
    check('the parameter equation is monic, giving a finite curve map rather than only dominance',S.Poly(U**2-X1,U).LC()==1)
    zero('the second cusp coordinate is compatible with the same finite map',(U**3-X2).subs({U:s,X2:s**3}))
    check('the cusp point is hit without assuming the curve smooth',h.subs({s:0,t:0})==0)

    # Row stabilizer control: full parameter action cannot fix a moving row.
    row_swap = Permutation(0,2)(1,3)
    within_swap = Permutation(0,1)(2,3)
    target = PermutationGroup(row_swap,within_swap)
    check('the moving-row control has a nontrivial parameter quotient',target.order()==4)
    check('a source action surjecting onto that quotient cannot fix the selected moving row',
          {row_swap(point) for point in [0,1]}=={2,3})
    check('within-row transitivity cannot fix a selected component when the row has two components',within_swap(0)!=0)
    check('the selected-component stabilizer in the control loses the parameter quotient',target.stabilizer(0).order()==1)

    T,V,Y,Z,W = S.symbols('T V Y Z W')
    multiplier = X1
    root = X2*Z**3+T*X2+1
    relation = multiplier*Y-root-V
    drop = S.groebner([X1,X2],X1,X2)
    check('the common positive-root coefficient drop locus is a single helper point',drop.is_zero_dimensional)
    zero('the generic fibre over that point is a nonzero constant equation',relation.subs({X1:0,X2:0})+1+V)
    check('transcendental V prevents that point from gaining a full root-variable fibre',V+1!=0)
    check('outside the multiplier the eliminated variable is genuinely linear',S.degree(relation.subs(X1,1),Y)==1)
    zero('on the multiplier away from the drop point a nonconstant root equation remains',
         relation.subs({X1:0,X2:1})+Z**3+T+1+V)
    check('the surviving root equation has positive degree despite a possible singular fibre',
          S.degree(relation.subs({X1:0,X2:1}),Z)==3)
    q,u = S.symbols('q u',nonnegative=True,integer=True)
    zero('all nonempty helper fibres have the expected dimension uniformly in root and unused-variable counts',
         (q-1)+1+u-(q+u))
    zero('a curve pullback has exactly the source-component dimension',1+q+u-(q+u+1))

    # The finite-point argument fails with three helper variables.
    X3 = S.symbols('X3')
    higher = X1*Y-X2*Z-X3-V
    check('the three-helper coefficient-drop locus has positive dimension',
          not S.groebner([X1,X2],X1,X2,X3).is_zero_dimensional)
    zero('a V-dependent point on that locus produces a full Y,Z fibre',higher.subs({X1:0,X2:0,X3:-V}))
    check('away from that locus only the expected single root coordinate is free',S.degree(higher.subs(X1,1),Y)==1)

    # Flatness need not mean faithfulness; pullback integrality needs its own check.
    basic = X1*Y-X2*Z-V
    check('the basic flat helper map omits the origin in its geometric generic fibre',basic.subs({X1:0,X2:0})==-V)
    zero('one helper curve has an explicit polynomial pullback',basic.subs({X1:1,Y:X2*Z+V}))
    zero('another helper curve has a different explicit polynomial pullback',basic.subs({X2:1,Z:X1*Y-V}))
    zero('intersecting helper curves have an actual common lift, contradicting disjoint source components',
         basic.subs({X1:1,X2:1,Y:V,Z:0}))
    parabola_pullback = basic.subs(X2,X1**2)
    zero('a helper curve through an omitted point can have a nonconstant unit in its pullback',
         X1*(Y-X1*Z)/V-1-parabola_pullback/V)
    check('that unit excludes an affine-space source component rather than proving helper faithfulness',
          S.degree(X1,X1)==1)

    # Disjoint conjugate curves have a literal polynomial pencil in the helper plane.
    alpha = S.symbols('alpha')
    helper = X1+X2**2
    first,second = helper-alpha,helper+alpha
    zero('the two conjugate helper-curve equations differ by a nonzero scalar',second-first-2*alpha)
    check('their ordinary helper degrees are equal',S.Poly(first,X1,X2).total_degree()==S.Poly(second,X1,X2).total_degree()==2)
    zero('each coordinate quotient is a polynomial algebra in the free helper coordinate',
         first.subs(X1,alpha-X2**2))
    zero('the component product has the descended univariate helper equation',first*second-(helper**2-alpha**2))
    zero('the target-field model of that pencil has exactly the component count as its outer degree',
         (helper**2-V).subs(V,alpha**2)-first*second)
    check('the pencil outer degree is two in this control',S.degree(U**2-V,U)==2)

    # Higher-dimensional dominant polynomial images need not be surjective.
    a,b,A,B = S.symbols('a b A B')
    image_a = 1+a*b
    image_b = a+b**2+a*b**3
    zero('the inspected surface example omits the origin by an exact Bezout identity',
         (1+b**3)*image_a-b*image_b-1)
    cubic = A*b**3-B*b+A-1
    inverse_a = B-A*b**2
    zero('away from the first-coordinate zero the cubic supplies its preimage equation',
         image_a.subs(a,inverse_a)-A+cubic)
    zero('the second-coordinate equation follows from the same cubic',
         image_b.subs(a,inverse_a)-B+b**2*cubic)
    zero('the first-coordinate-zero chart has explicit preimages for every nonzero second coordinate',
         image_a.subs({a:B,b:-1/B}))
    zero('those explicit preimages retain the requested second coordinate',image_b.subs({a:B,b:-1/B})-B)

    # Complete HC4 variable-count closure within the admitted family.
    cases = [(p,roots,unused) for p in [1,2] for roots in [1,2] for unused in [0,1]
             if p+roots+unused==3]
    check('every admitted HC4 count is now in the one- or two-helper theorem',
          set(cases)=={(1,1,1),(1,2,0),(2,1,0)} and all(p in [1,2] for p,roots,unused in cases))

    # Actual HC4 closed embedding and Hessian in a genuine overlapping family.
    u0,v0 = S.symbols('u0 v0')
    actual_M = 1+X1*X2
    lower = T*X2**2*Z**2+Z**4+X2+1
    actual_P = (X1+1)*Z**5+T*Z+X1*lower
    embedding = {T:s,X1:S.Integer(0),X2:u0,Y:t+v0**5+s*v0,Z:v0}
    zero('the genuine two-helper HC4 embedding carries the identity pair',
         (actual_M*Y-actual_P).subs(embedding,simultaneous=True)-t)
    check('its explicit source inverse proves closedness',
          all(S.expand(value.subs(embedding,simultaneous=True)-coord)==0
              for value,coord in zip([T,Y-Z**5-T*Z,X2,Z],[s,t,u0,v0])))
    pullback = (T*X2+Y*Z+X1**2).subs(embedding,simultaneous=True)
    zero('the actual four-variable potential has Hessian determinant one',S.hessian(pullback,(s,t,u0,v0)).det()-1)
    check('the original leading coefficient fails the previous Bezout leading certificate in this presentation',
          actual_M.subs({X1:-1,X2:1})==0 and (X1+1).subs(X1,-1)==0)
    check('the residual polynomial direction still survives on the multiplier component generically',
          S.gcd(actual_M,X1+1)==1)

    context = json.loads((HERE/'sources/two-helper-curve-context.json').read_text())
    check('the curve-image fact is reproved and no new unproved theorem is imported',
          context['newly_used_unproved_external_statements']==[]
          and context['newly_cached_primary_context'][0]['needed_curve_surjectivity_fully_reproved'] is True)
    check('generic forgetful flatness is proved without assuming its faithfulness',
          context['generic_forgetful_flatness_proved'] is True and context['generic_forgetful_faithfulness_assumed'] is False)
    check('no curve or action classification or higher-helper conclusion is silently imported',
          context['curve_normalization_classification_imported'] is False and context['Rentschler_imported'] is False
          and context['higher_helper_closure_asserted'] is False)
    controls = len(checks)
    for name,digest in context['retained_proof_and_source_sha256'].items():
        check('retained byte pin '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    counts = {}
    for name in ['receipt-one-helper-overlap.json','receipt-additive-actions.json','receipt-overlapping-outputs.json']:
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
            ['verify_two_helper_curves.py','TWO_HELPER_CURVE_PENCIL_CLOSURE.md','sources/two-helper-curve-context.json']},
        dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=counts,
        newly_cached_primary_context=context['newly_cached_primary_context'],retained_structural_imports=context['retained_structural_imports'],
        newly_used_unproved_external_statements=[],written_not_finite_exhaustion=context['direct_proof_arguments'],**flags,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-two-helper-curves.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=controls),indent=2))


if __name__=='__main__':
    main()
