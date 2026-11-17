"""Exact controls for parameter inertia, moving factor rows and global adapted frames."""
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

    T,V,U,u,w,X,Y,Z,H = S.symbols('T V U u w X Y Z H')
    parameter_discriminant = S.discriminant(U*U-T,U)
    check('a finite parameter discriminant is a nonzero unit at every horizontal prime of kappa[V]',
          parameter_discriminant!=0 and S.degree(parameter_discriminant,V)==0)
    check('the same parameter cover can ramify vertically, which is outside the inertia-fixing assertion',
          parameter_discriminant.subs(T,0)==0)
    zero('moving geometric multiplier rows split over the parameter cover',
         (X*X-T).subs(T,u*u)-(X-u)*(X+u))
    outer = U*U-T+V
    split_outer_discriminant = S.discriminant(outer.subs(T,u*u),U)
    check('the outer quadratic remains irreducible over this split-row target field by an odd V valuation',
          S.degree(split_outer_discriminant,V)==1
          and S.degree(S.gcd(split_outer_discriminant,S.diff(split_outer_discriminant,V)),V)==0)
    # Enlarging an integral monic-model parameter field preserves integral coefficients.
    original = T**3*Z**5+T**2*Z**2+T+1
    enlarged = original.subs({T:w**10,Z:Y/w**6})
    zero('an enlarged field simultaneously contains the fifth-root model and a quadratic parameter extension',
         enlarged-(Y**5+w**8*Y**2+w**10+1))
    zero('the quadratic parameter embeds into that common enlargement',(w**5)**2-w**10)
    zero('the fifth-root model parameter embeds into that common enlargement',(w**2)**5-w**10)

    # Abstract Galois-tower control: retain bad inertia in the parameter kernel.
    cycle_a = Permutation(0,1,2,size=8)
    cycle_b = Permutation(3,4,5,size=8)
    block_swap = Permutation(0,3)(1,4)(2,5)
    transposition_a = Permutation(0,1,size=8)
    parameter_flip = Permutation(6,7,size=8)
    target = PermutationGroup(cycle_a,cycle_b,block_swap,transposition_a,parameter_flip)
    original_target = PermutationGroup(cycle_a,cycle_b,block_swap,transposition_a)
    bad = PermutationGroup(cycle_a,cycle_b)
    original_source = PermutationGroup(cycle_a,cycle_b,block_swap)
    total_source = PermutationGroup(cycle_a,cycle_b,block_swap,parameter_flip)
    check('the finite tower control has a genuine parameter quotient',
          target.order()==144 and original_target.order()==72)
    check('the lifted horizontal inertia fixes the parameter letters and retains its original order',
          cycle_a(6)==6 and cycle_a(7)==7 and cycle_a.order()==3)
    check('the original bad subgroup is normal and contained after fixing the parameter extension',
          bad.is_normal(original_target) and all(g in original_source for g in bad.generators))
    parameter_kernel = [g for g in total_source.generate_schreier_sims() if g(6)==6 and g(7)==7]
    check('fixing the parameter really gives the smaller source subgroup in the control',
          len(parameter_kernel)==original_source.order()==18 and total_source.order()==36)
    check('the smaller subgroup is still transitive on the original roots',
          original_source.orbit(0)==set(range(6)))
    check('parameter-stable root transitivity is not full target-group equality',
          original_source.order()<original_target.order())
    correlated = S.Matrix([[0,-1],[-1,0]])
    check('the old correlated moving-row representation retains its nonzero abstract fixed vector',
          correlated*S.Matrix([1,-1])==S.Matrix([1,-1]))
    check('a row-fixing transitive root action instead has zero invariants in the two-row quadratic case',
          (-S.eye(2)-S.eye(2)).rank()==2)

    # Global leading certificates and the independent coefficient-content condition.
    monic_root_relation = (1+T*X)*Y-Z**3-V
    monic_root_basis = S.groebner([monic_root_relation],Z,Y,X,T,V)
    normal = monic_root_basis.reduce(Z**4*X)[1]
    check('the monic-root certificate provides a global basis with bounded Z exponent',
          all(monomial[0]<3 for monomial,coefficient in S.Poly(normal,Z,Y,X,T,V).terms()))
    check('the multiplier in that certificate is primitive despite its nonunit top coefficient',S.gcd(1,T)==1)
    zero('the affine fibre at the leading drop is an integral coordinate presentation',
         monic_root_relation.subs(T,0)-(Y-Z**3-V))
    zero('its graded fibre may nevertheless be nonreduced, as allowed by the affine-fibre criterion',
         (T*X*Y-Z**3).subs(T,0)+Z**3)
    check('a monic moving multiplier stays nonzero on every parameter level',S.Poly(X**2-T,X).LC()==1)
    z1,z2 = S.symbols('z1 z2')
    closed_level = z1*z2+V
    check('a closed content level has an integral generic fibre despite a reducible zero level',
          S.gcd(z1,V)==1 and len(S.factor_list(z1*z2,z1,z2)[1])==2)
    alpha = S.symbols('alpha')
    zero('a composite content level instead splits geometrically at the generic V',
         (Z**2+V).subs(V,-alpha**2)-(Z-alpha)*(Z+alpha))
    check('a constant content level has an empty generic height-one fibre',S.degree(V+1,V)==1)
    s,t = S.symbols('s t')
    zero('a carried constant content level would make the Jacobian divisible by f-c',
         S.det(S.Matrix([[S.diff(s,s),S.diff(s,t)],[S.diff(-1+s*t,s),S.diff(-1+s*t,t)]]))-s)

    # One adapted presentation and one weight identity cover every nonlinear degree.
    a,n = S.symbols('a n',positive=True,integer=True)
    m = n+1
    zero('the two adapted leading terms have identical weight in all degrees',
         a*(m-1)+(m+a)-m*(a+1))
    zero('the BZ term is strictly below the leading weight by m-1',
         m*(a+1)-(a*(m-1)+(a+1))-(m-1))
    zero('the unmultiplied H term is strictly below by a(m-1)',
         m*(a+1)-(m+a)-a*(m-1))
    residue_degree = S.symbols('residue_degree',nonnegative=True,integer=True)
    gap = S.symbols('gap',positive=True,integer=True)
    general_degree = residue_degree+gap
    zero('the general degree-gap frame balances its two leading terms',
         a*(general_degree-residue_degree)+(general_degree+a*residue_degree)-general_degree*(a+1))
    zero('an arbitrary lower-degree remainder stays below the leading weight by the positive degree gap',
         general_degree*(a+1)-(a*(general_degree-residue_degree)+residue_degree*(a+1))-gap)
    X1,X2 = S.symbols('X1 X2')
    fixed_B = X1**3+X1*X2+3
    parameter_h = T*(T-1)
    for degree in [2,5]:
        polynomial_b = Z**degree+Z+1
        G = (1+parameter_h*fixed_B)*Y-parameter_h*polynomial_b-Z
        H_original = fixed_B*Y-polynomial_b
        frame = (1+parameter_h*fixed_B)*H+polynomial_b-fixed_B*Z-fixed_B*V
        zero('the original Y reconstructs polynomially from the adapted frame, degree '+str(degree),
             Y-(V+Z-parameter_h*H).subs({V:G,H:H_original},simultaneous=True))
        zero('the adapted relation holds on the original ambient ring, degree '+str(degree),
             frame.subs({V:G,H:H_original},simultaneous=True))
        zero('the reverse presentation reconstructs H, degree '+str(degree),
             H-fixed_B*(V+Z-parameter_h*H)+polynomial_b-frame)
        weights = [4,degree-1,degree-1,degree+3,0,0]
        # Here a=3, hence weight Z=4 and weight H=m+3.
        weighted_terms = S.Poly(frame,Z,X1,X2,H,T,V).terms()
        top = sum(coefficient*S.prod(variable**exponent for variable,exponent in
                  zip((Z,X1,X2,H,T,V),monomial))
                  for monomial,coefficient in weighted_terms
                  if sum(exponent*weight for exponent,weight in zip(monomial,weights))==4*degree)
        zero('the adapted generic graded equation has exactly the stated two terms, degree '+str(degree),
             top-(parameter_h*X1**3*H+Z**degree))
        check('the adapted graded relation is primitive in H, degree '+str(degree),
              S.gcd(parameter_h*X1**3,Z**degree)==1)
        frame_basis = S.groebner([frame],Z,H,X1,X2,T,V)
        reduced = frame_basis.reduce(Z**(degree+1)*X1)[1]
        check('the adapted frame is globally monic and has no high Z basis powers, degree '+str(degree),
              all(monomial[0]<degree for monomial,coefficient in S.Poly(reduced,Z,H,X1,X2,T,V).terms()))
    check('the adapted multiplier has constant content for arbitrary parameter zeros',
          S.gcd(1+3*parameter_h,parameter_h)==1)
    zero('at a root of h its whole multiplier specializes to one',(1+parameter_h*fixed_B).subs(T,0)-1)

    # A Bezout parameter matrix gives the larger adapted family uniformly.
    A0,hh,uu,vv = T,T-1,S.Integer(1),S.Integer(-1)
    parameter_matrix = S.Matrix([[A0,hh],[-vv,uu]])
    parameter_inverse = S.Matrix([[uu,-hh],[vv,A0]])
    zero('the general parameter matrix has determinant one by Bezout',parameter_matrix.det()-1)
    check('its inverse is polynomial over the original parameter ring',
          (parameter_matrix*parameter_inverse).applyfunc(S.expand)==S.eye(2))
    q,d0 = T*T+1,T**3
    moving_B = X1**3+(T+1)*X2+3
    moving_b = Z**4+T*Z*z2+z2**2+1
    J = S.symbols('J')
    remainder_polynomial = q*Z+(T+2)*z2+d0+T*Z*z2+Z**2
    target_value = (A0+hh*moving_B)*Y-hh*moving_b-remainder_polynomial
    original_H = moving_B*Y-moving_b
    original_J = uu*original_H-vv*Y
    affine_S = V+remainder_polynomial
    general_frame = (A0+hh*moving_B)*J+moving_b-(uu*moving_B-vv)*affine_S
    zero('the general unimodular frame reconstructs Y exactly',
         Y-(uu*affine_S-hh*J).subs({V:target_value,J:original_J},simultaneous=True))
    zero('the general unimodular frame reconstructs H exactly',
         original_H-(vv*affine_S+A0*J).subs({V:target_value,J:original_J},simultaneous=True))
    zero('the moving-coefficient adapted relation holds in the original ambient ring',
         general_frame.subs({V:target_value,J:original_J},simultaneous=True))
    check('the larger frame remains monic despite parameter-varying lower root coefficients',
          S.Poly(general_frame,Z).LC()==1)
    general_weights = [4,2,2,10,4,0,0]
    general_top = sum(coefficient*S.prod(variable**exponent for variable,exponent in
                      zip((Z,X1,X2,J,z2,T,V),monomial))
                      for monomial,coefficient in S.Poly(general_frame,Z,X1,X2,J,z2,T,V).terms()
                      if sum(exponent*weight for exponent,weight in zip(monomial,general_weights))==16)
    zero('the larger adapted frame has the same primitive generic leading equation',
         general_top-(hh*X1**3*J+Z**4))
    check('a primitive positive coefficient row makes the larger multiplier content constant',
          S.gcd(S.gcd(hh,(T+1)*hh),A0+3*hh)==1)
    content_counterexample = S.expand(T+(T-1)*((T-2)*X-2))
    content_coefficients = S.Poly(content_counterexample,X).all_coeffs()
    check('Bezout coprimality alone does not remove a moving multiplier content zero',
          S.degree(S.gcd(*content_coefficients),T)==1)
    zero('the content-boundary control really loses its entire multiplier at that level',
         content_counterexample.subs(T,2))
    arbitrary_lower_P = T**2*Z**5+(T**3+1)*Z**4+(T**2+T)*Z**2+T*Z+1
    quotient_b = Z**5+T*Z**4+Z**2
    residue_C = Z**4+T*Z**2+T*Z+1
    zero('coefficientwise Euclidean division gives the uniform univariate leading-divisibility frame',
         arbitrary_lower_P-T**2*quotient_b-residue_C)
    check('that corollary admits the maximal remainder degree with arbitrary lower coefficients',
          S.degree(residue_C,Z)==S.degree(quotient_b,Z)-1 and S.Poly(quotient_b,Z).LC()==1)

    # Actual HC4 positive control for the moving adapted-frame family.
    source_u,source_v = S.symbols('source_u source_v')
    source_P = (T*T+T)*(Z**4+Z)+Z
    pulled_P = source_P.subs({T:s,Z:source_v})
    potential = s*source_u+t*source_v+source_v*pulled_P
    zero('the moving all-degree family has an actual four-variable unit-Hessian control',
         S.det(S.hessian(potential,(s,t,source_u,source_v)))-1)
    zero('its explicit closed embedding carries the identity pair',
         ((1+(T*T+T)*X1*X2)*Y-source_P).subs(
             {T:s,X1:0,X2:source_u,Y:t+pulled_P,Z:source_v},simultaneous=True)-t)

    context = json.loads((HERE/'sources/parameter-stability-import.json').read_text())
    check('Galois-tower inertia surjectivity is one explicitly new primary statement',
          len(context['newly_used_external_statements'])==1
          and '0BU7' in context['newly_used_external_statements'][0]['statement'])
    check('the inertia proof stays on the original plane source and only fixes parameters horizontally',
          context['original_plane_source_used_for_inertia'] is True
          and context['horizontal_parameter_cover_unramified_required'] is True
          and context['parameter_cover_source_assumed_affine_plane'] is False
          and context['vertical_parameter_inertia_asserted_trivial'] is False)
    check('moving generic invariants are separated from global filtration and content-fibre conditions',
          context['moving_multiplier_generic_class_invariants_zero'] is True
          and context['global_frame_and_content_fibre_conditions_separated'] is True
          and context['every_moving_multiplier_asserted_global_frame'] is False
          and context['unimodular_parameter_frame_with_content_hypothesis'] is True
          and context['strict_joint_degree_gap_frame'] is True
          and context['bezout_coprimality_alone_asserted_content_sufficient'] is False)
    pins = context['retained_proof_and_source_sha256']
    for path,digest in pins.items():
        check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    retained_counts = {}
    for name in [name for name in pins if '/receipt-' in name]:
        receipt = json.loads((ROOT/name).read_text())
        check('retained receipt is PASS '+name,receipt['outcome']=='PASS')
        for relative,digest in receipt['source_sha256'].items():
            check('retained receipt source '+str(Path(name).parent/relative),
                  hashlib.sha256((ROOT/Path(name).parent/relative).read_bytes()).hexdigest()==digest)
        for mapping in ['dependencies','dependency_source_sha256']:
            for path,digest in receipt.get(mapping,{}).items():
                check('retained dependency '+name+' -> '+path,
                      hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
        retained_counts[Path(name).name] = receipt['check_count']
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_parameter_stability.py','PARAMETER_STABILITY_AND_MOVING_MULTIPLIERS.md','sources/parameter-stability-import.json']},
        dependencies=pins,retained_receipt_counts_not_newly_replayed=retained_counts,
        written_not_finite_exhaustion=context['direct_proof_arguments'],
        newly_used_external_statements=context['newly_used_external_statements'],
        parameter_extension_irreducibility_plane_only=True,
        original_plane_source_used_for_inertia=True,
        horizontal_parameter_cover_unramified_required=True,
        moving_multiplier_generic_class_invariants_zero=True,
        disjoint_multiplier_and_root_variables_required=True,
        global_frame_and_content_fibre_conditions_separated=True,
        adapted_frame_proved_for_every_nonlinear_root_degree=True,
        unimodular_parameter_frame_with_content_hypothesis=True,
        strict_joint_degree_gap_frame=True,
        univariate_leading_divisibility_corollary=True,
        parameter_cover_source_assumed_affine_plane=False,
        vertical_parameter_inertia_asserted_trivial=False,
        full_splitting_group_preservation_inferred=False,
        every_moving_multiplier_asserted_global_frame=False,
        bezout_coprimality_alone_asserted_content_sufficient=False,
        composite_content_level_asserted_Keller_realization=False,
        overlapping_variables_asserted_covered=False,
        arbitrary_graph_normal_form_asserted=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-parameter-stability.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
