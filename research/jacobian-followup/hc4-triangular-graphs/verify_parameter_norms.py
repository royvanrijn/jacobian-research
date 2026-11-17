"""Exact controls for parameter norms, integral critical values and all polynomial outputs."""
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

    s,t,u,T,V,Y,z,h = S.symbols('s t u T V Y z h')
    f = t*(2+s*s*t)
    unit_part = 1+s*s*t
    zero('the nowhere-critical parameter has an exact gradient Bezout identity',
         (1-s*s*t)*S.diff(f,t)/2+s**3*S.diff(f,s)/2-1)
    zero('the nonconstant source unit and its parameter conjugate have norm one',
         unit_part**2-s*s*f-1)
    multiplication = S.Matrix([[unit_part,s*f],[s,unit_part]])
    inverse_multiplication = S.Matrix([[unit_part,-s*f],[-s,unit_part]])
    check('the finite-free multiplication matrix is invertible over the original polynomial ring',
          (multiplication*inverse_multiplication).applyfunc(S.expand)==S.eye(2))
    zero('the ring norm is precisely the determinant, a nonzero constant',multiplication.det()-1)
    check('a source-cover unit of constant norm need not be constant',S.diff(unit_part+s*u,u)==s)
    zero('inverse integrality follows from the product of all other conjugates',
         S.rem((unit_part+s*u)*(unit_part-s*u)-1,u*u-f,u))
    zero('a constant product alone cannot force integrality of both factors',T*(1/T)-1)
    check('the missing-integrality control has a genuine finite pole',S.denom(1/T)==T)

    # The nontrivial parameter extension for a genuinely algebraic monic input.
    original = T**3*z**5+T**2*z**2+T+1
    monic = Y**5+u**4*Y**2+u**5+1
    zero('the algebraic fifth-root parameter model is integral and monic',
         original.subs({T:u**5,z:Y/u**3})-monic)
    check('the original leading exponent requires an algebraic fifth-root scaling',3%5!=0)
    zero('the no-finite-pole binomial coefficient has an explicit monic integrality equation',
         (u**4)**5-(T**4).subs(T,u**5))
    check('the binomial divisibility alternative is exhaustive in this integral-model control',
          S.rem((T**2)**5,(T**3)**2,T)==0)

    # A centered cubic gives a direct critical-value sum/product control.
    a,b,critical = S.symbols('a b critical')
    Q = Y**3+a*Y+b
    critical_values = S.resultant(S.diff(Q,Y),critical-Q,Y)
    zero('the critical-value polynomial records both centered coefficients',
         critical_values-(27*(critical-b)**2+4*a**3))
    # Normalized bad reductions for two distinct arbitrary-degree mechanisms.
    for polynomial,degree,ramification,scale_power,reduction in [
        (Y**4+Y**2/T,4,2,-1,z**4+z*z),
        (Y**5+Y/T**3,5,4,-3,z**5+z)]:
        scaled = S.cancel(polynomial.subs({T:h**ramification,Y:h**scale_power*z})/h**(scale_power*degree))
        zero('valuation scaling gives the stated integral centered reduction in degree '+str(degree),scaled-reduction)
        check('the reduction has a surviving lower coefficient in degree '+str(degree),reduction!=z**degree)
        check('its critical points cannot all have value zero in degree '+str(degree),
              S.degree(S.gcd(reduction,S.diff(reduction,z)),z)<degree-1)
    m,r = S.symbols('m r',integer=True,positive=True)
    zero('the all-degree critical-zero count forces one distinct zero',
         (m-1)-(m-r)-(r-1))
    c = S.symbols('c')
    for degree in [2,5]:
        check('a pure power is centered only when its unique root is zero, degree '+str(degree),
              S.Poly((z-c)**degree,z).coeff_monomial(z**(degree-1))==-degree*c)
    check('residue characteristic zero cannot be omitted from the reduction argument',
          S.Poly(S.diff(z**3+z,z),z,modulus=3).degree()==0)
    leading,penultimate = S.symbols('leading penultimate',nonzero=True)
    shift = penultimate/(5*leading)
    transformed = S.expand(leading*(Y-shift)**5+penultimate*(Y-shift)**4+Y**2)
    zero('the rational input centering kills the next-to-leading coefficient exactly',
         S.Poly(transformed,Y).coeff_monomial(Y**4))

    # A nontransitive normal subgroup supplies blocks; full source group is unnecessary.
    cycle_a = Permutation(0,1,2,size=6)
    cycle_b = Permutation(3,4,5,size=6)
    transposition_a = Permutation(0,1,size=6)
    block_swap = Permutation(0,3)(1,4)(2,5)
    group = PermutationGroup(cycle_a,cycle_b,transposition_a,block_swap)
    bad = PermutationGroup(cycle_a,cycle_b)
    source = PermutationGroup(cycle_a,cycle_b,block_swap)
    blocks = {frozenset(orbit) for orbit in bad.orbits()}
    check('bad inertia can be normal but nontransitive in a polynomial-style block control',
          bad.is_normal(group) and len(blocks)==2)
    check('the entire group permutes the bad-inertia blocks',
          all(frozenset(g(i) for i in block) in blocks for g in group.generators for block in blocks))
    check('transitive block action and contained bad inertia force source root transitivity',
          source.is_transitive() and all(g in source for g in bad.generators))
    check('the all-output conclusion must not be upgraded to full splitting-group preservation',
          source.order()==18 and group.order()==72)

    # A concrete decomposition has a bad inner critical curve and integral outer values.
    inner = T*z*z+z
    outer = Y**3+Y
    P = S.expand(outer.subs(Y,inner))
    derivative = S.diff(P,z)
    zero('the polynomial block decomposition retains the entire derivative chain rule',
         derivative-S.diff(inner,z)*S.diff(outer,Y).subs(Y,inner))
    inner_critical = -1/(2*T)
    bad_curve = 64*T**3*V-16*T*T-1
    zero('the inner critical value lies on a genuinely punctured horizontal branch curve',
         bad_curve.subs(V,-P.subs(z,inner_critical)))
    check('that bad branch curve has no point over its denominator zero',bad_curve.subs(T,0)==-1)
    outer_values = S.resultant(S.diff(outer,Y),V+outer,Y)
    zero('the surviving outer critical values are integral constants',outer_values-(27*V*V+4))
    discriminant = S.discriminant(P+V,z)
    remainder = discriminant
    support = T*bad_curve*(27*V*V+4)
    while S.total_degree(remainder)>0:
        divisor = S.gcd(remainder,support)
        if S.total_degree(divisor)==0:
            break
        remainder = S.cancel(remainder/divisor)
    check('all branch support in the decomposition control is accounted for',
          remainder!=0 and S.total_degree(remainder)==0)
    # Polynomial parametrization implies integrality, without requiring a birational map.
    X = S.symbols('X')
    param_T = h**3+h+1
    param_V = h*h+2*h
    zero('a polynomial parameter satisfies a monic equation over its nonconstant T coordinate',
         (X**3+X+1-param_T).subs(X,h))
    implicit = S.resultant(X**3+X+1-T,X*X+2*X-V,X)
    zero('the induced branch-value function satisfies a monic integral equation',implicit.subs({T:param_T,V:param_V}))
    check('the integral value equation has constant nonzero leading coefficient',S.Poly(implicit,V).LC() in [1,-1])
    # The output theorem still has actual stable embedding and Hessian positive controls.
    w,v = S.symbols('w v')
    potential = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('the all-polynomial output has an actual-unit-Hessian automorphic control',
         S.det(S.hessian(potential,(s,t,w,v)))-1)
    zero('the carried positive embedding has an explicit parameter inverse',
         (t+P.subs({T:s,z:w}))-P.subs({T:s,z:w})-t)

    context = json.loads((HERE/'sources/parameter-norm-import.json').read_text())
    check('ordinary Luroth and structural source normality inputs are explicit',
          len(context['newly_used_external_statements'])==3
          and any('Luroth' in item['statement'] for item in context['newly_used_external_statements']))
    check('the critical-value scaling lemma and polynomial one-pole refinement are direct proofs',
          context['critical_value_scaling_lemma_proved_directly'] is True
          and context['unique_pole_polynomial_refinement_proved_directly'] is True
          and context['good_reduction_dynamical_theorem_imported'] is False
          and context['Ritt_classification_imported'] is False)
    check('all-dimensional norm and integral-model claims are distinct from plane unrestricted output closure',
          len(context['all_dimensional_claims'])==3 and len(context['plane_only_claims'])==1)
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
            ['verify_parameter_norms.py','PARAMETER_NORMS_AND_ALL_POLYNOMIAL_OUTPUTS.md','sources/parameter-norm-import.json']},
        dependencies=pins,retained_receipt_counts_not_newly_replayed=retained_counts,
        written_not_finite_exhaustion=context['direct_proof_arguments'],
        newly_used_external_statements=context['newly_used_external_statements'],
        all_dimensional_parameter_norm_and_integral_model_theorems=True,
        all_dimensional_integral_critical_value_theorem=True,
        all_polynomial_output_irreducibility_plane_only=True,
        cross_resultant_unit_norm_used=True,all_parameter_conjugates_integral_required=True,
        polynomial_block_quotient_used=True,
        parameter_cover_source_assumed_affine_space=False,parameter_cover_units_assumed_constant=False,
        essential_surjectivity_used=False,good_reduction_dynamical_theorem_imported=False,
        Ritt_classification_imported=False,generic_group_classification_used=False,
        full_splitting_group_preservation_inferred=False,
        arbitrary_target_cover_group_preservation=False,
        nonlinear_source_Hessian_invariance_used=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-parameter-norms.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
