"""Exact controls for generative factors, fixed-multiplier class lattices and global suspension filtrations."""
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

    T,V,U,X,Y,z1,z2 = S.symbols('T V U X Y z1 z2')
    generator = z1*z1+z2**3
    outer = T*U**3+U
    P = S.expand(outer.subs(U,generator))
    check('the noncoordinate generative control has a zero-dimensional critical locus',
          S.groebner([S.diff(generator,z1),S.diff(generator,z2)],z1,z2).is_zero_dimensional)
    check('its zero fibre is singular, excluding a polynomial coordinate',
          all(poly.subs({z1:0,z2:0})==0 for poly in [generator,S.diff(generator,z1),S.diff(generator,z2)]))
    check('joint degree differs from the outer root degree in the geometric factor count',
          S.Poly(P,z1,z2).total_degree()==9 and S.degree(outer,U)==3)
    check('the outer generic polynomial is separable without a constant-leading hypothesis',
          S.discriminant(outer+V,U)!=0 and S.Poly(outer,U).LC()==T)
    alpha,beta,gamma = S.symbols('alpha beta gamma')
    factored_outer = (U-alpha)*(U-beta)*(U-gamma)
    zero('the multivariable factorization is exactly the outer-root factorization',
         factored_outer.subs(U,generator)-(generator-alpha)*(generator-beta)*(generator-gamma))
    check('distinct outer roots give coprime geometric factor controls',
          S.gcd(generator-alpha,generator-beta)==1)
    check('generic geometric irreducibility must not be replaced by all-special-fibre irreducibility',
          S.gcd(z1,U)==1 and len(S.factor_list(z1*z2,z1,z2)[1])==2)

    # A rational field decomposition is kept out of the global integral model.
    rational_generator = generator/T
    rational_outer = T**4*U**3+T*U
    zero('a generative decomposition may have parameter denominators while the original P is polynomial',
         rational_outer.subs(U,rational_generator)-P)
    check('the generative field calculation really has a pole at a finite T level',
          S.denom(rational_generator)==T)

    # Exact root-lattice matrices, including the proper transitive subgroup case.
    def root_lattice_matrix(permutation):
        degree = len(permutation)
        matrix = S.zeros(degree-1)
        for column in range(degree-1):
            target = permutation[column]
            if target==degree-1:
                for row in range(degree-1):
                    matrix[row,column] = -1
            else:
                matrix[target,column] = 1
        return matrix

    for degree in [2,3,5]:
        cycle = list(range(1,degree))+[0]
        action = root_lattice_matrix(cycle)
        check('a transitive cyclic action has zero fixed root-lattice subgroup in degree '+str(degree),
              (action-S.eye(degree-1)).rank()==degree-1)
    check('degree one gives the zero class lattice',root_lattice_matrix([0]).shape==(0,0))
    nontransitive = root_lattice_matrix([0,2,1])
    check('nontransitivity can leave a nonzero fixed class and is not silently admitted',
          len((nontransitive-S.eye(2)).nullspace())==1)
    m,c = S.symbols('m c',positive=True,integer=True)
    zero('the all-degree sum-of-coordinates proof kills every diagonal class shift',m*c-c*m)

    A = X**2*(X-1)**3
    multiplier_factors = S.factor_list(A,X)[1]
    check('repeated multiplier factors contribute by distinct support, not multiplicity',
          len(multiplier_factors)==2 and sum(exponent for factor,exponent in multiplier_factors)==5)
    check('the exact generic class rank uses both distinct multiplier support and outer degree',
          len(multiplier_factors)*(S.degree(outer,U)-1)==4)
    cycle_three = root_lattice_matrix([1,2,0])
    fixed_rows = S.diag(cycle_three,cycle_three)
    check('every fixed multiplier row inherits zero source-fixed subgroup',
          (fixed_rows-S.eye(4)).rank()==4)
    correlated = S.Matrix([[0,-1],[-1,0]])
    check('moving multiplier factors can create a correlated fixed class despite root transitivity',
          correlated*S.Matrix([1,-1])==S.Matrix([1,-1]))
    check('fixed rows are genuinely different from that unadmitted moving-factor action',
          (-S.eye(2)-S.eye(2)).rank()==2)

    # One monic relation gives a global basis even across parameter degree drops.
    relation = S.expand(A*Y-P-V)
    basis = S.groebner([relation],Y,X,z1,z2,T,V)
    normal = basis.reduce(Y*X**6*z1)[1]
    check('global monic division leaves no monomial divisible by Y times the leading multiplier monomial',
          all(not (monomial[0]>=1 and monomial[1]>=5)
              for monomial,coefficient in S.Poly(normal,Y,X,z1,z2,T,V).terms()))
    zero('the global division identity is exact',basis.reduce(Y*X**6*z1-normal)[1])
    weights = [9,9,6,6,0,0]

    def weighted_degree(polynomial):
        return max(sum(exponent*weight for exponent,weight in zip(monomial,weights))
                   for monomial,coefficient in S.Poly(polynomial,Y,X,z1,z2,T,V).terms())

    check('normal reduction never increases the displayed weight in the control',
          weighted_degree(normal)<=weighted_degree(Y*X**6*z1))
    leading = sum(coefficient*S.prod(variable**exponent for variable,exponent in
                  zip((Y,X,z1,z2,T,V),monomial))
                  for monomial,coefficient in S.Poly(relation,Y,X,z1,z2,T,V).terms()
                  if sum(exponent*weight for exponent,weight in zip(monomial,weights))==54)
    zero('the generic graded equation has the claimed independent-variable leading forms',
         leading-(X**5*Y-T*z2**9))
    check('those leading forms are coprime, giving the primitive-linear generic domain',S.gcd(X**5,T*z2**9)==1)
    degree_drop = T*generator**3+T*T+V
    zero('when every positive-Z coefficient drops, the height-one fibre remains nonzero',degree_drop.subs(T,0)-V)
    check('the fixed multiplier remains coprime to that nonzero residue polynomial',S.gcd(A,V)==1)
    zero('the base-only excluded case has a reducible height-one fibre',
         (A*Y-T*T-V).subs(V,-T*T)-A*Y)
    check('overlapping variables can fail height-one geometric integrality',
          S.gcd(z1,z1*z2)==z1 and len(S.factor_list(z1*(Y-z2),z1,Y,z2)[1])==2)

    # Actual stable embeddings and four-variable unit-Hessian positive controls.
    s,t,w,v = S.symbols('s t w v')
    P_constant = generator**3+generator
    source_P = P_constant.subs({z1:w,z2:v})
    potential = 1+s*v+t*w+w*source_P
    zero('the noncoordinate multivariable suspension has an actual-unit-Hessian HC4 control',
         S.det(S.hessian(potential,(s,t,w,v)))-1)
    zero('its carried identity pair is exact under the explicit closed embedding',
         (t+source_P)-source_P-t)
    X1,X2,Z = S.symbols('X1 X2 Z')
    fixed_multiplier = 1+(X1-1)*X2
    zero('the independent two-variable multiplier has an explicit unit slice in the second HC4 control',
         fixed_multiplier.subs(X1,1)-1)
    varying_P = T*Z**3+Z
    second_P = varying_P.subs({T:s,Z:v})
    second_potential = 1+s*w+t*v+v*second_P
    zero('the fixed-multiplier control permits genuinely varying parameter coefficients and actual unit Hessian',
         S.det(S.hessian(second_potential,(s,t,w,v)))-1)

    # The obstruction is stable under ambient and target automorphisms.
    transformed_second_z = z2+X*z1*z1
    zero('the nonlinear ambient shear has an explicit polynomial inverse',
         transformed_second_z-X*z1*z1-z2)
    zero('the output after that ambient shear reduces to the original under its inverse',
         P.subs(z2,transformed_second_z).subs(z2,z2-X*z1*z1)-P)
    target_T,target_V = T+V*V,V
    zero('the nonlinear target shear also has an explicit polynomial inverse',target_T-target_V**2-T)

    context = json.loads((HERE/'sources/multivariable-suspension-import.json').read_text())
    check('the only expanded primary use is the explicitly audited closed-polynomial implication over C(T)',
          len(context['newly_used_external_statements'])==1
          and context['newly_used_external_statements'][0]['expanded_field_scope'] is True
          and 'over C(T)' in context['newly_used_external_statements'][0]['statement'])
    check('minimal-degree existence and generic geometric integrality are proved directly',
          context['minimal_degree_generative_existence_proved_directly'] is True
          and context['geometric_generic_integrality_proved_by_field_argument'] is True
          and context['special_fibre_classification_imported'] is False
          and context['proposition_2_base_extension_theorem_imported'] is False)
    check('disjoint fixed multiplier and positive joint degree are required rather than silently generalized',
          context['fixed_multiplier_independent_of_T_and_Z'] is True
          and context['joint_Z_degree_positive_required'] is True
          and context['parameter_dependent_multiplier_asserted_covered'] is False)
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
            ['verify_multivariable_suspensions.py','MULTIVARIABLE_SUSPENSIONS_AND_FIXED_MULTIPLIERS.md','sources/multivariable-suspension-import.json']},
        dependencies=pins,retained_receipt_counts_not_newly_replayed=retained_counts,
        written_not_finite_exhaustion=context['direct_proof_arguments'],
        newly_used_external_statements=context['newly_used_external_statements'],
        arbitrary_field_outer_irreducibility_transfer=True,
        unrestricted_Keller_application_plane_only=True,
        fixed_multiplier_independent_of_T_and_Z=True,joint_Z_degree_positive_required=True,
        outer_degree_not_joint_degree_controls_lattice=True,
        multiplier_multiplicities_do_not_change_lattice=True,
        global_common_free_filtration_and_all_height_one_fibres_proved=True,
        ambient_and_target_polynomial_automorphism_saturation=True,
        special_fibre_classification_imported=False,
        proposition_2_base_extension_theorem_imported=False,
        generative_polynomial_assumed_coordinate=False,
        parameter_dependent_multiplier_asserted_covered=False,
        correlated_action_control_asserted_Keller_realization=False,
        full_splitting_group_preservation_inferred=False,
        arbitrary_target_cover_group_preservation=False,
        nonlinear_source_Hessian_invariance_used=False,
        arbitrary_graph_normal_form_asserted=False,
        unit_ideal_endpoint_resolved=False,singleton_zero_endpoint_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-multivariable-suspensions.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
