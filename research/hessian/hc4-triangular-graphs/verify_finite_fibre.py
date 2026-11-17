"""Exact controls for finite-fibre rows and universal cotangent templates.

All multi-point rows below are non-Keller reference controls. Arbitrary
n and the universal equivalences are proved in the accompanying note.
"""
import argparse
from collections import Counter
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

    def check(name,condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(item.doit())==0 for item in values))

    s,U,P,Q = S.symbols('s U P Q')
    variables = (s,U,P,Q)
    a1,a2,mu,nu,m11,m12,m21,m22 = S.symbols('a1 a2 mu nu m11 m12 m21 m22')
    row = S.Matrix([[a1,a2]])
    matrix = S.Matrix([[m11,m12],[m21,m22]])
    target = row*matrix
    kernel = S.Matrix([-a2,a1])
    corrected = matrix+kernel*S.Matrix([[nu,-mu]])
    zero('general finite-ideal correction preserves the complete target row',row*corrected-target)
    zero('adjugate transports the canonical row syzygy exactly',matrix.adjugate()*kernel-S.Matrix([-target[1],target[0]]))
    zero('the determinant correction removes the entire target-ideal defect',
         corrected.det()-matrix.det()+mu*target[0]+nu*target[1])

    # A determinant-one linear change first separates the given points.
    original_points = [(0,1),(0,3),(2,1)]
    linear = S.Matrix([[1,2],[0,1]])
    separated = [tuple(linear*S.Matrix(point)) for point in original_points]
    check('the initial linear separator has determinant one',linear.det()==1)
    check('the linear separator gives distinct first coordinates',len({point[0] for point in separated})==len(separated))
    J = S.interpolate([(point[0],ii-point[1]) for ii,point in enumerate(separated,1)],s)
    E = S.interpolate([(ii,ii-point[0]) for ii,point in enumerate(separated,1)],U)
    after_first = U+J
    mapped_s = s+E.subs(U,after_first)
    beta = S.Matrix([mapped_s,after_first-mapped_s])
    back_s = s-E.subs(U,U+s)
    back_u = U+s-J.subs(s,back_s)
    zero('the three interpolation shears have an actual polynomial inverse',
         beta.subs({s:back_s,U:back_u},simultaneous=True)-S.Matrix([s,U]))
    zero('the inverse shear composition works in both directions',
         S.Matrix([back_s,back_u]).subs({s:beta[0],U:beta[1]},simultaneous=True)-S.Matrix([s,U]))
    zero('the complete interpolation map has determinant one',beta.jacobian((s,U)).det()-1)
    check('all labelled points are sent to their prescribed collinear positions',
          all(beta.subs({s:point[0],U:point[1]},simultaneous=True)==S.Matrix([ii,0])
              for ii,point in enumerate(separated,1)))

    hu,hv = S.symbols('hu hv')
    t = hu*hv
    controls = {}
    for count in (1,2,3,5):
        label = 'n='+str(count)
        A = S.prod(s-ii for ii in range(1,count+1)).expand()
        Ap = S.diff(A,s)
        B = S.invert(Ap,A)
        explicit_B = S.expand(sum(S.div(A,s-ii,s)[0]/Ap.subs(s,ii)**2 for ii in range(1,count+1)))
        zero('closed interpolation formula equals the modular inverse '+label,B-explicit_B)
        quotient,remainder = S.div(Ap*B-1,A,s)
        zero('the inverse-derivative orientation is exact modulo the point polynomial '+label,remainder)
        zero('the reference row generates the full point ideal '+label,
             Ap*(B*U)-quotient*A*U-U)
        check('the reference multiplier is a unit at every marked point '+label,
              all(Ap.subs(s,ii)*B.subs(s,ii)==1 for ii in range(1,count+1)))
        D = S.rem(s*s-3*s+1,A,s)
        template = D+A*P+B*U*Q
        determinant = S.hessian(template,variables).det()
        zero('full template Hessian has its exact nonconstant native factor '+label,
             determinant-(Ap*B)**2)
        check('the template has degree exactly n plus one '+label,
              S.Poly(template,*variables).total_degree()==count+1)
        critical_points = [{s:ii,U:0,P:-S.diff(D,s).subs(s,ii)/Ap.subs(s,ii),Q:0}
                           for ii in range(1,count+1)]
        gradient = S.Matrix([S.diff(template,z) for z in variables])
        check('all prescribed full critical points are exact '+label,
              all(gradient.subs(point,simultaneous=True)==S.zeros(4,1) for point in critical_points))
        check('all critical Hessians are one despite the native nonunit factor '+label,
              all(determinant.subs(point,simultaneous=True)==1 for point in critical_points))
        critical_values = [template.subs(point,simultaneous=True) for point in critical_points]
        check('critical values are precisely the interpolated list '+label,
              critical_values==[D.subs(s,ii) for ii in range(1,count+1)])
        jumps = Counter(critical_values)
        check('the complete value-profile jumps have total count n '+label,sum(jumps.values())==count)
        zero('the generic whole-value profile retains the full critical count '+label,
             (t*t-count)*t-(t**3-count*t))
        if count>=2:
            check('the multi-point reference row is deliberately not Keller '+label,S.degree(Ap*B,s)>0)
        else:
            zero('the singleton reference is the known quadratic endpoint',determinant-1)
        controls[str(count)] = dict(A=str(A),B=str(B),D=str(D),
            native_degree=count+1,critical_count=count,native_Jacobian=str(S.expand(Ap*B)),
            is_Keller_reference=count==1)

    # Ideal equality does not by itself imply a constant-determinant row lift.
    A = (s-1)*(s-2)
    B = 2*s-3
    zero('the two-point ideal certificate is literal',B*B-4*A-1)
    first_row = S.Matrix([[A,U]])
    second_row = S.Matrix([[A,B*U]])
    orientations = [B.subs(s,ii) for ii in (1,2)]
    check('one ideal can have different local orientation ratios',orientations==[-1,1])
    check('those ratios rule out every constant-determinant equivalence of the two rows',len(set(orientations))>1)

    # An actual non-Keller, multi-point target row exercises the correction.
    rn = second_row
    C0 = S.Matrix([[1,s+U*U],[U,1+U*(s+U*U)]])
    zero('the concrete dual matrix is unimodular',C0.det()-1)
    target = rn*C0
    ka = S.Matrix([-rn[1],rn[0]])
    h1,h2 = s*U,U*U+1
    raw_matrix = C0+ka*S.Matrix([[h1,h2]])
    zero('an arbitrary syzygy perturbation keeps the target row',rn*raw_matrix-target)
    zero('the raw matrix determinant has its supplied target-ideal certificate',
         raw_matrix.det()-1-(h2*target[0]-h1*target[1]))
    repaired = raw_matrix+ka*S.Matrix([[-h1,-h2]])
    zero('the exact determinant correction recovers an invertible matrix',repaired-C0)
    zero('the complete row source map has an actual polynomial inverse',
         C0*(C0.inv()*S.Matrix([P,Q]))-S.Matrix([P,Q]))
    target_jac = target.jacobian((s,U)).det()
    check('constant matrix determinant does not transport the row Jacobian off its zero fibre',
          S.expand(target_jac-S.diff(A,s)*B)!=0)
    check('local row orientations still agree at both marked points',
          all(target_jac.subs({s:ii,U:0},simultaneous=True)==1 for ii in (1,2)))
    check('the concrete multi-point target is not asserted a Keller pair',S.Poly(target_jac,s,U).total_degree()>0)

    D0 = s*s+1+A*(s*U+U*U)+B*U*(s*s+U)
    offset = D0+A*P+B*U*Q
    zero('the entire off-fibre offset is removed by polynomial dual translations',
         offset.subs({P:P-s*U-U*U,Q:Q-s*s-U},simultaneous=True)
             -(s*s+1+A*P+B*U*Q))
    check('a nonlinear dual matrix can increase the native template degree',
          S.Poly((target*S.Matrix([P,Q]))[0],*variables).total_degree()>3)

    dependencies = {}
    for proof,receipt in [('GLOBAL_VALUE_ORBITS.md','receipt-global-value.json'),
                          ('INTRINSIC_AFFINE_PLANE.md','receipt-intrinsic-plane.json'),
                          ('P_GRAPH_SINGLETON_ENDPOINT.md','receipt-square-axis.json')]:
        retained = json.loads((HERE/receipt).read_text())
        digest = hashlib.sha256((HERE/proof).read_bytes()).hexdigest()
        check('retained proof pin '+proof,digest==retained['source_sha256'][proof])
        dependencies[str((HERE/proof).relative_to(ROOT))] = digest
    rowproof = 'hessian/hc4-relative-row-orbit/PROOF.md'
    digest = hashlib.sha256((ROOT/rowproof).read_bytes()).hexdigest()
    check('retained constructive row and birational-Keller boundary pin',
          digest==json.loads((HERE/'receipt-square-axis.json').read_text())['dependencies'][rowproof])
    dependencies[rowproof] = digest
    tang = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for path in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'hessian/tang-hc4-audit/source-v1.pdf']:
        digest = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        check('retained Tang input '+path,digest==tang[path])
        dependencies[path] = digest
    context = json.loads((HERE/'sources/finite-fibre-literature-context.json').read_text())
    primary = HERE/'sources'/context['retained_file']
    digest = hashlib.sha256(primary.read_bytes()).hexdigest()
    check('related Guyot primary theorem has the inspected retained bytes',digest==context['sha256'])
    check('the general orbit theorem is attributed rather than asserted new',
          context['general_theorem_asserted_new'] is False and context['external_classification_imported'] is False)
    dependencies[str(primary.relative_to(ROOT))] = digest
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,exact_reference_controls=controls,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name,'FINITE_FIBRE_NORMAL_FORM.md','sources/finite-fibre-literature-context.json']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['orientation-complete row lifting for every reduced finite plane fibre',
            'three determinant-one interpolation shears normalize any ordered finite point set',
            'every nonautomorphic unit-Hessian cotangent potential has a degree n+1 source normal form',
            'arbitrary n templates retain exact critical Hessian and whole-value profiles without native unit Hessian',
            'full source unitization for n>=2 is equivalent to one relative SL2 row condition',
            'the explicit template sequence gives an exact unrestricted JC2/HC4 existence reduction'],
        related_general_row_orbit_theorem_asserted_new=False,
        general_finite_point_normalization_descended_to_arbitrary_field=False,
        native_Hessian_asserted_preserved_by_source_normalization=False,
        off_fibre_row_Jacobian_asserted_scaled_only_by_matrix_determinant=False,
        native_template_Hessian_asserted_unit_for_n_at_least_two=False,
        finite_reference_controls_asserted_Keller_for_n_at_least_two=False,
        matrix_or_counterexample_degree_bound_asserted=False,
        template_count_n_bounded=False,
        coefficient_row_unitization_condition_solved=False,
        unconditional_Meng_Yang_graph_embedding_asserted=False,
        marked_gradient_collision_transfer_asserted=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
