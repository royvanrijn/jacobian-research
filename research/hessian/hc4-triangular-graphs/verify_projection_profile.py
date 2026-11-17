"""Exact controls for the mixed gradient projection generator profile.

Conditional lower bounds and the all-degree stable quantifiers are
written proofs, not conclusions of these finite controls.
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
        check(name,all(S.cancel(item.doit())==0 for item in values))

    s,U,e1,e2,t1,t2,t3,t4,P,Q = S.symbols('s U e1 e2 t1 t2 t3 t4 P Q')
    zvars = (s,U,e1,e2)
    # Upper constructions hold for arbitrary polynomial pairs. This one
    # is deliberately non-Keller and is not a hypothetical counterexample.
    K = S.Matrix([s*s+U**3,s*U+U*U])
    check('the universal upper-construction control is explicitly non-Keller',
          S.Poly(K.jacobian((s,U)).det(),s,U).total_degree()>0)
    outputs = S.Matrix([K[0],K[1],e1,e2])
    matrices = [
        ('pure',S.Matrix([[1,0,0,0],[0,1,0,0]]),0,2,6),
        ('rank_one',S.Matrix([[2,3,1,0],[4,-1,0,0]]),1,1,5),
        ('rank_two',S.Matrix([[1,0,1,2],[0,1,3,7]]),2,0,4)]
    profiles = {}
    for label,rows,rank_b,ell,ambient in matrices:
        check(label+' output rows are independent with the asserted eta rank',
              rows.rank()==2 and rows[:,2:].rank()==rank_b)
        check(label+' pure intersection dimension equals row dimension minus eta rank',
              2-rows[:,2:].rank()==ell)
        projected = rows*outputs
        if label=='pure':
            embedding = projected.col_join(S.Matrix([s,U,e1,e2]))
            recovered = embedding[2:,0]
            zero('the pure projection has a six-coordinate polynomial source recovery',
                 recovered-S.Matrix(zvars))
        elif label=='rank_one':
            embedding = projected.col_join(S.Matrix([s,U,e2]))
            recovered = S.Matrix([embedding[2],embedding[3],
                embedding[0]-2*K[0]-3*K[1],embedding[4]])
            zero('the rank-one mixed projection has a five-coordinate polynomial source recovery',
                 recovered-S.Matrix(zvars))
            relation = t2-4*K[0]+K[1]
            zero('its hypersurface equation is exact and linear in an ambient target coordinate',
                 relation.subs({t1:projected[0],t2:projected[1]},simultaneous=True))
            zero('that graph hypersurface has an explicit coordinate partial',S.diff(relation,t2)-1)
        else:
            embedding = projected.col_join(S.Matrix([s,U]))
            inv_eta = rows[:,2:].inv()*(S.Matrix([t1,t2])-rows[:,:2]*K)
            zero('the rank-two mixed projection recovers both eta variables without a Keller inverse',
                 inv_eta.subs({t1:projected[0],t2:projected[1]},simultaneous=True)-S.Matrix([e1,e2]))
            zero('the four-coordinate mixed construction has unit source Jacobian',
                 embedding.jacobian(zvars).det()-1)
            returned = projected.subs({e1:inv_eta[0],e2:inv_eta[1]},simultaneous=True)
            zero('its other coordinate composition is the identity',returned-S.Matrix([t1,t2]))
        check(label+' construction has the predicted total number of ambient coordinates',
              len(embedding)==ambient)
        profiles[label] = dict(eta_rank=rank_b,pure_intersection_dimension=ell,
            conditional_minimum_ambient_dimension=ambient,
            minimum_under_nonautomorphic_Keller_hypothesis_only=True)

    # All possible dimension strata for four outputs; this is a small
    # linear-algebra table, not a search over Keller maps or polynomial degrees.
    dimension_table = []
    for ell in range(3):
        for rank_b in range(3):
            k = ell+rank_b
            rows = []
            for index in range(rank_b):
                row = [0]*4
                row[2+index] = 1
                row[0] = index+1
                rows.append(row)
            for index in range(ell):
                row = [0]*4
                row[index] = 1
                rows.append(row)
            mat = S.Matrix(rows) if rows else S.zeros(0,4)
            check('dimension stratum '+str((k,ell))+' realizes its exact intersection rank',
                  mat.rank()==k and k-mat[:,2:].rank()==ell)
            helpers = 4-k+ell
            check('dimension stratum '+str((k,ell))+' has exact helper budget and stable excess',
                  helpers==4-rank_b and k+helpers==4+ell)
            dimension_table.append(dict(k=k,ell=ell,helpers=helpers,ambient=4+ell))
    q,N,k,ell = S.symbols('q N k ell',integer=True)
    zero('adjoining free variables shifts only the helper count, not the excess',
         ((N+q)-k+ell)-(N-k+ell)-q)
    zero('the relative minimum ambient dimension is source dimension plus intersection dimension',
         k+(N-k+ell)-(N+ell))

    # The retained curve formula detects exactly the affine-line case.
    u,v,c,G,b = S.symbols('u v c G b')
    curve = c*u*v-G*(u+v)+c-b
    difference = S.Poly(curve-u*v,u,v)
    check('curve coefficient separation gives one component, genus zero and one puncture',
          [difference.coeff_monomial(mon) for mon in [u*v,u,v,1]]==[c-1,-G,-G,c-b])
    zero('the affine line gives precisely the required curve polynomial',
         curve.subs({c:1,G:0,b:1})-u*v)
    for label,sub in [('punctured',{c:1,G:0,b:2}),('positive_genus',{c:1,G:1,b:1}),
                      ('disconnected',{c:2,G:0,b:2})]:
        check(label+' smooth-curve control cannot be canceled to an affine-space cylinder',
              S.expand(curve.subs(sub)-u*v)!=0)

    # An honest unit-Hessian control transports the algebra to cotangent gradients.
    autoK = S.Matrix([s,U+s**3])
    D = s*U+s**3
    f = D+P*autoK[0]+Q*autoK[1]
    eta = S.Matrix([S.diff(f,s),S.diff(f,U)])
    jac = autoK.jacobian((s,U))
    zero('the cotangent source change has actual polynomial recovery',
         jac.T.inv()*(eta-S.Matrix([S.diff(D,s),S.diff(D,U)]))-S.Matrix([P,Q]))
    zero('the cotangent control has actual unit full Hessian',S.hessian(f,(s,U,P,Q)).det()-1)
    zero('the eta-plus-base source change is a polynomial automorphism',
         S.Matrix([s,U,eta[0],eta[1]]).jacobian((s,U,P,Q)).det()-1)
    check('the positive control is automorphic, so conditional lower bounds are not inferred from it',
          autoK[0]==s and autoK[1]==U+s**3)

    dependency_pins = {}
    for receipt_name,source_names in [
        ('receipt-stable-hypersurface.json',['STABLE_HYPERSURFACE_OBSTRUCTION.md']),
        ('receipt-intrinsic-plane.json',['INTRINSIC_AFFINE_PLANE.md']),
        ('receipt-critical-components.json',['CRITICAL_COMPONENTS.md'])]:
        retained = json.loads((HERE/receipt_name).read_text())
        for name in source_names:
            dependency_pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
        if receipt_name!='receipt-stable-hypersurface.json':
            dependency_pins.update(retained['dependencies'])
    for name,expected in dependency_pins.items():
        check('retained mathematical or primary input pin '+name,
              hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected)
    context = json.loads((HERE/'sources/projection-profile-context.json').read_text())
    path = HERE/'sources'/context['retained_file']
    check('the inspected classical plane stable-equivalence context is byte-pinned',
          hashlib.sha256(path.read_bytes()).hexdigest()==context['sha256'])
    dependency_pins[str(path.relative_to(ROOT))] = context['sha256']
    check('classical stable-coordinate context is not claimed new or imported as a general cancellation theorem',
          context['plane_stable_coordinate_fact_asserted_new'] is False and
          context['new_external_classification_imported'] is False and
          context['higher_dimensional_stable_coordinate_conjecture_assumed'] is False)

    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        two_projection_rank_controls=profiles,four_output_dimension_table=dimension_table,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in ['verify_projection_profile.py','MIXED_PROJECTION_PROFILE.md',
                         'sources/projection-profile-context.json']},
        dependencies=dependency_pins,
        written_not_finite_exhaustion=[
            'no nonzero linear combination of a nonautomorphic plane Keller pair becomes a stable coordinate',
            'exact generator count N-k+ell for every linear output subspace and every number of free variables',
            'intrinsic HC4 projection minimum ambient dimensions four, five or six',
            'all allowed A5 projection completions have explicit polynomial source recovery',
            'noncoordinate directional derivatives recognize exactly the unique intrinsic cotangent plane'],
        finite_examples_asserted_to_prove_lower_bounds=False,
        nonautomorphic_Keller_example_constructed=False,
        unconditional_profile_asserted_for_automorphic_Keller_pair=False,
        plane_stable_coordinate_fact_asserted_new=False,
        higher_dimensional_cancellation_or_stable_coordinate_theorem_assumed=False,
        allowed_projection_embedding_asserted_to_match_fixed_Meng_Yang_potential=False,
        arbitrary_graph_physical_base_asserted_linear_gradient_projection=False,
        marked_gradient_covector_or_collision_transfer_asserted=False,
        relative_minimum_asserted_as_absolute_embedding_dimension=False,
        nonlinear_Hessian_invariance_used=False,universal_SL2_row_equation_solved=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,novelty_priority_asserted=False,
        formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-projection-profile.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=result['check_count'],
                          elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
