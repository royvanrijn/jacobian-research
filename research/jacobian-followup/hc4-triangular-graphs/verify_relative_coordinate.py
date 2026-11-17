"""Exact controls for relative coordinates and the omitted-point endpoint.

These controls verify algebra and actual polynomial inverses. They neither
exhibit a unit-ideal Keller pair nor prove the imported Tang alternative.
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
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    checks = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def zero(name, expression):
        values = list(expression) if isinstance(expression, S.MatrixBase) else [expression]
        check(name, all(S.cancel(value.doit()) == 0 for value in values))

    x,y,p,q,r,s,U,P,Q,Z,V,W = S.symbols('x y p q r s U P Q Z V W')
    L,N = S.symbols('L N',nonzero=True)
    M = S.symbols('M')
    u = 1+x*y
    rowA = S.Matrix([[u**3,3*x*u*u,-x**3]])
    rowB = S.Matrix([[y*y*u*(4+3*x*y),y*(9*x*x*y*y+12*x*y+1),x*(2-3*x*y)]])
    dual = S.Matrix([p,q,r])
    A,B = (rowA*dual)[0],(rowB*dual)[0]
    psi = L*A*A+M*A+N*B
    zero('whole old base-origin slice is the p composite',psi.subs({x:0,y:0})-L*p*p-M*p)

    # Universal bordered composite-trace Hessian, with arbitrary normal jet.
    h,k,n0 = S.symbols('h k n0')
    aa = S.Matrix(S.symbols('a1 a2 a3'))
    nn = S.Matrix(S.symbols('n1 n2 n3'))
    t11,t12,t13,t22,t23,t33 = S.symbols('t11 t12 t13 t22 t23 t33')
    tangent = S.Matrix([[t11,t12,t13],[t12,t22,t23],[t13,t23,t33]])
    block = h*tangent+k*aa*aa.T
    adj = block.adjugate()
    zero('dual adjugate has the H-prime factor entrywise',adj.subs(h,0))
    block_det = block.det()
    zero('dual determinant has at least two H-prime factors, constant term',block_det.subs(h,0))
    zero('dual determinant has at least two H-prime factors, linear term',S.diff(block_det,h).subs(h,0))
    bordered = n0*block_det-(nn.T*adj*nn)[0]
    zero('full bordered determinant vanishes on the critical composite level',bordered.subs(h,0))
    constant_trace_data = dict.fromkeys([t11,t12,t13,t22,t23,t33]+list(aa),0)
    zero('constant boundary trace also makes the full Hessian singular',bordered.subs(constant_trace_data))
    composite = P+Q*Q+Z**3
    trace = L*composite**2+M*composite
    trace_gradient = S.Matrix([S.diff(composite,z) for z in (P,Q,Z)])
    zero('actual nonlinear three-dual composite has the required Hessian block',
         S.hessian(trace,(P,Q,Z))-(2*L*composite+M)*S.hessian(composite,(P,Q,Z))
         -2*L*trace_gradient*trace_gradient.T)
    zero('critical composite level is nonempty by an explicit control',
         (2*L*composite+M).subs(P,-M/(2*L)-Q*Q-Z**3))

    # Literal two-by-two minors and the complete origin-ideal certificate.
    minors = S.Matrix([rowA[0]*rowB[1]-rowA[1]*rowB[0],
        rowA[0]*rowB[2]-rowA[2]*rowB[0],rowA[1]*rowB[2]-rowA[2]*rowB[1]])
    expected = S.Matrix([y*u**3,x*u*(u+1),2*x*x*(2*u+1)])
    zero('all three original coefficient-row minors',minors-expected)
    squared_x = (2*u+1)*minors[2]/2-4*x*minors[1]
    zero('minor ideal certificate for x squared',squared_x-x*x)
    recovered_x = minors[1]/2-(3*y+x*y*y)*squared_x/2
    zero('minor ideal certificate for x',recovered_x-x)
    recovered_y = minors[0]-(3*y*y+3*x*y**3+x*x*y**4)*recovered_x
    zero('minor ideal certificate for y',recovered_y-y)
    ex = S.Matrix([[0,S.Rational(1,2)+2*x*(3*y+x*y*y),-(3*y+x*y*y)*(2*u+1)/4]])
    ey = S.Matrix([[1,0,0]])-(3*y*y+3*x*y**3+x*x*y**4)*ex
    zero('explicit coefficient vector recovers x from the minors',(ex*minors)[0]-x)
    zero('explicit coefficient vector recovers y from the minors',(ey*minors)[0]-y)
    bez_x,bez_y = S.symbols('bez_x bez_y')
    certificate = bez_x*ex+bez_y*ey
    third_row = S.Matrix([[certificate[2],-certificate[1],certificate[0]]])
    general_frame = S.Matrix.vstack(rowA,rowB,third_row)
    zero('universal relative-frame determinant is the base-curve Bezout certificate',
         general_frame.det()-(bez_x*x+bez_y*y))
    zero('nonconstant graph y=1+s cubed avoids the origin with a polynomial certificate',
         (-s*s*x+y).subs({x:s,y:1+s**3})-1)
    zero('opposite nonconstant graph x=1+s cubed has its own certificate',
         (x-s*s*y).subs({x:1+s**3,y:s})-1)

    # Actual simple frame and its polynomial inverse on one fixed graph.
    uc = 1+s
    Ac = uc**3*p+3*uc**2*q-r
    Tc = uc*(uc+1)*p+2*(2*uc+1)*q
    Wc = -(2*uc+1)*p/2-4*q
    Bc = B.subs({x:1,y:s})
    frame = S.Matrix([Ac,Tc,Wc])
    zero('fixed x=1 frame has constant determinant minus one',frame.jacobian((p,q,r)).det()+1)
    zero('fixed x=1 original B is the declared frame combination',Bc-Tc+(2-3*s)*Ac)
    original_fixed_potential = psi.subs({x:1,y:s})
    zero('native fixed-base potential remains Hessian-degenerate',S.hessian(original_fixed_potential,(s,p,q,r)).det())
    rec_T = (V-L*U*U-M*U)/N+(2-3*s)*U
    rec_p = -4*rec_T-2*(2*uc+1)*W
    rec_q = (2*uc+1)*rec_T/2+uc*(uc+1)*W
    rec_r = uc**3*rec_p+3*uc**2*rec_q-U
    physical_inverse = S.Matrix([rec_p,rec_q,rec_r])
    physical_coordinates = S.Matrix([Ac,original_fixed_potential,Wc])
    zero('fixed potential-coordinate inverse left',
         physical_inverse.subs(dict(zip((U,V,W),physical_coordinates)),simultaneous=True)-dual)
    zero('fixed potential-coordinate inverse right',
         physical_coordinates.subs(dict(zip((p,q,r),physical_inverse)),simultaneous=True)-S.Matrix([U,V,W]))
    zero('fixed potential-coordinate inverse is a genuine relative automorphism',
         physical_inverse.jacobian((U,V,W)).det()+1/N)

    # Plane orientations of a hypothetical unit-Hessian relative coordinate.
    F,G,D = S.Function('F')(s,U),S.Function('G')(s,U),S.Function('D')(s,U)
    cotangent = D+P*F+Q*G
    J = S.Matrix([F,G]).jacobian((s,U)).det()
    zero('rank-zero relative-coordinate normal form Hessian',S.hessian(cotangent,(s,U,P,Q)).det()-J*J)
    A1,C1,D1 = S.Function('A1')(U,Z),S.Function('C1')(U,Z),S.Function('D1')(U,Z)
    rank_one = D1+s*A1+V*C1
    J1 = S.Matrix([A1,C1]).jacobian((U,Z)).det()
    zero('rank-one relative-coordinate normal form Hessian',S.hessian(rank_one,(s,V,U,Z)).det()-J1*J1)
    solution_s = (-S.diff(D1,U)*S.diff(C1,Z)+S.diff(D1,Z)*S.diff(C1,U))/J1
    solution_V = (-S.diff(A1,U)*S.diff(D1,Z)+S.diff(A1,Z)*S.diff(D1,U))/J1
    zero('rank-one dual U derivative can always be killed on a zero of C',
         S.diff(D1,U)+solution_s*S.diff(A1,U)+solution_V*S.diff(C1,U))
    zero('rank-one dual Z derivative can always be killed at the same point',
         S.diff(D1,Z)+solution_s*S.diff(A1,Z)+solution_V*S.diff(C1,Z))
    zero('rank-one remaining dual derivative is exactly C',S.diff(rank_one,V)-C1)
    du,fu,gu,fs,gs = S.symbols('du fu gu fs gs')
    zero('unit plane Jacobian excludes simultaneous zero U derivatives',(fs*gu-fu*gs).subs({fu:0,gu:0}))
    zero('rank-zero dual critical point when F_U is nonzero',(du+P*fu+Q*gu).subs({P:-du/fu,Q:0}))
    zero('rank-zero dual critical point when G_U is nonzero',(du+P*fu+Q*gu).subs({P:0,Q:-du/gu}))

    # Complete the unimodular row explicitly; check formulas before imposing
    # the Bezout relation, so no symbolic inverse of the plane pair is hidden.
    ab,bb = S.Function('ab')(s,U),S.Function('bb')(s,U)
    bezout = ab*F+bb*G
    row_W = -bb*P+ab*Q
    coordinate_map = S.Matrix([U,cotangent,row_W])
    zero('relative completion determinant is precisely the Bezout polynomial',
         coordinate_map.jacobian((U,P,Q)).det()-bezout)
    back_P,back_Q = ab*(V-D)-G*W,bb*(V-D)+F*W
    zero('completion inverse left P modulo the exact Bezout relation',
         back_P.subs({V:cotangent,W:row_W},simultaneous=True)-P-(bezout-1)*P)
    zero('completion inverse left Q modulo the exact Bezout relation',
         back_Q.subs({V:cotangent,W:row_W},simultaneous=True)-Q-(bezout-1)*Q)
    zero('completion inverse right coordinate modulo the exact Bezout relation',
         cotangent.subs({P:back_P,Q:back_Q},simultaneous=True)-V-(bezout-1)*(V-D))
    zero('completion inverse right companion modulo the exact Bezout relation',
         row_W.subs({P:back_P,Q:back_Q},simultaneous=True)-W-(bezout-1)*W)
    zero('dual gradient has a certified unit combination',ab*S.diff(cotangent,P)+bb*S.diff(cotangent,Q)-bezout)

    # A genuine fixed-graph unit-row control that deliberately fails Keller.
    control_F,control_G = s,s*U-1
    control_a,control_b = U,-1
    control_f = P*control_F+Q*control_G
    control_W = P+U*Q
    zero('unimodular non-Keller control row has literal Bezout certificate',control_a*control_F+control_b*control_G-1)
    zero('unimodular control has nonconstant plane Jacobian',S.Matrix([control_F,control_G]).jacobian((s,U)).det()-s)
    zero('actual control cotangent Hessian is s squared',S.hessian(control_f,(s,U,P,Q)).det()-s*s)
    control_phi = physical_inverse.subs({V:control_f,W:control_W},simultaneous=True)
    zero('actual fixed-graph retained-dual control has constant nonzero Jacobian',control_phi.jacobian((U,P,Q)).det()+1/N)
    zero('original potential on the actual embedded control is the declared cotangent',
         original_fixed_potential.subs(dict(zip((p,q,r),control_phi)),simultaneous=True)-control_f)
    parameter_inverse = S.Matrix([Ac,Ac*original_fixed_potential-(s*Ac-1)*Wc,-original_fixed_potential+s*Wc])
    zero('actual graph parameter inverse left',
         parameter_inverse.subs(dict(zip((p,q,r),control_phi)),simultaneous=True)-S.Matrix([U,P,Q]))
    zero('actual graph parameter inverse right',
         control_phi.subs(dict(zip((U,P,Q),parameter_inverse)),simultaneous=True)-dual)
    zero('unit-ideal control gradient really omits zero',U*S.diff(control_f,P)-S.diff(control_f,Q)-1)

    # The opposite control has unit Hessian but fails relative coordinate.
    unit_f = s*P+U*Q
    zero('common-zero control has unit full Hessian',S.hessian(unit_f,(s,U,P,Q)).det()-1)
    zero('common-zero control has a dual critical point and is not a relative coordinate',
         S.Matrix([S.diff(unit_f,z) for z in (U,P,Q)]).subs({s:0,U:0,Q:0}))
    one_control = U*U+s*U+V*Z
    zero('rank-one excluded orientation can have unit Hessian without relative coordinate',S.hessian(one_control,(s,V,U,Z)).det()-1)
    zero('rank-one control has the reconstructed dual critical points',
         S.Matrix([S.diff(one_control,z) for z in (U,Z,V)]).subs({s:-2*U,V:0,Z:0}))

    pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    dependencies = {}
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md','jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in [Path(__file__).name,'RELATIVE_COORDINATE_ENDPOINT.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['bordered composite-hyperplane trace rank obstruction',
            'all base-only graph curves via universal minor ideal',
            'both plane-projection ranks for every relative coordinate',
            'Nullstellensatz and global unit-ideal necessity',
            'full equivalence to omitted targets of plane Keller maps'],
        zero_boundary_trace_has_no_unit_Hessian_after_relative_dual_changes=True,
        origin_avoiding_base_only_graphs_are_relative_coordinates=True,
        unit_Hessian_relative_coordinate_exists_iff_unit_ideal_plane_Keller_pair_exists=True,
        one_fixed_x_equals_one_graph_realizes_entire_endpoint=True,
        genuine_graph_embedding_and_relative_dual_inverses_explicit=True,
        unit_ideal_plane_Keller_pair_exhibited=False,
        plane_Keller_surjectivity_proved=False,
        equivalence_to_unrestricted_JC2_asserted=False,
        singleton_zero_properness_resolved=False,
        fixed_marked_collision_transferred=False,
        imported_theorems=['HCTANG2 for necessity in the relative-coordinate classification'],
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
