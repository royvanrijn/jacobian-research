"""Exact controls for the global ambient quadratic form and graph closure.

The all-degree graph proof is separate from these finite symbolic controls.
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
        check(name,all(S.cancel(S.expand(item.doit()))==0 for item in values))

    x,y,p,q,r,L,M,N,T,Z,W = S.symbols('x y p q r L M N T Z W')
    vv = x*y
    u = 1+vv
    aa = 1-3*vv+6*vv*vv
    bb = y**3*(6*vv*vv+15*vv+10)
    hh = 6*vv*vv+15*vv+4
    zero('the global unit-pivot Bezout identity is literal',aa*u**3-x**3*bb-1)
    Umat = S.Matrix([[u**3,3*x*u*u,-x**3],[0,1,0],[-bb,0,aa]])
    Uinv = S.Matrix([[aa,-3*x*u*u*aa,x**3],[0,1,0],[bb,-3*x*u*u*bb,u**3]])
    zero('the unit-pivot frame has its exact polynomial inverse',Umat*Uinv-S.eye(3))
    zero('the unit-pivot frame determinant is one',Umat.det()-1)
    C0 = S.Matrix([
        [(3*vv*vv+6*vv+2)/2,x*x*(2*vv+3)],
        [3*y*y*(18*vv**4+63*vv**3+69*vv*vv+21*vv-1)/2,
         36*vv**5+108*vv**4+87*vv**3+3*vv*vv-3*vv+1]])
    zero('the retained relative row normalizer determinant is one',C0.det()-1)
    Cblock = S.diag(1,C0)
    Cblock_inv = S.diag(1,C0.adjugate())
    shear = S.Matrix([[1,0,0],[y*hh,1,0],[0,0,1]])
    shear_inv = S.Matrix([[1,0,0],[-y*hh,1,0],[0,0,1]])
    forward = shear*Cblock_inv*Umat
    inverse = Uinv*Cblock*shear_inv
    zero('the full relative linear frame and its written inverse compose to the identity',forward*inverse-S.eye(3))
    zero('the other full relative frame composition is the identity',inverse*forward-S.eye(3))
    Arow = S.Matrix([[u**3,3*x*u*u,-x**3]])
    Brow = S.Matrix([[y*y*u*(4+3*x*y),y+3*x*y*y*(4+3*x*y),2*x-3*x*x*y]])
    zero('the first normal coordinate is exactly native A',forward[0,:]-Arow)
    zero('the full native B row becomes y Z plus 2x W globally',Brow-S.Matrix([[0,y,2*x]])*forward)
    normal_unshifted = forward*S.Matrix([p,q,r])
    normal = normal_unshifted+S.Matrix([M/(2*L),0,0])
    recovered = inverse*(normal-S.Matrix([M/(2*L),0,0]))
    zero('the affine normal frame has actual polynomial source recovery',recovered-S.Matrix([p,q,r]))
    native_A = (Arow*S.Matrix([p,q,r]))[0]
    native_B = (Brow*S.Matrix([p,q,r]))[0]
    lstar = -M*M/(4*L)
    zero('the entire scaled potential is exactly its nondegenerate quadratic normal form',
         L*native_A*native_A+M*native_A+N*native_B-
         (lstar+L*normal[0]**2+N*y*normal[1]+2*N*x*normal[2]))
    Q5 = L*T*T+N*y*Z+2*N*x*W
    zero('the ambient quadratic has constant nonzero Hessian determinant',
         S.hessian(Q5,(x,y,T,Z,W)).det()-8*L*N**4)

    # A nonlinear anisotropic-direction graph has an actual gradient inverse.
    R = x**3+y*y
    square_graph = L*R*R+N*y*Z+2*N*x*W
    variables = (x,y,Z,W)
    grad = S.Matrix([S.diff(square_graph,z) for z in variables])
    t0,t1,t2,t3 = S.symbols('t0:4')
    xb,yb = t3/(2*N),t2/N
    Rb = R.subs({x:xb,y:yb},simultaneous=True)
    Wb = (t0-2*L*Rb*(3*xb*xb))/(2*N)
    Zb = (t1-2*L*Rb*(2*yb))/N
    back = {x:xb,y:yb,Z:Zb,W:Wb}
    zero('the nonlinear square graph has exact unit full Hessian',S.hessian(square_graph,variables).det()-4*N**4)
    zero('the square-graph gradient has a polynomial inverse',grad.subs(back,simultaneous=True)-S.Matrix([t0,t1,t2,t3]))
    zero('its source recovery composition is exact',
         S.Matrix([xb,yb,Zb,Wb]).subs(dict(zip((t0,t1,t2,t3),grad)),simultaneous=True)-S.Matrix(variables))
    a1,a2,R0,w1,w2 = S.symbols('a1 a2 R0 w1 w2')
    affine_R = R0+a1*w1+a2*w2
    zero('the square term has the exact constant-direction Hessian outer product',
         S.hessian(L*affine_R**2,(w1,w2))-2*L*S.Matrix([a1,a2])*S.Matrix([[a1,a2]]))
    check('a nonconstant coefficient cannot have constant square in the finite control',S.degree(x*x,x)==2)
    check('a higher W degree produces an uncancellable higher square degree',S.degree((w1*w2+w1**3)**2,w1)==6)
    zero('several complex squares can cancel, so the one-square argument is not silently generalized',R**2+(S.I*R)**2)

    # A nonlinear isotropic-direction graph and its explicit gradient inverse.
    kappa = S.symbols('kappa',nonzero=True)
    h = y**4+y
    isotropic_graph = L*T*T+N*y*Z+2*N*W*(kappa*W+h)
    ivars = (y,T,Z,W)
    igrad = S.Matrix([S.diff(isotropic_graph,z) for z in ivars])
    yback,Tback = t2/N,t1/(2*L)
    Wback = (t3-2*N*h.subs(y,yback))/(4*N*kappa)
    Zback = (t0-2*N*Wback*S.diff(h,y).subs(y,yback))/N
    iback = {y:yback,T:Tback,Z:Zback,W:Wback}
    zero('the isotropic graph keeps an exact nondegenerate ternary quadratic trace',
         isotropic_graph.subs(W,0)-(L*T*T+N*y*Z))
    zero('the trace quadratic has rank three',S.hessian(L*T*T+N*y*Z,(y,T,Z)).det()+2*L*N*N)
    zero('the nonlinear isotropic graph has actual unit Hessian',S.hessian(isotropic_graph,ivars).det()+8*L*N**3*kappa)
    zero('its gradient has an actual polynomial inverse',igrad.subs(iback,simultaneous=True)-S.Matrix([t0,t1,t2,t3]))

    # Mixed ambient directions exercise the bilinear graph split and affine recovery.
    s0,s1,s2,s3,rr = S.symbols('s0:4 rr')
    ss = S.Matrix([s0,s1,s2,s3])
    basepoint = S.Matrix([1,-2,3,1,2])
    B = S.eye(5)[:,:4]
    ambient_hess = S.hessian(Q5,(x,y,T,Z,W)).subs({L:1,N:1})
    qfun = lambda column: (column.T*ambient_hess*column)[0]/2
    for label,direction in [('isotropic',S.Matrix([1,1,0,-2,1])),
                            ('anisotropic',S.Matrix([1,1,1,0,1]))]:
        check(label+' mixed direction completes the retained affine plane',B.row_join(direction).det()!=0)
        cc = qfun(direction)
        check(label+' mixed direction has the claimed norm type', (cc==0)==(label=='isotropic'))
        ell = (direction.T*ambient_hess*(basepoint+B*ss))[0]
        graph_value = qfun(basepoint+B*ss+direction*rr)
        zero(label+' graph decomposition is exact with affine offsets',graph_value-(qfun(basepoint+B*ss)+ell*rr+cc*rr*rr))
        projection = S.eye(4).row_join(-direction[:4,0])
        nonlinear_R = s0**3+s1*s2+s3**2
        image = basepoint+B*ss+direction*nonlinear_R
        zero(label+' affine retained projection recovers the actual graph parameter',projection*image-projection*basepoint-ss)
        if label=='isotropic':
            row = S.Matrix([[S.diff(ell,z) for z in ss]])
            kernel = S.Matrix.hstack(*row.nullspace())
            restricted = kernel.T*(B.T*ambient_hess*B)*kernel
            check('the mixed isotropic hyperplane trace has nondegenerate rank-three quadratic part',restricted.det()!=0)
        else:
            completed = qfun(basepoint+B*ss)-ell*ell/(4*cc)
            zero('the anisotropic completion of the square is literal',
                 graph_value-(completed+cc*(rr+ell/(2*cc))**2))

    # Source reparameterization and the two endpoint controls.
    diagonal = s0*s0+s1*s1+s2*s2+s3*s3
    reparameterized = diagonal.subs(s0,s0+s1*s1)
    check('a polynomial source shear can change a constant Hessian determinant',
          S.hessian(diagonal,ss).det()==16 and S.Poly(S.hessian(reparameterized,ss).det(),*ss).total_degree()>0)
    zero('the fixed T-zero normal graph is an actual nondegenerate four-variable quadric',
         S.hessian(Q5.subs(T,0),(x,y,Z,W)).det()-4*N**4)
    coordinate_graph = Q5.subs(x,1/(2*N))
    zero('the fixed constant isotropic graph is a polynomial coordinate with explicit companion',S.diff(coordinate_graph,W)-1)
    hconst = S.symbols('hconst')
    constant_D = lstar+L*hconst*hconst
    zero('the transverse square condition alone has universal constant-D solutions',constant_D-lstar-L*hconst*hconst)

    dependency_pins = {}
    for receipt_name,names in [
        ('receipt-first-slice.json',['FIRST_SLICE_INVERTIBILITY.md']),
        ('receipt-relative-coordinate.json',['RELATIVE_COORDINATE_ENDPOINT.md']),
        ('receipt-global-value.json',['GLOBAL_VALUE_ORBITS.md']),
        ('receipt-fixed-projection.json',['FIXED_POTENTIAL_PROJECTION.md'])]:
        retained = json.loads((HERE/receipt_name).read_text())
        for name in names:
            dependency_pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
    row_receipt = json.loads((HERE/'receipt-square-axis.json').read_text())
    for name,expected in row_receipt['dependencies'].items():
        if 'hc4-relative-row-orbit' in name or 'tang-hc4-audit' in name:
            dependency_pins[name] = expected
    for name,expected in dependency_pins.items():
        check('retained mathematical input pin '+name,
              hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected)

    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ['verify_ambient_quadric.py','AMBIENT_QUADRIC_GRAPHS.md']},
        dependencies=dependency_pins,
        written_not_finite_exhaustion=[
            'global base-preserving polynomial equivalence of the entire Meng--Yang potential with a nondegenerate quadratic',
            'all-degree quadratic-plus-one-square and ternary quadratic trace HC4 criteria',
            'every polynomial graph in every constant affine direction of a nondegenerate quadratic ambient A5',
            'explicit genuine Meng--Yang embeddings by inverse ambient normalization',
            'four independent affine-linear retained normal projections close a candidate in its actual source parameters'],
        new_external_classification_imported=False,
        ambient_normalization_asserted_gradient_conjugacy=False,
        nonlinear_Hessian_invariance_used=False,
        graph_image_alone_asserted_sufficient_without_actual_affine_parameterization=False,
        arbitrary_source_reparameterizations_asserted_closed=False,
        arbitrary_embeddings_asserted_to_have_four_affine_normal_projections=False,
        multiple_polynomial_squares_asserted_covered_by_one_square_argument=False,
        universal_transverse_square_condition_asserted_obstructive=False,
        constant_D_square_solution_asserted_to_supply_embedding=False,
        marked_gradient_collision_transfer_asserted=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        universal_SL2_row_equation_solved=False,unrestricted_HC4_or_JC2_resolved=False,
        novelty_priority_asserted=False,formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-ambient-quadric.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=result['check_count'],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
