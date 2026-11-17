"""Exact algebra for all-degree native base-eliminated HC4 graphs.

Trace degree arguments and the complete two-plane classification are written
proofs. The finite leading-coefficient controls do not certify all degrees.
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

    x,y,p,q,r,U,V,Z = S.symbols('x y p q r U V Z')
    L,N,kappa = S.symbols('L N kappa',nonzero=True)
    M,alpha,beta,gamma,c,delta = S.symbols('M alpha beta gamma c delta')
    u = 1+x*y
    A = u**3*p+3*x*u**2*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = L*A*A+M*A+N*B
    R0 = S.Function('R0')(p,q,r)
    hy = L*p*p+M*p+4*N*p*R0**2+N*q*R0
    ax = p+3*R0*q-R0**3*r
    hx = L*ax*ax+M*ax+2*N*R0*r
    zero('native y graph exact x=0 trace',psi.subs(x,0).subs(y,R0)-hy)
    zero('native x graph exact y=0 trace',psi.subs(y,0).subs(x,R0)-hx)

    # Controls on the edge cases in the unbounded affine-direction proof.
    p0,q0,r0,ap,bq,cr,lead,trans = S.symbols('p0 q0 r0 ap bq cr lead trans')
    pv,qv,rv = p0+ap*V,q0+bq*V,r0+cr*V
    Rlinear = trans+lead*V
    y_control = L*pv**2+M*pv+4*N*pv*Rlinear**2+N*qv*Rlinear
    zero('y trace moving-p degree-three coefficient',S.Poly(y_control,V).coeff_monomial(V**3)-4*N*ap*lead**2)
    y_constant = L*pv**2+M*pv+4*N*pv*trans**2+N*qv*trans
    zero('y trace V-independent graph leaves a moving-p quadratic',S.Poly(y_constant,V).coeff_monomial(V**2)-L*ap**2)
    zero('y trace remaining quadratic factors exactly',
         S.Poly(y_control.subs(ap,0),V).coeff_monomial(V**2)-N*lead*(4*p0*lead+bq))
    y_invariant = L*p0*p0+M*p0+4*N*p0*trans**2+N*qv*trans
    zero('y trace annihilating p and graph is sufficient for an affine direction',S.diff(y_invariant,V,2))
    zero('y trace affine-direction coefficient is N times v(q) times graph',S.diff(y_invariant,V)-N*bq*trans)
    for degree in [1,2]:
        graph = trans+lead*V**degree
        argument = pv+3*graph*qv-graph**3*rv
        trace = L*argument**2+M*argument+2*N*graph*rv
        zero('x trace nonconstant graph leading coefficient with moving r, degree '+str(degree),
             S.Poly(trace,V).coeff_monomial(V**(6*degree+2))-L*lead**6*cr**2)
        zero('x trace nonconstant graph leading coefficient with transverse r, degree '+str(degree),
             S.Poly(trace.subs(cr,0),V).coeff_monomial(V**(6*degree))-L*lead**6*r0**2)
    argument = pv+3*c*qv-c**3*rv
    trace_constant_x = L*argument**2+M*argument+2*N*c*rv
    zero('x constant graph affine-direction quadratic coefficient',
         S.Poly(trace_constant_x,V).coeff_monomial(V**2)-L*(ap+3*c*bq-c**3*cr)**2)
    tangent_trace = trace_constant_x.subs(ap,-3*c*bq+c**3*cr)
    zero('x constant graph allowed direction has constant first derivative',S.diff(tangent_trace,V)-2*N*c*cr)

    # Both possible projections of the affine two-plane to the free base.
    FF,CC,DD = S.Function('FF')(U,Z),S.Function('CC')(U,Z),S.Function('DD')(U,Z)
    normal_one = DD+x*FF+V*CC
    jac_one = S.Matrix([FF,CC]).jacobian((U,Z)).det()
    zero('one-dimensional base projection full Hessian determinant',S.hessian(normal_one,(x,V,U,Z)).det()-jac_one**2)
    const_coefficient = S.symbols('const_coefficient')
    zero('constant pure-dual coefficient makes rank-one projection Hessian singular',
         S.hessian(DD+x*FF+V*const_coefficient,(x,V,U,Z)).det())
    F0,G0,D0 = S.Function('F0')(x,p),S.Function('G0')(x,p),S.Function('D0')(x,p)
    normal_zero = D0+q*F0+r*G0
    jac_zero = S.Matrix([F0,G0]).jacobian((x,p)).det()
    zero('zero-dimensional base projection full Hessian determinant',S.hessian(normal_zero,(q,r,x,p)).det()-jac_zero**2)
    c0,c1 = S.symbols('c0 c1')
    f1,g1 = S.Function('f1')(x,p),S.Function('g1')(x,p)
    constant_axis_pair = S.Matrix([c0+x*f1,c1+x*g1])
    zero('x-chart rank-zero coefficient pair is singular on the entire boundary line',
         constant_axis_pair.jacobian((x,p)).det().subs(x,0))
    line_pair = S.Matrix([N*(kappa*p+delta)+x*f1,x*g1])
    zero('y-chart rank-zero boundary Jacobian has the declared unit factors',
         line_pair.jacobian((x,p)).det().subs(x,0)+N*kappa*g1.subs(x,0))
    line_output = S.symbols('line_output')
    zero('y-chart boundary restriction is an explicitly invertible affine line map',
         ((line_output/N-delta)/kappa).subs(line_output,line_pair[0].subs(x,0))-p)

    # Derive the first normal coefficient literally from the original family.
    Hp = 2*L*p+M
    original_xjet = 3*Hp*(y*p+q)+N*(7*y**3*p+12*y*y*q+2*r)
    zero('full original first x jet on x=0',S.diff(psi,x).subs(x,0)-original_xjet)
    zero('original p derivative on x=0',S.diff(psi,p).subs(x,0)-Hp-4*N*y*y)
    zero('original q derivative on x=0',S.diff(psi,q).subs(x,0)-N*y)
    zero('original r derivative on x=0',S.diff(psi,r).subs(x,0))
    zero('original y derivative on x=0',S.diff(psi,y).subs(x,0)-N*(8*y*p+q))
    transformation = S.Matrix([x,U+alpha*x,V+beta*x,Z+c*V+gamma*x])
    derivative = transformation.jacobian((x,V,U,Z))
    check('adapted source change is constant affine and invertible',S.cancel(derivative.det())!=0)
    inverse = S.Matrix([x,q-beta*x,p-alpha*x,r-c*q-(gamma-c*beta)*x])
    zero('adapted source change inverse left',
         inverse.subs(dict(zip((x,p,q,r),transformation)),simultaneous=True)-S.Matrix([x,V,U,Z]))
    zero('adapted source change inverse right',
         transformation.subs(dict(zip((x,V,U,Z),inverse)),simultaneous=True)-S.Matrix([x,p,q,r]))
    R = S.Function('R')(U,Z)
    E = S.Function('E')(U,V,Z)
    first = (original_xjet+alpha*(Hp+4*N*y*y)+beta*N*y+N*(8*y*p+q)*E).subs(
        {p:U,q:V,r:Z+c*V,y:R},simultaneous=True)
    expected_first = V*(3*(2*L*U+M)+12*N*R*R+2*N*c)+N*(V+8*U*R)*E \
        +2*N*Z+3*U*(2*L*U+M)*R+7*N*U*R**3+alpha*(2*L*U+M+4*N*R*R)+beta*N*R
    zero('adapted normal first jet retains every chain-rule contribution',first-expected_first)
    e0,e1 = S.symbols('e0 e1')
    zero('positive V degree of graph normal jet raises the coefficient degree',
         S.Poly(expected_first.subs(E,e0+e1*V),V).coeff_monomial(V**2)-N*e1)
    forced_E = -3*(2*L*U+M)/N-12*R*R-2*c
    forced_F = 2*N*Z-21*U*(2*L*U+M)*R-89*N*U*R**3-16*N*c*U*R \
        +alpha*(2*L*U+M+4*N*R*R)+beta*N*R
    zero('forced graph normal jet cancels the entire V coefficient',S.diff(expected_first.subs(E,forced_E),V))
    zero('forced normal coefficient has exact -21 and -89 terms',expected_first.subs(E,forced_E)-forced_F)
    forced_C = N*R
    J = S.Matrix([forced_F,forced_C]).jacobian((U,Z)).det()
    expected_J = N*((2*L*alpha-(84*L*U+21*M+16*N*c)*R-89*N*R**3)*S.diff(R,Z)-2*N*S.diff(R,U))
    zero('residual plane Jacobian is the declared first-order polynomial equation',J-expected_J)
    h,k = S.Function('h')(U),S.Function('k')(U)
    for degree in [1,2,3]:
        control = h+k*Z**degree
        control_J = N*((2*L*alpha-(84*L*U+21*M+16*N*c)*control-89*N*control**3)*S.diff(control,Z)-2*N*S.diff(control,U))
        zero('degree-control leading coefficient at Z degree '+str(degree),
             S.Poly(control_J,Z).coeff_monomial(Z**(4*degree-1))+89*N*N*degree*k**4)
    affine_R = kappa*U+delta
    F_tri = forced_F.subs(R,affine_R).doit()
    C_tri = N*affine_R
    zero('forced residual pair is triangular with unit plane Jacobian',
         S.Matrix([F_tri,C_tri]).jacobian((U,Z)).det()+2*N*N*kappa)
    zero('triangular first coefficient has constant Z derivative',S.diff(F_tri,Z)-2*N)
    D_tri = L*U*U+M*U+4*N*U*affine_R**2
    triangular_normal = D_tri+x*F_tri+V*C_tri
    zero('full triangular normal form unit Hessian',S.hessian(triangular_normal,(x,V,U,Z)).det()-4*N**4*kappa**2)
    # Inverse controls specialize constants only to keep their polynomial
    # expansions short. The inverse formula in the proof is general.
    concrete = triangular_normal.subs({L:1,M:13,N:2,kappa:1,delta:0,alpha:0,beta:0,c:0})
    concrete_F0 = F_tri.subs({L:1,M:13,N:2,kappa:1,delta:0,alpha:0,beta:0,c:0,Z:0})
    concrete_D = D_tri.subs({L:1,M:13,N:2,kappa:1,delta:0})
    z1,z2,z3,z4 = S.symbols('z1 z2 z3 z4')
    output = S.Matrix([z1,z2,z3,z4])
    rec_U,rec_x = z2/2,z4/4
    rec_Z = (z1-concrete_F0.subs(U,rec_U))/4
    rec_V = (z3-S.diff(concrete_D,U).subs(U,rec_U)-rec_x*S.diff(concrete_F0,U).subs(U,rec_U))/2
    inverse_gradient = S.Matrix([rec_x,rec_V,rec_U,rec_Z])
    variables = (x,V,U,Z)
    gradient = S.Matrix([S.diff(concrete,z) for z in variables])
    zero('triangular normal-form gradient inverse left',
         inverse_gradient.subs(dict(zip(output,gradient)),simultaneous=True)-S.Matrix(variables))
    zero('triangular normal-form gradient inverse right',
         gradient.subs(dict(zip(variables,inverse_gradient)),simultaneous=True)-output)

    # Explicitly retain the actual external inputs, not just bibliography.
    manifest_name = 'sources/critical-jet-import.json'
    manifest = json.loads((HERE/manifest_name).read_text())
    source_name = 'sources/'+manifest['retained_file']
    check('primary Gwozdziewicz theorem source byte pin',
          hashlib.sha256((HERE/source_name).read_bytes()).hexdigest()==manifest['sha256']
          =='4cce6a7334640ea47036b1fb8c676e28da7946dbd15d0a48868b30d94c68585a')
    pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    dependencies = {}
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md','hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in [Path(__file__).name,'BASE_ELIMINATED_GRAPHS.md',manifest_name,source_name]},
        dependencies=dependencies,
        written_not_finite_exhaustion=['all-degree affine-direction classification of both traces',
            'both projection ranks for both native base graphs','arbitrary-degree 4m-1 residual equation',
            'polynomial-unit argument and line-injectivity import','characteristic-zero descent'],
        native_coordinate_graph_charts_closed=['x','y','p','q','r'],
        x_chart_unit_Hessian_has_no_jointly_affine_two_plane=True,
        y_chart_surviving_plane_gradient_is_invertible=True,
        nonlinear_reparameterizations_of_base_eliminated_charts_closed=False,
        singleton_zero_properness_resolved=False,
        arbitrary_nonlinear_base_dual_mixing_closed=False,
        general_HC4_or_JC2_resolved=False,
        imported_theorems=['HCTANG2','HC4MYGL0 for the y-chart pure-dual-plane case'],
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
