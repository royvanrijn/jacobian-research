"""Exact algebra for the q square-axis theorem and fixed p-graph orbit.

The degree lemma, all-plane exhaustion, singleton ideal and complete row
orbit are written arguments. No bounded search certifies those arguments.
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

    x, y, p, q, r, P, Q, U, V, T, W = S.symbols('x y p q r P Q U V T W')
    L, N, alpha = S.symbols('L N alpha', nonzero=True)
    M, eta, tau, sigma = S.symbols('M eta tau sigma')
    v, u = x*y, 1+x*y
    A = u**3*p+3*x*u*u*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = L*A*A+M*A+N*B
    zero('q-chart y=0 axis identity retains the entire graph function',
         psi.subs(y, 0)-(L*(p+3*x*q-x**3*r)**2+M*(p+3*x*q-x**3*r)+2*N*x*r))
    zero('base origin is the retained-p nonlinear square', psi.subs({x: 0, y: 0})-(L*p*p+M*p))
    # Affine dual normalization and independence of the kernel coefficient.
    a,b,c,d,e,g = S.symbols('a b c d e g')
    p0, r0 = a*P+b*Q+c, d*P+e*Q+g
    kernel = S.Matrix([b, -a])
    zero('initial pure-dual kernel annihilates p0', (S.Matrix([[a, b]])*kernel)[0])
    zero('initial r0 kernel derivative is the nonzero affine determinant',
         (S.Matrix([[d, e]])*kernel)[0]+a*e-b*d)
    k = S.Function('k')(x)
    C0 = 2*N*x*(1+x*k)
    zero('rank-one boundary coefficient has the forced nonzero initial slope', S.diff(C0, x).subs(x, 0)-2*N)
    transported = C0.subs(x, x-eta*y)
    zero('rank-one transport solves the entire constant-flow equation',
         eta*S.diff(transported, x)+S.diff(transported, y)+tau*S.diff(transported, P))
    transform = S.Matrix([y, x-eta*y, P-tau*y, Q-sigma*y])
    zero('rank-one adapted source change is constant and invertible',
         transform.jacobian((x,y,P,Q)).det()+1)
    zero('rank-one base generator becomes the first fibre direction',
         transform.jacobian((x,y,P,Q))*S.Matrix([eta,1,tau,sigma])-S.Matrix([1,0,0,0]))
    F1, B1 = S.Function('F1')(T,U), S.Function('B1')(T,U)
    C1 = S.Function('C1')(T)
    normal_one = C1*V+y*F1+B1
    zero('rank-one full Hessian determinant',
         S.hessian(normal_one,(y,V,T,U)).det()-(S.diff(C1,T)*S.diff(F1,U))**2)
    F2,G2 = S.Function('F2')(U,V),S.Function('G2')(U,V)
    normal_two = L*U*U+M*U+x*F2+y*G2
    J2 = S.Matrix([F2,G2]).jacobian((U,V))
    zero('rank-two full Hessian determinant', S.hessian(normal_two,(x,y,U,V)).det()-J2.det()**2)
    ax = S.Function('ax')(x,U,V)
    rx = S.Function('ra')(x)*U+S.Function('rb')(x)*V+S.Function('rc')(x)
    axis = L*ax*ax+M*ax+2*N*x*rx
    target = L*U*U+M*U+x*F2
    zero('rank-two comparison is the exact polynomial-square identity',
         4*L*(axis-target)-((2*L*ax+M)**2-(2*L*U+M)**2-4*L*x*F2+8*L*N*x*rx))

    # A negative control on the essential affine-residual hypothesis.
    square_control = U+x*V*V
    Hcontrol, Econtrol = 2*U*V*V, x*V**4
    zero('nonlinear residual permits a nonlinear H in the square identity',
         square_control**2-U*U-x*Hcontrol-x*Econtrol)
    check('square lemma negative control really violates the affine residual bound', S.total_degree(Econtrol,U,V)>1)
    outside_axis = L*U*U+M*U+2*N*x*(V+U*U)+y*U
    zero('non-affine-r axis control can still have unit Hessian',
         S.hessian(outside_axis,(x,y,U,V)).det()-4*N*N)
    check('non-affine-r axis control x coefficient is not affine', S.Poly(S.diff(outside_axis,x),U,V).total_degree()==2)

    # Actual relative affine and triangular-base automorphisms permitted
    # in the graph corollary, without truncating the q graph itself.
    p_affine=(1+x*y)*P+x*Q
    r_affine=y*P+Q
    zero('relative affine retained-dual map determinant',
         S.Matrix([p_affine,r_affine]).jacobian((P,Q)).det()-1)
    Rgraph = S.Function('Rgraph')(x,y,P,Q)
    full_q = psi.subs({p:p_affine,q:Rgraph,r:r_affine},simultaneous=True)
    axis_A = p_affine.subs(y,0)+3*x*Rgraph.subs(y,0)-x**3*r_affine.subs(y,0)
    zero('arbitrary q graph with relative affine dual changes satisfies the square axis',
         full_q.subs(y,0)-(L*axis_A**2+M*axis_A+2*N*x*Q))
    s,t=S.symbols('s t')
    base=S.Matrix([s+t**3+2,t])
    zero('permitted nonlinear triangular base map determinant',base.jacobian((s,t)).det()-1)
    zero('old y=0 has an affine parameter axis',base.subs({s:x-2,t:0})-S.Matrix([x,0]))

    # A literal square-axis unit-Hessian control and both gradient inverses.
    f_axis=L*P*P+M*P+2*N*x*Q+y*P
    grad_axis=S.Matrix([S.diff(f_axis,z) for z in (x,y,P,Q)])
    z1,z2,z3,z4=S.symbols('z1 z2 z3 z4')
    axis_inverse=S.Matrix([z4/(2*N),z3-2*L*z2-M,z2,z1/(2*N)])
    zero('square-axis explicit gradient inverse left composition',
         axis_inverse.subs(dict(zip((z1,z2,z3,z4),grad_axis)),simultaneous=True)-S.Matrix([x,y,P,Q]))
    zero('square-axis explicit gradient inverse right composition',
         grad_axis.subs(dict(zip((x,y,P,Q),axis_inverse)),simultaneous=True)-S.Matrix([z1,z2,z3,z4]))

    # p-chart constancy identity, valid before imposing any graph.
    zero('p-chart fibre constancy completed-square identity',
         u*u*psi-L*u*u*A*A-(M*u*u+N*y*y*(3*u+1))*A-N*y*u*u*q-N*x*(u+1)*r)
    zero('p-chart hyperbola fibre has a nonzero q coefficient',
         S.diff(psi.subs(y,-1/x),q)-2*N/x)
    zero('p-chart linear constancy coefficients vanish together only at the origin on u nonzero',
         (N*x*(u+1)/(u*u)).subs(y,0)-2*N*x)

    # One fixed degree-six graph realizes the entire oriented row orbit.
    pstar=-2*x*(v*v+1)*q-2*x**3*r
    check('fixed physical p graph has total degree six',S.Poly(pstar,x,y,q,r).total_degree()==6)
    D=S.Matrix([
        [u*u*(1-2*v-2*v*v-2*v**3),-x*x*(2*v**3+6*v*v+6*v+3)],
        [-y*y*(2*v**4+4*v**3+3*v*v-2),1-v*v-4*v**3-2*v**4]])
    zero('fixed p-graph frame has determinant one',D.det()-1)
    zero('fixed frame inverse left and right',D*D.adjugate()-S.eye(2))
    td,wd=D*S.Matrix([q,r])
    Ap=A.subs(p,pstar)
    Bp=B.subs(p,pstar)
    zero('fixed p-graph original A equals x times the frame coordinate',Ap-x*td)
    zero('fixed p-graph original B is the normalized row',Bp-y*td-2*x*wd)
    qr=D.adjugate()*S.Matrix([T,W-L*x*T*T/(2*N)])
    frame_after=D*qr
    zero('polynomial frame and shear recover t',frame_after[0]-T)
    zero('polynomial frame and inverse shear recover W',frame_after[1]+L*x*frame_after[0]**2/(2*N)-W)
    old_frame=D*S.Matrix([q,r])
    recovered_new=S.Matrix([old_frame[0],old_frame[1]+L*x*old_frame[0]**2/(2*N)])
    zero('inverse shear followed by frame inverse recovers q,r',
         qr.subs({T:recovered_new[0],W:recovered_new[1]},simultaneous=True)-S.Matrix([q,r]))
    zero('actual fixed-graph potential cancels its entire quadratic part',
         (L*Ap**2+M*Ap+N*Bp).subs({q:qr[0],r:qr[1]},simultaneous=True)
         -((M*x+N*y)*T+2*N*x*W))

    c11,c12,c21,c22=S.symbols('c11 c12 c21 c22')
    C=S.Matrix([[c11,c12],[c21,c22]])
    linear=C*S.Matrix([P,Q])
    general_qr=qr.subs({T:linear[0],W:linear[1]},simultaneous=True)
    zero('full retained-dual Jacobian equals the input matrix determinant',
         general_qr.jacobian((P,Q)).det()-C.det())
    row0=S.Matrix([[M*x+N*y,2*N*x]])
    zero('every matrix gives the declared row potential on the same graph',
         ((M*x+N*y)*T+2*N*x*W).subs({T:linear[0],W:linear[1]},simultaneous=True)
         -(row0*C*S.Matrix([P,Q]))[0])
    H0=S.Matrix([[1,0],[M/(2*N),1]])
    zero('row orientation normalization is determinant one',H0.det()-1)
    zero('fixed graph row is the retained oriented row orbit',
         row0-N*S.Matrix([[y,2*x]])*H0)
    zero('row0 ideal explicitly contains x',row0[1]/(2*N)-x)
    zero('row0 ideal explicitly contains y',row0[0]/N-M*row0[1]/(2*N*N)-y)
    zero('fixed graph row Jacobian orientation',S.Matrix(list(row0)).jacobian((x,y)).det()+2*N*N)
    G1,G2=S.Function('G1')(x,y),S.Function('G2')(x,y)
    generic_cotangent=G1*P+G2*Q
    JG=S.Matrix([G1,G2]).jacobian((x,y))
    zero('realized row cotangent Hessian is exactly its plane Jacobian squared',
         S.hessian(generic_cotangent,(x,y,P,Q)).det()-JG.det()**2)

    # An unbounded triangular-row family has exact inverses; x^3 is a
    # polynomial specialization, not a coefficient-box search.
    triangular_row=(M*x+N*y+2*N*x**4,2*N*x)
    f_tri=triangular_row[0]*P+triangular_row[1]*Q
    grad_tri=S.Matrix([S.diff(f_tri,z) for z in (x,y,P,Q)])
    xr=z4/(2*N)
    Pr=z2/N
    yr=(z3-M*xr-2*N*xr**4)/N
    Qr=(z1-(M+8*N*xr**3)*Pr)/(2*N)
    inv_tri=S.Matrix([xr,yr,Pr,Qr])
    zero('fixed graph row-family Hessian constant',S.hessian(f_tri,(x,y,P,Q)).det()-4*N**4)
    zero('fixed graph nonlinear row gradient inverse left',
         inv_tri.subs(dict(zip((z1,z2,z3,z4),grad_tri)),simultaneous=True)-S.Matrix([x,y,P,Q]))
    zero('fixed graph nonlinear row gradient inverse right',
         grad_tri.subs(dict(zip((x,y,P,Q),inv_tri)),simultaneous=True)-S.Matrix([z1,z2,z3,z4]))

    # The fixed graph contains the marked points, but C=I loses their
    # collision. In general transfer requires row-value matching.
    common=S.Matrix([0,0,-S.Rational(1,2),0,0])
    embedded_p=pstar.subs({q:qr[0],r:qr[1]},simultaneous=True)
    embedding=S.Matrix([x,y,embedded_p,qr[0],qr[1]])
    pull=embedding.jacobian((x,y,T,W)).T*common
    ambient_points=[{x:1,y:-S.Rational(3,2),T:0,W:0},
                    {x:-1,y:S.Rational(3,2),T:0,W:0}]
    for index,point in enumerate(ambient_points):
        zero('fixed graph contains marked point '+str(index),
             embedding.subs(point).subs({L:1,M:13,N:2})-S.Matrix([point[x],point[y],0,0,0]))
        expected=S.Matrix([0,0,row0[0],row0[1]]).subs(point).subs({L:1,M:13,N:2})
        zero('marked covector agrees with the coefficient row '+str(index),
             pull.subs(point).subs({L:1,M:13,N:2})-expected)

    pins=json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    dependencies={}
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('pinned Tang dependency '+name,dependencies[name]==pins[name])
    dependencies['hessian/hc4-relative-row-orbit/PROOF.md']=hashlib.sha256((ROOT/'hessian/hc4-relative-row-orbit/PROOF.md').read_bytes()).hexdigest()
    result=dict(outcome='PASS',check_count=len(checks),checks=checks,
        q_chart_scope='Every standard q graph; arbitrary relative affine retained-dual changes; polynomial base automorphisms with affine old y coordinate',
        native_dual_graph_charts_closed=['p','q','r'],
        all_nonlinear_q_base_or_dual_changes_closed=False,
        p_chart_scope='One fixed degree-six graph with genuine polynomial dual inverse realizes every oriented singleton-zero row; full nonlinear-dual p invertibility is equivalent to the existing singleton-zero Keller problem',
        p_graph_embedding_and_dual_inverse_explicit=True,
        singleton_zero_properness_resolved=False,
        elementary_q_plane_proof_uses_Rentschler=False,
        written_not_finite_exhaustion=['square-degree lemma','three projection ranks','singleton ideal and descent','complete row orbit equivalence'],
        dependencies=dependencies,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name,'Q_GRAPH_SQUARE_AXIS.md','P_GRAPH_SINGLETON_ENDPOINT.md']},
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start,3),formal_verification=False,external_review=False)
    rendered=json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
