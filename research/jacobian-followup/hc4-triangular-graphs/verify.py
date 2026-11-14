"""Universal identities supporting the written triangular-graph theorem.

The degree arguments and two-plane exhaustion are written in PROOF.md.
No coefficient search or sampled-directions census is used.
"""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
checks = []


def check(name, condition):
    assert condition, name
    checks.append(name)


def zero(name, expression):
    values = list(expression) if isinstance(expression,S.MatrixBase) else [expression]
    check(name,all(S.cancel(value)==0 for value in values))


def main():
    start = time.monotonic()
    x,y,P,Q,p,q,r,L,M,N,a,b,s,t,c,d,h,eta,tau,T,V = S.symbols(
        'x y P Q p q r L M N a b s t c d h eta tau T V')
    variables = (x,y,P,Q)
    u = 1+x*y
    A = u**3*p+3*x*u**2*q-x**3*r
    B = y**2*u*(4+3*x*y)*p+(y+3*x*y**2*(4+3*x*y))*q+(2*x-3*x**2*y)*r
    psi = L*A*A+M*A+N*B
    K0 = M-N*(2-3*x*y)/x**2
    alpha = N*u*(u+1)/x**2
    gamma = 2*N*(2*u+1)/x
    U = A+K0/(2*L)
    zero('complete universal elimination identity',x*x*B+(2-3*x*y)*A-u*(u+1)*p-2*x*(2*u+1)*q)
    zero('universal completed square',psi-(L*U*U+alpha*p+gamma*q-K0*K0/(4*L)))
    zero('Gamma is affine in y',S.diff(gamma,y,2))
    zero('alpha is quadratic in y',S.diff(alpha,y,3))
    check('alpha y-second derivative is 2N',S.cancel(S.diff(alpha,y,2))==2*N)
    G = S.Function('G')(x,y)
    H1 = S.Function('H1')(x,y,P)
    H2 = S.Function('H2')(x,y,Q)
    R = S.Function('R')(x,y,P,Q)
    D = lambda expression:S.diff(expression,y)+s*S.diff(expression,P)+t*S.diff(expression,Q)
    p1,q1 = a*P+G,b*Q+H1
    p2,q2 = a*P+H2,b*Q+G
    phi1,phi2 = S.Matrix([x,y,p1,q1]),S.Matrix([x,y,p2,q2])
    zero('first triangular dual Jacobian',phi1.jacobian(variables).det()-a*b)
    zero('opposite triangular dual Jacobian',phi2.jacobian(variables).det()-a*b)
    A1 = A.subs({p:p1,q:q1,r:R},simultaneous=True)
    A2 = A.subs({p:p2,q:q2,r:R},simultaneous=True)
    U1,U2 = A1+K0/(2*L),A2+K0/(2*L)
    E1 = alpha*p1+gamma*H1-K0*K0/(4*L)
    E2 = alpha*H2+gamma*q2-K0*K0/(4*L)
    zero('first completed-square remainder independent of Q',S.diff(E1,Q))
    zero('opposite completed-square remainder independent of P',S.diff(E2,P))
    zero('first immutable Q derivative',S.diff(U1,Q)-(3*b*x*u*u-x**3*S.diff(R,Q)))
    zero('first immutable mixed jet under every s,t',D(S.diff(U1,Q))-(6*b*x*x*u-x**3*D(S.diff(R,Q))))
    check('first mixed-jet obstruction residue is 6b',(6*b*u).subs(x,0)==6*b)
    zero('opposite immutable P derivative',S.diff(U2,P)-(a*u**3-x**3*S.diff(R,P)))
    zero('opposite immutable mixed jet under every s,t',D(D(S.diff(U2,P)))-(6*a*x*x*u-x**3*D(D(S.diff(R,P)))))
    check('opposite mixed-jet obstruction residue is 6a',(6*a*u).subs(x,0)==6*a)
    G00 = S.Symbol('G00')
    first_slice = L*(a*P+G00)**2+M*(a*P+G00)
    dual_D = lambda expression:c*S.diff(expression,P)+d*S.diff(expression,Q)
    zero('first pure-dual slice forces c zero',dual_D(dual_D(first_slice))-2*L*a*a*c*c)
    zero('first Q-degree derivative obstruction',S.diff(A1,Q)/x-(3*b*u*u-x*x*S.diff(R,Q)))
    h0,g0 = S.Function('h0')(y),S.Function('g0')(y)
    p0,q0 = a*P+h*Q+h0,b*Q+g0
    f0 = L*p0*p0+(M+4*N*y*y)*p0+N*y*q0
    trace = f0.subs({P:(T-h*V)/a,Q:V},simultaneous=True)
    zero('opposite exceptional-trace coefficient is Nby',S.diff(trace,V)-N*b*y)
    zero('opposite exceptional-trace affine direction',S.diff(trace,V,2))
    transported = N*b*(y-eta*x)
    zero('transported coefficient satisfies constant-flow equation',
         S.diff(transported,x)+eta*S.diff(transported,y)+tau*S.diff(transported,T))
    zero('transported coefficient has exact x=0 data',transported.subs(x,0)-N*b*y)
    q_special = b*V+S.Function('g')(x)
    restricted = psi.subs({y:-1/x,q:q_special},simultaneous=True)
    A_special = -x**3*r
    zero('second-divisor complete potential identity',restricted-(L*A_special**2+(M-5*N/x**2)*A_special+2*N*q_special/x))
    zero('second-divisor coefficient discrepancy',
         x*(transported.subs(y,-1/x)-2*N*b/x)/(N*b)+(eta*x*x+3))
    check('discrepancy cannot be an identity for constant eta',S.Poly(eta*x*x+3,x).coeff_monomial(1)==3)
    # An exact positive control: the opposite order can have one affine
    # direction. The theorem must exclude planes, not all such directions.
    Hcontrol = 3*b*x*(x*y-1)*Q
    Rcontrol = 3*b*y*y*u*u*Q
    Acontrol = S.expand(A.subs({p:a*P+Hcontrol,q:b*Q,r:Rcontrol},simultaneous=True))
    psicontrol = psi.subs({p:a*P+Hcontrol,q:b*Q,r:Rcontrol},simultaneous=True)
    zero('positive control A independent of Q',Acontrol-a*u**3*P)
    zero('positive control has the claimed affine coefficient',S.diff(psicontrol,Q)-N*b*y*(3*u*u-2))
    zero('positive control genuinely affine in Q',S.diff(psicontrol,Q,2))
    check('positive control coefficient is not a transported linear form',
          S.Poly(S.expand(N*b*y*(3*u*u-2)),x,y).coeff_monomial(x*x*y**3)==3*N*b)
    # Determinant-one reparameterization does not preserve a unit Hessian.
    phi = S.Matrix([x,y,P,Q+P*P])
    J = phi.jacobian(variables)
    quadratic_pullback = sum(z*z for z in phi)/2
    chain = J.T*J+sum((phi[i]*S.hessian(phi[i],variables) for i in range(4)),S.zeros(4))
    zero('full Hessian chain rule in nonlinear control',S.hessian(quadratic_pullback,variables)-chain)
    check('nonlinear control parameterization has determinant one',J.det()==1)
    zero('nonlinear control Hessian is nonconstant',S.hessian(quadratic_pullback,variables).det()-(1+2*Q+2*P*P))
    ambient = S.Matrix([0,0,-S.Rational(1,2),0,0])
    for number,oldp,oldq,expected in [
        (1,p1,q1,S.Matrix([S.diff(G,x),S.diff(G,y),a,0])),
        (2,p2,q2,S.Matrix([S.diff(H2,x),S.diff(H2,y),a,S.diff(H2,Q)]))]:
        embedding = S.Matrix([x,y,oldp,oldq,R])
        zero(f'order {number} marked gradient pullback formula',embedding.jacobian(variables).T*ambient+expected/2)
    # All base-independent Keller substitutions: universal directional jets.
    ell = S.Symbol('ell')
    pa,qa = S.Function('pa')(ell),S.Function('qa')(ell)
    Djet = lambda expression:S.diff(expression,y)+S.diff(expression,ell)
    slice_jet = L*pa*pa+M*pa+4*N*y*y*pa+N*y*qa
    coefficients = S.Poly(S.expand(Djet(Djet(slice_jet))),y)
    zero('base-independent y^2 directional coefficient',
         coefficients.coeff_monomial(y*y)-4*N*S.diff(pa,ell,2))
    zero('base-independent y directional coefficient',
         coefficients.coeff_monomial(y)-16*N*S.diff(pa,ell)-N*S.diff(qa,ell,2))
    reduced_constant = 2*L*S.diff(pa,ell)**2+8*N*pa+2*N*S.diff(qa,ell)
    zero('base-independent constant coefficient after delta_second_p zero',
         coefficients.coeff_monomial(1).subs(S.diff(pa,ell,2),0)-reduced_constant)
    zero('base-independent compatibility forces -24N delta_p',
         S.diff(reduced_constant,ell).subs({S.diff(pa,ell,2):0,
              S.diff(qa,ell,2):-16*S.diff(pa,ell)},simultaneous=True)+24*N*S.diff(pa,ell))
    j11,j12,j21,j22 = S.symbols('j11 j12 j21 j22')
    determinant = j11*j22-j12*j21
    dp,dq = j11*c+j12*d,j21*c+j22*d
    zero('polynomial Jacobian ideal identity for c',j22*dp-j12*dq-determinant*c)
    zero('polynomial Jacobian ideal identity for d',-j21*dp+j11*dq-determinant*d)
    # Exact whole-chart regression for a nontriangular alternating word.
    q3 = Q+P*P
    p3 = P+q3*q3
    check('base-independent alternating word has Jacobian one',
          S.Matrix([p3,q3]).jacobian((P,Q)).det().expand()==1)
    vy = S.Symbol('v_y')
    Dcontrol = lambda expression:vy*S.diff(expression,y)+c*S.diff(expression,P)+d*S.diff(expression,Q)
    slice3 = (L*p3*p3+M*p3+4*N*y*y*p3+N*y*q3).subs({L:1,M:13,N:2})
    necessary = S.Poly(S.expand(Dcontrol(Dcontrol(slice3))),y,P,Q).coeffs()
    for coordinate,others in [(vy,(c,d)),(c,(vy,d)),(d,(vy,c))]:
        basis = S.groebner([f.subs(coordinate,1) for f in necessary],*others)
        check('alternating-word full projective affine-direction chart '+str(coordinate),
              len(basis.polys)==1 and basis.polys[0].as_expr()==1)
    embedding3 = S.Matrix([x,y,p3,q3,R])
    zero('base-independent marked covector formula',
         embedding3.jacobian(variables).T*ambient+
         S.Matrix([0,0,S.diff(p3,P),S.diff(p3,Q)])/2)
    # Pin the external theorem separately from the elementary algebra.
    pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    source_hashes = {}
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        source_hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('pinned structural dependency: '+name,source_hashes[name]==pins[name])
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),checks=checks,
        theorem_scope='One base-dependent triangular dual change in either order, or any base-independent plane Keller substitution; arbitrary polynomial degrees',
        elementary_proof='Written polynomial-degree, constant-plane and two-divisor arguments',
        global_invertibility_import='HCTANG2 only; no plane degree bound or assumption of JC2',
        base_independent_words_covered=True,base_dependent_alternating_words_covered=False,
        nonlinear_base_changes_covered=False,
        marked_collision='Requires equality of pulled-back p covectors, not point containment alone',
        independent_controls=['an opposite-order affine direction exists','a determinant-one coordinate change can destroy unit Hessian determinant'],
        dependency_source_sha256=source_hashes,elapsed_seconds=round(time.monotonic()-start,3),
        runtime={'python':platform.python_version(),'sympy':S.__version__},
        formal_verification=False,external_review=False),indent=2))


if __name__=='__main__':
    main()
