"""Exact algebra for FIRST_SLICE_INVERTIBILITY.md; no coefficient search.

Direction exhaustion, polynomial degree and arbitrary word lengths are
written arguments. Tang's structural theorem is a separate pinned import.
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

    x, y, P, Q, ell = S.symbols('x y P Q ell')
    L, N, a, k, nu, alpha, kappa = S.symbols(
        'L N a k nu alpha kappa', nonzero=True)
    M, d0, eta, tau, sigma = S.symbols('M d0 eta tau sigma')
    u, v = 1+x*y, x*y
    p0, q0 = S.Function('p0')(P, Q), S.Function('q0')(P, Q)
    E, H, R = (S.Function(name)(x, y, P, Q) for name in ('E', 'H', 'R'))
    p, q = p0+x*E, q0+x*H
    A = u**3*p+3*x*u*u*q-x**3*R
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*R
    psi = L*A*A+M*A+N*B
    first_slice = L*p0*p0+M*p0+4*N*y*y*p0+N*y*q0
    zero('arbitrary x-divisible corrections have the required complete trace',
         psi.subs(x, 0)-first_slice)

    pa, qa = S.Function('pa')(ell), S.Function('qa')(ell)
    slice_jet = L*pa*pa+M*pa+4*N*y*y*pa+N*y*qa
    D = lambda f: S.diff(f, y)+S.diff(f, ell)
    coefficients = S.Poly(S.expand(D(D(slice_jet))), y)
    zero('moving-y direction y^2 coefficient',
         coefficients.coeff_monomial(y*y)-4*N*S.diff(pa, ell, 2))
    zero('moving-y direction y coefficient',
         coefficients.coeff_monomial(y)-16*N*S.diff(pa, ell)-N*S.diff(qa, ell, 2))
    reduced = 2*L*S.diff(pa, ell)**2+8*N*pa+2*N*S.diff(qa, ell)
    zero('moving-y direction reduced constant coefficient',
         coefficients.coeff_monomial(1).subs(S.diff(pa, ell, 2), 0)-reduced)
    zero('moving-y direction compatibility',
         S.diff(reduced, ell).subs({S.diff(pa, ell, 2): 0,
             S.diff(qa, ell, 2): -16*S.diff(pa, ell)}, simultaneous=True)
         +24*N*S.diff(pa, ell))
    j11, j12, j21, j22, c, d = S.symbols('j11 j12 j21 j22 c d')
    J0 = S.Matrix([[j11, j12], [j21, j22]])
    zero('polynomial Jacobian matrix preserves directional component ideal',
         J0.adjugate()*J0*S.Matrix([c, d])-J0.det()*S.Matrix([c, d]))
    X, Y, T, Z = S.symbols('X Y T Z')
    tracep, traceq = a*T+d0, k*Z+S.Function('h')(T)
    zero('adapted trace determinant',
         S.Matrix([tracep, traceq]).jacobian((T, Z)).det()-a*k)
    adapted_slice = L*tracep**2+M*tracep+4*N*y*y*tracep+N*y*traceq
    zero('adapted fibre coefficient initial value', S.diff(adapted_slice, Z)-N*k*y)
    C = N*k*(y-eta*x)
    zero('linear transported coefficient satisfies flow equation',
         S.diff(C, x)+eta*S.diff(C, y)+tau*S.diff(C, T))
    zero('linear transported coefficient has exact initial value', C.subs(x, 0)-N*k*y)
    transformation = S.Matrix([x, y-eta*x, P-tau*x, Q-sigma*x])
    zero('adapted linear source change determinant',
         transformation.jacobian((x, y, P, Q)).det()-1)
    zero('adapted source change sends plane generator to X direction',
         transformation.jacobian((x, y, P, Q))*S.Matrix([1, eta, tau, sigma])
         -S.Matrix([1, 0, 0, 0]))

    # A completely arbitrary cotangent remainder and fibre coefficient.
    U, V = S.Function('U')(Y, T), S.Function('V')(Y, T)
    normal = nu*Y*Z+X*U+V
    variables = (X, Z, Y, T)
    hess = S.hessian(normal, variables)
    zero('whole affine fibre block vanishes', hess[:2, :2])
    zero('complete off-diagonal block',
         hess[:2, 2:]-S.Matrix([[S.diff(U, Y), S.diff(U, T)], [nu, 0]]))
    zero('unit-Hessian triangular gate with arbitrary remainder',
         hess.det()-nu**2*S.diff(U, T)**2)
    beta = S.Function('beta')(Y)
    triangular = nu*Y*Z+X*(alpha*T+beta)+V
    grad = S.Matrix([S.diff(triangular, z) for z in variables])
    outputs = S.symbols('z_X z_S z_Y z_T')
    zx, zs, zy, zt = outputs
    Yr = zs/nu
    Tr = (zx-beta.subs(Y, Yr))/alpha
    evaluate = lambda value: value.subs({Y: Yr, T: Tr}, simultaneous=True).doit()
    Xr = (zt-evaluate(S.diff(V, T)))/alpha
    Zr = (zy-Xr*evaluate(S.diff(beta, Y))-evaluate(S.diff(V, Y)))/nu
    inverse = S.Matrix([Xr, Zr, Yr, Tr])
    zero('gradient after inverse is the identity with arbitrary beta,V',
         grad.subs(dict(zip(variables, inverse)), simultaneous=True)-S.Matrix(outputs))
    zero('inverse after gradient is the identity with arbitrary beta,V',
         inverse.subs(dict(zip(outputs, grad)), simultaneous=True)-S.Matrix(variables))
    zero('triangular normal form has nonzero constant determinant',
         S.hessian(triangular, variables).det()-nu**2*alpha**2)
    # A nonlinear polynomial specialization independently replays the two
    # inverse compositions without any unevaluated-function substitutions.
    fixed = triangular.subs({beta: Y**5-Y**2+2,
        V: Y**3*T**4+Y*T+T**6}, simultaneous=True).doit()
    fixed_grad = S.Matrix([S.diff(fixed, z) for z in variables])
    Yf = zs/nu
    Tf = (zx-(Yf**5-Yf**2+2))/alpha
    Xf = (zt-(4*Yf**3*Tf**3+Yf+6*Tf**5))/alpha
    Zf = (zy-Xf*(5*Yf**4-2*Yf)-(3*Yf**2*Tf**4+Tf))/nu
    fixed_inverse = S.Matrix([Xf, Zf, Yf, Tf])
    zero('nonlinear specialization inverse after gradient',
         fixed_inverse.subs(dict(zip(outputs, fixed_grad)), simultaneous=True)-S.Matrix(variables))

    # Arbitrary-length closure is the written ideal argument; this actual
    # nontriangular word checks the larger x-only sufficient condition.
    qw = Q+x*P**2
    pw = P+x*qw**2
    zero('x-only alternating word dual Jacobian',
         S.Matrix([pw, qw]).jacobian((P, Q)).det()-1)
    zero('x-only alternating word first trace identity',
         S.Matrix([pw, qw]).subs(x, 0)-S.Matrix([P, Q]))
    check('x-only word can fail both matching second traces',
          S.expand(pw.subs(y, -1/x)-P) != 0 and S.expand(qw.subs(y, -1/x)-Q) != 0)
    badp, badq = P, Q-4*y*P
    zero('y-dependent trace failure control is still Keller over K[y]',
         S.Matrix([badp, badq]).jacobian((P, Q)).det()-1)
    badslice = L*badp**2+M*badp+4*N*y*y*badp+N*y*badq
    zero('y-dependent trace has a moving-y affine direction', S.diff(badslice, y, 2))
    check('y-dependent trace is outside the theorem', S.diff(badq, y) != 0)
    paux, qaux, faux = S.symbols('paux qaux faux')
    trace_shear = L*paux*paux/(N*y)+(M+4*N*y*y)*paux/(N*y)+qaux
    zero('genuine first-trace automorphism gives a generic-coordinate shear',
         S.Matrix([paux, trace_shear]).jacobian((paux, qaux)).det()-1)
    zero('generic-coordinate trace shear has an explicit inverse over K(y)',
         (faux/(N*y)-L*paux*paux/(N*y)-(M+4*N*y*y)*paux/(N*y))
         .subs(faux, N*y*trace_shear)-qaux)

    # Construct an actual polynomial unit-Hessian member of the enlarged
    # Meng--Yang class using the retained HC5T1 unit-pivot frame.
    pbez = 1-3*v+6*v*v
    qbez = y**3*(6*v*v+15*v+10)
    zero('HC5T1 polynomial Bezout identity', pbez*u**3-x**3*qbez-1)
    C0 = S.Matrix([
        [(3*v*v+6*v+2)/2, x*x*(2*v+3)],
        [S.Rational(3, 2)*y*y*(18*v**4+63*v**3+69*v*v+21*v-1),
         36*v**5+108*v**4+87*v**3+3*v*v-3*v+1]])
    zero('retained C0 normalizer determinant', C0.det()-1)
    gamma = S.Matrix([[-y*(18*v**5+81*v**4+120*v**3+60*v*v-1), x*(v+1)*(v+2)]])
    zero('retained C0 normalizes the entire residual B row',
         gamma*C0-S.Matrix([[y, 2*x]]))
    z = Q-y*(6*v*v+15*v)*P
    w = kappa*P/(2*N)
    rho, sig = C0*S.Matrix([z, w])
    pc = pbez*(P-3*x*u*u*rho)+x**3*sig
    qc = rho
    rc = qbez*(P-3*x*u*u*rho)+u**3*sig
    Ac = u**3*pc+3*x*u*u*qc-x**3*rc
    Bc = y*y*u*(4+3*x*y)*pc+(y+3*x*y*y*(4+3*x*y))*qc+(2*x-3*x*x*y)*rc
    zero('actual polynomial member has A=P', Ac-P)
    zero('actual polynomial member has the target B',
         Bc-(4*y*y*P+y*Q+kappa*x*P/N))
    zero('actual polynomial member first p,q traces',
         S.Matrix([pc, qc]).subs(x, 0)-S.Matrix([P, Q]))
    for name, value, trace in [('p', pc, P), ('q', qc, Q)]:
        check('actual '+name+' correction is polynomial and divisible by x',
              S.Poly(S.cancel((value-trace)/x), x, y, P, Q).is_zero is not True)
    zero('actual polynomial member fails the old second q trace',
         S.diff(qc.subs(y, -1/x), Q)+S.Rational(1, 2))
    F = L*P**2+M*P+4*N*y*y*P+N*y*Q+kappa*x*P
    zero('actual member Hessian determinant',
         S.hessian(F, (x, Q, y, P)).det()-(N*kappa)**2)
    zero('actual member has a genuine affine fibre plane',
         S.Matrix([S.diff(F, x, 2), S.diff(F, x, Q), S.diff(F, Q, 2)]))
    # Recover the two parameter coordinates through the invertible frames.
    recovered_t = Ac
    recovered_sig = -qbez*pc+pbez*rc
    recovered_z, recovered_w = C0.adjugate()*S.Matrix([qc, recovered_sig])
    zero('actual embedding recovers P polynomially', recovered_t-P)
    zero('actual embedding recovers Q polynomially',
         recovered_z+y*(6*v*v+15*v)*recovered_t-Q)
    zero('actual embedding image satisfies its transformed graph equation',
         recovered_w-kappa*recovered_t/(2*N))

    # Sharp boundary: allow a y-dependent first-trace Keller pair, and
    # every plane Keller cotangent lift has an explicit polynomial lift.
    Fmap, Gmap = S.Function('Fmap')(P, Q), S.Function('Gmap')(P, Q)
    zmap = Gmap-y*(6*v*v+15*v+4)*Fmap
    rhomap, sigmap = C0*S.Matrix([zmap, Fmap/2])
    pmap = pbez*(Fmap-3*x*u*u*rhomap)+x**3*sigmap
    qmap = rhomap
    rmap = qbez*(Fmap-3*x*u*u*rhomap)+u**3*sigmap
    Amap = u**3*pmap+3*x*u*u*qmap-x**3*rmap
    Bmap = y*y*u*(4+3*x*y)*pmap+(y+3*x*y*y*(4+3*x*y))*qmap+(2*x-3*x*x*y)*rmap
    zero('arbitrary plane-map encoding A identity', Amap-Fmap)
    zero('arbitrary plane-map encoding B identity', Bmap-x*Fmap-y*Gmap)
    zero('arbitrary plane-map encoding first traces',
         S.Matrix([pmap, qmap]).subs(x, 0)-S.Matrix([Fmap, Gmap-4*y*Fmap]))
    Jplane = S.Matrix([Fmap, Gmap]).jacobian((P, Q))
    zero('encoded y-dependent trace has the same plane Jacobian',
         S.Matrix([Fmap, Gmap-4*y*Fmap]).jacobian((P, Q)).det()-Jplane.det())
    zero('encoded full dual Jacobian is the declared nonunit multiple',
         S.Matrix([pmap, qmap]).jacobian((P, Q)).det()
         -Jplane.det()*(pbez*C0[0, 0]+x**3/2))
    check('encoded full dual Jacobian multiplier is nonconstant',
          S.diff(pbez*C0[0, 0]+x**3/2, x) != 0)
    hmap = L*Fmap**2+M*Fmap
    fmap = hmap+N*(x*Fmap+y*Gmap)
    zero('encoded Hessian determinant is N^4 times plane Jacobian squared',
         S.hessian(fmap, (x, y, P, Q)).det()-N**4*Jplane.det()**2)
    gradmap = S.Matrix([S.diff(fmap, z) for z in (x, y, P, Q)])
    source_fibres = S.Matrix([S.diff(hmap, P), S.diff(hmap, Q)])+N*Jplane.T*S.Matrix([x, y])
    zero('encoded gradient factors as plane map times identity after fibre change',
         gradmap-S.Matrix([N*Fmap, N*Gmap, source_fibres[0], source_fibres[1]]))
    zero('fibre source change inverse uses only the plane Jacobian denominator',
         Jplane.T.adjugate()*(source_fibres-S.Matrix([S.diff(hmap, P), S.diff(hmap, Q)]))
         -N*Jplane.det()*S.Matrix([x, y]))

    pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    dependencies = {}
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('pinned structural dependency '+name, dependencies[name] == pins[name])
    result = dict(outcome='PASS', check_count=len(checks), checks=checks,
        elementary_scope='Any polynomial with the stated complete first slice: every affine fibre plane yields a triangular cotangent endpoint and an explicit gradient inverse if the Hessian is a unit',
        global_scope='p=p0+x*E, q=q0+x*H, r=R with base-independent trace Keller pair; no full dual Jacobian or second-divisor hypothesis',
        all_first_trace_preserving_words_covered=True,
        arbitrary_y_dependent_first_traces_covered=False,
        second_trace_gate_needed_for_gradient_invertibility=False,
        no_plane_claim_for_enlarged_class=False,
        actual_polynomial_unit_Hessian_control=True,
        y_dependent_trace_class_contains_every_plane_Keller_cotangent_lift=True,
        encoding_is_a_full_dual_automorphism=False,
        encoding_asserted_to_be_a_closed_embedding=False,
        global_invertibility_import='HCTANG2; no plane degree bound or assumption of JC2',
        source_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name, 'FIRST_SLICE_INVERTIBILITY.md']},
        dependency_source_sha256=dependencies,
        runtime=dict(python=platform.python_version(), sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start, 3),
        formal_verification=False, external_review=False)
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
