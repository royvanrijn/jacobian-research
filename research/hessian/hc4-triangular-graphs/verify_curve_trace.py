"""Exact identities and controls for CURVE_TRACE.md.

Rational-pole nilpotence and all projection ranks are written proofs.
Finite examples do not establish the unbounded theorem or its imports.
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
    started = time.monotonic()
    checks = []

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def zero(name, expression):
        values = list(expression) if isinstance(expression, S.MatrixBase) else [expression]
        check(name, all(S.cancel(value.doit()) == 0 for value in values))

    x, y, p, q, r, z, eps, s, t, P, Q, U, V, T = S.symbols('x y p q r z eps s t P Q U V T')
    L, N = S.symbols('L N', nonzero=True)
    M = S.Symbol('M')
    u = 1+x*y
    A = u**3*p+3*x*u*u*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = L*A*A+M*A+N*B
    zero('universal elimination identity in all five ambient variables',
         x*x*B+(2-3*x*y)*A-u*(u+1)*p-2*x*(2*u+1)*q)
    Aq, Bq = S.diff(A, q), S.diff(B, q)
    v = x*y
    Cq, Fq, Gq = 9*v*v+12*v+1, -27*v/2, (3*v+2)/2
    zero('omitted-q quadratic and linear coefficients',
         S.Matrix([Aq-3*x*u*u, Bq-y*Cq]))
    zero('omitted-q coefficient Bezout identity', Fq*u*u+Gq*Cq-1)
    v_certificate = Fq*y*Aq/3+Gq*x*Bq
    zero('omitted-q coefficient ideal contains xy', v_certificate-v)
    zero('omitted-q coefficient ideal contains x', Aq/3-x*(2+v)*v_certificate-x)
    zero('omitted-q coefficient ideal contains y', Bq-y*(12+9*v)*v_certificate-y)
    trace_p = L*z**6*V*V+(5*N*z-M*z**3)*V+2*N*Q/z
    zero('omitted-p hyperbola trace ignores the arbitrary p graph',
         psi.subs({x: z, y: -1/z, q: Q, r: V}, simultaneous=True)-trace_p)
    zero('hyperbola coordinate trace has nonzero quadratic dual Hessian',
         S.diff(trace_p, V, 2)-2*L*z**6)
    zero('hyperbola trace is a coordinate over K(z)',
         S.Matrix([V, trace_p*z/(2*N)]).jacobian((Q, V)).det()+1)
    hout = S.Symbol('hout')
    inverse_q = z*(hout-L*z**6*V*V-(5*N*z-M*z**3)*V)/(2*N)
    zero('hyperbola trace shear has an explicit polynomial-in-duals inverse',
         inverse_q.subs(hout, trace_p)-Q)

    # A base automorphism with both outputs nonlinear and an actual inverse.
    Xbase = s+t*t
    Ybase = t+Xbase*Xbase
    Sinv = x-(y-x*x)**2
    Tinv = y-x*x
    zero('nonlinear base automorphism determinant',
         S.Matrix([Xbase, Ybase]).jacobian((s, t)).det()-1)
    zero('nonlinear base inverse left composition',
         S.Matrix([Sinv, Tinv]).subs({x: Xbase, y: Ybase}, simultaneous=True)-S.Matrix([s, t]))
    zero('nonlinear base inverse right composition',
         S.Matrix([Xbase, Ybase]).subs({s: Sinv, t: Tinv}, simultaneous=True)-S.Matrix([x, y]))
    gamma_r = S.Matrix([Sinv, Tinv]).subs({x: 0, y: z}, simultaneous=True)
    zero('standard-r curve after nonlinear base change', gamma_r-S.Matrix([-z*z, z]))
    zero('standard-r curve has the required old base coordinates',
         S.Matrix([Xbase, Ybase]).subs(dict(zip((s, t), gamma_r)), simultaneous=True)-S.Matrix([0, z]))
    gamma_p = S.Matrix([Sinv, Tinv]).subs({x: z, y: -1/z}, simultaneous=True)
    zero('omitted-p rational curve after nonlinear base change',
         S.Matrix([Xbase, Ybase]).subs(dict(zip((s, t), gamma_p)), simultaneous=True)-S.Matrix([z, -1/z]))
    ca, cb, cc = S.symbols('ca cb cc')
    line = S.Poly(S.expand(z*z*(ca*gamma_p[0]+cb*gamma_p[1]+cc)), z)
    zero('nonlinear-base hyperbola control first affine-line obstruction', line.coeff_monomial(z**6)+ca)
    zero('nonlinear-base hyperbola control second affine-line obstruction',
         line.coeff_monomial(z**4)+cb)

    # Genuinely nonlinear dual coordinates along the standard-r trace.
    qz = Q+z*P*P
    pz = P+z*qz*qz
    hz = L*pz*pz+(M+4*N*z*z)*pz+N*z*qz
    zero('specialized nonlinear dual map Jacobian',
         S.Matrix([pz, qz]).jacobian((P, Q)).det()-1)
    zero('r-chart origin still has the nonlinear square trace', hz.subs(z, 0)-(L*P*P+M*P))
    po, ho = S.symbols('po ho')
    qr = ho-L*po*po/(N*z)-(M+4*N*z*z)*po/(N*z)
    Pr = po-z*qr*qr
    Qr = qr-z*Pr*Pr
    zero('nonlinear trace coordinate pair inverse right composition',
         S.Matrix([Pr, Qr]).subs({po: pz, ho: hz/(N*z)}, simultaneous=True)-S.Matrix([P, Q]))

    # All three base-projection normal forms retain the exact cotangent
    # determinant, regardless of the free bottom-right Hessian block.
    br, ar, cr = (S.Function(name)(T, U) for name in ('br', 'ar', 'cr'))
    rank_one = br+s*ar+V*cr
    J_one = S.Matrix([ar, cr]).jacobian((T, U))
    zero('rank-one full Hessian determinant',
         S.hessian(rank_one, (s, V, T, U)).det()-J_one.det()**2)
    Cscalar = S.Function('Cscalar')(T)
    rank_one_triangular = br+s*ar+V*Cscalar
    zero('rank-one triangular endpoint determinant',
         S.hessian(rank_one_triangular, (s, V, T, U)).det()
         -(S.diff(ar, U)*S.diff(Cscalar, T))**2)
    Bdual = S.Function('Bdual')(z, U)
    Hdual = Bdual*V+S.Function('Hdual')(z, U)
    zero('coordinate Hamiltonian restricts to the scalar-coefficient derivation',
         S.diff(Hdual, U)*S.diff(U, V)-S.diff(Hdual, V)*S.diff(U, U)+Bdual)
    b2, a2, c2 = (S.Function(name)(U, V) for name in ('b2', 'a2', 'c2'))
    rank_two = b2+s*a2+t*c2
    J_two = S.Matrix([a2, c2]).jacobian((U, V))
    zero('rank-two full Hessian determinant',
         S.hessian(rank_two, (s, t, U, V)).det()-J_two.det()**2)
    b0, a0, c0 = (S.Function(name)(s, t) for name in ('b0', 'a0', 'c0'))
    rank_zero = b0+P*a0+Q*c0
    J_zero = S.Matrix([a0, c0]).jacobian((s, t))
    zero('rank-zero exact plane-Keller residual determinant',
         S.hessian(rank_zero, (P, Q, s, t)).det()-J_zero.det()**2)
    grad_base = S.Matrix([S.diff(rank_zero, s), S.diff(rank_zero, t)])
    fibre = J_zero.T.adjugate()*(grad_base-S.Matrix([S.diff(b0, s), S.diff(b0, t)]))
    zero('rank-zero fibre recovery uses only the plane Jacobian denominator',
         fibre-J_zero.det()*S.Matrix([P, Q]))

    # A coordinate family with rational poles at zero and infinity.
    # Verify specialization of iterates at each pole, including cancellation
    # of different pole orders. The arbitrary-order proof is in the note.
    endpoint_A = U*U+V
    endpoint_B = U
    sz, tz = z*z+1/z, z*z+z
    family = sz*endpoint_A+tz*endpoint_B
    bracket = lambda f, g: S.diff(f, U)*S.diff(g, V)-S.diff(f, V)*S.diff(g, U)
    zero('rational trace family is a coordinate over C(z)', bracket(U, family/sz)-1)
    zero('rational trace endpoint has a constant-Jacobian mate', bracket(endpoint_A, endpoint_B)+1)
    for pole, scaled, limit in [
        ('zero', eps*family.subs(z, eps), endpoint_A),
        ('infinity', eps*eps*family.subs(z, 1/eps), endpoint_A+endpoint_B)]:
        scaled = S.cancel(scaled)
        zero(pole+' scaled trace is regular with the correct pole endpoint', scaled.subs(eps, 0)-limit)
        for generator in (U, V):
            iterate, endpoint_iterate = generator, generator
            for order in range(1, 4):
                iterate = S.expand(bracket(scaled, iterate))
                endpoint_iterate = S.expand(bracket(limit, endpoint_iterate))
                zero(f'{pole} pole specialization of iterate {generator} order {order}',
                     iterate.subs(eps, 0)-endpoint_iterate)
            zero(pole+' pole coordinate Hamiltonian kills '+str(generator)+' at order three', iterate)

    # A nonlinear dual change can flatten just the hyperbola trace, while
    # leaving a nonzero quartic term in the actual full potential.
    shifted_q = Q-L*x**7*V*V/(2*N)
    pchart_control = psi.subs({p: 0, q: shifted_q, r: V}, simultaneous=True)
    zero('nonlinear p-chart dual change is a polynomial automorphism',
         S.Matrix([shifted_q, V]).jacobian((Q, V)).det()-1)
    zero('nonlinear dual change makes the hyperbola trace affine',
         pchart_control.subs({x: z, y: -1/z}, simultaneous=True)
         -((5*N*z-M*z**3)*V+2*N*Q/z))
    zero('affine curve trace is not a globally affine-dual potential',
         S.Poly(pchart_control, V).coeff_monomial(V**4)-9*L**3*x**16*u**4/(4*N*N))
    check('the full quartic dual coefficient is nonzero', 9*L**3*x**16*u**4/(4*N*N) != 0)

    # Marked points in the p chart need a graph-jet match as well as point
    # containment. This polynomial graph contains both but loses the collision.
    graph_p = x*x-1
    embedding = S.Matrix([x, y, graph_p, Q, V])
    common_ambient_covector = S.Matrix([0, 0, -S.Rational(1, 2), 0, 0])
    pulled = embedding.jacobian((x, y, Q, V)).T*common_ambient_covector
    plus = {x: 1, y: -S.Rational(3, 2), Q: 0, V: 0}
    minus = {x: -1, y: S.Rational(3, 2), Q: 0, V: 0}
    ambient_gradient = S.Matrix([S.diff(psi, v) for v in (x, y, p, q, r)])
    for label, point in [('plus', plus), ('minus', minus)]:
        zero('actual v2 ambient common gradient at '+label+' marked point',
             ambient_gradient.subs(point).subs({L: 1, M: 13, N: 2, p: 0, q: 0, r: 0})
             -common_ambient_covector)
    zero('p-chart marked containment control at both points',
         S.Matrix([graph_p.subs(plus), graph_p.subs(minus)]))
    zero('p-chart marked pullback is minus one half the graph gradient',
         pulled+S.Matrix([S.diff(graph_p, v) for v in (x, y, Q, V)])/2)
    check('p-chart point containment alone loses the marked gradient collision', pulled.subs(plus) != pulled.subs(minus))

    manifest = json.loads((HERE/'sources/automorphic-trace-import.json').read_text())
    source_pdf = HERE/'sources'/manifest['retained_file']
    check('retained Rentschler source matches its reviewed byte pin',
          hashlib.sha256(source_pdf.read_bytes()).hexdigest()==manifest['sha256'])
    source_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'hessian/tang-hc4-audit/source-v1.pdf']:
        check('pinned Tang structural input '+name,
              hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==source_pins[name])
    result = dict(outcome='PASS', check_count=len(checks), checks=checks,
        curve_trace_criterion='Generic coordinate on a nonconstant rational base curve; special nonlinear square gives full closure, while an affinely spanning curve gives invertibility or the exact globally affine-dual endpoint',
        standard_r_chart_all_base_and_dual_automorphisms_covered=True,
        standard_p_chart_all_graph_degrees_covered=True,
        p_chart_all_base_and_relative_affine_dual_automorphisms_covered=True,
        p_chart_nonlinear_dual_changes='Invertible gradient or globally affine-dual plane-Keller endpoint; not unrestricted closure',
        other_graph_chart_q_closed=False,
        omitted_q_graph_independent_curve_locus='Exactly the origin: (A_q,B_q)=(x,y)',
        arbitrary_nonlinear_base_dual_mixing_covered=False,
        constant_affine_parameter_changes_covered=True,
        word_length_or_degree_bound=None,
        written_proof_not_finite_exhaustion=['all projection ranks', 'pole specialization at arbitrary nilpotence orders',
            'affine-span argument for every base automorphism', 'characteristic-zero descent'],
        external_inputs=['HCTANG2 Tang alternative', 'Rentschler endpoint retained in HC4MYGV0'],
        source_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name, 'CURVE_TRACE.md']},
        runtime=dict(python=platform.python_version(), sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started, 3),
        formal_verification=False, external_review=False)
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__=='__main__':
    main()
