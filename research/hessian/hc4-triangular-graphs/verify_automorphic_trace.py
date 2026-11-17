"""Exact algebra and controls for AUTOMORPHIC_TRACE.md.

Written plane exhaustion and all-degree LND arguments are not certified by
finite controls. Tang and the Rentschler criterion are pinned imports.
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

    x, y, P, Q, U, V, X, Y, T, Z, lam = S.symbols('x y P Q U V X Y T Z lam')
    L, N, alpha, nu = S.symbols('L N alpha nu', nonzero=True)
    M, eta, tau, sigma, c, d, e, g = S.symbols('M eta tau sigma c d e g')
    bracket = lambda f, h: S.diff(f, U)*S.diff(h, V)-S.diff(f, V)*S.diff(h, U)

    # The first trace is a coordinate over the parameter field whenever
    # its dual pair has an actual inverse over K[y].
    pa, qa, haux = S.symbols('pa qa haux')
    htrace = L*pa*pa+(M+4*N*y*y)*pa+N*y*qa
    zero('trace shear has unit plane Jacobian over K(y)',
         S.Matrix([pa, htrace/(N*y)]).jacobian((pa, qa)).det()-1)
    zero('trace shear inverse is polynomial in the dual variables over K(y)',
         (haux/(N*y)-L*pa*pa/(N*y)-(M+4*N*y*y)*pa/(N*y))
         .subs(haux, htrace)-qa)

    # Noncommuting derivations: the top parameter coefficient is always
    # D_a^n, not an assumed commutative binomial formula.
    aa = U+(V+U*U)**3
    bb = U*U
    Da = lambda f: S.expand(bracket(aa, f))
    Db = lambda f: S.expand(bracket(bb, f))
    check('Hamiltonian control derivations do not commute', Db(Da(U))-Da(Db(U)) != 0)
    for generator in (U, V):
        combined, top = generator, generator
        for order in range(1, 5):
            combined = S.expand(Db(combined)+lam*Da(combined))
            top = Da(top)
            zero(f'noncommuting pencil top coefficient {generator} order {order}',
                 S.Poly(combined, lam).coeff_monomial(lam**order)-top)

    # A generic coordinate need not have a coordinate leading member.
    # The square has an LND Hamiltonian, but no constant-Jacobian mate.
    degenerate = U+lam*V*V
    zero('degenerating pencil is a literal coordinate with mate V',
         bracket(degenerate, V)-1)
    Ddeg = lambda f: bracket(degenerate, f)
    zero('degenerating coordinate pencil kills U at order three', Ddeg(Ddeg(Ddeg(U))))
    zero('degenerating coordinate pencil kills V at order two', Ddeg(Ddeg(V)))
    Hmate = S.Function('Hmate')(U, V)
    zero('square leading member cannot have a unit Jacobian mate', bracket(Hmate, V*V).subs(V, 0))
    H = S.Function('H')(T)
    Amate = S.Function('Amate')(T, Z)
    zero('Rentschler univariate factor splits the Jacobian-mate equation',
         S.Matrix([Amate, H]).jacobian((T, Z)).det()+S.diff(Amate, Z)*S.diff(H, T))

    C = S.Function('C')(y, T)
    B = S.Function('B')(y, T)
    pure_trace = C*Z+B
    zero('pure-dual trace Hamiltonian restricts to minus C times d/dT',
         S.diff(pure_trace, T)*S.diff(T, Z)-S.diff(pure_trace, Z)*S.diff(T, T)+C)
    # Leading-degree regressions for the written one-variable argument.
    for degree in (1, 2, 3):
        Ct = T**degree+1
        unary = lambda f: -Ct*S.diff(f, T)
        iterate = T
        for order in range(1, 5):
            iterate = S.expand(unary(iterate))
            check(f'nonconstant pure coefficient survives degree {degree} iteration {order}', iterate != 0)
            check(f'pure coefficient iterate degree {degree} iteration {order}',
                  S.degree(iterate, T)==1+order*(degree-1))
    scalar_C = S.Function('scalar_C')(y)
    zero('tangent whole-plane mixed derivative is C prime',
         S.diff(scalar_C*Z+B, Z, y)-S.diff(scalar_C, y))
    scalar_T = S.Function('scalar_T')(T)
    square_trace = L*scalar_T**2+M*scalar_T
    zero('special square trace is independent of the annihilated dual direction', S.diff(square_trace, Z))
    tangent_control = (x+1)*Q+y*P
    zero('linear-special-slice control has unit Hessian',
         S.hessian(tangent_control, (x, y, P, Q)).det()-1)
    zero('linear-special-slice control has a generic-coordinate trace',
         S.Matrix([tangent_control.subs(x, 0), P]).jacobian((P, Q)).det()+1)
    zero('linear-special-slice control has a genuine tangent fibre plane',
         S.Matrix([S.diff(tangent_control, y, 2), S.diff(tangent_control, y, Q),
                   S.diff(tangent_control, Q, 2)]))
    zero('tangent control fails the nonlinear-square special-trace hypothesis',
         tangent_control.subs({x: 0, y: 0})-Q)

    C_Y = S.Function('C_Y')(Y)
    transported = S.Function('scalar_C')(y-eta*x)
    zero('pure-dual coefficient constant-flow transport',
         S.diff(transported, x)+eta*S.diff(transported, y))
    Ufun, Vfun = S.Function('Ufun')(Y, T), S.Function('Vfun')(Y, T)
    transverse_pure = C_Y*Z+X*Ufun+Vfun
    zero('pure-dual whole-plane Hessian determinant',
         S.hessian(transverse_pure, (X, Z, Y, T)).det()
         -(S.diff(C_Y, Y)*S.diff(Ufun, T))**2)
    transform = S.Matrix([x, y, P-c*y-e*x, Q-d*y-g*x])
    Jtransform = transform.jacobian((x, y, P, Q))
    zero('moving-y adapted source determinant', Jtransform.det()-1)
    zero('moving-y plane generator becomes Y',
         Jtransform*S.Matrix([0, 1, c, d])-S.Matrix([0, 1, 0, 0]))
    zero('transverse generator with y component removed becomes X',
         Jtransform*S.Matrix([1, 0, e, g])-S.Matrix([1, 0, 0, 0]))
    Afun, afun, bfun = (S.Function(name)(U, V) for name in ('Afun', 'afun', 'bfun'))
    transverse_moving = bfun+X*Afun+Y*afun
    zero('moving-y whole-plane fibre block vanishes',
         S.hessian(transverse_moving, (X, Y, U, V))[:2, :2])
    zero('moving-y whole-plane Hessian determinant',
         S.hessian(transverse_moving, (X, Y, U, V)).det()-bracket(Afun, afun)**2)

    # A genuinely nonlinear coordinate endpoint with an exact plane
    # inverse; lift that inverse to the entire cotangent gradient.
    a = U+V*V
    A = V+a*a
    b = U*U+U*V+V**3
    Jbase = S.Matrix([A, a]).jacobian((U, V))
    zero('nonlinear LND endpoint plane Jacobian', Jbase.det()+1)
    za, zA, zu, zv = S.symbols('za zA zu zv')
    Vr = zA-za*za
    Ur = za-Vr*Vr
    base_inverse = {U: Ur, V: Vr}
    zero('nonlinear endpoint plane inverse right composition',
         S.Matrix([A, a]).subs(base_inverse, simultaneous=True)-S.Matrix([zA, za]))
    zero('nonlinear endpoint plane inverse left composition',
         S.Matrix([Ur, Vr]).subs({zA: A, za: a}, simultaneous=True)-S.Matrix([U, V]))
    f = b+X*A+Y*a
    grad = S.Matrix([S.diff(f, t) for t in (X, Y, U, V)])
    base_gradient = S.Matrix([S.diff(b, U), S.diff(b, V)])
    recovered_fibres = (Jbase.T.adjugate()*(S.Matrix([zu, zv])-base_gradient)/Jbase.det())
    recovered_fibres = recovered_fibres.subs(base_inverse, simultaneous=True)
    inv = S.Matrix([recovered_fibres[0], recovered_fibres[1], Ur, Vr])
    zero('nonlinear cotangent gradient inverse left composition',
         inv.subs(dict(zip((zA, za, zu, zv), grad)), simultaneous=True)-S.Matrix([X, Y, U, V]))
    zero('nonlinear cotangent gradient inverse right composition',
         grad.subs(dict(zip((X, Y, U, V), inv)), simultaneous=True)-S.Matrix([zA, za, zu, zv]))

    # A y-dependent nontriangular dual automorphism over K[x,y].
    qw = Q+y*P*P
    pw = P+y*qw*qw
    Pi, Qi = P-y*Q*Q, Q-y*(P-y*Q*Q)**2
    zero('nontriangular y-dependent word dual Jacobian',
         S.Matrix([pw, qw]).jacobian((P, Q)).det()-1)
    zero('nontriangular y-dependent word inverse left',
         S.Matrix([Pi, Qi]).subs({P: pw, Q: qw}, simultaneous=True)-S.Matrix([P, Q]))
    zero('nontriangular y-dependent word inverse right',
         S.Matrix([pw, qw]).subs({P: Pi, Q: Qi}, simultaneous=True)-S.Matrix([P, Q]))
    check('word first trace genuinely depends on y in both entries',
          S.diff(pw, y) != 0 and S.diff(qw, y) != 0)

    # Nagata over K[y], followed by the q shear removing the 4*y^2 term.
    delta = P*P+y*Q
    pn = P+y*delta
    qn = Q-2*P*delta-y*delta**2
    qtrace = qn-4*y*pn
    zero('Nagata invariant identity', pn*pn+y*qn-delta)
    zero('Nagata trace pair dual Jacobian',
         S.Matrix([pn, qtrace]).jacobian((P, Q)).det()-1)
    trace_invariant = P*P+y*(Q+4*y*P)
    Pinv = P-y*trace_invariant
    Qinv = Q+4*y*P+2*P*trace_invariant-y*trace_invariant**2
    zero('Nagata trace inverse left composition',
         S.Matrix([Pinv, Qinv]).subs({P: pn, Q: qtrace}, simultaneous=True)-S.Matrix([P, Q]))
    zero('Nagata trace inverse right composition',
         S.Matrix([pn, qtrace]).subs({P: Pinv, Q: Qinv}, simultaneous=True)-S.Matrix([P, Q]))
    nagata_h = pn*pn+4*y*y*pn+y*qtrace
    zero('Nagata moving-y affine trace normal form', nagata_h-(P*P+y*Q))
    zero('Nagata moving-y trace second derivative vanishes', S.diff(S.expand(nagata_h), y, 2))

    # External mathematical imports are named and byte-pinned separately.
    manifest = json.loads((HERE/'sources/automorphic-trace-import.json').read_text())
    pdf = HERE/'sources'/manifest['retained_file']
    check('Rentschler criterion primary research source is a PDF', pdf.read_bytes().startswith(b'%PDF'))
    check('Rentschler criterion source matches reviewed byte pin',
          hashlib.sha256(pdf.read_bytes()).hexdigest()==manifest['sha256'])
    pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    dependencies = {}
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('pinned structural dependency '+name, dependencies[name]==pins[name])
    result = dict(outcome='PASS', check_count=len(checks), checks=checks,
        scope='All polynomial extensions of an automorphic K[y] dual trace, hence every dual automorphism over K[x,y] and every polynomial R',
        all_base_dependent_dual_automorphisms_covered=True,
        first_trace_y_dependence_allowed=True,
        word_length_or_degree_bound=None,
        full_off_slice_dual_inverse_required=False,
        merely_Keller_first_trace_sufficient=False,
        plane_cases=['contained in x=0: impossible', 'transverse with pure-dual intersection: triangular',
                     'transverse with moving-y intersection: LND with Jacobian mate'],
        elementary_pencil_argument='Top parameter coefficient in (D_b+t*D_a)^n is D_a^n without commutativity',
        external_imports=['Rentschler criterion via Furter Lemma 4.2', 'HCTANG2 structural alternative'],
        import_sha256={manifest['retained_file']: manifest['sha256']},
        dependency_source_sha256=dependencies,
        source_sha256={name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name, 'AUTOMORPHIC_TRACE.md', 'sources/automorphic-trace-import.json']},
        finite_control_is_plane_exhaustion=False,
        nonlinear_base_changes_covered=False, other_graph_charts_covered=False,
        runtime=dict(python=platform.python_version(), sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start, 3),
        formal_verification=False, external_review=False)
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__=='__main__':
    main()
