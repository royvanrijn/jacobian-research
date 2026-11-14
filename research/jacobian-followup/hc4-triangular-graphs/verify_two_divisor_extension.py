"""Exact identities and controls for TWO_DIVISOR_EXTENSION.md.

Unbounded degree, direction exhaustion and word length are written proofs.
Run from any directory; --output stores an actual replay receipt.
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
        check(name, all(S.cancel(value) == 0 for value in values))

    x, y, P, Q, T, V, ell = S.symbols('x y P Q T V ell')
    L, M, N, a, b, k, d0, eta, tau, c, d = S.symbols('L M N a b k d0 eta tau c d')
    u = 1+x*y
    p0, q0 = S.Function('p0')(P, Q), S.Function('q0')(P, Q)
    E, H, R = (S.Function(name)(x, y, P, Q) for name in ('E', 'H', 'R'))
    p, q = p0+x*u*E, q0+x*u*H
    A = u**3*p+3*x*u**2*q-x**3*R
    B = y**2*u*(4+3*x*y)*p+(y+3*x*y**2*(4+3*x*y))*q+(2*x-3*x**2*y)*R
    psi = L*A*A+M*A+N*B
    zero('two divisor ideals are comaximal', u-x*y-1)
    zero('universal elimination with arbitrary corrections',
         x*x*B+(2-3*x*y)*A-u*(u+1)*p-2*x*(2*u+1)*q)
    zero('first divisor p trace', p.subs(x, 0)-p0)
    zero('first divisor q trace', q.subs(x, 0)-q0)
    zero('second divisor p trace', p.subs(y, -1/x)-p0)
    zero('second divisor q trace', q.subs(y, -1/x)-q0)
    f0 = L*p0*p0+M*p0+4*N*y*y*p0+N*y*q0
    zero('complete first divisor potential', psi.subs(x, 0)-f0)
    Astar = -x**3*R.subs(y, -1/x)
    fstar = L*Astar*Astar+(M-5*N/x**2)*Astar+2*N*q0/x
    zero('complete second divisor potential', psi.subs(y, -1/x)-fstar)

    # ell parametrizes an arbitrary constant dual direction.
    pa, qa = S.Function('pa')(ell), S.Function('qa')(ell)
    slice_jet = L*pa*pa+M*pa+4*N*y*y*pa+N*y*qa
    D = lambda f: S.diff(f, y)+S.diff(f, ell)
    coeff = S.Poly(S.expand(D(D(slice_jet))), y)
    zero('moving-y coefficient y^2', coeff.coeff_monomial(y*y)-4*N*S.diff(pa, ell, 2))
    zero('moving-y coefficient y',
         coeff.coeff_monomial(y)-16*N*S.diff(pa, ell)-N*S.diff(qa, ell, 2))
    constant = 2*L*S.diff(pa, ell)**2+8*N*pa+2*N*S.diff(qa, ell)
    zero('moving-y reduced constant coefficient',
         coeff.coeff_monomial(1).subs(S.diff(pa, ell, 2), 0)-constant)
    zero('moving-y compatibility forces delta p zero',
         S.diff(constant, ell).subs({S.diff(pa, ell, 2): 0,
             S.diff(qa, ell, 2): -16*S.diff(pa, ell)}, simultaneous=True)
         +24*N*S.diff(pa, ell))
    j11, j12, j21, j22 = S.symbols('j11 j12 j21 j22')
    J = S.Matrix([[j11, j12], [j21, j22]])
    zero('Jacobian inverse preserves directional component ideal',
         J.adjugate()*J*S.Matrix([c, d])-J.det()*S.Matrix([c, d]))
    h = S.Function('h')(T)
    tracep, traceq = a*T+d0, k*V+h
    zero('adapted trace Keller determinant',
         S.Matrix([tracep, traceq]).jacobian((T, V)).det()-a*k)
    adapted = L*tracep**2+M*tracep+4*N*y*y*tracep+N*y*traceq
    zero('adapted first divisor coefficient', S.diff(adapted, V)-N*k*y)
    transported = N*k*(y-eta*x)
    zero('transported coefficient initial value', transported.subs(x, 0)-N*k*y)
    zero('transported coefficient constant-flow equation',
         S.diff(transported, x)+eta*S.diff(transported, y)+tau*S.diff(transported, T))
    initial = S.Function('initial')(y, T)
    flowed = initial.subs({y: y-eta*x, T: T-tau*x}, simultaneous=True)
    zero('general constant-flow solution identity',
         S.diff(flowed, x)+eta*S.diff(flowed, y)+tau*S.diff(flowed, T))
    zero('second divisor coefficient when A is S-independent',
         S.diff(L*ell**2+(M-5*N/x**2)*ell+2*N*traceq/x, V)-2*N*k/x)
    zero('two divisor discrepancy',
         x*(transported.subs(y, -1/x)-2*N*k/x)/(N*k)+eta*x*x+3)
    check('two divisor discrepancy has nonzero constant coefficient',
          S.Poly(eta*x*x+3, x).coeff_monomial(1) == 3)
    Kstar = M-5*N/x**2
    Vstar = S.Function('Vstar')(x, T, V)
    Hstar = S.Function('Hstar')(x, T)
    qbar_required = -L*x*Vstar**2/(2*N)-k*(1+eta*x*x)*V/2+Hstar
    Abar_required = Vstar-Kstar/(2*L)
    boundary_required = L*Abar_required**2+Kstar*Abar_required+2*N*qbar_required/x
    zero('general second-trace square-plus-affine necessary form',
         S.diff(boundary_required, V)+N*k*(1+eta*x*x)/x)
    slope = S.Function('slope')(x, T)
    zero('affine second-trace slope resonance',
         x*(2*N*slope/x-transported.subs(y, -1/x))/(2*N)
         -(slope+k*(1+eta*x*x)/2))
    # Symbolic leading coefficients for arbitrary positive degree n.
    n = S.Symbol('n', integer=True, positive=True)
    an = S.Symbol('a_n', nonzero=True)
    check('square has strictly larger degree for every n>=1', (2*n > n) is S.true)
    check('square leading coefficient is nonzero for nonzero L and a_n',
          S.Poly(L*(an*V**3+V+1)**2+(M-5*N/x**2)*(an*V**3+V+1)
                 +2*N*traceq/x, V).LC() == L*an**2)

    # A genuine two-shear base-dependent word, outside either single order.
    divisor = x*u
    qword = Q+divisor*P**2
    pword = P+divisor*qword**2
    Jword = S.Matrix([pword, qword]).jacobian((P, Q))
    zero('nontriangular base-dependent word Jacobian is one', Jword.det()-1)
    check('word is nonlinear in both dual directions',
          S.diff(pword, Q, 2) != 0 and S.diff(qword, P, 2) != 0)
    for label, sub in [('first', {x: 0}), ('second', {y: -1/x})]:
        zero(label+' divisor whole-word identity',
             S.Matrix([pword, qword]).subs(sub)-S.Matrix([P, Q]))
    # Compose with a nonlinear base-independent Keller pair in both orders.
    basep, baseq = a*P+d0, b*Q+P**2
    compositions = [
        (pword.subs({P: basep, Q: baseq}, simultaneous=True),
         qword.subs({P: basep, Q: baseq}, simultaneous=True)),
        (basep.subs({P: pword, Q: qword}, simultaneous=True),
         baseq.subs({P: pword, Q: qword}, simultaneous=True))]
    for order, pair in enumerate(compositions):
        for label, sub in [('first', {x: 0}), ('second', {y: -1/x})]:
            zero(f'composition order {order} {label} trace',
                 S.Matrix(pair).subs(sub)-S.Matrix([basep, baseq]))
    Aword = u**3*pword+3*x*u**2*qword
    Bword = y**2*u*(4+3*x*y)*pword+(y+3*x*y**2*(4+3*x*y))*qword
    word_potential = L*Aword**2+M*Aword+N*Bword
    zero('whole-potential nonlinear-word leading Q coefficient',
         S.Poly(word_potential, Q).coeff_monomial(Q**4)-L*u**6*divisor**2)
    # Full substitutions need not be Keller: only their two traces do.
    nonkeller = S.Matrix([P+divisor*P**2, Q]).jacobian((P, Q)).det()
    zero('non-Keller interior is allowed', nonkeller-(1+2*divisor*P))
    check('non-Keller control has nonconstant full Jacobian', S.diff(nonkeller, P) != 0)
    # Replacing xu by x loses the second trace, even for a dual automorphism.
    badq = Q+x*P**2
    zero('x-only shear is still a dual automorphism',
         S.Matrix([P, badq]).jacobian((P, Q)).det()-1)
    zero('x-only shear has the first required trace', badq.subs(x, 0)-Q)
    check('x-only shear fails the second required trace',
          S.expand(badq.subs(y, -1/x)-Q) == x*P**2)
    # A different-trace word passes the necessary boundary slope but fails
    # the actual affine-direction equation away from both divisors.
    q1 = Q+x*P
    p1 = P-S.Rational(3, 2)*x*y*y*q1
    q2 = q1+x*p1
    zero('three-shear boundary survivor Jacobian',
         S.Matrix([p1, q2]).jacobian((P, Q)).det()-1)
    zero('three-shear boundary survivor first trace',
         S.Matrix([p1, q2]).subs(x, 0)-S.Matrix([P, Q]))
    zero('three-shear boundary survivor second q trace',
         q2.subs(y, -1/x)-(-Q/2+x*P/2))
    Asurvive = u**3*p1+3*x*u**2*q2
    Bsurvive = y*y*u*(4+3*x*y)*p1+(y+3*x*y*y*(4+3*x*y))*q2
    fsurvive = L*Asurvive**2+M*Asurvive+N*Bsurvive
    zero('boundary survivor f_Q matches first divisor',
         S.diff(fsurvive.subs(x, 0), Q)-N*y)
    zero('boundary survivor f_Q matches transported second divisor',
         S.diff(fsurvive.subs(y, -1/x), Q)+N/x)
    zero('boundary survivor fails interior affine-direction equation',
         S.diff(fsurvive, Q, 2)-2*L*S.diff(Asurvive, Q)**2)
    check('boundary survivor interior obstruction is nonzero',
          S.expand(S.diff(Asurvive, Q)) != 0)
    ambient = S.Matrix([0, 0, -S.Rational(1, 2), 0, 0])
    variables = (x, y, P, Q)
    embedding = S.Matrix([x, y, p, q, R])
    zero('marked collision still requires matching p covectors',
         embedding.jacobian(variables).T*ambient+S.Matrix([S.diff(p, z) for z in variables])/2)
    for point in [(1, -S.Rational(3, 2)), (-1, S.Rational(3, 2))]:
        check('marked point is outside both proof divisors '+str(point),
              (x*u).subs(dict(zip((x, y), point))) != 0)

    pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    hashes = {}
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        hashes[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('pinned structural dependency '+name, hashes[name] == pins[name])
    local_hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                    for name in [Path(__file__).name, 'TWO_DIVISOR_EXTENSION.md']}
    result = dict(outcome='PASS', check_count=len(checks), checks=checks,
        theorem_scope='p=p0+x*(1+xy)*E, q=q0+x*(1+xy)*H; (p0,q0) any base-independent plane Keller pair; E,H,R arbitrary polynomials',
        arbitrary_word_length_proof='Every shear congruent to identity modulo x*(1+xy) preserves that ideal under composition',
        full_base_dependent_Jacobian_required=False,
        unrestricted_base_dependent_words_covered=False,
        nonlinear_base_changes_covered=False,
        elementary_proof='Written slice, component-ideal, polynomial-degree and two-divisor plane exhaustion',
        global_invertibility_import='HCTANG2 only; no plane degree bound or assumption of JC2',
        source_sha256=local_hashes, dependency_source_sha256=hashes,
        runtime=dict(python=platform.python_version(), sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started, 3),
        formal_verification=False, external_review=False)
    rendered = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end='')


if __name__ == '__main__':
    main()
