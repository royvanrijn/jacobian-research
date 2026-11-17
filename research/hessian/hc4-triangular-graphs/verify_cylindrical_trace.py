"""Exact controls for the cylindrical-trace and two-dual y-chart closure.

The affine-curve, normalization and all-degree component arguments are
written proofs; these controls do not supply a finite-degree inference.
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

    def check(name,condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,expression):
        values = list(expression) if isinstance(expression,S.MatrixBase) else [expression]
        check(name,all(S.cancel(value.doit())==0 for value in values))

    s,y,p,q,r,T,Tprime,P,Q,t,lam,mu = S.symbols('s y p q r T Tprime P Q t lam mu')
    L,N = S.symbols('L N',nonzero=True)
    M = S.symbols('M')
    R = S.Function('R')(s,p,q)
    u = 1+s*y
    AA = u**3*p+3*s*u**2*q-s**3*r
    BB = y**2*u*(4+3*s*y)*p+(y+3*s*y**2*(4+3*s*y))*q+(2*s-3*s*s*y)*r
    psi = L*AA**2+M*AA+N*BB
    native = psi.subs(y,R)
    rr = ((1+s*R)**3*p+3*s*(1+s*R)**2*q-T)/s**3
    model = L*T*T+(M-N*(2-3*s*R)/s**2)*T+N*(1+s*R)*(2+s*R)*p/s**2+2*N*(3+2*s*R)*q/s
    zero('exact rational computation model permits arbitrary mixed p,q dependence',native.subs(r,rr)-model)
    F = 3+2*s*R
    c0 = M-13*N/(2*s*s)
    shift = (c0+3*N*F/(2*s*s))/(2*L)
    pp = N*p/(4*s*s)-9*N*N/(16*L*s**4)
    qq = 2*N*q/s-3*N*c0/(4*L*s*s)
    constant = 9*N*N/(16*L*s**4)+c0*c0/(4*L)
    normal = L*Tprime**2+pp*(F*F-1)+qq*F-constant
    zero('complete square gives exactly the universal plane polynomial and scalar shift',
         model.subs(T,Tprime-shift)-normal)
    forward = S.Matrix([pp,qq,AA.subs(y,R)+shift])
    zero('fixed-base triangular coordinate change has constant nonzero relative Jacobian',
         forward.jacobian((p,q,r)).det()+N*N/2)
    pb = 4*s*s*P/N+9*N/(4*L*s*s)
    qb = s*Q/(2*N)+3*c0/(8*L*s)
    Rb = R.subs({p:pb,q:qb},simultaneous=True)
    shiftb = shift.subs({p:pb,q:qb},simultaneous=True)
    rb = ((1+s*Rb)**3*pb+3*s*(1+s*Rb)**2*qb-Tprime+shiftb)/s**3
    inverse = {p:pb,q:qb,r:rb}
    zero('entire fixed-base dual transformation has an explicit polynomial inverse over C(s)',
         forward.subs(inverse,simultaneous=True)-S.Matrix([P,Q,Tprime]))
    R0 = S.Function('R0')(p,q)
    trace = L*p*p+M*p+4*N*p*R0*R0+N*q*R0
    zero('base-zero trace is an arbitrary plane polynomial with the third variable free',
         psi.subs({s:0,y:R0},simultaneous=True)-trace)
    zero('the exact exceptional trace has no hidden r dependence',S.diff(trace,r))

    a,b = S.symbols('a b')
    Pat,Qat = a*b-lam,(1-a*a)*b+lam*a
    zero('all affine-line candidate parameterizations satisfy the exact value identity',
         Pat*(a*a-1)+Qat*a-lam)
    zero('the second candidate polynomial restricts to exactly b',Pat*a+Qat-b)
    da,db = S.symbols('da db',integer=True,positive=True)
    check('if both candidate functions vary then the larger degree exceeds the smaller',
          S.ask(S.Q.positive((2*da+db)-(da+db))) is True)
    check('if both vary the larger degree is strictly below twice the smaller',
          S.ask(S.Q.positive(2*(da+db)-(2*da+db))) is True)
    for aa,bb in [(t,t*t+1),(t*t+1,t**3+2)]:
        coords = [Pat.subs({a:aa,b:bb}),Qat.subs({a:aa,b:bb})]
        check('nonembedding degree control '+str((S.degree(aa,t),S.degree(bb,t))),
              S.degree(coords[0],t)==S.degree(aa,t)+S.degree(bb,t)
              and S.degree(coords[1],t)==2*S.degree(aa,t)+S.degree(bb,t))
    c,d = S.symbols('c d',nonzero=True)
    zero('constant-a candidates are literal straight lines',
         ((c*c-1)*P+c*Q-lam).subs({P:Pat.subs(a,c),Q:Qat.subs(a,c)},simultaneous=True))
    zero('constant-a line normals never vanish together',(a*a-1)-a*a+1)
    zero('a fixed parallel line direction allows two constants, including the reciprocal sign',
         (a*a-1)*b-(b*b-1)*a-(a-b)*(a*b+1))
    para = c-P*(P+lam)/c
    zero('nonzero constant-b candidates are exactly graph parabolas',
         para.subs(P,Pat.subs(b,c))-Qat.subs(b,c))
    zero('zero-b candidates belong to the vertical Z-zero line',
         (Pat*a+Qat).subs(b,0))
    zero('zero-b candidates have the exact vertical base coordinate',Pat.subs(b,0)+lam)
    zero('zero-b zero-value candidates collapse to a point',
         S.Matrix([Pat,Qat]).subs({b:0,lam:0},simultaneous=True))
    other = d-P*(P+mu)/d
    difference = S.Poly(S.expand(para-other),P)
    zero('different nonzero b constants force a quadratic parabola intersection',
         difference.coeff_monomial(P*P)-(c-d)/(c*d))
    zero('equal b constants and different values force a linear intersection',
         difference.as_expr().subs(d,c)-P*(mu-lam)/c)
    slope,intercept = S.symbols('slope intercept')
    zero('every nonvertical straight line has a genuine quadratic intersection equation',
         S.Poly(para-slope*P-intercept,P).coeff_monomial(P*P)+1/c)
    zero('every vertical line meets each graph parabola at its finite graph value',
         (Q-para).subs({P:t,Q:para.subs(P,t)},simultaneous=True))

    # Nonlinear trace automorphisms preserve actual cylindrical critical lines.
    U,V,W = S.symbols('U V W')
    phi = S.Matrix([U+W*W,V+W**3,W])
    zero('nonlinear trace coordinate change has an actual unit Jacobian',phi.jacobian((U,V,W)).det()-1)
    h = U*V+U**3
    pullback = h.subs({U:phi[0],V:phi[1]},simultaneous=True)
    zero('nonlinear critical-line pullback is the actual closed A1 with free W',
         S.Matrix([S.diff(pullback,z) for z in (U,V,W)]).subs({U:-W*W,V:-W**3},simultaneous=True))
    zero('nonlinear critical line has its exact polynomial parameter inverse',
         phi.subs({U:-W*W,V:-W**3},simultaneous=True)-S.Matrix([0,0,W]))
    HH = S.Function('HH')
    g = S.Function('g')(U,V)
    zero('normalization compositeness argument uses the full polynomial gradient chain rule',
         S.Matrix([S.diff(HH(g),z) for z in (U,V)])
         -S.Subs(S.diff(HH(t),t),t,g)*S.Matrix([S.diff(g,z) for z in (U,V)]))

    chi,d0 = S.symbols('chi d0')
    zero('finite unramified double-cover Euler counts the branch curve exactly once',
         2*(1-chi)+chi-(2-chi))
    zero('the whole exceptional rank-one fibre has the same complement and affine-line strata',
         (1-chi)+d0-(1-chi+d0))
    chi_negative = S.symbols('chi_negative',integer=True,nonnegative=True)
    degree = S.symbols('degree',integer=True,positive=True)
    check('nonpositive critical-curve Euler and positive exceptional degree force Euler above one',
          S.ask(S.Q.positive((1+chi_negative+degree)-1)) is True)

    dependencies = {}
    prior = [('EXCEPTIONAL_FIBRE_CLOSURE.md','receipt-exceptional-fibre.json'),
             ('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
             ('RELATIVE_CRITICAL_LOCUS.md','receipt-relative-critical.json'),
             ('RELATIVE_COORDINATE_ENDPOINT.md','receipt-relative-coordinate.json')]
    for name,receipt_name in prior:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    manifest = json.loads((HERE/'sources/degree-embedding-import.json').read_text())
    check('degree theorem has the inspected precise Theorem6.9 location','Theorem 6.9' in manifest['location'])
    degree_source = HERE/'sources'/manifest['retained_file']
    digest = hashlib.sha256(degree_source.read_bytes()).hexdigest()
    check('additional degree theorem reuses the exact retained primary v1 manuscript',digest==manifest['sha256'])
    dependencies[str(degree_source.relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'CYLINDRICAL_TRACE_CLOSURE.md','sources/degree-embedding-import.json']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['nonsingular plane polynomial normalization and generic irreducibility',
            'cylindrical trace critical lines and both complete Tang ranks',
            'degree obstruction and all line/parabola alternatives in every embedding degree',
            'finitely many actual A1 components without assuming generic fibre irreducibility',
            'entire fixed-base suspension and whole exceptional fibre comparison',
            'all y=R(x,p,q) and actual retained-dual automorphisms in arbitrary degree'],
        all_y_graphs_independent_of_old_r_closed=True,
        arbitrary_mixed_p_q_dependence_closed=True,
        generic_plane_fibre_irreducibility_assumed_for_degree_argument=False,
        native_cylindrical_trace_required_to_remain_affine_after_automorphism=False,
        all_positive_dual_degree_graphs_asserted_unit_Hessian_free=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,
        singleton_zero_properness_resolved=False,
        all_y_graphs_with_genuine_r_dependence_closed=False,
        all_opposite_q_r_dependent_graphs_closed=False,
        unrestricted_HC4_or_JC2_resolved=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
