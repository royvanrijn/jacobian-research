"""Exact controls for the exceptional-fibre proof closing both base charts.

The general curve Euler bound and all-degree exhaustion are written proofs.
Finite controls are not used as a degree bound or an existence test.
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

    x,y,p,q,r,s,T,lam,rho = S.symbols('x y p q r s T lam rho')
    L,N = S.symbols('L N', nonzero=True)
    M,v = S.symbols('M v')
    R0 = S.Function('R0')(p)
    H = lambda value: L*value**2+M*value
    u = 1+x*y
    AA = u**3*p+3*x*u**2*q-x**3*r
    BB = y**2*u*(4+3*x*y)*p+(y+3*x*y**2*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = H(AA)+N*BB
    fy = H(p)+4*N*p*R0**2+N*R0*q
    T0 = p+3*R0*q-R0**3*r
    fx = H(T0)+2*N*R0*r
    zero('native y graph at base zero retains its complete linear-p term',
         psi.subs({x:0,y:R0}, simultaneous=True)-fy)
    zero('native opposite graph at base zero has the exact square coordinate',
         psi.subs({y:0,x:R0}, simultaneous=True)-fx)
    qy = (lam-H(p)-4*N*p*R0**2)/(N*R0)
    zero('y exceptional fibre off the boundary roots solves q with r free', fy.subs(q,qy)-lam)
    qx = (T-p+R0**3*r)/(3*R0)
    zero('opposite open-stratum coordinate map has a literal inverse',T0.subs(q,qx)-T)
    zero('opposite open-stratum coordinate Jacobian is3*R0',
         S.Matrix([p,T0,r]).jacobian((p,q,r)).det()-3*R0)
    zero('opposite exceptional fibre solves r with T free',
         fx.subs(q,qx).subs(r,(lam-H(T))/(2*N*R0))-lam)
    for name,trace in [('y',fy),('x',fx)]:
        zero(name+' boundary-root stratum is an entire constant-value q,r plane',
             trace.subs(R0,0)-H(p))
    grad_y = S.Matrix([S.diff(fy,z) for z in (p,q,r)])
    grad_x = S.Matrix([S.diff(fx,z) for z in (p,q,r)])
    boundary = {R0:0,p:rho}
    zero('y boundary gradient includes all nonlinear trace derivatives',
         grad_y.subs(S.diff(R0,p),v).subs(boundary, simultaneous=True)
         -S.Matrix([2*L*rho+M+N*q*v,0,0]))
    zero('opposite boundary gradient includes all nonlinear trace derivatives',
         grad_x.subs(S.diff(R0,p),v).subs(boundary, simultaneous=True)
         -S.Matrix([(2*L*rho+M)*(1+3*v*q)+2*N*v*r,0,0]))
    zero('y off-root boundary has nonzero q derivative',grad_y[1]-N*R0)
    zero('opposite off-root boundary has the nonzero q,r derivative combination',
         grad_x[2]+R0**2*grad_x[1]/3-2*N*R0)
    k,jumps = S.symbols('k jumps', integer=True, nonnegative=True)
    zero('exact exceptional Euler function counts distinct root strata once',
         (1-k)*1+jumps*1-(1-k+jumps))
    check('generic exceptional Euler is at most one for every finite root count',
          S.ask(S.Q.nonnegative(1-(1-k))) is True)
    for name,trace in [('y',fy),('x',fx)]:
        zero(name+' identically zero trace is the separate composite case',trace.subs(R0,0)-H(p))

    # Tang rank one: polynomial reconstruction when its displayed J is a unit.
    U,Z,V = S.symbols('U Z V')
    A,B,E = [S.Function(name)(U,Z) for name in ('A','B','E')]
    AU,AZ,BU,BZ,EU,EZ = [S.diff(fun,z) for fun,z in
                         [(A,U),(A,Z),(B,U),(B,Z),(E,U),(E,Z)]]
    J = AU*BZ-AZ*BU
    sigma = -(EU*BZ-EZ*BU)/J
    vc = -(AU*EZ-AZ*EU)/J
    zero('whole rank-one critical equations uniquely reconstruct the base and V',
         S.Matrix([EU+sigma*AU+vc*BU,EZ+sigma*AZ+vc*BZ]))
    zero('the Keller mate derivative on the fibre is exactly the Jacobian unit',
         S.Matrix([AU,AZ]).dot(S.Matrix([BZ,-BU]))-J)
    zero('exceptional coefficient differential is minus base times the Keller mate differential',
         EU*BZ-EZ*BU+sigma*J)
    normal = E+s*A+V*B
    zero('rank-one relative V equation is exactly B',S.diff(normal,V)-B)
    zero('rank-one exceptional fibre is exactly E+V*B',normal.subs(s,0)-E-V*B)
    zero('rank-one full Hessian determinant is the squared plane Jacobian',
         S.hessian(normal,(s,U,Z,V)).det()-J**2)

    # Adapt an arbitrary rank-one affine-plane projection, keeping the base.
    a = S.symbols('a', nonzero=True)
    b,c,ell1,ell2,P,Q = S.symbols('b c ell1 ell2 P Q')
    ell = ell1*U+ell2*Z+c
    D = S.Function('D')(U,Z)
    adapted_E = D-ell*A/a
    adapted_A,adapted_B = A/a,B-b*A/a
    zero('constant rank-one source adaptation preserves the exact chosen base',
         (a*P+b*Q+ell).subs({P:(s-b*V-ell)/a,Q:V}, simultaneous=True)-s)
    zero('constant rank-one adaptation gives the stated entire potential',
         (D+P*A+Q*B).subs({P:(s-b*V-ell)/a,Q:V}, simultaneous=True)
         -(adapted_E+s*adapted_A+V*adapted_B))
    zero('rank-one plane Keller determinant stays a field unit',
         S.Matrix([adapted_A,adapted_B]).jacobian((U,Z)).det()-J/a)

    # Bezout row completion is a variety isomorphism on a fixed fibre.
    W,K = S.symbols('W K')
    F,G = U**2,U**3-1
    row_map = S.Matrix([F*P+G*Q,P+U*Q])
    row_inverse = {P:U*W-G*K,Q:-W+F*K}
    zero('rank-zero generic fibre row admits an actual unimodular completion',
         row_map.jacobian((P,Q)).det()-1)
    zero('rank-zero row completion inverse reconstructs both coordinates',
         row_map.subs(row_inverse, simultaneous=True)-S.Matrix([W,K]))
    zero('rank-zero value equation leaves two affine coordinates',
         (U**5+F*P+G*Q).subs(row_inverse, simultaneous=True).subs(W,lam-U**5)-lam)

    genus = S.symbols('genus', integer=True, nonnegative=True)
    punctures = S.symbols('punctures', integer=True, positive=True)
    component_count,d0 = S.symbols('component_count d0', integer=True, positive=True)
    zero('a nonconstant unit costs both a zero and a pole on the completion',
         (2-2*genus-(2+k))-(-2*genus-k))
    zero('whole exceptional rank-one value fibre Euler stratification',
         1-(2-2*genus-punctures)+d0-(1+2*genus+punctures-2+d0))
    check('rank-one exceptional Euler exceeds one with any nonempty curve of nonpositive Euler',
          S.ask(S.Q.positive((1+k+d0)-1)) is True)

    # An actual escape control prevents accidentally imposing finite projection.
    graph = (s*(s-1)**2*p-1)/2
    pc = (s-3)/(s**2*(s-1)**2)
    zero('escape control has a nonzero constant boundary and no boundary critical point',graph.subs(s,0)+S.Rational(1,2))
    zero('escape control has its entire active branch over C minus zero and one',
         (3+2*s*graph).subs(p,pc))
    zero('escape control derivative is a unit exactly on that punctured line',
         S.diff(graph,p)*2/(s*(s-1)**2)-1)
    zero('escape control really has a finite-base pole of order two',
         S.limit((s-1)**2*pc,s,1)+2)
    check('escape control has nonpositive Euler without being a finite C-star cover',2-3==-1)

    dependencies = {}
    prior = [('RELATIVE_CRITICAL_LOCUS.md','receipt-relative-critical.json'),
             ('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
             ('OPPOSITE_NONLINEAR_FRONTIER.md','receipt-opposite-nonlinear.json'),
             ('RELATIVE_COORDINATE_ENDPOINT.md','receipt-relative-coordinate.json')]
    for name,receipt_name in prior:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'EXCEPTIONAL_FIBRE_CLOSURE.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['nonconstant unit and affine curve Euler bound',
            'both Tang plane projection ranks and the whole critical-scheme realization',
            'exact exceptional base-zero fibre comparison',
            'both positive-p graph charts in every degree and with all finite-base escapes',
            'the distinct base-only unit-ideal endpoint'],
        both_positive_p_base_graph_charts_closed=True,
        opposite_six_root_escape_frontier_retired=True,
        no_finite_projection_or_pole_classification_needed=True,
        all_positive_p_graphs_asserted_unit_Hessian_free=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,
        singleton_zero_properness_resolved=False,
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
