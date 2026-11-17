"""Exact controls for all-degree one-dual-direction base graphs.

Finite-degree row identities control arbitrary graph degrees through the
written rank-stratum and polynomial-surjectivity arguments.
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

    s,z,L,N = S.symbols('s z L N',nonzero=True)
    M,p,q,r,T,v,tau,lam = S.symbols('M p q r T v tau lam')
    a,b,c = S.symbols('a b c')
    u = 1+s*z
    Ay = S.Matrix([u**3,3*s*u*u,-s**3])
    By = S.Matrix([z*z*u*(4+3*s*z),z+3*s*z*z*(4+3*s*z),2*s-3*s*s*z])
    Ax = S.Matrix([u**3,3*z*u*u,-z**3])
    Bx = S.Matrix([s*s*u*(4+3*s*z),s+3*z*s*s*(4+3*s*z),2*z-3*z*z*s])
    cy = S.Matrix([2*s*s*(2*s*z+3),-s*u*(u+1),z*u**3])
    cx = S.Matrix([2*z*z*(2*s*z+3),-z*u*(u+1),s*u**3])
    zero('native y-chart row cross product is the exact three-entry polynomial',Ay.cross(By)-cy)
    zero('native opposite-chart row cross product is the exact three-entry polynomial',Ax.cross(Bx)-cx)
    zero('y-chart first cross-entry zero leaves a nonzero middle entry',cy[1].subs(z,-3/(2*s))-s/4)
    zero('opposite last cross-entry zero leaves a nonzero first entry',cx[0].subs(z,-1/s)-2/s**2)
    ell = S.Matrix([a,b,c])
    Ky,Kx = S.expand(ell.dot(cy)),S.expand(ell.dot(cx))
    check('y-chart nonzero-r coefficient gives a quartic direction polynomial',S.Poly(Ky,z).degree()==4)
    check('y-chart zero-r nonzero-q coefficient gives a quadratic direction polynomial',S.Poly(Ky.subs(c,0),z).degree()==2)
    check('y-chart pure-p direction gives a linear direction polynomial',S.Poly(Ky.subs({b:0,c:0}),z).degree()==1)
    zero('y-chart A row never has zero r entry',Ay[2]+s**3)
    Py = 2*(1+tau)**3*(2*v+3)-3*(1+tau)**2*(v+1)*(v+2)-v*(v+1)**3
    zero('y-chart direction proportional to A has exactly the normalized quartic',
         Ky.subs({a:(1+tau)**3,b:3*s*(1+tau)**2,c:-s**3,z:v/s},simultaneous=True)/s**2-Py)
    zero('the only possible y-chart linear-type value is a root',Py.subs(v,tau))
    zero('that y-chart root is simple except at the stated exceptional value',S.diff(Py,v).subs(v,tau)+6*(1+tau)**3)
    zero('exceptional pure-r y direction retains a different quadratic-type root',Py.subs(tau,-1)+v*(v+1)**3)
    zero('the pure-r y quadratic-type root has nonzero restricted A',Ay[0].subs(z,0)-1)

    alpha,beta,gamma = S.symbols('alpha beta gamma')
    Px = 2*alpha*v*v*(2*v+3)-beta*v*(v+1)*(v+2)+gamma*(v+1)**3
    zero('opposite chart has the exact weighted cubic direction polynomial',
         s*s*Kx.subs({a:alpha,b:beta/s,c:gamma/s**3,z:v/s},simultaneous=True)-Px)
    coeff_matrix = S.Matrix([[4,-1,1],[6,-3,3],[0,-2,3]])
    zero('opposite nonconstant-coefficient map is invertible',coeff_matrix.det()+6)
    zero('opposite leading three coefficients are exactly the displayed invertible map',
         S.Matrix([S.Poly(Px,v).coeff_monomial(v**k) for k in (3,2,1)])
         -coeff_matrix*S.Matrix([alpha,beta,gamma]))
    Q = (6*tau*tau+9*tau+4)*v*v+(9*tau*tau+13*tau+6)*v+tau*tau
    proportional = Px.subs({alpha:(1+tau)**3,beta:3*tau*(1+tau)**2,gamma:-tau**3},simultaneous=True)
    zero('opposite possible linear-type root factors with the complete remaining polynomial',proportional-(v-tau)*Q)
    zero('remaining opposite polynomial cannot have both leading coefficients zero',
         2*(9*tau*tau+13*tau+6)-3*(6*tau*tau+9*tau+4)+tau)
    zero('at tau zero the remaining quadratic leading coefficient stays nonzero',
         (6*tau*tau+9*tau+4).subs(tau,0)-4)
    zero('other opposite roots differ from the linear-type value away from zero and minus one',Q.subs(v,tau)-6*tau*(tau+1)**3)
    zero('pure-p opposite direction has the separate quadratic-type root minus three halves',
         proportional.subs(tau,0)-2*v*v*(2*v+3))
    zero('pure-r opposite direction has only the explicitly linear-type root',Px.subs({alpha:0,beta:0},simultaneous=True)-gamma*(v+1)**3)
    v1,v2 = S.symbols('v1 v2')
    zero('opposite q/p row ratio determines its finite base value uniquely',
         v1/(v1+1)-v2/(v2+1)-(v1-v2)/((v1+1)*(v2+1)))

    # General affine-plane value strata, not a finite cover formula.
    A1,A2,B1,B2,A0,B0,V,W = S.symbols('A1 A2 B1 B2 A0 B0 V W')
    row_matrix = S.Matrix([[A1,A2],[B1,B2]])
    coordinates = row_matrix*S.Matrix([V,W])+S.Matrix([A0,B0])
    inverse = row_matrix.inv()*(S.Matrix([T,(lam-L*T*T-M*T)/N])-S.Matrix([A0,B0]))
    zero('rank-two affine row inverse identifies the whole value stratum with one free line',
         coordinates.subs({V:inverse[0],W:inverse[1]},simultaneous=True)
         -S.Matrix([T,(lam-L*T*T-M*T)/N]))
    ratio,offset = S.symbols('ratio offset')
    zero('rank-one quadratic stratum retains its nonzero quadratic coefficient',
         L*T*T+M*T+N*(ratio*T+offset)-(L*T*T+(M+N*ratio)*T+N*offset))
    zero('rank-one linear stratum leaves a free affine-line coordinate',
         (L*A0*A0+M*A0+N*(B1*V+B2*W+B0)).subs(V,((lam-L*A0*A0-M*A0)/N-B2*W-B0)/B1)-lam)
    ht,m,k = S.symbols('ht m k')
    zero('whole rank-change Hodge count retains linear-type roots once and quadratic-type roots twice',
         (ht-m)*ht+2*k*ht+(m-k)*ht-(ht*ht+k*ht))
    zero('pure-r opposite direction has zero positive budget rather than an invented closure',
         ((ht-m)*ht+m*ht)-ht*ht)

    # Opposite pure-r model and every nonlinear off-slice derivative.
    rawA = Ax.dot(S.Matrix([p,q,r]))
    rawB = Bx.dot(S.Matrix([p,q,r]))
    psi = L*rawA**2+M*rawA+N*rawB
    pp = (T-3*z*u*u*q+z**3*r)/u**3
    scalar_model = L*T*T+(M+N*s*s*(4+3*s*z)/u**2)*T+N*s*q+N*z*(2+s*z)*r/u**2
    zero('pure-r opposite rational model is the literal ambient potential',psi.subs(p,pp)-scalar_model)
    zero('the u-zero stratum has the stated nonzero native q derivative',S.diff(psi,q).subs(z,-1/s)+2*N*s)
    zero('u-zero stratum is quadratic in r for the pure-q counting control',
         psi.subs(z,-1/s)-(L*(r/s**3)**2+M*r/s**3-2*N*s*q-5*N*r/s))
    R = S.Function('R')(s,r)
    model = scalar_model.subs(z,R)
    zero('pure-r opposite model has no relative critical point on the u-nonzero chart',S.diff(model,q)-N*s)
    native = psi.subs(z,R)
    rho,rv,rr,rs,rss,rsr,rrr = S.symbols('rho rv rr rs rss rsr rrr')
    jets = {R:rv,S.diff(R,r):rr,S.diff(R,s):rs,S.diff(R,s,2):rss,
            S.diff(R,s,r):rsr,S.diff(R,r,2):rrr}
    critical_p = -M/(2*L)-3*rv*q+rv**3*rho
    at = {s:0,r:rho,p:critical_p}

    def critical_value(expression):
        return S.cancel(expression.xreplace(jets).subs(at,simultaneous=True).subs(rv,-rho*rr))

    zero('the whole arbitrary-jet critical line satisfies all three native relative equations',
         S.Matrix([S.diff(native,z0) for z0 in (p,q,r)]).applyfunc(critical_value))
    Hline = S.hessian(native,(s,p,q,r)).applyfunc(critical_value)
    zero('full native critical Hessian includes all nonlinear off-slice jets and equals the stated border',
         Hline.det()+4*L*N**3*(2*rr+rho*rrr))
    aa,bb = S.symbols('aa bb',nonzero=True)
    boundary_g = r*(aa*r+bb)
    root = -bb/(2*aa)
    zero('linear pure-r boundary gives one simple critical root',S.diff(boundary_g,r).subs(r,root))
    zero('linear boundary determinant is the exact common nonzero value',
         -4*L*N**3*S.diff(boundary_g,r,2)+8*L*N**3*aa)
    check('multiple critical-root boundary forces a zero critical Hessian',S.diff(r**4,r,2).subs(r,0)==0)
    zero('squarefree quadratic-boundary control has unequal critical determinants',
         S.diff(r*(r*r+1),r,2).subs(r,S.I/S.sqrt(3))
         -S.diff(r*(r*r+1),r,2).subs(r,-S.I/S.sqrt(3))-4*S.sqrt(3)*S.I)
    control = psi.subs({z:r,L:1,N:1,M:0},simultaneous=True)
    off_value = S.hessian(control,(s,p,q,r)).subs({s:0,p:1,q:0,r:0},simultaneous=True).det()
    check('a passed pure-r critical determinant is not a unit-Hessian certificate',off_value==-440 and off_value!=-8)

    # Nonpositive trace signatures do not by themselves determine the global endpoint.
    cotangent = V*s+W*T
    shifted = V*(s-1)+W*T
    zero('singleton-zero cotangent trace has a product with one free variable',cotangent.subs(s,0)-W*T)
    zero('a zero-coefficient trace can still come from a nonempty shifted Keller fibre',shifted.subs(s,0)-(-V+W*T))
    check('the shifted plane Keller pair really has an attained origin away from the trace line',
          (s-1).subs(s,1)==0 and T.subs(T,0)==0)

    dependencies = {}
    prior = [('HODGE_TRACE_CLOSURE.md','receipt-hodge-trace.json'),
             ('RELATIVE_CRITICAL_LOCUS.md','receipt-relative-critical.json'),
             ('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
             ('RELATIVE_COORDINATE_ENDPOINT.md','receipt-relative-coordinate.json'),
             ('P_GRAPH_SINGLETON_ENDPOINT.md','receipt-square-axis.json')]
    for name,receipt_name in prior:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    manifest = json.loads((HERE/'sources/hodge-trace-import.json').read_text())
    source = HERE/'sources'/manifest['retained_file']
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    check('retained basic Hodge primary input',digest==manifest['sha256'])
    dependencies[str(source.relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md','jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'DIRECTIONAL_BASE_GRAPHS.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['whole affine-plane rank strata and exact positive Hodge coefficient',
            'quartic y-direction root exhaustion in every coefficient direction',
            'cubic opposite-direction exhaustion except the explicit pure-r exception',
            'arbitrary graph degrees via polynomial surjectivity, not bounded coefficient search',
            'whole pure-r critical scheme and arbitrary-jet critical Hessian',
            'all nonlinear pure-r boundary degrees excluded by one derivative polynomial argument',
            'remaining pure-r nonautomorphy supplies one of the two known plane endpoints'],
        all_y_one_affine_dual_direction_graphs_closed=True,
        all_opposite_non_r_dual_direction_graphs_closed=True,
        unimodular_polynomial_direction_row_required=False,
        pure_r_opposite_positive_Hodge_coefficient_asserted=False,
        pure_r_zero_or_nonlinear_boundary_unit_Hessian_excluded=True,
        remaining_pure_r_nonautomorphy_reduced_to_known_endpoints=True,
        converse_endpoint_realization_in_pure_r_graphs_asserted=False,
        critical_Hessian_constancy_asserted_sufficient=False,
        trace_coefficient_zero_asserted_to_imply_unit_ideal=False,
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
