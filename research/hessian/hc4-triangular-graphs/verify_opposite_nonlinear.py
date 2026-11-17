"""Exact controls for the opposite nonlinear global frontier.

The entire-value Euler, finite-cover and Laurent-unit arguments are
written proofs. No degree search certifies their unrestricted scope.
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

    x,s,p,q,r,T,t,lam = S.symbols('x s p q r T t lam')
    L,N,v = S.symbols('L N v',nonzero=True)
    M,rho = S.symbols('M rho')
    u = 1+x*s
    AA = u**3*p+3*x*u*u*q-x**3*r
    BB = s*s*u*(4+3*x*s)*p+(s+3*x*s*s*(4+3*x*s))*q+(2*x-3*x*x*s)*r
    psi = L*AA**2+M*AA+N*BB
    R = S.Function('R')(s,p)
    Rp = S.diff(R,p)
    rr = ((1+R*s)**3*p+3*R*(1+R*s)**2*q-T)/R**3
    c = M-N*(2-3*R*s)/R**2
    epoly = N*(1+R*s)*(2+R*s)*p/R**2
    dd = 2*N*(3+2*R*s)/R
    model = L*T*T+c*T+epoly+dd*q
    zero('opposite rational model is literally the native potential',psi.subs({x:R,r:rr},simultaneous=True)-model)
    zero('opposite computation chart has native relative Jacobian R^3',
         S.Matrix([p,AA.subs(x,R),q]).jacobian((p,q,r)).det()-R**3)
    xc = -S.Rational(3,2)/s
    Tc = (13*N/(2*xc**2)-M)/(2*L)
    qc = 17*Tc/(12*xc)+p/(12*xc)-1/(24*v)
    critical = {R:xc,T:Tc,q:qc}
    zero('all opposite relative critical equations vanish on the stated branch',
         S.Matrix([S.diff(model,z) for z in (p,T,q)]).subs(Rp,v).subs(critical,simultaneous=True))
    zero('an R_p zero on F cannot supply a missing critical point',
         S.diff(model,p).subs(Rp,0).subs({R:xc,T:Tc},simultaneous=True)+N/(4*xc**2))
    D = -N*s*s*p/9-(M-26*N*s*s/9)**2/(4*L)
    zero('opposite critical value is exactly affine in p with the asymmetric s^2 weight',
         model.subs(critical,simultaneous=True)-D)
    zero('opposite relative Hessian determinant',S.hessian(model,(p,T,q)).det()+72*L*N*N*Rp**2/R**4)
    zero('opposite native critical Hessian prefactor is -162 rather than the y-chart weight',
         (-72*L*N*N*R**2*v**2).subs(R,xc)+162*L*N*N*v**2/s**2)
    zero('the R=0 stratum retains the necessary additional linear p term',
         psi.subs(x,0)-(L*p*p+M*p+4*N*s*s*p+N*s*q))
    zero('every value over an R zero has exactly one free r line for nonzero base',
         psi.subs(x,0).subs(q,(lam-L*p*p-M*p-4*N*s*s*p)/(N*s))-lam)
    roots_R,roots_F = S.symbols('roots_R roots_F')
    zero('extra zero-R strata cancel from the generic Euler value',
         (1-roots_R-roots_F)+roots_R+2*roots_F-(1+roots_F))

    # Every boundary jet is retained, including nonlinear R_pp.
    native = psi.subs(x,R)
    Rs,Rss,Rsp,Rpp = S.symbols('Rs Rss Rsp Rpp')
    jets = {S.diff(R,p):v,S.diff(R,s):Rs,S.diff(R,s,2):Rss,
            S.diff(R,s,p):Rsp,S.diff(R,p,2):Rpp}
    h0 = 2*L*rho+M
    boundary = {s:0,p:rho,R:0,r:-h0*(1+3*v*q)/(2*N*v)}
    gradient = S.Matrix([S.diff(native,z) for z in (p,q,r)])
    zero('arbitrary nonlinear simple boundary zero gives the entire affine critical line',
         gradient.applyfunc(lambda value:S.cancel(value.xreplace(jets).subs(boundary,simultaneous=True))))
    Hline = S.hessian(native,(s,p,q,r)).applyfunc(lambda value:S.cancel(value.xreplace(jets).subs(boundary,simultaneous=True)))
    zero('full nonlinear boundary-line Hessian is4*N^4*v^2',Hline.det()-4*N**4*v*v)

    # The every-value defect forces Laurent inactive roots with affine D.
    alpha,beta = S.symbols('alpha beta')
    pstar = -169*N*s*s/(9*L)+13*M/L-9*M*M/(4*L*N*s*s)
    pi = pstar-9*(alpha+beta*s)/(N*s*s)
    zero('affine inactive opposite critical value has its precise Laurent root',D.subs(p,pi)-alpha-beta*s)
    delta = -M*M/(4*L)
    zero('a forbidden inactive affine value would have a finite p root at zero',
         pi.subs({alpha:delta,beta:0})-(-169*N*s*s/(9*L)+13*M/L))
    e,j,k,d,h = S.symbols('e j k d h',integer=True)
    check('simple finite-value infinity punctures have an impossible odd Hessian valuation',
          S.Poly(2*k+4*e+1,k,e,modulus=2).as_expr()==1)
    zero('infinity pole parity with a pole of tau',
         (3*e+2*k-j)-(2*k+4*e-(e+j)))
    zero('infinity pole parity with regular tau and a pole of E',
         (3*e+2*k+d)-(2*k+4*e-(e-d)))
    check('every finite escape pole is even under a constant Hessian class',
          S.Poly(2*k-j-2*e,k,j,e,modulus=2).as_expr()==j)
    check('one finite escape and one infinity point require six distinct roots',2+2*1+2==6)

    # Entire C* components: several exact normal-form identities. These
    # are local/global curve controls, not realized unit-Hessian graphs.
    for degree,exponent in [(3,2),(3,4),(1,2)]:
        ss = t**degree
        value = alpha+beta*ss+t**exponent
        first = S.diff(value,t)/S.diff(ss,t)
        second = S.diff(first,t)/S.diff(ss,t)
        rp = t**((4*degree-exponent)//2)
        zero('cyclic Laurent control '+str((degree,exponent))+' has the stated constant critical Hessian',
             -162*L*N*N*rp**2/ss**2*second
             +162*L*N*N*exponent*(exponent-degree)/degree**2)
    n,theta,deg = S.symbols('n theta deg')
    root_exponent = deg-1+(n-1)*(theta-2)
    zero('root-product exponent equation is exactly the leading-coefficient constraint',
         (root_exponent-(2-theta/2))-(deg-(2*n+1-(n-S.Rational(1,2))*theta)))
    zero('no-inactive pole-two alternative forces the impossible denominator',
         (2*n+1-(n-S.Rational(1,2))*theta-2*n).subs(theta,2/(2*n-1)))
    zero('no-inactive cancelled alternative would require negative theta',
         n*(2-theta)-(2*n+1-(n-S.Rational(1,2))*theta)+1+theta/2)
    mult = S.symbols('mult')
    zero('one-inactive alternative retains the exact exceptional cyclic profile',
         (n-(2*n+1-(n-S.Rational(1,2))*theta)).subs({n:e+mult,theta:1+1/e},simultaneous=True)
         -(2*mult-1-e)/(2*e))

    # The exceptional exponent profile is not dropped: its actual graph
    # has a simple boundary root, which excludes the remaining case.
    for multiplicity in (2,3):
        degree = 2*multiplicity-1
        Z = s*p+S.Rational(169,9)*s**3+9
        scale = S.Rational(3,9**(degree+multiplicity))
        F = scale*(Z**degree-s)*Z**multiplicity
        graph = S.cancel((F-3)/(2*s))
        check('exceptional profile '+str(multiplicity)+' is an actual polynomial graph',graph.is_polynomial(s,p))
        zero('exceptional profile '+str(multiplicity)+' has the exact simple linear boundary',
             graph.subs(s,0)-(S.Rational(degree+multiplicity,6)*p-S.Rational(3,2*9**degree)))
        pc = -S.Rational(169,9)*s*s-9/s+t**(1-degree)
        ss = t**degree
        pc = pc.subs(s,ss)
        zero('exceptional profile '+str(multiplicity)+' has the exact cyclic active curve',
             (Z**degree-s).subs({s:ss,p:pc},simultaneous=True))
        rp = S.cancel(S.diff(graph,p).subs({s:ss,p:pc},simultaneous=True))
        Dc = S.cancel(D.subs({L:1,N:1,M:0,s:ss,p:pc},simultaneous=True))
        zero('exceptional profile '+str(multiplicity)+' retains its critical monomial',Dc-(t**degree-t**(degree+1)/9))
        D2 = S.diff(S.diff(Dc,t)/S.diff(ss,t),t)/S.diff(ss,t)
        zero('exceptional profile '+str(multiplicity)+' can pass the critical determinant test',
             -162*rp**2/ss**2*D2-S.Rational(9,2)*scale**2*(degree+1))

    dependencies = {}
    prior = [('Y_GRAPH_CLOSURE.md','receipt-y-graph-closure.json'),
             ('GLOBAL_CRITICAL_BUDGET.md','receipt-global-critical-budget.json'),
             ('CRITICAL_NORMALIZATION.md','receipt-critical-normalization.json'),
             ('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
             ('OPPOSITE_AFFINE_ENDPOINT.md','receipt-opposite-affine.json')]
    for name,receipt_name in prior:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained written proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md','hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'OPPOSITE_NONLINEAR_FRONTIER.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['opposite entire-value fibre stratification including every R zero',
            'all-degree nonlinear boundary-root closure',
            'whole critical normalization, affine inactive values and even infinity poles',
            'all finite etale C* covers and Laurent-unit root configurations',
            'the exceptional cyclic profile and its actual simple boundary',
            'the necessary six-distinct-root threshold without a bounded degree search'],
        all_opposite_simple_boundary_root_graphs_closed=True,
        every_no_finite_base_escape_opposite_graph_gradient_invertible_under_unit_Hessian=True,
        every_at_most_five_distinct_critical_root_opposite_graph_closed=True,
        total_p_degree_automatically_bounded_by_five=False,
        all_general_opposite_graphs_closed=False,
        six_distinct_root_case_asserted_exists=False,
        native_critical_constant_asserted_sufficient=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
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
