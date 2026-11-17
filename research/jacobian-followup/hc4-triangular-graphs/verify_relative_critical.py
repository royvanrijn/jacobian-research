"""Exact algebra for the relative critical-locus Hessian obstruction.

The arbitrary-degree rational-function argument is a written proof. The
controls do not assert that passing a critical branch gives a unit Hessian.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent


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
        values = list(expression) if isinstance(expression,S.MatrixBase) else [expression]
        check(name,all(S.cancel(value.doit())==0 for value in values))

    s,P,Q,Z,z1,z2,z3,t = S.symbols('s P Q Z z1 z2 z3 t')
    outer = s*z1+z2*z3+z1**3
    phi = S.Matrix([P+Q*Q+s*Z,Q+Z*Z,Z])
    phi_inverse = S.Matrix([z1-(z2-z3*z3)**2-s*z3,z2-z3*z3,z3])
    source = (s,P,Q,Z)
    target = (s,z1,z2,z3)
    full_phi = S.Matrix([s,*phi])
    jac = full_phi.jacobian(source)
    zero('relative automorphism has a constant unit full Jacobian',jac.det()-1)
    zero('relative automorphism polynomial inverse left',
         phi_inverse.subs(dict(zip((z1,z2,z3),phi)),simultaneous=True)-S.Matrix([P,Q,Z]))
    zero('relative automorphism polynomial inverse right',
         phi.subs(dict(zip((P,Q,Z),phi_inverse)),simultaneous=True)-S.Matrix([z1,z2,z3]))
    inner = outer.subs(dict(zip((z1,z2,z3),phi)),simultaneous=True)
    grad_outer = S.Matrix([S.diff(outer,z) for z in (z1,z2,z3)]).subs(dict(zip((z1,z2,z3),phi)),simultaneous=True)
    grad_inner = S.Matrix([S.diff(inner,z) for z in (P,Q,Z)])
    zero('relative gradient ideal transforms by the invertible dual Jacobian',
         grad_inner-phi.jacobian((P,Q,Z)).T*grad_outer)
    Houter = S.hessian(outer,target).subs(dict(zip((z1,z2,z3),phi)),simultaneous=True)
    Hinner = S.hessian(inner,source)
    remainder = S.zeros(4)
    for index in range(3):
        remainder += grad_outer[index]*S.hessian(phi[index],source)
    zero('actual full Hessian chain rule has only relative-gradient corrections',
         Hinner-jac.T*Houter*jac-remainder)
    zero('outer control has unit Hessian',S.hessian(outer,target).det()-1)
    critical_control = {s:-3*t*t,P:t,Q:0,Z:0}
    zero('transported critical curve has zero relative gradient',grad_inner.subs(critical_control))
    zero('Hessian determinant agrees on the entire critical curve',Hinner.subs(critical_control).det()-1)
    actual_det = S.expand(Hinner.det())
    ideal = S.groebner(list(grad_inner),Z,Q,s,P)
    zero('determinant congruence holds modulo the actual gradient ideal as a scheme',ideal.reduce(actual_det-1)[1])
    zero('off-critical nonlinear change really changes the Hessian determinant',Hinner.subs({s:0,P:1,Q:0,Z:0}).det()+323)

    x,y,p,q,r,T = S.symbols('x y p q r T')
    L,N = S.symbols('L N',nonzero=True)
    M = S.symbols('M')
    u = 1+x*y
    A = u**3*p+3*x*u*u*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = L*A*A+M*A+N*B
    R = S.Function('R')(x,p)
    ur = 1+x*R
    r_from_T = (ur**3*p+3*x*ur*ur*q-T)/x**3
    rational = L*T*T+(M-N*(2-3*x*R)/x**2)*T+N*ur*(ur+1)*p/x**2+2*N*(2*ur+1)*q/x
    zero('general y=R(x,p) rational computation chart is the original potential',
         psi.subs({y:R,r:r_from_T},simultaneous=True)-rational)
    graph_A = A.subs(y,R)
    zero('general native-to-computation relative Jacobian is x cubed',
         S.Matrix([p,graph_A,q]).jacobian((p,q,r)).det()-x**3)
    curve = 2*x*R+3
    zero('critical q equation is precisely the whole displayed curve',S.diff(rational,q)-2*N*curve/x)
    Rp,Rx,Rpp,Rxp,Rxx = S.symbols('Rp Rx Rpp Rxp Rxx')
    jets = {S.diff(R,p):Rp,S.diff(R,x):Rx,S.diff(R,p,2):Rpp,S.diff(R,x,p):Rxp,S.diff(R,x,2):Rxx}
    Tc = (13*N/(2*x*x)-M)/(2*L)
    qc = -3*Tc/(4*x)+1/(16*x*x*Rp)

    def at_critical(expression, include_q=True):
        substitutions = {R:-S.Rational(3,2)/x,T:Tc}
        if include_q:
            substitutions[q] = qc
        return S.cancel(expression.xreplace(jets).subs(substitutions,simultaneous=True))

    for variable in (p,T,q):
        zero('actual critical '+str(variable)+' equation',at_critical(S.diff(rational,variable)))
    zero('R_p zero has no critical point on the curve in this chart',
         at_critical(S.diff(rational,p),include_q=False).subs(Rp,0)+N/(4*x*x))
    D = -N*p/(4*x*x)-(M-13*N/(2*x*x))**2/(4*L)
    zero('actual critical value on every general curve component',at_critical(rational)-D)
    p1 = (3/(2*x*x)-Rx)/Rp
    p2 = (-3/x**3-Rxx-2*Rxp*p1-Rpp*p1*p1)/Rp
    zero('curve derivation preserves its localized defining ideal',
         (S.diff(curve,x)+S.diff(curve,p)*(3/(2*x*x)-S.diff(R,x))/S.diff(R,p))-curve/x)
    Hd = S.hessian(rational,(p,T,q)).applyfunc(at_critical)
    zero('critical dual Hessian determinant includes the exact R_p square',Hd.det()+32*L*N*N*Rp*Rp)
    Hfull = S.hessian(rational,(x,p,T,q)).applyfunc(at_critical)
    schur = Hfull[0,0]-(Hfull[0,1:]*Hd.inv()*Hfull[1:,0])[0]
    along = S.diff(D,x,2)+2*S.diff(D,x,p)*p1+S.diff(D,p)*p2
    zero('Schur complement is the second derivative along the whole curve',schur-along)
    zero('native critical determinant has the exact -32 L N squared x sixth factor',
         x**6*Hd.det()*schur+32*L*N*N*x**6*Rp*Rp*along)

    # Specialize to arbitrary polynomial coefficient functions a and b.
    # The polynomial-degree hypothesis is used in the written infinity proof.
    af,bf = S.Function('a')(x),S.Function('b')(x)
    pc = (-3/(2*x)-bf)/af
    Dab = N*(3+2*x*bf)/(8*af*x**3)-(M-13*N/(2*x*x))**2/(4*L)
    zero('affine-p graph has the declared rational critical curve',2*x*(af*pc+bf)+3)
    zero('affine-p graph critical value is exactly the stated rational function',D.subs(p,pc)-Dab)
    affine_jets = {Rp:af,Rx:S.diff(af,x)*p+S.diff(bf,x),Rpp:0,
        Rxp:S.diff(af,x),Rxx:S.diff(af,x,2)*p+S.diff(bf,x,2)}
    zero('general curve derivative specializes to the ordinary rational second derivative',
         along.subs(affine_jets,simultaneous=True).subs(p,pc)-S.diff(Dab,x,2))
    numerator = 2*L*N*x*(3+2*x*bf)-af*(2*M*x*x-13*N)**2
    zero('critical value has the universal degree d+4 denominator certificate',Dab-numerator/(16*L*af*x**4))
    order = S.symbols('order',integer=True,positive=True)
    zero('formal leading infinity order raises by exactly two under two derivatives',
         S.simplify(S.diff(x**(-order),x,2)*x**(order+2)-order*(order+1)))
    a0 = S.symbols('a0',nonzero=True)
    non_b_part = 3*N/(8*a0*x**3)-(M-13*N/(2*x*x))**2/(4*L)
    zero('constant-a impossible x-minus-five coefficient is nonzero and fixed',
         S.Poly(S.expand(x**6*S.diff(non_b_part,x,2)),x).coeff_monomial(x)-9*N/(2*a0))
    degree = S.symbols('degree',integer=True,nonnegative=True)
    zero('arbitrary polynomial b monomial contributes only powers degree minus four',
         S.simplify(S.diff(x**degree/x**2,x,2)*x**(4-degree)-(degree-2)*(degree-3)))

    # Opposite-chart control: a passed critical determinant is not a unit
    # determinant on the full source, nor a new reparameterization witness.
    control_coeff = S.Rational(1013,54)
    opposite = L*T*T+(M-N*(2-3*x*y)/x**2)*T+N*u*(u+1)*(x-control_coeff*y*y)/x**2+2*N*(2*u+1)*q/x
    opposite = opposite.subs({L:1,M:0,N:1})
    xc = -S.Rational(3,2)/y
    opposite_Tc = 13/(4*xc*xc)
    derivative_x = S.cancel(S.diff(opposite,x).subs({x:xc,T:opposite_Tc},simultaneous=True))
    opposite_qc = S.cancel(-derivative_x.subs(q,0)/S.diff(derivative_x,q))
    opposite_point = {x:xc,T:opposite_Tc,q:opposite_qc}
    for variable in (x,T,q):
        zero('opposite graph critical '+str(variable)+' equation',S.diff(opposite,variable).subs(opposite_point,simultaneous=True))
    opposite_D = y/6-y**4/486
    zero('opposite critical value has the tuned quartic form',opposite.subs(opposite_point,simultaneous=True)-opposite_D)
    opposite_dual = S.hessian(opposite,(x,T,q)).subs(opposite_point,simultaneous=True)
    zero('opposite native dual Hessian is minus 162 over y squared',xc**6*opposite_dual.det()+162/y**2)
    opposite_full = S.hessian(opposite,(y,x,T,q)).subs(opposite_point,simultaneous=True)
    opposite_schur = opposite_full[0,0]-(opposite_full[0,1:]*opposite_dual.inv()*opposite_full[1:,0])[0]
    zero('opposite Schur complement equals the second critical-value derivative',opposite_schur-S.diff(opposite_D,y,2))
    zero('opposite computed critical branch has full determinant four',xc**6*opposite_dual.det()*opposite_schur-4)
    opposite_native = psi.subs({L:1,M:0,N:1}).subs(x,p+control_coeff*y*y)
    opposite_H = S.hessian(opposite_native,(y,p,q,r))
    zero('opposite boundary line is genuinely relatively critical',
         S.Matrix([S.diff(opposite_native,z) for z in (p,q,r)]).subs({y:0,p:0,r:0}))
    zero('opposite boundary critical line also has determinant four',opposite_H.subs({y:0,p:0,r:0}).det()-4)
    zero('opposite full source is nevertheless not unit-Hessian',opposite_H.subs({y:0,p:1,q:0,r:0}).det()-1556)

    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in [Path(__file__).name,'RELATIVE_CRITICAL_LOCUS.md']},
        dependencies={},
        written_not_finite_exhaustion=['relative critical-scheme determinant congruence',
            'complete critical curve on x nonzero for every y=R(x,p)',
            'all-degree infinity-order contradiction for nonconstant a',
            'uncancellable x-minus-five coefficient for constant a and arbitrary polynomial b'],
        all_y_equals_a_x_p_plus_b_x_with_nonzero_a_closed_under_relative_dual_automorphisms=True,
        coefficient_degree_or_word_length_bound=None,
        global_invertibility_of_p_to_a_p_plus_b_assumed=False,
        rational_computation_chart_asserted_global=False,
        general_y_R_x_p_graphs_closed=False,
        all_opposite_x_chart_affine_p_graphs_closed=False,
        passed_critical_branch_implies_unit_Hessian=False,
        base_only_unit_ideal_endpoint_resolved=False,
        singleton_zero_properness_resolved=False,
        new_external_classification_inputs=[],
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
