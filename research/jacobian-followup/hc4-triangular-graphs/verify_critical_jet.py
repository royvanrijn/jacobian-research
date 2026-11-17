"""Exact controls for the critical-coordinate first-jet obstruction.

The plane-rank exhaustion and imported line theorem are written mathematics,
not conclusions from these finite algebraic controls.
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

    x,y,p,q,r,s,t,P,Q,U,V,T = S.symbols('x y p q r s t P Q U V T')
    L,N = S.symbols('L N', nonzero=True)
    M = S.symbols('M')
    u = 1+x*y
    A = u**3*p+3*x*u**2*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = L*A*A+M*A+N*B
    origin = {x:0, y:0}
    Hp = 2*L*p+M
    zero('ambient origin trace is the retained-p square', psi.subs(origin)-(L*p*p+M*p))
    for variable, expected in [(x,3*Hp*q+2*N*r),(y,N*q),(p,Hp),(q,0),(r,0)]:
        zero('ambient first jet in '+str(variable), S.diff(psi,variable).subs(origin)-expected)

    # Arbitrary base and substitution first jets, without truncating them away.
    aa,bb,cc,dd,ps,pt = S.symbols('aa bb cc dd ps pt')
    base_jac = S.Matrix([[aa,bb],[cc,dd]])
    delta = base_jac.det()
    ex = S.Matrix([dd,-cc])/delta
    ey = S.Matrix([-bb,aa])/delta
    zero('constant parameter direction realizing old x direction',base_jac*ex-S.Matrix([1,0]))
    zero('constant parameter direction realizing old y direction',base_jac*ey-S.Matrix([0,1]))
    ds = Hp*(ps+3*q*aa)+N*q*cc+2*N*r*aa
    dt = Hp*(pt+3*q*bb)+N*q*dd+2*N*r*bb
    zero('full q-graph derivative formula after arbitrary base first jet',
         ds*ex[0]+dt*ex[1]-Hp*(ps*ex[0]+pt*ex[1]+3*q)-2*N*r)
    zero('full r-graph derivative formula after arbitrary base first jet',
         ds*ey[0]+dt*ey[1]-Hp*(ps*ey[0]+pt*ey[1])-N*q)

    F,G,B0 = S.Function('F')(T,U),S.Function('G')(T,U),S.Function('B0')(T,U)
    normal_one = B0+s*F+V*G
    pair_one = S.Matrix([F,G])
    zero('rank-one full Hessian determinant is the actual residual Jacobian squared',
         S.hessian(normal_one,(s,V,T,U)).det()-pair_one.jacobian((T,U)).det()**2)
    k = S.Function('k')(T,U)
    C = T*k
    zero('rank-one zero-line Jacobian factors into two polynomial units',
         S.Matrix([F,C]).jacobian((T,U)).det().subs(T,0)
         +S.diff(F,U).subs(T,0)*k.subs(T,0))
    F2,G2,H0 = S.Function('F2')(P,Q),S.Function('G2')(P,Q),S.Function('H0')(P,Q)
    normal_two = H0+s*F2+t*G2
    zero('rank-two full Hessian determinant does not depend on the dual trace',
         S.hessian(normal_two,(s,t,P,Q)).det()-S.Matrix([F2,G2]).jacobian((P,Q)).det()**2)

    # Constant base-linear dual translations preserve the first-jet class
    # modulo H'(a), without a nonlinear four-variable Hessian change.
    a = P+Q**2
    b = Q+a**2
    H = L*a*a+M*a
    ca,cb,cd,ce = S.symbols('ca cb cd ce')
    direct = (H+s*F2+t*G2).subs({P:P+ca*s+cb*t,Q:Q+cd*s+ce*t},simultaneous=True)
    for variable, cP,cQ,coefficient in [(s,ca,cd,F2),(t,cb,ce,G2)]:
        zero('constant dual translation changes '+str(variable)+' jet by H-prime times a polynomial',
             S.diff(direct,variable).subs({s:0,t:0})-coefficient
             -(2*L*a+M)*(cP*S.diff(a,P)+cQ*S.diff(a,Q)))
    coordinate_pair = S.Matrix([a,b])
    inverse_pair = S.Matrix([U-(V-U**2)**2,V-U**2])
    zero('nonlinear initial coordinate pair has unit Jacobian',coordinate_pair.jacobian((P,Q)).det()-1)
    zero('initial coordinate pair explicit inverse left',
         inverse_pair.subs({U:a,V:b},simultaneous=True)-S.Matrix([P,Q]))
    zero('initial coordinate pair explicit inverse right',
         coordinate_pair.subs(dict(zip((P,Q),inverse_pair)),simultaneous=True)-S.Matrix([U,V]))
    cjet = S.Function('cjet')(P,Q)
    mu = S.Function('mu')(a)
    lam = S.symbols('lam', nonzero=True)
    jet = (2*L*a+M)*cjet+lam*b+mu
    normalized_jet = jet.subs(dict(zip((P,Q),inverse_pair)),simultaneous=True)
    critical = -M/(2*L)
    zero('critical coordinate line makes the tested coefficient affine injective',
         normalized_jet.subs(U,critical)-(lam*V+S.Function('mu')(critical)))

    # Actual nonlinear base and dual controls outside the old square-axis
    # class. No unit-Hessian claim is made about this arbitrary graph control.
    old_x = s+t**2
    old_y = t+old_x**2
    beta = S.Matrix([old_x,old_y])
    beta_inverse = S.Matrix([x-(y-x*x)**2,y-x*x])
    zero('nonlinear base map has unit Jacobian',beta.jacobian((s,t)).det()-1)
    check('old y is genuinely nonlinear in the parameter base',S.Poly(old_y,s,t).total_degree()>1)
    zero('nonlinear base inverse left',beta_inverse.subs({x:old_x,y:old_y},simultaneous=True)-S.Matrix([s,t]))
    zero('nonlinear base inverse right',beta.subs(dict(zip((s,t),beta_inverse)),simultaneous=True)-S.Matrix([x,y]))
    pp = P+Q*Q+s*Q**3
    rr = Q+pp*pp+t*pp**3
    phi = S.Matrix([pp,rr])
    zero('relative nonlinear retained p,r map has unit dual Jacobian',phi.jacobian((P,Q)).det()-1)
    zero('relative nonlinear dual specialization is exactly the initial coordinate pair',phi.subs({s:0,t:0})-coordinate_pair)
    R = S.Function('R')(s,t,P,Q)
    full = psi.subs({x:old_x,y:old_y,p:pp,q:R,r:rr},simultaneous=True)
    full_zero = {s:0,t:0}
    zero('actual q graph on nonlinear base and dual maps has the declared origin trace',full.subs(full_zero)-H)
    actual_jet = S.diff(full,s).subs(full_zero)
    zero('actual full q-graph first jet retains every graph term',
         actual_jet-(2*L*a+M)*(Q**3+3*R.subs(full_zero))-2*N*b)
    zero('actual q-graph critical coordinate jet is injective on the full line',
         actual_jet.subs(dict(zip((P,Q),inverse_pair)),simultaneous=True).subs(U,critical)-2*N*V)
    weak_p = a+s*P*P
    weak_r = Q+weak_p**2+t*Q**3
    weak_pair = S.Matrix([weak_p,weak_r])
    zero('single-slice specialization remains an automorphism without a global dual inverse',
         weak_pair.subs(full_zero)-coordinate_pair)
    check('single-slice control has nonconstant full retained-dual Jacobian',
          S.Poly(weak_pair.jacobian((P,Q)).det(),s,t,P,Q).total_degree()>0)
    weak_base = S.Matrix([s+s*s,t])
    zero('base substitution need only have an invertible Jacobian at the chosen point',
         weak_base.jacobian((s,t)).det().subs(full_zero)-1)
    check('weaker base substitution is not a polynomial automorphism',
          S.Poly(weak_base.jacobian((s,t)).det(),s,t).total_degree()>0)

    z1,z2,z3,z4 = S.symbols('z1 z2 z3 z4')
    outputs = S.Matrix([z1,z2,z3,z4])
    # General cotangent reconstruction checked on a nonlinear coordinate
    # initial trace and an explicitly automorphic residual pair.
    f_control = H+s*b+t*a
    gradient = S.Matrix([S.diff(f_control,z) for z in (s,t,P,Q)])
    rec_a = z2
    rec_Q = z1-z2*z2
    rec_s = z4-2*rec_Q*z3
    rec_t = z3-2*L*z2-M-2*z2*rec_s
    inverse = S.Matrix([rec_s,rec_t,rec_a-rec_Q**2,rec_Q])
    zero('nonlinear-coordinate critical-jet unit-Hessian control',S.hessian(f_control,(s,t,P,Q)).det()-1)
    zero('nonlinear-coordinate gradient inverse left',inverse.subs(dict(zip(outputs,gradient)),simultaneous=True)-S.Matrix([s,t,P,Q]))
    zero('nonlinear-coordinate gradient inverse right',gradient.subs(dict(zip((s,t,P,Q),inverse)),simultaneous=True)-outputs)
    # A higher-degree H checks that the criterion is not a quadratic-only
    # calculation. Its generality still comes from the written argument.
    quartic = a**4+a+s*b+t*a
    quartic_gradient = S.Matrix([S.diff(quartic,z) for z in (s,t,P,Q)])
    quartic_inverse = S.Matrix([rec_s,z3-4*z2**3-1-2*z2*rec_s,rec_a-rec_Q**2,rec_Q])
    zero('quartic initial trace still gives unit Hessian',S.hessian(quartic,(s,t,P,Q)).det()-1)
    zero('quartic-trace gradient inverse left',quartic_inverse.subs(dict(zip(outputs,quartic_gradient)),simultaneous=True)-S.Matrix([s,t,P,Q]))
    zero('quartic-trace gradient inverse right',quartic_gradient.subs(dict(zip((s,t,P,Q),quartic_inverse)),simultaneous=True)-outputs)

    # Rank-one reconstruction is checked in both directions as well.
    f_rank_one = s*(U+T*T)+V*T+L*U*U+M*U+T*U**3
    variables = (s,V,T,U)
    gradient_one = S.Matrix([S.diff(f_rank_one,z) for z in variables])
    rec_T,rec_U = z2,z1-z2*z2
    rec_S = z4-2*L*rec_U-M-3*rec_T*rec_U**2
    rec_V = z3-2*rec_S*rec_T-rec_U**3
    inverse_one = S.Matrix([rec_S,rec_V,rec_T,rec_U])
    zero('rank-one explicit control has unit Hessian',S.hessian(f_rank_one,variables).det()-1)
    zero('rank-one explicit gradient inverse left',inverse_one.subs(dict(zip(outputs,gradient_one)),simultaneous=True)-S.Matrix(variables))
    zero('rank-one explicit gradient inverse right',gradient_one.subs(dict(zip(variables,inverse_one)),simultaneous=True)-outputs)

    # The first-jet criterion is sufficient, not automatic from unit Hessian
    # and a coordinate square. The critical-line coefficient image is curved.
    control_F = Q+(P+Q*Q)**2
    control_G = P+Q*Q
    control_pair = S.Matrix([control_F,control_G])
    zero('negative first-jet control is still a plane automorphism',control_pair.jacobian((P,Q)).det()+1)
    control_inv = S.Matrix([V-(U-V*V)**2,U-V*V])
    zero('negative control plane inverse left',control_inv.subs({U:control_F,V:control_G},simultaneous=True)-S.Matrix([P,Q]))
    zero('negative control plane inverse right',control_pair.subs(dict(zip((P,Q),control_inv)),simultaneous=True)-S.Matrix([U,V]))
    f_negative = P*P+s*control_F+t*control_G
    zero('negative first-jet control really has unit Hessian',S.hessian(f_negative,(s,t,P,Q)).det()-1)
    direction_1,direction_2 = S.symbols('direction_1 direction_2')
    critical_row = S.Poly((direction_1*control_F+direction_2*control_G).subs(P,0),Q)
    zero('critical-line quartic coefficient forces first direction coefficient zero',critical_row.coeff_monomial(Q**4)-direction_1)
    zero('critical-line quadratic coefficient forces second direction coefficient zero',critical_row.coeff_monomial(Q**2)-direction_2)
    check('no nonzero base direction has the prescribed affine critical jet in the control',
          S.solve([critical_row.coeff_monomial(Q**4),critical_row.coeff_monomial(Q**2)],(direction_1,direction_2))=={direction_1:0,direction_2:0})

    # The primary theorem is retained byte-for-byte; these checks verify
    # provenance, not the external theorem's proof.
    manifest_name = 'sources/critical-jet-import.json'
    manifest = json.loads((HERE/manifest_name).read_text())
    source_name = 'sources/'+manifest['retained_file']
    source_hash = hashlib.sha256((HERE/source_name).read_bytes()).hexdigest()
    expected = '4cce6a7334640ea47036b1fb8c676e28da7946dbd15d0a48868b30d94c68585a'
    check('primary line theorem manuscript byte pin',source_hash==manifest['sha256']==expected)
    check('primary line theorem location and characteristic-zero hypotheses retained',
          manifest['locations']==['Theorem 1.1','Lemma 2.3'] and 'characteristic-zero' in manifest['used_implication'])
    dependencies = {}
    pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md','jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang dependency '+name,dependencies[name]==pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in [Path(__file__).name,'CRITICAL_COORDINATE_JET.md',manifest_name,source_name]},
        dependencies=dependencies,
        imported_theorems=['HC4MYGL0: Gwozdziewicz injectivity on one line','HCTANG2: reproduced Tang alternative'],
        written_not_finite_exhaustion=['three plane-projection ranks','coordinate-critical-line injectivity','rank-one zero-line argument','characteristic-zero descent'],
        q_chart_all_base_and_relative_dual_automorphisms_closed=True,
        r_chart_all_base_and_relative_dual_automorphisms_closed=True,
        single_specialized_retained_pair_inverse_sufficient=True,
        global_dual_inverse_required=False,
        merely_Keller_specialized_pair_covered=False,
        p_chart_singleton_zero_properness_resolved=False,
        base_eliminated_charts_closed=False,
        arbitrary_nonlinear_base_dual_mixing_closed=False,
        general_HC4_or_JC2_resolved=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
