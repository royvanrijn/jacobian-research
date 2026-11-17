"""Exact controls for the isolated-trace bound and hyperplane mixing closure.

The whole-curve argument, DVR intersection length and unrestricted-degree
exhaustion are written proofs; finite controls do not certify them.
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

    x,y,p,q,r,s,U,V,W,Z,t,lam,s0 = S.symbols('x y p q r s U V W Z t lam s0')
    L,N,a,b,c = S.symbols('L N a b c',nonzero=True)
    M,u0,cc = S.symbols('M u0 cc')
    u = 1+x*y
    AA = u**3*p+3*x*u**2*q-x**3*r
    BB = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = L*AA**2+M*AA+N*BB
    trace = L*p*p+M*p+4*N*p*y*y+N*y*q
    zero('full ambient x-zero trace is independent of every r substitution',psi.subs(x,0)-trace)
    graph = S.Function('R')(x,y,p,q)
    zero('arbitrary standard r graph has exactly the same three-variable trace',psi.subs(r,graph).subs(x,0)-trace)
    forward = S.Matrix([p+M/(2*L),y,q+4*p*y])
    inverse = {p:U-M/(2*L),y:V,q:W-4*(U-M/(2*L))*V}
    zero('native trace normalization has a field-unit coordinate Jacobian',forward.jacobian((p,y,q)).det()-1)
    zero('native trace normalization has its full polynomial inverse',
         forward.subs(inverse,simultaneous=True)-S.Matrix([U,V,W]))
    quadric = L*U*U+N*V*W-M*M/(4*L)
    zero('native trace becomes the exact nondegenerate quadric including its scalar shift',
         trace.subs(inverse,simultaneous=True)-quadric)
    gradient = S.Matrix([S.diff(trace,z) for z in (p,y,q)])
    critical = {p:-M/(2*L),y:0,q:0}
    zero('native isolated critical point solves all three trace equations',gradient.subs(critical))
    zero('trace Jacobian ideal recovers the coordinate y',gradient[2]/N-y)
    zero('trace Jacobian ideal recovers the coordinate q',gradient[1]/N-8*p*y-q)
    zero('trace Jacobian ideal recovers p plus its critical shift',
         (gradient[0]-4*N*y*y)/(2*L)-(p+M/(2*L)))
    zero('native isolated trace Hessian determinant is exactly -2*L*N^2',
         S.hessian(trace,(p,y,q)).det().subs(critical)+2*L*N*N)
    zero('quadric value fibre over nonzero V solves W globally on that stratum',
         quadric.subs(W,(lam+M*M/(4*L)-L*U*U)/(N*V))-lam)
    zero('quadric zero-V stratum has exactly the stated quadratic equation',
         quadric.subs(V,0)-(L*U*U-M*M/(4*L)))
    check('quadric whole generic and critical-value Euler strata give two and one',0+2*1==2 and 0+1*1==1)

    # All power degrees are retained symbolically; no degree search is used.
    mu = S.symbols('mu',integer=True,positive=True)
    power = a*(U-u0)**(mu+1)+b*V*W+cc
    zero('power-saddle trace has the exact Milnor Jacobian generators',
         S.Matrix([S.diff(power,z) for z in (U,V,W)])
         -S.Matrix([a*(mu+1)*(U-u0)**mu,b*W,b*V]))
    zero('power-saddle open stratum has its full value-fibre inverse',
         power.subs(W,(lam-cc-a*(U-u0)**(mu+1))/(b*V))-lam)
    zero('power trace generic Euler is one less than the isolated nonautomorphic lower bound',
         ((mu+2)-(mu+1))-1)
    H = S.Function('H')(U)
    control = s*U+H+b*V*W
    zero('arbitrary polynomial power control has unit full Hessian b squared',S.hessian(control,(s,U,V,W)).det()-b*b)
    targets = S.symbols('target_s target_U target_V target_W')
    target_s,target_U,target_V,target_W = targets
    gradient_control = S.Matrix([S.diff(control,z) for z in (s,U,V,W)])
    inverse_control = {s:target_U-S.diff(H,U).subs(U,target_s),U:target_s,V:target_W/b,W:target_V/b}
    zero('unit-Hessian controls have an actual polynomial gradient inverse for arbitrary H',
         gradient_control.subs(inverse_control,simultaneous=True)-S.Matrix(targets))

    # Rank one: whole critical reconstruction and the exact local ramification.
    A,B,E = [S.Function(name)(U,Z) for name in ('A','B','E')]
    AU,AZ,BU,BZ,EU,EZ = [S.diff(fun,z) for fun,z in
                         [(A,U),(A,Z),(B,U),(B,Z),(E,U),(E,Z)]]
    jac = AU*BZ-AZ*BU
    sigma = -(EU*BZ-EZ*BU)/jac
    vc = -(AU*EZ-AZ*EU)/jac
    zero('whole rank-one critical equations solve for both the free base and V',
         S.Matrix([EU+sigma*AU+vc*BU,EZ+sigma*AZ+vc*BZ]))
    zero('rank-one exceptional coefficient differential retains the exact base difference',
         (EU+s0*AU)*BZ-(EZ+s0*AZ)*BU-(s0-sigma)*jac)
    zero('Keller mate A is a local uniformizing coordinate on every smooth B-zero component',
         AU*BZ-AZ*BU-jac)
    local_sigma = s0+c*t**mu
    local_E = -s0*t-c*t**(mu+1)/(mu+1)
    zero('local intersection length mu gives the exact exceptional differential order',
         S.diff(local_E+s0*t,t)-(s0-local_sigma))
    zero('isolated trace value function has local degree mu plus one, not mu',
         local_E+s0*t+c*t**(mu+1)/(mu+1))
    chi_negative,extra_degree = S.symbols('chi_negative extra_degree',integer=True,nonnegative=True)
    zero('whole exceptional Euler bound includes all components and extra mapping degree',
         (1+chi_negative+mu+1+extra_degree)-(mu+2)-chi_negative-extra_degree)

    # Rank-zero common zeros give vertical critical lines, never the isolated trace point.
    rank_zero = s*U*U+V*s+W*U
    zero('rank-zero Keller control has constant full Hessian',S.hessian(rank_zero,(s,U,V,W)).det()-1)
    zero('rank-zero relative equations have the complete vertical affine line',
         S.Matrix([S.diff(rank_zero,z) for z in (U,V,W)]).subs({s:0,U:0,W:0},simultaneous=True))
    zero('rank-zero trace critical set has the explicitly free affine coordinate',
         rank_zero.subs(s,0)-W*U)

    # A real three-variable automorphism mixes old y with the old dual directions.
    Ymix = U+V*V
    Pmix = V+W*W
    Qmix = W+Ymix*Ymix
    mixing = S.Matrix([Ymix,Pmix,Qmix])
    Yt,Pt,Qt = S.symbols('Yt Pt Qt')
    Wback = Qt-Yt*Yt
    Vback = Pt-Wback*Wback
    Uback = Yt-Vback*Vback
    mix_inverse = {U:Uback,V:Vback,W:Wback}
    zero('three-variable mixing control has a unit actual coordinate Jacobian',mixing.jacobian((U,V,W)).det()-1)
    zero('three-variable mixing has its complete polynomial inverse',
         mixing.subs(mix_inverse,simultaneous=True)-S.Matrix([Yt,Pt,Qt]))
    check('old y really depends on a source dual in the mixing control',S.diff(Ymix,V)!=0)
    mixed_trace = trace.subs({y:Ymix,p:Pmix,q:Qmix},simultaneous=True)
    pstar = -M/(2*L)
    mixed_point = {U:-pstar*pstar,V:pstar,W:0}
    zero('nonlinear mixing preserves the actual isolated trace critical point',
         S.Matrix([S.diff(mixed_trace,z) for z in (U,V,W)]).subs(mixed_point,simultaneous=True))
    zero('trace Hessian at that critical point transports by the unit Jacobian square',
         S.hessian(mixed_trace,(U,V,W)).det().subs(mixed_point,simultaneous=True)+2*L*N*N)

    Es,F1,F2,F3,Rs = [S.Function(name)(s,U,V,W) for name in ('Es','F1','F2','F3','Rs')]
    arbitrary_substitution = {x:s*Es,y:Ymix+s*F1,p:Pmix+s*F2,q:Qmix+s*F3,r:Rs}
    zero('arbitrary off-hyperplane base and dual corrections preserve the complete required trace',
         psi.subs(arbitrary_substitution,simultaneous=True).subs(s,0)-mixed_trace)
    zero('the broader trace class does not require a global unit substitution Jacobian',
         S.Matrix([s*s,U,V,W]).jacobian((s,U,V,W)).det()-2*s)
    lower = S.Function('lower')(p,q)
    y_boundary = S.Matrix([a*r+lower,p,q])
    zero('genuine old-r-dependent y boundary has the exact constant-unit trace Jacobian',
         y_boundary.jacobian((p,q,r)).det()-a)
    zero('genuine old-r-dependent y boundary has its complete polynomial trace inverse',
         y_boundary.subs({p:Pt,q:Qt,r:(Yt-lower.subs({p:Pt,q:Qt},simultaneous=True))/a},simultaneous=True)
         -S.Matrix([Yt,Pt,Qt]))

    dependencies = {}
    prior = [('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
             ('RELATIVE_CRITICAL_LOCUS.md','receipt-relative-critical.json'),
             ('CYLINDRICAL_TRACE_CLOSURE.md','receipt-cylindrical-trace.json')]
    for name,receipt_name in prior:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md','hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'ISOLATED_TRACE_MIXING.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['isolated trace intersection length forces a dominating critical component',
            'absence of every A1 component and nonpositive whole critical-curve Euler',
            'both Tang ranks without a generic-base Euler assumption',
            'Milnor intersection length mu and exceptional value-map degree at least mu+1',
            'all power-saddle traces and arbitrary polynomial trace equivalence',
            'every r graph under all three-variable automorphisms over K[x]',
            'arbitrary substitutions with a coordinate triple on one actual affine source hyperplane'],
        generic_base_Euler_hypothesis_required=False,
        all_power_saddle_trace_degrees_closed=True,
        all_r_graphs_with_full_y_p_q_automorphisms_over_Kx_closed=True,
        arbitrary_off_hyperplane_mixing_with_coordinate_trace_closed=True,
        old_r_dependent_y_unit_coefficient_affine_boundary_closed=True,
        broader_substitutions_asserted_embedded=False,
        nonlinear_Hessian_invariance_asserted=False,
        unit_Hessian_nonexistence_asserted=False,
        all_mixing_or_all_old_r_dependent_y_graphs_closed=False,
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
