"""Exact controls for the Hodge-polynomial single-trace criterion.

The imported Hodge invariant and all-degree curve/rank arguments are
written mathematical inputs, not consequences of these finite controls.
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

    hu,hv = S.symbols('hu hv')
    tt = hu*hv
    comps,genus,punctures,d0,m,kappa = S.symbols('comps genus punctures d0 m kappa',integer=True)
    curve_E = comps*tt-genus*(hu+hv)+comps-punctures
    trace_E = tt**2-curve_E+d0*tt
    zero('whole rank-one Hodge formula retains genus and every missing point',
         trace_E-(tt**2+(d0-comps)*tt+genus*(hu+hv)+punctures-comps))
    poly = S.Poly(trace_E,hu,hv)
    zero('the hu coefficient is exactly total genus',poly.coeff_monomial(hu)-genus)
    zero('the hv coefficient is exactly total genus',poly.coeff_monomial(hv)-genus)
    zero('the constant coefficient is exactly punctures minus components',poly.coeff_monomial(1)-(punctures-comps))
    zero('zero genus and one puncture per component give the actual affine-line polynomial',
         curve_E.subs({genus:0,punctures:comps},simultaneous=True)-comps*tt)
    zero('rank-zero fixed-base generic value fibre has the opposite nonpositive sign',
         (tt-m)*tt-(tt**2-m*tt))
    positive_k = S.Symbol('positive_k',integer=True,positive=True)
    nonnegative_m = S.Symbol('nonnegative_m',integer=True,nonnegative=True)
    check('a positive coefficient cannot equal the rank-zero nonpositive coefficient',
          S.ask(S.Q.positive(positive_k+nonnegative_m)) is True)
    zero('the kappa-zero coordinate boundary is deliberately outside the positive test',(tt**2+kappa*tt).subs(kappa,0)-tt**2)

    # Euler alone loses precisely the curve information used in the new proof.
    punctured_rational = tt-2
    punctured_elliptic = tt-hu-hv
    check('two legitimate affine curve types have the same Euler but different Hodge polynomials',
          punctured_rational.subs({hu:1,hv:1})==punctured_elliptic.subs({hu:1,hv:1})
          and S.expand(punctured_rational-punctured_elliptic)!=0)
    zero('a rational three-puncture curve leaves a forbidden constant trace coefficient',
         tt**2-punctured_rational+d0*tt-(tt**2+(d0-1)*tt+2))
    zero('a once-punctured genus-one curve leaves forbidden linear Hodge coefficients',
         tt**2-punctured_elliptic+d0*tt-(tt**2+(d0-1)*tt+hu+hv))

    U,V,W,P,Q,s,s0,lam = S.symbols('U V W P Q s s0 lam')
    # Complete the primitive rank-zero row over C[U], keeping the gcd.
    gcd = U**2*(U-1)**3
    row = S.Matrix([U*P+(U*U-1)*Q,P+U*Q])
    row_inverse = {P:U*V-(U*U-1)*W,Q:-V+U*W}
    zero('rank-zero primitive row completion has a unit determinant',row.jacobian((P,Q)).det()-1)
    zero('rank-zero row completion has an actual global trace inverse',
         row.subs(row_inverse,simultaneous=True)-S.Matrix([V,W]))
    D = U**4+U
    raw = D+gcd*(U*P+(U*U-1)*Q)
    zero('rank-zero trace normal form keeps the full common factor',
         raw.subs(row_inverse,simultaneous=True)-(D+gcd*V))
    zero('rank-zero generic trace value leaves the second affine coordinate free',
         (D+gcd*V).subs(V,(lam-D)/gcd)-lam)
    check('rank-zero repeated gcd roots count once in the missing-base set',S.Poly(gcd,U).sqf_part().degree()==2)
    ff,gg = S.Function('ff')(s,U),S.Function('gg')(s,U)
    Fboth,Gboth = (s-s0)*ff,(s-s0)*gg
    zero('both rank-zero rows cannot vanish identically at a base of a Keller pair',
         S.Matrix([Fboth,Gboth]).jacobian((s,U)).det().subs(s,s0))

    # The full rank-one conditional gradient reconstruction uses an actual plane inverse.
    A,B = U+W*W,W
    Dfun = S.Function('Dfun')(U,W)
    potential = Dfun+s*A+V*B
    outputs = S.Matrix([S.diff(potential,z) for z in (s,V,U,W)])
    ta,tb,tu,tw = S.symbols('ta tb tu tw')
    Ub,Wb = ta-tb*tb,tb
    Du = S.diff(Dfun,U).subs({U:Ub,W:Wb},simultaneous=True)
    Dw = S.diff(Dfun,W).subs({U:Ub,W:Wb},simultaneous=True)
    sb = tu-Du
    vb = tw-Dw-2*Wb*sb
    zero('an actual plane inverse recovers every full cotangent gradient output',
         outputs.subs({U:Ub,W:Wb,s:sb,V:vb},simultaneous=True)-S.Matrix([ta,tb,tu,tw]))

    coeff = V*V*(V-1)**3*(V-2)
    H = (V-2)*(V*U**4+(V-1)*U*U)+V
    linear_trace = H+coeff*W
    zero('linear-coefficient trace open stratum has an exact whole-fibre inverse',
         linear_trace.subs(W,(lam-H)/coeff)-lam)
    degrees = [S.Poly(H.subs(V,rho),U).degree() for rho in (0,1,2)]
    check('root budget retains degree drops and constant specializations',degrees==[2,4,0])
    check('root multiplicities are excluded from the distinct-root budget',S.Poly(coeff,V).sqf_part().degree()==3)
    budget = sum(value-1 for value in degrees)
    check('the mixed degree-drop control has a positive exact budget',budget==3)
    zero('linear-coefficient strata give the complete positive Hodge polynomial',
         (tt-3)*tt+sum(degrees)*tt-(tt**2+budget*tt))
    cubic = U**3-3*U
    check('arbitrary univariate traces include multiple distinct critical values',
          [cubic.subs(U,value) for value in (-1,1)]==[2,-2]
          and all(S.diff(cubic,U).subs(U,value)==0 for value in (-1,1)))
    check('the cubic example was not excluded by the earlier one-point Euler bound',S.degree(cubic,U)==1+2)

    # Native y boundary and the two nonlinear coordinate families.
    p,q,r,w = S.symbols('p q r w')
    L,N,a = S.symbols('L N a',nonzero=True)
    M = S.symbols('M')
    J = S.Function('J')
    R0 = J(r)
    native_trace = L*p*p+M*p+4*N*p*R0*R0+N*q*R0
    zero('univariate old-r boundary gives exactly the required linear-coefficient trace',
         native_trace-(L*p*p+(M+4*N*J(r)**2)*p+N*J(r)*q))
    zero('every zero of the boundary coefficient retains the nonzero quadratic p coefficient',
         native_trace.subs(J(r),0)-(L*p*p+M*p))
    repeated_trace = native_trace.subs(J(r),r*r)
    zero('a repeated old-r boundary may have no isolated trace critical point',
         S.Matrix([S.diff(repeated_trace,z) for z in (p,q,r)]).subs({p:-M/(2*L),r:0},simultaneous=True))
    zero('the repeated-root trace nevertheless has the positive distinct-root Hodge budget',
         (tt-1)*tt+2*tt-(tt**2+tt))
    lower = S.Function('lower')(p,q)
    zcoord = a*r+lower
    zero('unit-affine old-r coordinate has its exact polynomial inverse',
         zcoord.subs(r,(w-lower)/a)-w)
    transformed = native_trace.subs(r,zcoord)
    zero('a nonlinear polynomial of the unit-affine r coordinate reduces by the actual trace inverse',
         transformed.subs(r,(w-lower)/a)-native_trace.subs(r,w))
    Ccoord,Dcoord = r+(p+r*r)**2,p+r*r
    alpha,T = S.symbols('alpha T')
    pb,rb = T-(alpha-T*T)**2,alpha-T*T
    zero('nonvertical nonlinear plane coordinate has unit Jacobian',
         S.Matrix([Ccoord,Dcoord]).jacobian((p,r)).det()+1)
    zero('nonvertical nonlinear plane coordinate has its actual inverse',
         S.Matrix([Ccoord,Dcoord]).subs({p:pb,r:rb},simultaneous=True)-S.Matrix([alpha,T]))
    check('nonvertical coordinate restriction has exact degree four on every level',S.Poly(pb,T).degree()==4)
    check('quadratic critical values have exact degree eight on those levels',S.Poly(L*pb*pb+M*pb,T).degree()==8)
    zero('two distinct nonlinear-coordinate levels have the full positive budget fourteen',
         (tt**2-2*tt)+2*8*tt-(tt**2+14*tt))
    zero('the excluded vertical plane coordinate has a nonpositive budget',
         tt**2+(2*0-m)*tt-(tt**2-m*tt))

    # Branched triples are not passed off as coordinate or Keller triples.
    zero('a nonlinear branched trace triple really has a nonunit Jacobian',
         S.Matrix([r*r,p,q]).jacobian((p,q,r)).det()-2*r)
    Hfun = S.Function('Hfun')(U)
    unit_control = s*U+Hfun+V*W
    zero('arbitrary univariate trace controls have unit full Hessian',S.hessian(unit_control,(s,U,V,W)).det()-1)

    dependencies = {}
    prior = [('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
             ('RELATIVE_CRITICAL_LOCUS.md','receipt-relative-critical.json'),
             ('ISOLATED_TRACE_MIXING.md','receipt-isolated-trace.json')]
    for name,receipt_name in prior:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    manifest = json.loads((HERE/'sources/hodge-trace-import.json').read_text())
    hodge_source = HERE/'sources'/manifest['retained_file']
    digest = hashlib.sha256(hodge_source.read_bytes()).hexdigest()
    check('basic Hodge input has the exact retained primary v1 bytes',digest==manifest['sha256'])
    check('Hodge input records the ring statement and precise projective curve formula',
          any('3.6' in value for value in manifest['locations'])
          and any('ring homomorphism' in value for value in manifest['locations']))
    dependencies[str(hodge_source.relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md','hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'HODGE_TRACE_CLOSURE.md','sources/hodge-trace-import.json']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['both Tang ranks from one entire Hodge trace value fibre',
            'curve genus and puncture coefficients force actual closed A1 components',
            'all traces H(U,V)+a(V)*W with positive distinct-root degree budget',
            'every univariate H(U)+VW, with arbitrary critical-value pattern',
            'all polynomial unit-affine-r boundaries and nonvertical plane-coordinate boundaries',
            'branched trace triples with arbitrary off-hyperplane mixing'],
        generic_base_Euler_hypothesis_required=False,
        isolated_trace_point_required=False,
        all_univariate_H_plus_VW_degrees_closed=True,
        nonlinear_univariate_polynomials_of_unit_affine_r_boundary_closed=True,
        nonvertical_plane_coordinate_boundaries_closed=True,
        branched_trace_triples_with_positive_budget_closed=True,
        arbitrary_finite_cover_Hodge_multiplicativity_used=False,
        zero_or_negative_Hodge_coefficient_asserted_closed=False,
        broader_substitutions_asserted_embedded=False,
        unit_Hessian_nonexistence_asserted=False,
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
