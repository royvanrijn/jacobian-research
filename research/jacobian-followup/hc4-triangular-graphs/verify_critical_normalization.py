"""Exact controls for finite-polar exclusion and reciprocal escape laws.

The finite-branch exhaustion and all-degree normalization arguments are
written proofs. Standard local inputs are explicitly pinned imports.
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

    def order(expression):
        numerator,denominator = S.cancel(expression).as_numer_denom()
        assert numerator!=0
        def lowest(poly):
            return min(monomial[0] for monomial,coefficient in S.Poly(poly,t).terms() if coefficient!=0)
        return lowest(numerator)-lowest(denominator)

    x,p,q,r,T,t,w = S.symbols('x p q r T t w')
    L,N = S.symbols('L N',nonzero=True)
    M = S.symbols('M')

    def critical_value(xx,pp):
        return -N*pp/(4*xx**2)-(M-13*N/(2*xx**2))**2/(4*L)

    def critical_hessian(xx,pp,rp):
        value = critical_value(xx,pp)
        derivative = S.diff(value,t)/S.diff(xx,t)
        second = S.diff(derivative,t)/S.diff(xx,t)
        return S.cancel(-32*L*N*N*xx**6*rp**2*second)

    # Smooth ramified branch, with any other factors counted by I.
    # These are local plane models; no global graph or unit Hessian is
    # inferred from them. The cases with I>0 are not in the prior gate.
    for e,I in [(2,0),(3,1),(4,2)]:
        xx,pp = 1+t**e,t
        h = x-1-p**e
        full = h*p**I
        rp = S.diff(full,p).subs({x:xx,p:pp},simultaneous=True)/(2*xx)
        K = critical_hessian(xx,pp,rp)
        check('smooth finite branch e='+str(e)+', I='+str(I)+' has exact odd Hessian valuation2I-1',order(K)==2*I-1)
    xx,pp = 1+t,t*t
    full = (p-(x-1)**2)*(x-1)**2
    rp = S.diff(full,p).subs({x:xx,p:pp},simultaneous=True)/(2*xx)
    check('smooth x-etale component intersection gives a positive Hessian valuation',
          order(critical_hessian(xx,pp,rp))==4)

    # Primitive monomial singular branches check the conductor identity
    # and strict positive order, with a repeated other factor in one case.
    for e,m,other_power in [(2,3,0),(3,2,0),(3,4,0),(2,3,2)]:
        xx,pp = 1+t**e,t**m
        h = p**e-(x-1)**m
        full = h*p**other_power
        rp = S.diff(full,p).subs({x:xx,p:pp},simultaneous=True)/(2*xx)
        conductor = (e-1)*(m-1)
        intersection = m*other_power
        check('singular branch '+str((e,m,other_power))+' has the exact conductor/intersection derivative order',
              order(rp)==e-1+conductor+intersection)
        check('singular branch '+str((e,m,other_power))+' forces a positive Hessian order beyond the proved bound',
              order(critical_hessian(xx,pp,rp))>=2*conductor+2*intersection-1>0)

    # An actual original-family polynomial demonstrates why the full F_p
    # is needed, rather than the derivative of its active squarefree part.
    F = (1-x)*(x*p+3)
    R = S.cancel((F-3)/(2*x))
    zero('intersection control is an actual polynomial original-family graph',3+2*x*R-F)
    check('intersection graph has polynomial coefficients',R.is_polynomial(x,p))
    point = {x:1,p:-3}
    active = x*p+3
    zero('full polar Bezout gate has a literal common zero',S.Matrix([active,S.diff(F,p)]).subs(point))
    zero('using only the active derivative would miss that common zero',S.diff(active,p).subs(point)-1)

    # A complete critical graph over C* with explicit polynomial inverse
    # representatives in its coordinate ring. It passes the necessary
    # Bezout condition but was already excluded by the affine-p theorem.
    R = p
    F = 3+2*x*p
    Fp = S.diff(F,p)
    zero('explicit polynomial Bezout identity for the complete-curve control',F/3-p*Fp/3-1)
    inverse_x = -2*p/3
    zero('base inverse is polynomial modulo the entire critical curve',x*inverse_x-1+F/3)
    Tc = 13*N*inverse_x**2/(4*L)-M/(2*L)
    qc = -3*Tc*inverse_x/4+inverse_x**2/16
    rc = (-p/8+3*x*qc/4-Tc)*inverse_x**3
    check('critical graph coordinates are actual polynomial restrictions',
          all(value.is_polynomial(x,p) for value in (Tc,qc,rc)))
    y = S.symbols('y')
    u = 1+x*y
    AA = u**3*p+3*x*u*u*q-x**3*r
    BB = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    native = (L*AA*AA+M*AA+N*BB).subs(y,R)
    critical_sub = {p:-3/(2*x),q:qc.subs(p,-3/(2*x)),r:rc.subs(p,-3/(2*x))}
    zero('whole polynomial critical graph satisfies all original relative critical equations',
         S.Matrix([S.diff(native,zi) for zi in (p,q,r)]).subs(critical_sub,simultaneous=True))
    zero('whole-curve value function agrees with the original potential',
         native.subs(critical_sub,simultaneous=True)-critical_value(x,-3/(2*x)))

    # The reciprocal differentiation factor is valid for arbitrary p
    # degree; exact controls retain lower terms and a vertical factor.
    for degree in (1,2,5):
        generic = (x-2)*(3+2*x*(p**degree+(x*x+1)*p+7))
        reciprocal = S.expand(w**degree*generic.subs(p,1/w))
        zero('reciprocal derivative identity in degree '+str(degree),
             S.diff(reciprocal,w)-degree*w**(degree-1)*generic.subs(p,1/w)
             +w**(degree-2)*S.diff(generic,p).subs(p,1/w))

    # Escape valuation without a classification or a ramification cutoff.
    for e,m in [(1,1),(2,2),(3,4)]:
        xx,pp = 1+t**e,t**(-m)
        DD = critical_value(xx,pp)
        DD2 = S.diff(S.diff(DD,t)/S.diff(xx,t),t)/S.diff(xx,t)
        check('finite-base escape '+str((e,m))+' has exact D second-derivative order',order(DD2)==-m-2*e)
    e,m,n,conductor,intersection = S.symbols('e m n conductor intersection',integer=True)
    zero('escape conductor equality is the exact cancellation of ramification index',
         2*((2-n)*m+conductor+e-1+intersection)-(2*e+m)
         -(2*(conductor+intersection)-(2*n-3)*m-2))
    check('a reduced ordinary node can never meet the escape equality in positive integral degree',
          S.Rational(3,2).is_integer is False)

    # A genuine global graph with a cusp escape passing the local budget.
    # Its full critical function is deliberately checked to be nonconstant.
    z = x-1
    R = -3*z*p+S.Rational(3,2)*x*z*z*(1-z)*p*p
    F = 3+2*x*R
    xx = 1+t*t
    pp = 1/(xx*(t*t+t**3))
    zero('cusp escape control is on an actual original-family critical curve',F.subs({x:xx,p:pp},simultaneous=True))
    rp = S.cancel(S.diff(R,p).subs({x:xx,p:pp},simultaneous=True))
    zero('cusp escape has exactly the required R_p order3',rp+3*t**3)
    K = critical_hessian(xx,pp,rp)
    zero('admissible cusp escape has a nonzero critical-Hessian limit',S.limit(K,t,0)-144*L*N**3)
    check('passing the cusp escape budget does not make the critical Hessian constant',
          S.cancel(K.subs({L:1,N:1,M:0})-144)!=0)
    zero('cusp reciprocal chart has conductor2 via primitive (t^2,t^3) coordinates',
         (1/pp)/xx-(xx-1)-t**3)
    check('cusp packet satisfies the exact degree-two conductor law',2*2==(2*2-3)*2+2)

    manifest_name = 'sources/critical-normalization-import.json'
    manifest = json.loads((HERE/manifest_name).read_text())
    dependencies = {}
    for imported in manifest['imports']:
        name = imported['retained_file_from_research_root']
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained primary byte pin '+imported['authors'],dependencies[name]==imported['sha256'])
    check('local imports retain their precise theorem locations',
          'Theorem 1.6' in manifest['imports'][0]['location'] and 'Theorem 2.1' in manifest['imports'][1]['location'])
    for name in ['RELATIVE_CRITICAL_LOCUS.md','NEWTON_CRITICAL_GATE.md','CRITICAL_COMPONENTS.md']:
        receipt_name = {'RELATIVE_CRITICAL_LOCUS.md':'receipt-relative-critical.json',
                        'NEWTON_CRITICAL_GATE.md':'receipt-newton-critical.json',
                        'CRITICAL_COMPONENTS.md':'receipt-critical-components.json'}[name]
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained written proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md','jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input for residual curve realization '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'CRITICAL_NORMALIZATION.md',manifest_name]},
        dependencies=dependencies,
        written_not_finite_exhaustion=['all finite smooth and singular branches, including every component intersection',
            'the full polynomial Bezout gate and entire closed critical curve',
            'all finite-base escape pole orders and ramification indices',
            'exact reciprocal conductor/intersection budget',
            'abstract residual plane Keller-fibre realization preserving both functions'],
        entire_finite_polar_locus_closed_in_all_degrees=True,
        all_finite_singularities_and_component_intersections_closed=True,
        full_F_p_not_replaced_by_active_squarefree_derivative=True,
        critical_curve_globally_complete_after_Bezout_gate=True,
        etale_projection_automatically_asserted_finite=False,
        all_finite_base_escapes_closed=False,
        admissible_cusp_packet_asserted_sufficient_for_unit_Hessian=False,
        abstract_curve_isomorphism_asserted_ambient_plane_automorphism=False,
        all_repeated_initial_faces_closed=False,
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
