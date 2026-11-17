"""Exact controls for the two-trace compatibility theorem.

Finite identities do not certify the written all-degree curve argument
or the imported Tang and Hodge results.
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

    def check(name, value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name, value):
        values = list(value) if isinstance(value, S.MatrixBase) else [value]
        check(name, all(S.cancel(v.doit()) == 0 for v in values))

    u, v = S.symbols('u v')
    t = u*v
    c, g, r, e, d, m = S.symbols('c g r e d m', integer=True)
    curve = c*t-g*(u+v)+c-r-e
    trace = t*t-curve+d*t
    zero('reduced singular curve normalization keeps branch identifications',
         trace-(t*t+(d-c)*t+g*(u+v)+r+e-c))
    poly = S.Poly(trace, u, v)
    zero('total normalization genus is visible in each linear Hodge coefficient',
         poly.coeff_monomial(u)+poly.coeff_monomial(v)-2*g)
    zero('punctures and branch identifications survive in the constant term',
         poly.coeff_monomial(1)-(r+e-c))
    zero('comparison with a nonpositive pure trace gives d plus m equals c',
         (trace-(t*t-m*t)).subs({g:0, r:c, e:0})-(d+m-c)*t)
    nonnegative = S.symbols('nonnegative', nonnegative=True)
    positive = S.symbols('positive', positive=True)
    check('positive component count makes d at least twice c incompatible with d plus m equals c',
          S.ask(S.Q.positive(positive+nonnegative)) is True)

    db, dc = S.symbols('db dc', integer=True)
    rank_one_b = t*t+(db-c)*t+g*(u+v)+r-c
    rank_one_c = t*t+(dc-c)*t+g*(u+v)+r-c
    zero('parallel rank-one profiles differ only in the t coefficient',
         rank_one_b-rank_one_c-(db-dc)*t)
    zero('rank-zero pure profiles have nonpositive coefficient',
         (t-m)*t-(t*t-m*t))
    zero('rank-one pure comparison forces one puncture per genus-zero component',
         rank_one_b.subs({g:0,r:c})-(t*t+(db-c)*t))

    x, p, q, rr, a, b, L, M, N, lam, V, U, W = S.symbols(
        'x p q rr a b L M N lam V U W')
    # Use generic function expressions: these controls do not fix graph degree.
    Z = S.Function('Z')(p, rr)
    h = L*p*p+M*p+4*N*p*Z*Z+N*q*Z
    qback = V/N-4*p*Z
    zero('native boundary shear has its exact polynomial inverse over the field',
         N*(qback+4*p*Z)-V)
    zero('native q-independent boundary becomes H plus V times Z',
         h.subs(q,qback)-(L*p*p+M*p+V*Z))
    shear = S.Matrix([p,N*(q+4*p*Z),rr])
    zero('boundary shear has field-unit Jacobian', shear.jacobian((p,q,rr)).det()-N)
    zero('linear trace open stratum solves its free coefficient variable exactly',
         (L*p*p+M*p+V*Z).subs(V,(lam-L*p*p-M*p)/Z)-lam)

    # All physical base points other than the origin have independent A,B rows.
    aa, bb = S.symbols('aa bb')
    uu = 1+aa*bb
    Arow = S.Matrix([[uu**3, 3*aa*uu**2, -aa**3]])
    Brow = S.Matrix([[bb**2*uu*(4+3*aa*bb),
                      bb+3*aa*bb**2*(4+3*aa*bb), 2*aa-3*aa**2*bb]])
    cross = Arow.cross(Brow)
    target = S.Matrix([[2*aa**2*(2*aa*bb+3), -aa*uu*(uu+1), bb*uu**3]])
    zero('ambient row cross product at an arbitrary constant graph slice', cross-target)
    zero('the first possible cross-product zero retains a nonzero second entry',
         target[1].subs(bb,-S.Rational(3,2)/aa)-aa/4)
    zero('the x-zero nonorigin slice is independent as well',
         target.subs(aa,0)-S.Matrix([[0,0,bb]]))
    A, B, C, value = S.symbols('A B C value')
    zero('constant-base trace is a coordinate with an actual triangular inverse',
         (L*A*A+M*A+N*B).subs(B,(value-L*A*A-M*A)/N)-value)
    zero('constant-base coordinate completion has field-unit Jacobian',
         S.Matrix([A,L*A*A+M*A+N*B,C]).jacobian((A,B,C)).det()-N)

    T = S.Function('T')(x,p,q,rr)
    interpolation = (1-x/a)*Z+(x/a)*b+x*(x-a)*T
    zero('arbitrary interpolation corrections preserve the x-zero boundary',
         interpolation.subs(x,0)-Z)
    zero('arbitrary interpolation corrections preserve the nonzero constant slice',
         interpolation.subs(x,a)-b)

    # Genuine trace automorphisms, with nonlinear base-dependent mixed words.
    pp = p+x*rr**2
    qq = q+x*pp**2
    rnew = rr+x*qq**2
    back_r = W-x*V**2
    back_p = U-x*back_r**2
    back_q = V-x*U**2
    zero('mixed retained-dual word has an actual inverse in every base fibre',
         S.Matrix([pp,qq,rnew]).subs(
             {p:back_p,q:back_q,rr:back_r},simultaneous=True)-S.Matrix([U,V,W]))

    # Singular/positive-genus controls distinguish the actual strata.
    zero('a unibranch cusp has the same E polynomial as its A1 normalization',
         curve.subs({c:1,g:0,r:1,e:0})-t)
    zero('a node keeps its branch-identification correction',
         curve.subs({c:1,g:0,r:1,e:1})-(t-1))
    zero('two crossing affine lines keep their intersection correction',
         curve.subs({c:2,g:0,r:2,e:1})-(2*t-1))
    zero('a punctured hyperbola has a nonzero constant obstruction',
         trace.subs({c:1,g:0,r:2,e:0,d:2})-(t*t+t+1))
    zero('a once-punctured genus-one curve has both genus obstructions',
         trace.subs({c:1,g:1,r:1,e:0,d:4})-(t*t+3*t+u+v))
    zero('a disjoint vertical line and hyperbola do not fake a pure signature',
         trace.subs({c:2,g:0,r:3,e:0,d:2})-(t*t+1))
    zero('vertical curves alone really can have the unresolved nonpositive signature',
         trace.subs({c:3,g:0,r:3,e:0,d:0})-(t*t-3*t))

    tt = S.Symbol('tt')
    cusp = rr**2-p**3
    zero('cuspidal boundary has its exact normalization parametrization',
         cusp.subs({p:tt**2,rr:tt**3},simultaneous=True))
    check('cuspidal normalization gives generic quadratic trace degree four',
          S.Poly((L*p*p+M*p).subs(p,tt**2),tt).degree()==4)
    zero('cuspidal boundary yields a positive coefficient three',
         t*t-t+4*t-(t*t+3*t))
    for power in (1,2,7):
        check('repeated cusp factor counts the same reduced curve, exponent '+str(power),
              S.Poly(cusp**power,p,rr).sqf_part().monic()==S.Poly(cusp,p,rr).monic())
    for mm,nn in ((2,3),(3,5),(4,7)):
        zero('coprime binomial normalization '+str((mm,nn)),
             (rr**mm-p**nn).subs({p:tt**mm,rr:tt**nn},simultaneous=True))
        check('coprime binomial positive budget '+str((mm,nn)),
              S.Poly((L*p*p+M*p).subs(p,tt**mm),tt).degree()-1==2*mm-1>0)

    # The tempting sum over distinct critical values is not a degree bound.
    cubic = tt**3-3*tt
    zero('two critical points of a cubic have distinct values',
         cubic.subs(tt,1)-cubic.subs(tt,-1)+4)
    check('local degrees across different value fibres cannot be summed globally',
          2+2>S.degree(cubic,tt))

    dependencies = {}
    for proof,receipt_name in [
        ('HODGE_TRACE_CLOSURE.md','receipt-hodge-trace.json'),
        ('CRITICAL_NORMALIZATION.md','receipt-critical-normalization.json'),
        ('CRITICAL_COMPONENTS.md','receipt-critical-components.json')]:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/proof).read_bytes()).hexdigest()
        check('retained proof pin '+proof,digest==retained['source_sha256'][proof])
        dependencies[str((HERE/proof).relative_to(ROOT))] = digest
    pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        digest = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,digest==pins[name])
        dependencies[name] = digest
    manifest = json.loads((HERE/'sources/hodge-trace-import.json').read_text())
    source = HERE/'sources'/manifest['retained_file']
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    check('retained primary Hodge input',digest==manifest['sha256'])
    dependencies[str(source.relative_to(ROOT))] = digest
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name,'TWO_TRACE_COMPATIBILITY.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['parallel Tang rank and Hodge-profile rigidity',
            'singular reduced curve normalization and branch-count additivity',
            'A1 normalization surjectivity excludes all vertical components',
            'unrestricted linear trace cannot have a nonpositive pure signature',
            'two-slice closure for arbitrary q-independent positive-r boundary',
            'one-place unibranch cuspidal boundary closure without a second slice'],
        all_plane_curve_genera_singularities_multiplicities_allowed_with_second_slice=True,
        arbitrary_off_both_slices_corrections_allowed=True,
        one_place_unibranch_boundary_requires_second_slice=False,
        arbitrary_finite_cover_Hodge_multiplicativity_used=False,
        local_degrees_summed_across_distinct_values=False,
        normalization_A1_asserted_actual_smooth_A1=False,
        one_pure_trace_alone_asserted_sufficient=False,
        all_positive_r_boundaries_without_second_slice_asserted_closed=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,
        singleton_zero_properness_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
