"""Exact controls for the full positive-p y-graph closure.

Every-value Euler comparison, normalization root products and all-degree
exhaustion are written proofs. These controls are not a formalization.
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

    x,p,q,r,t,s,Y,lam,zeta = S.symbols('x p q r t s Y lam zeta')
    L,N = S.symbols('L N',nonzero=True)
    M,alpha,beta,alpha2,beta2 = S.symbols('M alpha beta alpha2 beta2')
    D = -N*p/(4*x*x)-(M-13*N/(2*x*x))**2/(4*L)
    pstar = -169*N/(4*L*x*x)+13*M/L-M*M*x*x/(L*N)
    zero('the full critical value has the exact universal-shift form',D+N*(p-pstar)/(4*x*x))
    p1,p2 = S.symbols('p1 p2')
    zero('distinct native p roots have distinct exceptional values at fixed nonzero base',
         D.subs(p,p1)-D.subs(p,p2)+N*(p1-p2)/(4*x*x))
    pinactive = pstar-4*x*x*(alpha+beta*x)/N
    zero('an affine exceptional value has exactly the stated rational inactive root',D.subs(p,pinactive)-alpha-beta*x)
    zero('every such inactive root has the same nonzero quadratic pole',S.limit(x*x*pinactive,x,0)+169*N/(4*L))

    # Exact every-value defects for a curve pencil, including one missing
    # boundary point. This is an abstract C* control, not a unit graph.
    Hs = t*t+s*t
    zero('abstract jump control has the generic two-sheet fibre discriminant',S.discriminant(Hs-lam,t)-s*s-4*lam)
    zero('its only affine critical point is simple',S.diff(Hs,t).subs(t,-s/2))
    zero('its active critical value and exact double-root fibre',Hs+s*s/4-(t+s/2)**2)
    zero('its boundary-defect value has one excluded and one affine root',Hs-t*(t+s))
    zero('its boundary local degree is exactly one generically',S.diff(Hs,t).subs(t,0)-s)
    zero('critical Legendre identity on the abstract pencil',S.diff(t*t,t)+(-2*t)*S.diff(t,t))
    check('a repeated boundary defect is distinguishable from a distinct native inactive value',S.discriminant((Y-alpha-beta*x)**2,Y)==0)
    check('two distinct affine exceptional values form a squarefree polynomial generically',
          S.discriminant((Y-alpha-beta*x)*(Y-alpha2-beta2*x),Y)!=0)

    # Boundary valuation bookkeeping in the every-value comparison.
    e,j = S.symbols('e j',integer=True)
    zero('a simple missing infinity value gives the odd critical-Hessian order',
         -6*e+2*j+(2*e+1)-(2*j-4*e+1))
    check('that infinity order cannot be zero for integers e,j',
          S.Poly(2*j-4*e+1,j,e,modulus=2).as_expr()==1)
    zero('the zero-base simple missing-value order',6*e+2*j+(1-e)-(5*e+2*j+1))
    # Both regular critical values give reciprocal-root differences of
    # order at least six in x, including all universal base terms.
    pc1 = pstar-4*x*x*(alpha+beta*x)/N
    pc2 = pstar-4*x*x*(alpha2+beta2*x)/N
    zero('exact reciprocal-root difference identity',
         1/pc1-1/pc2-4*x*x*((alpha+beta*x)-(alpha2+beta2*x))/(N*pc1*pc2))
    c = -169*N/(4*L)
    zero('regular-value reciprocal roots have the exact leading order-six coefficient',
         S.limit((1/pc1-1/pc2)/x**6,x,0)-4*(alpha-alpha2)/(N*c*c))
    # A primitive cubic normalization control with d=2; constant and
    # linear-in-x terms cancel between its conjugates.
    xx = t**3
    Dc1 = alpha+beta*xx+t*t
    Dc2 = alpha+beta*xx+zeta**2*t*t
    pc1 = pstar.subs(x,xx)-4*xx*xx*Dc1/N
    pc2 = pstar.subs(x,xx)-4*xx*xx*Dc2/N
    zero('self-conjugate reciprocal-root contact includes the extra d=2 order',
         S.limit((1/pc1-1/pc2)/t**20,t,0)-4*(1-zeta**2)/(N*c*c))

    # Universal root-product and linear contradiction, with no cutoffs.
    a,A,ell,P,mu,v,b,d = S.symbols('a A ell P mu v b d')
    n = a+v
    upper = (2*n-5)*e-d/2
    lower = (e-1)*(6*e+d)+6*e*(A-e)+6*e*v+2*e*(a-A-ell)+e*mu
    desired_slack = 2*ell-mu-4*(A+v)+1-(1-1/(2*e))*d
    zero('the contact table counts all n-1 other roots with inactive multiplicity',
         (e-1)+(A-e)+v+(a-A-ell)+ell-(n-1))
    zero('the root-product inequality rearranges exactly to the stated all-degree bound',
         (upper-lower)/e-desired_slack)
    zero('the even-d improvement is a sum of nonnegative terms',
         -1+(1-1/(2*e))*d-((d-2)*(1-1/(2*e))+1-1/e))
    zero('the global pole budget gives the stated sheet-count slack',
         a+b-5*ell-(a-A-ell)-(A+b-4*ell))
    zero('the final linear contradiction accounts for every inactive multiplicity',
         (2*mu+7*A+8*v-b)-(2*mu+7*A+15*b+8*(v-2*b)))
    zero('the exact reciprocal derivative order at the zero-base regular-D branch',
         (3-2*n)*e+upper-(-2*e-d/2))
    # The first non-affine power t^d determines D'', while the constant
    # and x-linear terms disappear, at arbitrary formal integer e,d.
    de,dd = S.symbols('de dd',integer=True,positive=True)
    coefficient = dd*(dd-de)/de**2
    zero('formal second-derivative leading coefficient',
         (dd/de)*((dd-de)/de)-coefficient)
    check('odd leading exponent cannot yield an integral critical derivative order',
          S.Poly(4*e+2*j+1,j,e,modulus=2).as_expr()==1)

    # The repeated simple-boundary example is retained: the closure is
    # gradient invertibility, not universal nonexistence of unit Hessian.
    repeated_R = 3*p+S.Rational(3,2)*x*p*p
    zero('repeated whole curve has no active component',3+2*x*repeated_R-3*(1+x*p)**2)
    zero('the repeated example has an actual simple boundary root',repeated_R.subs(x,0)-3*p)
    zero('relative coordinate control has no relative critical points',S.diff(p*p+q,q)-1)

    dependencies = {}
    prior = [('GLOBAL_CRITICAL_BUDGET.md','receipt-global-critical-budget.json'),
             ('CRITICAL_NORMALIZATION.md','receipt-critical-normalization.json'),
             ('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
             ('RELATIVE_COORDINATE_ENDPOINT.md','receipt-relative-coordinate.json'),
             ('RELATIVE_CRITICAL_LOCUS.md','receipt-relative-critical.json')]
    for name,receipt_name in prior:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained written proof byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md','jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang source '+name,dependencies[name]==tang_pins[name])
    local_manifest = json.loads((HERE/'sources/critical-normalization-import.json').read_text())
    for item in local_manifest['imports']:
        name = item['retained_file_from_research_root']
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained local normalization source '+item['authors'],dependencies[name]==item['sha256'])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'Y_GRAPH_CLOSURE.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['every-value Euler jump and inactive factor polynomial identity',
            'all simple missing-value boundary points occur over zero',
            'all primitive Puiseux-root contacts including repeated inactive multiplicity',
            'the universal root-product/global-budget linear contradiction',
            'all positive p degrees in the y graph, every retained-dual automorphism',
            'the entire chart has only the existing base-only nonautomorphic endpoint'],
        all_positive_p_y_graphs_gradient_invertible_under_unit_Hessian=True,
        no_simple_root_remaining_y_boundary_forbids_unit_Hessian=True,
        active_positive_p_factor_forbids_unit_Hessian=True,
        whole_y_chart_nonautomorphic_endpoint_is_unit_ideal_plane_Keller=True,
        every_positive_p_graph_asserted_no_unit_Hessian=False,
        repeated_simple_boundary_case_asserted_unit_Hessian_sufficient=False,
        base_only_unit_Hessian_existence_asserted_resolved=False,
        general_nonlinear_opposite_x_graphs_closed=False,
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
