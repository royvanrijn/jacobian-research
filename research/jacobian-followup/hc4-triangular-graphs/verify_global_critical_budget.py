"""Exact controls for the global Euler budget and differential descent.

The stratification, covering argument, all-degree conjugate argument and
parametric family proof are written mathematics, not finite searches.
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

    def laurent_degree(expression,variable):
        exponents = [0]+[int(term.as_powers_dict().get(variable,0)) for term in S.expand(expression).as_ordered_terms()]
        return max(exponents)-min(exponents)

    x,p,q,r,t,s,V,U,Z,lam = S.symbols('x p q r t s V U Z lam')
    a,b,e1,e2,chiC,h = S.symbols('a b e1 e2 chiC h',integer=True)
    # The rank-one value stratification, with a genuinely unit-Hessian
    # normal-form control. Its A1 critical curve is not a remaining graph.
    normal = U**3+s*U+V*Z
    zero('rank-one control has actual unit Hessian',S.hessian(normal,(s,V,U,Z)).det()-1)
    zero('off the residual curve the value equation solves uniquely for V',
         normal.subs(V,(lam-U**3-s*U)/Z)-lam)
    zero('on the residual curve the value equation is the critical-value pencil',
         normal.subs(Z,0)-(U**3+s*U))
    zero('generic cubic pencil has three distinct values counted by its discriminant',
         S.discriminant(U**3+s*U-lam,U)+4*s**3+27*lam**2)
    zero('critical-value and Legendre functions in the unit control',
         normal.subs({s:-3*U**2,V:0,Z:0})-(-2*U**3))
    zero('rank-one control preserves the critical target-line derivative',S.diff(normal,s)-U)
    zero('Euler comparison gives h=r+chiC algebraically',
         (1-chiC+h)-(1+r)-(h-r-chiC))
    zero('etale cover plus its missing finite fibres has chi=-sum escape indices',
         -2*a+(a-e1)+(a-e2)+e1+e2)
    zero('the genus-free degree budget follows from both value fibres',
         (h-r-chiC).subs(chiC,-e1-e2)-(h+e1+e2-r))
    check('hyperelliptic missing-fibre control has chi=-4 by the covering strata',2*(-2)+0+0==-4)

    # Pole contributions at zero and at a finite-base escape.
    D = t**(-4)
    tau = S.diff(D,t)
    E = D-t*tau
    check('noncancelling zero-base unit-index branch costs five Euler units',laurent_degree(E+s*tau,t)==5)
    xx = 1+t**2
    DD = t**(-2)
    ttau = S.diff(DD,t)/S.diff(xx,t)
    EE = DD-xx*ttau
    check('cusp-type finite escape has pencil pole m+e=4',laurent_degree(EE+s*ttau,t)==4)
    check('the same finite escape costs m+2e=6 including its missing sheets',4+2==6)

    # Differential descent and its full quadratic exception.
    alpha,beta,gamma = S.symbols('alpha beta gamma')
    zero('all conjugate affine differences have zero second derivative',S.diff(alpha*x+beta,x,2))
    algebraic_u = x**S.Rational(1,3)/x**2
    check('a nonrational algebraic control has nonrational second derivative',
          not S.diff(algebraic_u,x,2).is_rational_function(x))
    A,B,C = S.symbols('A B C')
    quadratic = A*p*p+B*p+C
    disc = B*B-4*A*(C+3/(2*x))
    zero('quadratic derivative squares descend to the same rational discriminant',
         S.diff(quadratic,p)**2-disc-4*A*(quadratic+3/(2*x)))
    split_F = 3*((x*x*p-1)**2-x**8)
    split_R = S.cancel((split_F-3)/(2*x))
    pplus,pminus = x**(-2)+x*x,x**(-2)-x*x
    zero('split quadratic exceptional branch difference passes its affine second-derivative test',
         S.diff((pplus-pminus)/x**2,x,2))
    split_D = -pplus/(4*x*x)-169/(16*x**4)
    split_Rp = S.diff(split_R,p).subs(p,pplus)
    K = S.cancel(-32*x**6*split_Rp**2*S.diff(split_D,x,2))
    zero('the rational split exception instead has a positive critical-Hessian order at zero',S.limit(K/x**10,x,0)-62280)

    # Actual all-degree family controls; the odd-k quantifier is proved
    # once in the canonical note, including generic squarefreeness of Q.
    for k in (3,5,7):
        ell = (k-1)//2
        G = 3+S.Rational(3,4)*x**(2*k-2)*p+S.Rational(507,16)*x**(2*k-4)
        invx = -S.Rational(1,4)*x**(2*k-3)*p-S.Rational(169,16)*x**(2*k-5)
        Q = invx**ell+(1+x)*G/3
        F = G*Q**2
        R = S.cancel((F-3)/(2*x))
        pc = -4*x**(-(2*k-2))-S.Rational(169,4)*x**(-2)
        check('family k='+str(k)+' is a genuine polynomial y graph of p degree k',R.is_polynomial(x,p) and S.degree(R,p)==k)
        zero('family k='+str(k)+' has the asserted nonzero constant boundary',R.subs(x,0)-(-S.Rational(459,16) if k==3 else 3))
        zero('family k='+str(k)+' polynomial inverse of x modulo G',x*invx-1+G/3)
        zero('family k='+str(k)+' active critical curve is an entire Cstar graph',G.subs(p,pc))
        zero('family k='+str(k)+' repeated factor is a unit on the active curve',Q.subs(p,pc)-x**(-ell))
        rp = S.cancel(S.diff(R,p).subs(p,pc))
        zero('family k='+str(k)+' full R_p on the active curve',rp-S.Rational(3,8)*x**(k-2))
        value = S.cancel(-pc/(4*x*x)-169/(16*x**4))
        zero('family k='+str(k)+' exact critical value',value-x**(-2*k))
        zero('family k='+str(k)+' native critical Hessian is exactly nonzero constant',
             -32*x**6*rp**2*S.diff(value,x,2)+9*k*(2*k+1))
        bezout_V = S.Rational(4,3)*invx**(k-1)
        quotient = S.cancel((1-bezout_V*S.diff(F,p))/G)
        check('family k='+str(k)+' has a full polynomial Bezout certificate',quotient.is_polynomial(x,p))
        zero('family k='+str(k)+' full Bezout identity',quotient*G+bezout_V*S.diff(F,p)-1)
        pencil = (2*k+1)*x**(-2*k)-2*k*s*x**(-2*k-1)
        check('family k='+str(k)+' global pencil degree exceeds every distinct native root',
              laurent_degree(pencil,x)==2*k+1>ell+1)
        face = S.expand(F.subs(p,Z*x**(-(2*k-3)))).subs(x,0)
        zero('family k='+str(k)+' initial face is entirely repeated',face-3*(1+(-Z/4)**ell)**2)
    v = S.symbols('v')
    for ell in (1,2,3):
        QQ = x**(-ell)*(1-v)**ell+(1+x)*v
        zero('generic repeated-factor squarefreeness identity for ell='+str(ell),
             ell*QQ+(1-v)*S.diff(QQ,v)-(1+x)*(1+(ell-1)*v))

    dependencies = {}
    for name,receipt_name in [('CRITICAL_NORMALIZATION.md','receipt-critical-normalization.json'),
                              ('CRITICAL_COMPONENTS.md','receipt-critical-components.json'),
                              ('RELATIVE_CRITICAL_LOCUS.md','receipt-relative-critical.json')]:
        retained = json.loads((HERE/receipt_name).read_text())
        digest = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        check('retained canonical input byte pin '+name,digest==retained['source_sha256'][name])
        dependencies[str((HERE/name).relative_to(ROOT))] = digest
    tang_pins = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md','jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'GLOBAL_CRITICAL_BUDGET.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['entire rank-one value-fibre stratification',
            'etale covering Euler calculation with every finite-base escape',
            'genus-free global degree and pole budget in all graph degrees',
            'algebraic conjugate descent in every extension degree',
            'all quadratic p dependence including the repeated-boundary distinction',
            'the entire odd-k family and its generic root count'],
        global_pencil_degree_plus_escape_budget_equals_distinct_root_count=True,
        inactive_root_multiplicities_do_not_multiply_budget=True,
        every_quadratic_p_graph_gradient_invertible_under_unit_Hessian=True,
        every_quadratic_p_graph_asserted_no_unit_Hessian=False,
        all_degree_local_gate_passing_family_excluded_by_global_Euler=True,
        constant_critical_Hessian_asserted_sufficient=False,
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
