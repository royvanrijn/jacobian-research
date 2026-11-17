"""Exact algebra for the all-degree Meng--Yang graph obstruction.

The universal polynomial-degree and constant-plane arguments, and the
all-order binomial coefficient formula, are written in PROOF.md.
The finite series below is a regression, not the nontermination proof.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time
import flint
import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
checks = []


def check(name, condition):
    assert condition, name
    checks.append(name)


def zero(name, expression):
    if isinstance(expression, S.MatrixBase):
        check(name, all(S.cancel(e) == 0 for e in expression))
    else:
        check(name, S.cancel(expression) == 0)


def coefficient_polynomials(expression, variables, ctx):
    polynomial = S.Poly(expression, *variables)
    result = {}
    for exponent, coefficient in polynomial.terms():
        degree, tail = exponent[0], exponent[1:]
        result.setdefault(degree, {})[tail] = flint.fmpq(str(coefficient))
    return {degree: ctx.from_dict(terms) for degree, terms in result.items()}


def normal_series(C, D, variables, delta, c_formula):
    """Actual symbolic series at L=1,M=13,N=2,kappa=1; y,p,q remain free."""
    x,y,p,q = variables
    ctx = flint.fmpq_mpoly_ctx.get(['y','p','q'], 'lex')
    cc = coefficient_polynomials(C, variables, ctx)
    dd = coefficient_polynomials(D, variables, ctx)
    z = ctx.constant(0)
    check('implicit-series constant coefficient is a unit', cc[0] == 4)
    result = []
    for n in range(9):
        residual = dd.get(n, z)
        for i, coefficient in cc.items():
            if 1 <= i <= n:
                residual += coefficient*result[n-i]
        if n >= 5:
            residual += sum((result[i]*result[n-5-i] for i in range(n-4)), z)
        result.append(-residual*flint.fmpq(1,4))
    jet_degrees = [int(max(map(sum,value.to_dict()))) for value in result[:2]]
    check('formal example escapes the old first-jet degree-four gate',jet_degrees==[4,5])
    # Replay the relation as a whole, including the next generated coefficient.
    for n in range(9):
        residual = dd.get(n,z)
        residual += sum((coefficient*result[n-i] for i,coefficient in cc.items() if i<=n),z)
        if n>=5:
            residual += sum((result[i]*result[n-5-i] for i in range(n-4)),z)
        check(f'full y,p,q formal relation through x order {n}',residual == 0)
    for j in range(2,5):
        order = 3*j-5
        expected = S.cancel(c_formula(j).subs({S.Symbol('L'):1,S.Symbol('N'):2}))
        coefficient = result[order].to_dict().get((0,0,j),flint.fmpq(0))
        check(f'first q^{j} coefficient at x order {order}', coefficient == flint.fmpq(str(expected)))
        check(f'no earlier q^{j} term',all(not any(e[2]>=j for e in value.to_dict()) for value in result[:order]))
    # A direct square-root expansion after specialization is a different
    # implementation of the retained recurrence on that smaller slice.
    specialized_delta = delta.subs({y:0,p:0})
    root = S.sqrt(specialized_delta).series(x,0,14).removeO().expand()
    specialized_C = C.subs({y:0,p:0})
    radical = S.Poly(S.expand(-specialized_C+root),x)
    for n in range(9):
        actual = S.Poly(str(result[n]),y,p,q).as_expr().subs({y:0,p:0})
        expected = radical.coeff_monomial(x**(n+5))/2
        zero(f'independent radical expansion on y=p=0, x order {n}',actual-expected)
    return dict(parameters={'L':1,'M':13,'N':2,'kappa':1},free_coefficients=['y','p','q'],
                trace_and_first_normal_coefficient_degrees=jet_degrees,
                normal_orders=list(range(9)),coefficients=[str(f) for f in result],
                all_order_nontermination='Written binomial proof, not inferred from these nine coefficients')


def main():
    start = time.monotonic()
    x,y,p,q,r,L,M,N,kappa,b,c,d = S.symbols('x y p q r L M N kappa b c d')
    u = 1+x*y
    A = u**3*p+3*x*u**2*q-x**3*r
    B = y**2*u*(4+3*x*y)*p+(y+3*x*y**2*(4+3*x*y))*q+(2*x-3*x**2*y)*r
    psi = L*A*A+M*A+N*B
    zero('universal r-elimination identity',x*x*B+(2-3*x*y)*A-u*(u+1)*p-2*x*(2*u+1)*q)
    slice_psi = S.expand(psi.subs(x,0))
    zero('universal graph-independent slice',slice_psi-(L*p*p+M*p+N*(4*y*y*p+y*q)))
    direction = lambda expression:b*S.diff(expression,y)+c*S.diff(expression,p)+d*S.diff(expression,q)
    zero('complete second directional derivative on x=0',
         direction(direction(slice_psi))-(8*N*b*b*p+16*N*b*c*y+2*L*c*c+2*N*b*d))
    R = S.Function('R')(x,y,p,q)
    aq = S.diff(A.subs(r,R),q)
    zero('q derivative has immutable x coefficient',aq-(3*x*u*u-x**3*S.diff(R,q)))
    check('q derivative coefficient survives polynomial localization', (3*u*u).subs(x,0)==3)
    degree = S.Symbol('n',integer=True,positive=True)
    check('all-degree second derivative coefficient nonzero in characteristic zero',
          (2*degree*(2*degree-1)).is_positive)
    target = L*p*p+M*p+4*N*y*y*p+N*y*q+kappa*x*p
    zero('target cotangent Hessian determinant',S.hessian(target,(x,y,p,q)).det()-N*N*kappa*kappa)
    A0,psi0 = S.expand(A.subs(r,0)),S.expand(psi.subs(r,0))
    C = S.expand(N*(2-3*x*y)-M*x*x-2*L*x*x*A0)
    D = S.cancel((psi0-target)/x)
    check('divided target equation is polynomial',S.denom(D)==1)
    relation = L*x**5*r*r+C*r+D
    zero('complete graph-target quadratic relation',psi-target-x*relation)
    check('formal implicit derivative is the unit 2N',C.subs(x,0)==2*N)
    expected_trace = (kappa*p-3*(2*L*p+M)*(y*p+q)-N*(7*y**3*p+12*y*y*q))/(2*N)
    zero('exact formal trace',-D.subs(x,0)/(2*N)-expected_trace)
    delta = S.expand(C*C-4*L*x**5*D)
    check('discriminant is exactly linear in q',S.Poly(delta,q).degree()==1)
    zero('discriminant q coefficient',S.Poly(delta,q).coeff_monomial(q)+12*L*N*x**3*(x*y+2))
    check('discriminant constant x coefficient is 4N^2',delta.subs(x,0)==4*N*N)
    c_formula = lambda j:(N/L)*S.binomial(S.Rational(1,2),j)*(-6*L/N)**j
    leading = {j:str(S.factor(c_formula(j))) for j in range(2,7)}
    v2 = {L:1,M:13,N:2,kappa:1}
    points = [(1,-S.Rational(3,2),0,0,0),(-1,S.Rational(3,2),0,0,0)]
    grad = S.Matrix([S.diff(psi,v) for v in (x,y,p,q,r)]).subs(v2)
    g = S.Matrix([0,0,-S.Rational(1,2),0,0])
    for i,point in enumerate(points):
        values = dict(zip((x,y,p,q,r),point))
        zero(f'marked ambient collision at point {i}',grad.subs(values)-g)
        check(f'algebraic graph contains marked point {i}',relation.subs(v2).subs(values)==0)
        check(f'double root over marked point {i}',S.diff(relation,r).subs(v2).subs(values)==0)
        check(f'no regular graph through marked point {i}',S.diff(relation,q).subs(v2).subs(values)==3)
        check(f'simple discriminant ramification at point {i}',S.diff(delta,q).subs(v2).subs(values)==(-12 if i==0 else 12))
    series = normal_series(C.subs(v2),D.subs(v2),(x,y,p,q),delta.subs(v2),c_formula)
    sources = {}
    pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'hessian/tang-hc4-audit/source-v1.pdf']:
        sources[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('pinned structural dependency: '+name,sources[name]==pins[name])
    print(json.dumps(dict(outcome='PASS',checks=checks,check_count=len(checks),
        all_degree_arguments='Written no-plane, polynomial q-degree, nonsquare discriminant and binomial proofs',
        graph_degree_bound=None,formal_trace=str(S.factor(expected_trace)),
        leading_higher_q_coefficients=leading,finite_series=series,
        dependency_source_sha256=sources,elapsed_seconds=round(time.monotonic()-start,3),
        runtime={'python':platform.python_version(),'sympy':S.__version__,'flint':flint.__version__},
        formal_verification=False,external_review=False),indent=2))


if __name__=='__main__':
    main()
