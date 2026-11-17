"""Exact bounded diagonal/unit controls, not a general holonomic source solver."""
from pathlib import Path
from math import prod
import hashlib,json
import sympy as s


def symbol_jet(expr,xis,L):
    if not isinstance(L,int) or L<0 or L>12:raise ValueError('Control jet order0..12 required')
    expr=s.sympify(expr)
    if expr.has(s.Float):raise ValueError('Exact coefficients required')
    eps=s.Dummy('eps')
    shifted=s.expand(expr).subs({v:eps*v for v in xis},simultaneous=True)
    return s.expand(s.cancel(s.expand(s.series(shifted,eps,0,L+1).removeO().subs(eps,1))))


def apply_symbol(expr,f,xis,xs):
    f=s.Poly(s.expand(f),*xs)
    if f.is_zero:return s.Integer(0)
    L=int(f.total_degree());jet=s.Poly(symbol_jet(expr,xis,L),*xis)
    total=0
    for alpha,c in jet.terms():
        derivative=f.as_expr()
        for x,a in zip(xs,alpha):
            if a:derivative=s.diff(derivative,x,a)
        total+=c*derivative
    return s.expand(total)


def iterate_symbol(mu,f,m,xis,xs):
    for _ in range(m):f=apply_symbol(mu,f,xis,xs)
    return s.expand(f)


def diagonal_coefficient(mu,P,alpha,m,xis,xs):
    f=s.expand(prod(x**a for x,a in zip(xs,alpha))*P**m)
    L=int(s.Poly(f,*xs).total_degree());am=s.Poly(symbol_jet(mu**m,xis,L),*xis)
    ws=s.symbols('w:'+str(len(xs)));shifted=s.Poly(f.subs({x:x+w for x,w in zip(xs,ws)},simultaneous=True),*ws)
    rho=s.symbols('r:'+str(len(xs)));polynomial=0
    for beta,c in am.terms():
        polynomial+=c*shifted.coeff_monomial(prod(w**b for w,b in zip(ws,beta)))*prod(s.factorial(b) for b in beta)*prod(r**b for r,b in zip(rho,beta))
    direct=iterate_symbol(mu,f,m,xis,xs)
    actual=s.expand(polynomial.subs({r:1 for r in rho}))
    assert s.expand(actual-direct)==0
    assert s.Poly(polynomial,*rho).total_degree()<=int(s.Poly(P,*xs).total_degree())*m+sum(alpha) if polynomial!=0 else True
    qs=s.symbols('q:'+str(len(xs)));translated=s.Poly(polynomial.subs({r:1+q for r,q in zip(rho,qs)},simultaneous=True),*qs)
    assert s.expand(translated.coeff_monomial(1)-actual)==0
    return actual


def run():
    x,y,z=s.symbols('x y z');X,Y,Z=s.symbols('X Y Z');a,t,r=s.symbols('a t r')
    models=[([x],[X],s.sqrt(1+X)-1,x+x*x),([x,y],[X,Y],X+Y*s.sqrt(1+X),x*y+y*y),([x,y,z],[X,Y,Z],X+Y*s.sqrt(1-Z),z*z)]
    diagonals=[];units=[]
    for xs,xis,mu,P in models:
        for alpha in [tuple(0 for _ in xs),tuple(1 if i==0 else 0 for i in range(len(xs)))]:
            for m in [1,2]:
                value=diagonal_coefficient(mu,P,alpha,m,xis,xs)
                diagonals.append(dict(dimension=len(xs),marker_powers=list(alpha),index=m,output=str(value)))
        U=s.exp(sum(xis))/(1-xis[-1]);Uinv=s.exp(-sum(xis))*(1-xis[-1])
        f=s.expand(xs[0]**2+xs[-1]**3+xs[0]*xs[-1])
        lhs=iterate_symbol(U*mu,f,2,xis,xs)
        rhs=iterate_symbol(U,iterate_symbol(mu,f,2,xis,xs),2,xis,xs)
        assert s.expand(lhs-rhs)==0
        recovered=iterate_symbol(Uinv,lhs,2,xis,xs)
        assert s.expand(recovered-iterate_symbol(mu,f,2,xis,xs))==0
        units.append(dict(dimension=len(xs),nonzero_unit_constant=True,inverse_recovered=True))
    # Square-zero coefficient model U=1+aD_x, inverse1-aD_x moduloa^2.
    f=x**3+2*x*x+x
    truncate=lambda p:s.Poly(s.expand(p),a).nth(0)+a*s.Poly(s.expand(p),a).nth(1)
    transformed=truncate(f+a*s.diff(f,x));recovered=truncate(transformed-a*s.diff(transformed,x));assert s.expand(recovered-f)==0
    # The degree-triangular unit matrix has invertible diagonal; determinantu0^3.
    u0,u1,u2=s.symbols('u0 u1 u2');matrix=s.Matrix([[u0,u1,2*u2],[0,u0,2*u1],[0,0,u0]])
    assert matrix.det()==u0**3 and s.simplify(matrix.inv()*matrix)==s.eye(3)
    J=sum(s.factorial(k)*r**k for k in range(9))
    ode=r*r*s.diff(J,r,2)+(3*r-1)*s.diff(J,r)+J
    assert s.series(ode,r,0,7).removeO()==0
    # Illustrative b(E_s)+lower-source-grade relation for1/(1-ts).
    q=s.symbols('q');F=1/(1-t*q)
    assert s.simplify(q*s.diff(F,q)-t*q*(q*s.diff(F,q)+F))==0
    # A genuinely algebraic nonpolynomial branch with a fixed transverse cutoff.
    mu=X+Y*s.sqrt(1-Z);P=z**3;safe=[]
    for m in [2,3]:
        out=iterate_symbol(mu,(x+y)*P**m,m,[X,Y,Z],[x,y,z]);assert out==0;safe.append(m)
    refusals=0
    for action in [lambda:symbol_jet(s.Float('0.5')*X,[X],2),lambda:symbol_jet(X,[X],-1),lambda:symbol_jet(X,[X],13)]:
        try:action()
        except ValueError:refusals+=1
        else:raise AssertionError('Expected guard refusal')
    return dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                diagonal_controls=diagonals,formal_unit_controls=units,square_zero_unit_inverse=True,unit_matrix_determinant=str(matrix.det()),
                factorial_kernel_ode_zero_through_degree=6,illustrative_taylor_recurrence=True,algebraic_transverse_zero_indices=safe,
                guard_refusals=refusals,supporting_checker_hashes={},bounded_jet_action_and_diagonal_calculators_implemented=True,
                universal_holonomicity_proved_by_finite_controls=False,general_transcendental_symbol_source_claimed=False,
                full_GVC_for_algebraic_symbols_claimed=False,all_ring_algebraic_linearity_claimed=False,
                general_effective_source_solver_implemented=False,full_exceptional_MZ_decided=False,
                independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)


if __name__=='__main__':
    receipt=run();output=Path(__file__).with_name('algebraic_formal_symbols_verification.json')
    if output.exists():assert json.loads(output.read_text())==receipt,'Preserve immutable receipt'
    else:output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(completed=True,checker_sha256=receipt['checker_sha256'],diagonal_controls=len(receipt['diagonal_controls']),unit_models=len(receipt['formal_unit_controls']),guard_refusals=receipt['guard_refusals'])))
