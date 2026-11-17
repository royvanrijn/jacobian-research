"""Exact plane-map certificates for affine-linear coefficient graphs."""
import sympy as sp

def pencil_data(H,J,L,T,parameters=()):
    domain=sp.QQ.frac_field(*parameters) if parameters else sp.QQ
    pol=[sp.Poly(f,T,domain=domain) for f in (H,J,L)]
    D=pol[0].gcd(pol[1]).gcd(pol[2]).monic()
    reduced=[f.exquo(D).as_expr() for f in pol]
    h,j,l=reduced;S=sp.Symbol("S")
    A=sp.expand(h.subs(T,S)*l-h*l.subs(T,S))
    C=sp.expand(j.subs(T,S)*l-j*l.subs(T,S))
    N=max(sp.degree(f,T) for f in reduced)
    fiber_domain=sp.QQ.frac_field(*(tuple(parameters)+(T,)))
    distinct_degrees=len({sp.degree(f,T) for f in reduced})==3
    if sp.isprime(N) and distinct_degrees:
        # Basepoint-free nondegenerate plane series: N=e*deg(image),
        # and deg(image)>=2. This is an exact all-coefficient certificate,
        # avoiding an unnecessary large rational-function Euclidean run.
        fiber=sp.Poly(S-T,S,domain=fiber_domain)
        method="prime degree and three independent sections"
    else:
        fiber=sp.Poly(A,S,domain=fiber_domain).gcd(sp.Poly(C,S,domain=fiber_domain)).monic()
        method="exact fiber polynomial gcd"
    e=fiber.degree()
    if N%e:raise AssertionError("plane map degree does not divide section degree")
    return dict(H=H,J=J,L=L,D=D.as_expr(),base_degree=D.degree(),
      h=h,j=j,l=l,moving_degree=int(N),cover_degree=e,image_degree=int(N//e),
      fiber_polynomial=sp.factor(fiber.as_expr()),degree_certificate=method,moving_root_group=(
        "S_"+str(N) if e==1 else f"imprimitive with {N//e} blocks of size {e}"))

def graph_pencil(n,graphs):
    Q,R,B,T=sp.symbols("Q R B T")
    H=T**n;J=T**(n-1);L=T-B
    for k in range(2,n-1):
        h=sp.Poly(sp.sympify(graphs.get(k,0)),Q,R)
        if h.total_degree()>1:raise ValueError("graph not affine-linear in Q,R")
        H+=h.coeff_monomial(Q)*T**k
        J+=h.coeff_monomial(R)*T**k
        L+=h.coeff_monomial(1)*T**k
    return tuple(map(sp.expand,(H,J,L)))

def homogenize(f,T,X,Z,degree=None):
    f=sp.Poly(f,T);degree=f.degree() if degree is None else degree
    return sp.expand(sum(c*X**i*Z**(degree-i) for (i,),c in f.terms()))

def moving_field(H,J,L,T,D):
    Q,R,A=sp.symbols("Q R A")
    h=sp.cancel(H/D);j=sp.cancel(J/D);l=sp.cancel(L/D)
    W=sp.cancel(sp.diff(j,T)-j*sp.diff(h,T)/h)
    V=sp.cancel(sp.diff(l,T)-l*sp.diff(h,T)/h)
    q=sp.cancel(-(R*j+l)/h)
    derivative=sp.cancel(D*(R*W+V))
    if W==0:raise AssertionError("consecutive section degrees must have nonzero Wronskian")
    r=sp.cancel((A**sp.Symbol("m")/D-V)/W)
    return dict(q=q,W=W,V=V,derivative=derivative,r_formal=r)

def monic_cofactor_family(f,T):
    """Polynomial graph P=(f(T)-B)*(1+R*T^2+Q*(T^3-c*T^2)).

    f monic with f(0)=0,f'(0)=1. Its top two coefficients are Q,R.
    """
    Q,R,B=sp.symbols("Q R B");f=sp.Poly(f,T)
    if f.LC()!=1 or f.nth(0)!=0 or f.nth(1)!=1:raise ValueError("normalized monic f required")
    c=f.nth(f.degree()-1)
    D=f.as_expr()-B;C=1+R*T*T+Q*(T**3-c*T*T)
    return sp.expand(D*C),D,C
