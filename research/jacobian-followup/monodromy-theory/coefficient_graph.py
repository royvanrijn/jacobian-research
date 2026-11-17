"""Polynomial completion of lower marked-root coefficients minus h_j(A0).

The all-degree induction and monodromy argument are in COEFFICIENT_GRAPHS.md.
This module constructs honest polynomial coordinates by exact division.
"""
import sympy as sp
from marked_wreath import polynomial_map

def exact_x2(num,x):
    return sp.Poly(sp.expand(num),x).exquo(sp.Poly(x*x,x)).as_expr()

def completion(n,graphs=None):
    """graphs maps j to a polynomial in B, 2<=j<=n-2.

    Returns the full inverse source chart and its three-variable zero graph.
    The current-coordinate first jets are used, never rational elimination.
    """
    if n<3:raise ValueError("n>=3 required")
    graphs={} if graphs is None else graphs
    B=sp.Symbol("B");orig,T,phi,F=polynomial_map(n);x,y=orig[:2]
    Cs={j:phi.nth(j) for j in range(2,n+1)};b=-phi.nth(0)
    mapping={v:v for v in orig};slots={2:y,**{j:orig[j-1] for j in range(3,n+1)}}
    stages=[];vs=[]
    for k in range(2,n-1):
        u,z=slots[k],slots[k+2]
        current_b=sp.expand(b.xreplace(mapping))
        h=sp.sympify(graphs.get(k,0))
        if h.free_symbols-{B}:raise ValueError("h_j must be polynomial in B only")
        sp.Poly(h,B)
        f=sp.expand(Cs[k].xreplace(mapping)-h.subs(B,current_b))
        if sp.diff(f,z)!=x*x:raise AssertionError(("highest coordinate coefficient",n,k))
        f0=f.subs(z,0)
        leading=sp.expand(f0.subs(x,0));c=sp.diff(leading,u)
        if c.free_symbols or c==0:raise AssertionError(("nonconstant slope",n,k,c))
        beta=sp.expand(leading-c*u)
        jet=sp.expand(sp.diff(f0,x).subs(x,0))
        v=sp.Symbol(f"v{k}");w=sp.Symbol(f"s{k+2}")
        psi=sp.expand((v-beta)/c-x/c*jet.subs(u,(v-beta)/c))
        newu=sp.expand(psi+x*x*w)
        newz=exact_x2(v-f0.subs(u,newu),x)
        forward_w=exact_x2(u-psi.subs(v,f),x)
        stages.append(dict(k=k,u=u,z=z,v=v,w=w,f=f,f0=f0,c=c,beta=beta,jet=jet,
          psi=psi,newu=newu,newz=newz,forward_w=forward_w))
        repl={u:newu,z:newz}
        mapping={old:sp.expand(value.xreplace(repl)) for old,value in mapping.items()}
        slots[k]=v;slots[k+2]=w;vs.append(v)
    endpoint_vars=(x,slots[n-1],slots[n])
    selected={v:0 for v in vs}
    zero_mapping={old:sp.expand(value.subs(selected)) for old,value in mapping.items()}
    # Output convention is (r=C_(n-1),q=C_n,b=A0).
    endpoint=[sp.expand(expr.xreplace(zero_mapping)) for expr in (Cs[n-1],Cs[n],b)]
    return dict(n=n,m=n-2,original=orig,T=T,phi=phi,coefficients=Cs,b=b,
      graphs=graphs,mapping=mapping,stages=stages,parameters=vs,endpoint_vars=endpoint_vars,
      zero_mapping=zero_mapping,endpoint=endpoint,
      expected_jacobian=(-1)**(n-2)*sp.Rational(n-2,2))

def morse_weights(n,graphs):
    """Sufficient graph gate after r=0: c+j-1 > a*(n-1) for q^a*b^c."""
    Q,R,B=sp.symbols("Q R B");rows=[]
    for j in range(2,n-1):
        h=sp.expand(sp.sympify(graphs.get(j,0)).subs(R,0))
        for (a,c),coeff in sp.Poly(h,Q,B).terms():
            if coeff:
                rows.append(dict(j=j,q_power=a,b_power=c,coefficient=str(coeff),
                  epsilon_power=c+j-1-a*(n-1)))
    return rows

def passes_morse_gate(n,graphs):
    return all(r["epsilon_power"]>0 for r in morse_weights(n,graphs))

def scaled_root_polynomial(n,graphs):
    Q,R,B,T,e,Z,V=sp.symbols("Q R B T e Z V")
    P=Q*T**n+R*T**(n-1)+sum(sp.sympify(graphs.get(j,0))*T**j for j in range(2,n-1))+T-B
    scaled=sp.expand(P.subs({Q:e**(-(n-1)),R:0,B:e*V,T:e*Z},simultaneous=True)/e)
    return P,scaled,(Q,R,B,T,e,Z,V)
