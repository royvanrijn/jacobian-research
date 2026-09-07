from sage.all import *

# Exact geometric I4 classifier on the verified GF(31) E6 surface.
# We only need the distinction:
#   0       = identity component (section misses the node)
#   outer   = component 1 or 3
#   middle  = component 2
#
# For an A3 (= I4) surface singularity, a section through the node:
#   - resolves after the first blow-up -> outer component (1/3)
#   - still passes through the residual singularity and resolves after
#     the second blow-up -> middle component (2)
#
# Work over GF(31^2) so all exceptional geometry is visible.

K=GF(31)
L.<w>=GF(31^2)

TU=PolynomialRing(L,"u")
u=TU.gen()
FU=FractionField(TU)

R.<U,X,Y>=PolynomialRing(L,3)

# Original t-polynomials, re-expanded locally using t=a+u.
T=PolynomialRing(K,"t")
t=T.gen()

A0=T(8*t^5 + 8*t^4 + 10*t^3 + 11*t^2 + 23*t + 20)
B0=T(t^8 + 28*t^7 + 6*t^6 + 27*t^5 + 19*t^4 + 17*t^2 + 20*t + 8)

def to_local_poly(poly,a):
    # K[t] -> L[u], t |-> a+u
    return TU(sum(L(poly[i])*(L(a)+u)^i for i in range(poly.degree()+1)))

def rf_to_local(f,a):
    f=factorless = f
    num=T(f.numerator())
    den=T(f.denominator())
    return FU(to_local_poly(num,a))/FU(to_local_poly(den,a))

def at0(f):
    f=FU(f)
    if f.denominator()(0)==0:
        raise ValueError("section coordinate has a pole in current blow-up chart")
    return L(f.numerator()(0))/L(f.denominator()(0))

def divide_max_U(F):
    if F==0:
        return F,Infinity
    m=min(mon[0] for mon,coeff in F.dict().items() if coeff)
    if m:
        F=R(sum(coeff*U^(mon[0]-m)*X^mon[1]*Y^mon[2]
                for mon,coeff in F.dict().items()))
    return F,ZZ(m)

def blowup_u_chart(F,xsec,ysec,x0,y0):
    # Blow up the center (U,X,Y)=(0,x0,y0), taking the U-chart:
    # X = x0 + U*X1, Y = y0 + U*Y1.
    G=R(F(U, L(x0)+U*X, L(y0)+U*Y))
    G,m=divide_max_U(G)

    xnew=(xsec-L(x0))/u
    ynew=(ysec-L(y0))/u
    return G,FU(xnew),FU(ynew),m

def is_singular(F,x0,y0):
    pt={U:L(0),X:L(x0),Y:L(y0)}
    vals=[
        F.subs(pt),
        F.derivative(U).subs(pt),
        F.derivative(X).subs(pt),
        F.derivative(Y).subs(pt),
    ]
    return all(v==0 for v in vals),vals

def classify_i4(name,a,s,x_t,y_t,verbose=True):
    a=L(a); s=L(s)

    Al=to_local_poly(A0,a)
    Bl=to_local_poly(B0,a)

    # Center x at the nodal x-coordinate s.
    # Local surface:
    #   Y^2 = (s+X)^3 + A(a+u)(s+X) + B(a+u)
    F=R(Y^2-(s+X)^3-R(Al)*(s+X)-R(Bl))

    xsec=rf_to_local(x_t,a)-s
    ysec=rf_to_local(y_t,a)

    x0=at0(xsec)
    y0=at0(ysec)

    if verbose:
        print(f"I4BLOW|section={name}|a={a}|stage=0|landing=({x0},{y0})",flush=True)

    if x0!=0 or y0!=0:
        print(f"I4BLOW|section={name}|a={a}|class=0|depth=0",flush=True)
        return 0

    # Resolve along the section. For I4 at most two point blow-ups are needed
    # to distinguish outer vs middle.
    for depth in [1,2,3]:
        F,xsec,ysec,m=blowup_u_chart(F,xsec,ysec,x0,y0)
        x0=at0(xsec)
        y0=at0(ysec)
        sing,derivs=is_singular(F,x0,y0)

        if verbose:
            print(
                f"I4BLOW|section={name}|a={a}|stage={depth}"
                f"|divideU={m}|landing=({x0},{y0})|singular={sing}"
                f"|derivs={derivs}",
                flush=True
            )

        if not sing:
            if depth==1:
                cls="outer(1/3)"
            elif depth==2:
                cls="middle(2)"
            else:
                cls=f"resolved_depth_{depth}"
            print(f"I4BLOW|section={name}|a={a}|class={cls}|depth={depth}",flush=True)
            return cls

        # next blow-up is centered at the residual singular point hit by section
        # in this U-chart
        # Verify it really is on U=0 strict transform before proceeding.
        x0=L(x0); y0=L(y0)

    print(f"I4BLOW|section={name}|a={a}|class=UNRESOLVED|depth=3",flush=True)
    return None


# Build section coordinates over K(t).
FK=FractionField(T)

P1x=FK(29*t^2+27*t+11)
P1y=FK(t^4+6*t^3+19*t^2+2*t+3)

P2x=FK(22*t^4+27*t^3+15*t^2+21*t+5)/FK(t-8)^2
P2y=FK(t^7+16*t^6+22*t^4+11*t^3+12*t^2)/FK(t-8)^3

# Use EC group law to get combinations exactly.
EK=EllipticCurve(FK,[0,0,0,FK(A0),FK(B0)])
P1=EK(P1x,P1y)
P2=EK(P2x,P2y)

sections=[
    ("P1",P1),
    ("-P1",-P1),
    ("P2",P2),
    ("-P2",-P2),
    ("P1+P2",P1+P2),
    ("-(P1+P2)",-(P1+P2)),
    ("P1-P2",P1-P2),
]

for name,P in sections:
    for fn,a,s in [
        ("I4_0",0,18),
        ("I4_1",1,5),
    ]:
        print(f"I4BLOW|BEGIN|section={name}|fiber={fn}",flush=True)
        classify_i4(name+"@"+fn,a,s,P[0],P[1])
