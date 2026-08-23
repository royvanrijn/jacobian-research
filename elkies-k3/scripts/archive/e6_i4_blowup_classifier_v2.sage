from sage.all import *

K=GF(31)
L.<w>=GF(31^2)

T=PolynomialRing(K,"t")
t=T.gen()
FK=FractionField(T)

TU=PolynomialRing(L,"u")
u=TU.gen()
FU=FractionField(TU)

R.<U,X,Y>=PolynomialRing(L,3)

A0=T(8*t^5 + 8*t^4 + 10*t^3 + 11*t^2 + 23*t + 20)
B0=T(t^8 + 28*t^7 + 6*t^6 + 27*t^5 + 19*t^4 + 17*t^2 + 20*t + 8)

def to_local_uni(poly,a):
    return TU(sum(L(poly[i])*(L(a)+u)^i for i in range(poly.degree()+1)))

def uni_to_R(poly):
    poly=TU(poly)
    return R(sum(L(poly[i])*U^i for i in range(poly.degree()+1)))

def rf_to_local(f,a):
    f=FK(f)
    num=T(f.numerator())
    den=T(f.denominator())
    return FU(to_local_uni(num,a))/FU(to_local_uni(den,a))

def at0(f):
    f=FU(f)
    if f.denominator()(0)==0:
        raise ValueError("pole in chosen U-chart")
    return L(f.numerator()(0))/L(f.denominator()(0))

def divide_max_U(F):
    d=F.dict()
    if not d:
        return F, Infinity
    m=min(mon[0] for mon,c in d.items() if c)
    G=R(sum(c*U^(mon[0]-m)*X^mon[1]*Y^mon[2]
            for mon,c in d.items()))
    return G,ZZ(m)

def blowup_u_chart(F,xsec,ysec,x0,y0):
    G=R(F(U, L(x0)+U*X, L(y0)+U*Y))
    G,m=divide_max_U(G)
    return G,(xsec-L(x0))/u,(ysec-L(y0))/u,m

def eval3(F,x0,y0):
    return F(U=L(0),X=L(x0),Y=L(y0))

def singular_at(F,x0,y0):
    vals=[
        eval3(F,x0,y0),
        eval3(F.derivative(U),x0,y0),
        eval3(F.derivative(X),x0,y0),
        eval3(F.derivative(Y),x0,y0),
    ]
    return all(v==0 for v in vals),vals

def classify_i4(name,a,s,x_t,y_t):
    a=L(a); s=L(s)

    Al=to_local_uni(A0,K(a))
    Bl=to_local_uni(B0,K(a))

    AR=uni_to_R(Al)
    BR=uni_to_R(Bl)

    F=R(Y^2-(s+X)^3-AR*(s+X)-BR)

    xsec=rf_to_local(x_t,K(a))-s
    ysec=rf_to_local(y_t,K(a))

    x0=at0(xsec)
    y0=at0(ysec)

    print(f"I4BLOW|section={name}|stage=0|landing=({x0},{y0})",flush=True)

    if x0!=0 or y0!=0:
        print(f"I4BLOW|section={name}|class=0|depth=0",flush=True)
        return 0

    for depth in [1,2,3]:
        F,xsec,ysec,m=blowup_u_chart(F,xsec,ysec,x0,y0)
        x0=at0(xsec)
        y0=at0(ysec)

        on_surface=eval3(F,x0,y0)
        sing,derivs=singular_at(F,x0,y0)

        print(
            f"I4BLOW|section={name}|stage={depth}"
            f"|divideU={m}|landing=({x0},{y0})"
            f"|on_surface={on_surface}|singular={sing}|derivs={derivs}",
            flush=True
        )

        if on_surface != 0:
            print(f"I4BLOW|section={name}|ERROR=section_not_on_strict_transform",flush=True)
            return None

        if not sing:
            # For A3 resolution:
            # first exceptional hit => outer component
            # second residual blow-up => middle component
            cls = "outer(1/3)" if depth==1 else ("middle(2)" if depth==2 else f"resolved_depth_{depth}")
            print(f"I4BLOW|section={name}|class={cls}|depth={depth}",flush=True)
            return cls

    print(f"I4BLOW|section={name}|class=UNRESOLVED",flush=True)
    return None

E=EllipticCurve(FK,[0,0,0,FK(A0),FK(B0)])

P1=E(
    FK(29*t^2+27*t+11),
    FK(t^4+6*t^3+19*t^2+2*t+3)
)

P2=E(
    FK(22*t^4+27*t^3+15*t^2+21*t+5)/FK(t-8)^2,
    FK(t^7+16*t^6+22*t^4+11*t^3+12*t^2)/FK(t-8)^3
)

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
    for fn,a,s in [("I4_0",0,18),("I4_1",1,5)]:
        print(f"I4BLOW|BEGIN|section={name}|fiber={fn}",flush=True)
        classify_i4(name+"@"+fn,a,s,P[0],P[1])
