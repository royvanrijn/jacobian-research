from sage.all import *
from itertools import product

K=GF(31)
T=PolynomialRing(K,"t"); t=T.gen()
F=FractionField(T)

# Verified rank-2 candidate surface.
A=T(8*t^5 + 8*t^4 + 10*t^3 + 11*t^2 + 23*t + 20)
B=T(t^8 + 28*t^7 + 6*t^6 + 27*t^5 + 19*t^4 + 17*t^2 + 20*t + 8)
E=EllipticCurve(F,[0,0,0,F(A),F(B)])

P1=E(
    F(29*t^2+27*t+11),
    F(t^4+6*t^3+19*t^2+2*t+3)
)
P2=E(
    F(22*t^4+27*t^3+15*t^2+21*t+5)/F(t-8)^2,
    F(t^7+16*t^6+22*t^4+11*t^3+12*t^2)/F(t-8)^3
)
Q=P1+P2

fibers={
    "I4_0":(K(0),K(18),4),
    "I4_1":(K(1),K(5),4),
    "I2_lam":(K(23),K(14),2),
    "I2_mu":(K(10),K(16),2),
}

def ord_at_poly(poly,a):
    poly=T(poly)
    if poly==0:
        return Infinity
    q=t-a
    n=0
    while poly(a)==0:
        poly=poly//q
        n+=1
    return n

def ord_at(fr,a):
    fr=F(fr)
    return ord_at_poly(fr.numerator(),a)-ord_at_poly(fr.denominator(),a)

def value_at(fr,a):
    fr=F(fr)
    if fr.denominator()(a)==0:
        return None
    return K(fr.numerator()(a))/K(fr.denominator()(a))

def component(section,a,s,n):
    x,y=section[0],section[1]
    xv=value_at(x,a); yv=value_at(y,a)
    if xv is not None and yv is not None and not (xv==s and yv==0):
        return {"class":"0","oriented":0,"valuations":None}
    if n==2:
        return {"class":"1","oriented":1,"valuations":None}
    if n!=4:
        raise NotImplementedError(n)

    tangent=3*s
    if not tangent.is_square():
        return {"class":"nonsplit?","oriented":None,"valuations":None}
    c=tangent.sqrt()
    X=x-F(s)
    vp=ord_at(y+F(c)*X,a)
    vm=ord_at(y-F(c)*X,a)
    # For uv=u_base^4, valuations are j and 4-j, up to branch orientation.
    vals=sorted([int(vp),int(vm)])
    if vals==[2,2]:
        cls="2"
    elif vals==[1,3]:
        cls="1/3"
    elif vals==[0,4]:
        # Section should then be on identity component; useful diagnostic.
        cls="0"
    else:
        cls=f"unexpected:{vals}"
    return {"class":cls,"oriented":(int(vp)%4,int(vm)%4),"valuations":vals}

print("E6CHECK|surface|A="+str(A)+"|B="+str(B),flush=True)
for name,P in [("P1",P1),("P2",P2),("P1+P2",Q),("-(P1+P2)",-Q)]:
    print(f"E6CHECK|section={name}|point={P}",flush=True)
    prof=[]
    for fn,(a,s,n) in fibers.items():
        c=component(P,a,s,n)
        prof.append(c["class"])
        print(f"E6CHECK|section={name}|fiber={fn}|class={c['class']}|oriented={c['oriented']}|vals={c['valuations']}",flush=True)
    print(f"E6CHECK|section={name}|profile={prof}",flush=True)

rels=[]
for m in range(-16,17):
    for n in range(-16,17):
        if m==n==0: continue
        if m*P1+n*P2==E(0):
            rels.append((m,n))
print(f"E6CHECK|relations16={rels}",flush=True)
