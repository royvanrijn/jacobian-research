from sage.all import *
from itertools import product

K=GF(31)
T=PolynomialRing(K,"t"); t=T.gen()
F=FractionField(T)

A=T(8*t^5 + 8*t^4 + 10*t^3 + 11*t^2 + 23*t + 20)
B=T(t^8 + 28*t^7 + 6*t^6 + 27*t^5 + 19*t^4 + 17*t^2 + 20*t + 8)
E=EllipticCurve(F,[0,0,0,F(A),F(B)])

P1=E(F(29*t^2+27*t+11),F(t^4+6*t^3+19*t^2+2*t+3))
P2=E(
    F(22*t^4+27*t^3+15*t^2+21*t+5)/F(t-8)^2,
    F(t^7+16*t^6+22*t^4+11*t^3+12*t^2)/F(t-8)^3
)

def solve_Y(X,top):
    # Y=y0+y1 t+y2 t^2+y3 t^3+top t^4.
    ys=[K(0)]*4
    def mk():
        return sum(ys[i]*t^i for i in range(4))+top*t^4
    # Coefficients 7,6,5,4 solve y3,y2,y1,y0 successively.
    for k,idx in [(7,3),(6,2),(5,1),(4,0)]:
        Y=mk()
        e=Y^2-X^3-A*X-B
        c=e[k]
        # coefficient in y_idx is 2*top
        ys[idx]-=c/(2*top)
    Y=mk()
    return Y if Y^2==X^3+A*X+B else None

# Precompute bounded subgroup <P1,P2>.
known={}
for m in range(-12,13):
    for n in range(-12,13):
        P=m*P1+n*P2
        known[(P[0],P[1])]=(m,n)

hits=0
new=[]
tested=0

for top in (K(1),K(-1)):
    for x0 in K:
        for x1 in K:
            for x2 in K:
                tested+=1
                X=T(x0+x1*t+x2*t^2)
                Y=solve_Y(X,top)
                if Y is None:
                    continue
                P=E(F(X),F(Y))
                hits+=1
                rel=known.get((P[0],P[1]))
                tag="known" if rel is not None else "NEW"
                print(f"E6POLY|{tag}|X={X}|Y={Y}|relation={rel}",flush=True)
                if rel is None:
                    new.append(P)

print(f"E6POLY|done|tested={tested}|hits={hits}|new={len(new)}",flush=True)
