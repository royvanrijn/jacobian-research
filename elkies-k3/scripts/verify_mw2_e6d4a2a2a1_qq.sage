from sage.all import *


R=PolynomialRing(QQ,"t")
t=R.gen()
lam=QQ(9)/25
mu=QQ(49)/25
pole=QQ(16)/25
s=QQ(100)/1323

A=R([0,0,-QQ(32447500)/583443,-QQ(906250)/194481,
     QQ(31250000)/194481,-QQ(19531250)/194481])
B=R([0,0,0,QQ(300827000000)/2315685267,QQ(340001171875)/1029193452,
     -QQ(498857421875)/257298363,QQ(29541015625)/10501974,
     -QQ(152587890625)/85766121,QQ(152587890625)/343064484])
P1X=R([1,-QQ(800)/1323,QQ(625)/147])
P1Y=R([1,-QQ(400)/441,-QQ(394375)/18522,QQ(484375)/9261,
       -QQ(390625)/18522])
N=R([0,QQ(77824)/33075,QQ(12400)/1323,-QQ(30500)/1323,QQ(5000)/441])
M=R([0,0,QQ(4096)/343,QQ(281408)/9261,-QQ(1517000)/9261,
     QQ(1296875)/6174,-QQ(1015625)/9261,QQ(390625)/18522])
q=t-pole

assert P1Y**2==P1X**3+A*P1X+B
assert M**2==N**3+A*N*q**4+B*q**6
assert N(pole)!=0 and M(pole)!=0

Delta=-16*(4*A**3+27*B**2)


def valuation_at(poly,point):
    factor=t-point
    value=0
    while poly and poly(point)==0:
        poly//=factor
        value+=1
    return value


valuations=(
    valuation_at(Delta,0),valuation_at(Delta,1),valuation_at(Delta,lam),
    valuation_at(Delta,mu),24-Delta.degree(),
)
assert valuations==(6,3,3,2,8)
residual=Delta/(t**6*(t-1)**3*(t-lam)**3*(t-mu)**2)
residual=R(residual)
assert residual.degree()==2 and gcd(residual,residual.derivative())==1
assert all(residual(point)!=0 for point in (0,1,lam,mu))


def multiplicative_steps(A0,B0,X0,Y0,fiber_point,node):
    P=PolynomialRing(QQ,("u","xx","yy"))
    u,xx,yy=P.gens()
    At=R(A0(t+fiber_point)); Bt=R(B0(t+fiber_point))
    Ap=sum(P(c)*u**i for i,c in enumerate(At.list()))
    Bp=sum(P(c)*u**i for i,c in enumerate(Bt.list()))
    surface=yy**2-(node+xx)**3-Ap*(node+xx)-Bp
    Xt=R(X0(t+fiber_point)); Yt=R(Y0(t+fiber_point))
    if Xt(0)!=node or Yt(0)!=0:
        return 0
    sx=R((Xt-node)//t); sy=R(Yt//t)
    surface=P(surface(u,u*xx,u*yy)//u**2)
    steps=1
    while True:
        cx,cy=sx(0),sy(0)
        point={u:0,xx:cx,yy:cy}
        if any(surface.derivative(v).subs(point) for v in (u,xx,yy)):
            return steps
        sx=R((sx-cx)//t); sy=R((sy-cy)//t)
        surface=P(surface(u,cx+u*xx,cy+u*yy)//u**2)
        steps+=1
        if steps>5:
            raise RuntimeError("multiplicative resolution did not terminate")


def local_rational(numerator,power,point,precision=7):
    S=PowerSeriesRing(QQ,"u",default_prec=precision)
    u=S.gen()
    num=sum(S(numerator[i])*(S(point)+u)**i for i in range(numerator.degree()+1))
    den=(S(point)+u-pole)**power
    expansion=num/den
    return R([expansion[i] for i in range(precision)])


node_lam=P1X(lam)
node_mu=P1X(mu)
assert multiplicative_steps(A,B,P1X,P1Y,lam,node_lam)==1
assert multiplicative_steps(A,B,P1X,P1Y,mu,node_mu)==1
assert (P1X(1),P1Y(1))!=(s,0)

P2X1=local_rational(N,2,1)
P2Y1=local_rational(M,3,1)
A1=R(A(t+1)); B1=R(B(t+1))
assert multiplicative_steps(A1,B1,P2X1,P2Y1,0,s)==1
for point,node in ((lam,node_lam),(mu,node_mu)):
    assert (N(point),M(point))!=(node*(point-pole)**2,0)

# At I0*, P1 specializes to the smooth point (1,1).  P2 follows the simple
# root c of the first exceptional cubic, hence a nonzero D4 triality class.
assert (P1X(0),P1Y(0))==(1,1)
c=N[1]/pole**2
d4_cubic=c**3+A[2]*c+B[3]
assert d4_cubic==0 and 3*c**2+A[2]!=0
assert local_rational(M,3,0).valuation()==2

# At IV* infinity, both sections have (v_u(xbar),v_u(ybar))=(2,2).
# The second blowup has exceptional equation y2^2=B8; opposite leading signs
# are the two inverse nonzero E6 component classes.
assert P1X.degree()==2 and P1Y.degree()==4
assert N.degree()-2==2 and M.degree()-3==4
assert M[7]==-P1Y[4] and P1Y[4]**2==B[8]

pair_gcd=gcd(P1X*q**2-N,P1Y*q**3-M)
assert pair_gcd.monic()==t**2-QQ(67)/25*t+QQ(1008)/625
assert pair_gcd.degree()==2 and gcd(pair_gcd,pair_gcd.derivative())==1

# Shioda replay for profiles
# P1=(1,0,0,1,1;0), P2=(2,d1,1,0,0;1).
local11=QQ(4)/3+QQ(2)/3+QQ(1)/2
local22=QQ(4)/3+1+QQ(2)/3
local12=QQ(2)/3
h11=4-2*0-local11
h22=4+2*1-local22
h12=2+0+1-2-local12
H=matrix(QQ,[[h11,h12],[h12,h22]])
assert H==matrix(QQ,[[QQ(3)/2,QQ(1)/3],[QQ(1)/3,3]])
assert H.det()==QQ(79)/18

print("MW2QQ|fibers=IV*,I0*,I3,I3,I2,2I1|valuations="+",".join(map(str,valuations)),flush=True)
print("MW2QQ|profiles=P1:1,0,0,1,1;P2:2,d1,1,0,0|O=0,1|pair=2",flush=True)
print("MW2QQ|height_gram=(1/6)*[9,2;2,18]|det=79/18",flush=True)
print(f"MW2QQ|residual_delta={residual}|pair_gcd={pair_gcd.monic()}",flush=True)
print("MW2QQ|PASS",flush=True)
