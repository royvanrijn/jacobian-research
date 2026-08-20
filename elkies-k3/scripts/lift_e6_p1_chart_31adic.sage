from sage.all import *
import argparse

ap=argparse.ArgumentParser(description='Lift the verified E6/P1 GF(31) point to a local 31-adic chart.')
ap.add_argument('--prec',type=int,default=20,help='31-adic precision in base-31 digits')
ap.add_argument('--r0',type=int,default=4)
ap.add_argument('--s0',type=int,default=18)
ap.add_argument('--x1',type=int,default=27)
args=ap.parse_args()

p=31
# Solved coordinates on the verified residue-class point.
sol_names=['a1','a2','a4','lam','mu','s1','sl','sm']
seed=[4,16,6,24,18,23,23,4]

# Build the characteristic-zero equations BEFORE reduction mod p.
R=PolynomialRing(QQ,sol_names,order='degrevlex')
d=R.gens_dict(); RF=FractionField(R)
a1,a2,a4,lam,mu,s1,sl,sm=[RF(d[n]) for n in sol_names]
r0=QQ(args.r0); s0=QQ(args.s0); x1=QQ(args.x1)
x0=r0**2-2*s0
y0=r0*(r0**2-3*s0)
y1=(a1+3*(r0**2-s0)*x1)/(2*r0)

# Reintroduce only the four triangularly eliminated variables locally.
U=PolynomialRing(RF,['a3','b4','b5','y2'])
a3,b4,b5,y2=U.gens()
Ut=PolynomialRing(U,'t'); t=Ut.gen()

def L(z): return U(RF(z))
a0=-3*s0**2
b0=2*s0**3
b1=-s0*a1
b2=a1**2/(12*s0)-s0*a2
b3=(a1**3+36*a1*a2*s0**2-216*RF(a3)*s0**4)/(216*s0**3)
a5=-3*s1**2-(a0+a1+RF(a3)+a2+a4)

# Solve b6,b7 linearly from the t=1 singular-root equations.
V=PolynomialRing(FractionField(U),['B6']); B6=V.gen()
Vt=PolynomialRing(V,'t'); tt=Vt.gen()
AA=sum(V(L(c))*tt**i for i,c in enumerate([a0,a1,a2,a3,a4,a5]))
B7=2*V(L(s1))**3-(V(L(b0))+V(L(b1))+V(L(b2))+V(L(b3))+V(b4)+V(b5)+B6)-1
BB=V(L(b0))+V(L(b1))*tt+V(L(b2))*tt**2+V(L(b3))*tt**3+V(b4)*tt**4+V(b5)*tt**5+B6*tt**6+B7*tt**7+tt**8
EE=BB.derivative(tt)(1)+V(L(s1))*AA.derivative(tt)(1)
c6=EE.derivative(B6)
b6=FractionField(U)(-EE.subs({B6:0})/c6)
b7=FractionField(U)(B7.subs({B6:V(b6)}))

FU=FractionField(U)
T=PolynomialRing(FU,'t'); t=T.gen()
A=T(FU(L(a0)))+T(FU(L(a1)))*t+T(FU(L(a2)))*t**2+T(FU(a3))*t**3+T(FU(L(a4)))*t**4+T(FU(L(a5)))*t**5
B=T(FU(L(b0)))+T(FU(L(b1)))*t+T(FU(L(b2)))*t**2+T(FU(L(b3)))*t**3+T(FU(b4))*t**4+T(FU(b5))*t**5+T(b6)*t**6+T(b7)*t**7+t**8
X=T(FU(L(x0)))+T(FU(L(x1)))*t+T(FU(L(s1-x0-x1)))*t**2
Y=T(FU(L(y0)))+T(FU(L(y1)))*t+T(FU(y2))*t**2+T(FU(L(-1-y0-y1))-FU(y2))*t**3+t**4
Sec=Y**2-X**3-A*X-B

# Exact triangular elimination in the already-known order.
subs={}
def se(e):
    e=FU(e)
    for _ in range(5):
        old=e
        if subs: e=FU(e.subs(subs))
        if e==old: break
    return e

def solve_linear(e,v,label):
    e=se(e); vv=FU(v)
    c=e.derivative(vv)
    if c==0 or c.derivative(vv)!=0: raise RuntimeError(label+' not affine-linear')
    out=FU(-e.subs({vv:0})/c)
    subs[vv]=out
    print(f'E6PADIC|elim={label}',flush=True)
    return out

b4e=solve_linear(Sec[7],b4,'b4')
b5e=solve_linear(Sec[6],b5,'b5')
y2e=solve_linear(Sec[5],y2,'y2')
a3e=solve_linear(Sec[4],a3,'a3')

# Substitute the four eliminations into A,B and build the eight independent fiber equations.
def subpoly(P):
    return sum(se(P[i])*t**i for i in range(P.degree()+1))
A4=subpoly(A); B4=subpoly(B)
Ap=A4.derivative(t); Bp=B4.derivative(t)
Delta=-16*(4*A4**3+27*B4**2)
D2=Delta.derivative(t,2); D3=Delta.derivative(t,3)
eq_fr=[D2(1),D3(1),
       A4(lam)+3*sl**2, B4(lam)-2*sl**3, Bp(lam)+sl*Ap(lam),
       A4(mu)+3*sm**2, B4(mu)-2*sm**3, Bp(mu)+sm*Ap(mu)]

# Convert each rational equation to its numerator. Denominators are units at the seed.
eqs=[]
for e in eq_fr:
    e=se(e)
    eqs.append(R(e.numerator()))

# Verify exact equations reduce to the already-verified GF(31) point.
K=GF(p); RK=PolynomialRing(K,sol_names); kd=RK.gens_dict()
phi=R.hom([kd[n] for n in sol_names],RK)
P0={kd[n]:K(v) for n,v in zip(sol_names,seed)}
mods=[phi(e).subs(P0) for e in eqs]
print('E6PADIC|mod31_zero='+str(all(v==0 for v in mods)),flush=True)
J=Matrix(R,[[e.derivative(d[n]) for n in sol_names] for e in eqs])
JK=Matrix(K,[[phi(x).subs(P0) for x in row] for row in J.rows()])
print(f'E6PADIC|jac_rank={JK.rank()}|jac_det={int(JK.det())}',flush=True)
if JK.rank()!=8: raise RuntimeError('singular Jacobian')

# Newton iteration in Q_31. The fixed r0,s0,x1 values define this local slice.
Q=Qp(p,prec=args.prec+8,type='capped-rel')
RQ=PolynomialRing(Q,sol_names); qd=RQ.gens_dict()
phiq=R.hom([qd[n] for n in sol_names],RQ)
Fq=[phiq(e) for e in eqs]
Jq=Matrix(RQ,[[f.derivative(qd[n]) for n in sol_names] for f in Fq])
x=vector(Q,[Q(v) for v in seed])

def ev(poly,xv): return poly(**{n:xv[i] for i,n in enumerate(sol_names)})

for it in range(12):
    fv=vector(Q,[ev(f,x) for f in Fq])
    vals=[v.valuation() if v!=0 else Infinity for v in fv]
    minv=min(vals)
    print(f'E6PADIC|iter={it}|min_v={minv}',flush=True)
    if minv>=args.prec: break
    jv=Matrix(Q,[[ev(Jq[i,j],x) for j in range(8)] for i in range(8)])
    dx=jv.solve_right(fv)
    x=x-dx

M=ZZ(p)**args.prec
print(f'E6PADIC|modulus=31^{args.prec}',flush=True)
for n,v in zip(sol_names,x):
    # Integral branch: print canonical residue modulo p^prec.
    if v.valuation()<0: raise RuntimeError(f'nonintegral lift for {n}: {v}')
    r=ZZ(v.lift()) % M
    bal=r if r<=M//2 else r-M
    print(f'E6PADIC|coord={n}|residue={r}|balanced={bal}',flush=True)

print('E6PADIC|params|r0=%s|s0=%s|x1=%s'%(r0,s0,x1),flush=True)
