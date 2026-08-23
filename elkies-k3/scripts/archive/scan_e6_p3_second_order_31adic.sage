from sage.all import *
import argparse, time

ap=argparse.ArgumentParser(description="Second-order 31-adic obstruction scan for P3 on the smooth E6/P1 chart.")
ap.add_argument("--progress",type=int,default=1000)
ap.add_argument("--max-examples",type=int,default=30)
args=ap.parse_args()

K=GF(31)
MOD=31
inv2=K(2)**-1

surf_names=["a1","a2","a4","lam","mu","s1","sl","sm"]
param_names=["r0","s0","x1"]
p3_names=["u1","u2","v1","v2","v3"]
all_names=surf_names+param_names+p3_names

R=PolynomialRing(K,all_names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)

a1,a2,a4,lam,mu,s1,sl,sm=[RF(d[n]) for n in surf_names]
r0,s0,x1=[RF(d[n]) for n in param_names]
u1,u2,v1,v2,v3=[RF(d[n]) for n in p3_names]

# Auxiliary variables used only for the exact triangular P1 elimination.
U=PolynomialRing(RF,["a3","b4","b5","y2"])
a3,b4,b5,y2=U.gens()
FU=FractionField(U)
Ut=PolynomialRing(FU,"t"); t=Ut.gen()

def L(z): return FU(U(RF(z)))
def LU(z): return FU(U(z))

x0=r0**2-2*s0
y0=r0*(r0**2-3*s0)
y1=(a1+3*(r0**2-s0)*x1)/(2*r0)

a0=-3*s0**2
b0=2*s0**3
b1=-s0*a1
b2=a1**2/(12*s0)-s0*a2
b3=(L(a1**3+36*a1*a2*s0**2)-216*FU(a3)*s0**4)/(216*s0**3)
a5=L(-3*s1**2-(a0+a1+a2+a4))-FU(a3)

V=PolynomialRing(FU,["B6"]); B6=V.gen()
Vt=PolynomialRing(V,"q"); q=Vt.gen()
AA=sum(V(LU(c))*q**i for i,c in enumerate([a0,a1,a2,a3,a4,a5]))
B7=2*V(L(s1))**3-(V(L(b0))+V(L(b1))+V(L(b2))+V(b3)+V(b4)+V(b5)+B6)-1
BB=(V(L(b0))+V(L(b1))*q+V(L(b2))*q**2+V(b3)*q**3+
    V(b4)*q**4+V(b5)*q**5+B6*q**6+B7*q**7+q**8)
ee=BB.derivative(q)(1)+V(L(s1))*AA.derivative(q)(1)
b6=FU(-ee.subs({B6:0})/ee.derivative(B6))
b7=FU(B7.subs({B6:V(b6)}))

A=L(a0)+L(a1)*t+L(a2)*t**2+FU(a3)*t**3+L(a4)*t**4+FU(a5)*t**5
B=L(b0)+L(b1)*t+L(b2)*t**2+FU(b3)*t**3+FU(b4)*t**4+FU(b5)*t**5+b6*t**6+b7*t**7+t**8
X1=L(x0)+L(x1)*t+L(s1-x0-x1)*t**2
Y1=L(y0)+L(y1)*t+FU(y2)*t**2+(L(-1-y0-y1)-FU(y2))*t**3+t**4
Sec=Y1**2-X1**3-A*X1-B

subs={}
def settle(e):
    e=FU(e)
    for _ in range(8):
        old=e
        if subs:
            e=FU(e.subs(subs))
        if e==old:
            break
    return e

def solve_linear(e,v):
    e=settle(e); vv=FU(v)
    c=e.derivative(vv)
    if c==0 or c.derivative(vv)!=0:
        raise RuntimeError(f"{v} not affine-linear")
    out=FU(-e.subs({vv:0})/c)
    subs[vv]=out
    return out

print("E6P3ORD2|stage=triangular",flush=True)
solve_linear(Sec[7],b4)
solve_linear(Sec[6],b5)
solve_linear(Sec[5],y2)
solve_linear(Sec[4],a3)

def subpoly(P):
    return sum(settle(P[i])*t**i for i in range(P.degree()+1))

A=subpoly(A); B=subpoly(B)
Ap=A.derivative(t); Bp=B.derivative(t)
Delta=-16*(4*A**3+27*B**2)

fiber_fu=[
    Delta.derivative(t,2)(1),
    Delta.derivative(t,3)(1),
    A(L(lam))+3*L(sl)**2,
    B(L(lam))-2*L(sl)**3,
    Bp(L(lam))+L(sl)*Ap(L(lam)),
    A(L(mu))+3*L(sm)**2,
    B(L(mu))-2*L(sm)**3,
    Bp(L(mu))+L(sm)*Ap(L(mu)),
]

X3=L(s0)+L(u1)*t+L(u2)*t**2
Y3=L(v1)*t+L(v2)*t**2+L(v3)*t**3+t**4
S3=Y3**2-X3**3-A*X3-B
p3_fu=[S3[k] for k in range(8)]

# Collapse Frac(U) -> Frac(R), with no cross-parent substitution.
def collapse(e):
    e=settle(e)
    nu=U(e.numerator()); de=U(e.denominator())
    zero=(0,)*U.ngens()
    for mon in nu.dict():
        if tuple(mon)!=zero:
            raise RuntimeError("aux variable survived in numerator")
    for mon in de.dict():
        if tuple(mon)!=zero:
            raise RuntimeError("aux variable survived in denominator")
    n=RF(next((c for mon,c in nu.dict().items() if tuple(mon)==zero),0))
    q0=RF(next((c for mon,c in de.dict().items() if tuple(mon)==zero),0))
    return n/q0

fiber=[collapse(f) for f in fiber_fu]
p3=[collapse(f) for f in p3_fu]

base={
    "a1":4,"a2":16,"a4":6,"lam":24,"mu":18,"s1":23,"sl":23,"sm":4,
    "r0":4,"s0":18,"x1":27,
    "u1":8,"u2":28,"v1":0,"v2":1,"v3":29
}
P={d[n]:K(v) for n,v in base.items()}

def eval_rf(f):
    f=RF(f)
    num=R(f.numerator()).subs(P)
    den=R(f.denominator()).subs(P)
    if den==0:
        raise ZeroDivisionError("denominator vanished at base point")
    return K(num)/K(den)

print("E6P3ORD2|stage=first_derivatives",flush=True)

Js=Matrix(K,[[eval_rf(f.derivative(d[n])) for n in surf_names] for f in fiber])
Jp=Matrix(K,[[eval_rf(f.derivative(d[n])) for n in param_names] for f in fiber])
Gs=Matrix(K,[[eval_rf(f.derivative(d[n])) for n in surf_names] for f in p3])
Gp=Matrix(K,[[eval_rf(f.derivative(d[n])) for n in param_names] for f in p3])
Gu=Matrix(K,[[eval_rf(f.derivative(d[n])) for n in p3_names] for f in p3])

if Js.rank()!=8:
    raise RuntimeError(f"Js rank={Js.rank()} != 8")

JsInv=Js.inverse()
Spar=-JsInv*Jp                    # ds1 = Spar * dp
Meff=Gp+Gs*Spar                   # == Gp-Gs Js^-1 Jp

left_basis=Gu.left_kernel().basis()
ker_basis=Gu.right_kernel().basis()

print(f"E6P3ORD2|Js_rank={Js.rank()}|Gu_rank={Gu.rank()}|left_dim={len(left_basis)}|ker_dim={len(ker_basis)}|Meff_rank={Meff.rank()}",flush=True)

# Particular first-order P3 correction linear in dp.
# Solve Gu * Upar[:,j] = -Meff[:,j].
def particular(A,b):
    # Sage's tiny solve may call FFLAS on this Mac, so use explicit RREF over GF(31).
    m=A.nrows(); n=A.ncols()
    aug=[[int(A[i,j])%MOD for j in range(n)]+[int(b[i])%MOD] for i in range(m)]
    piv=[]; r=0
    for c in range(n):
        q=next((i for i in range(r,m) if aug[i][c]),None)
        if q is None: continue
        aug[r],aug[q]=aug[q],aug[r]
        inv=pow(aug[r][c],-1,MOD)
        aug[r]=[(x*inv)%MOD for x in aug[r]]
        for i in range(m):
            if i!=r and aug[i][c]:
                z=aug[i][c]
                aug[i]=[(aug[i][j]-z*aug[r][j])%MOD for j in range(n+1)]
        piv.append(c); r+=1
    for i in range(r,m):
        if all(aug[i][j]==0 for j in range(n)) and aug[i][n]!=0:
            raise RuntimeError("inconsistent first-order P3 correction")
    x=[0]*n
    for row,c in enumerate(piv):
        x[c]=aug[row][n]
    return vector(K,x)

Ucols=[]
for j in range(3):
    Ucols.append(particular(Gu,-Meff.column(j)))
Upar=Matrix(K,5,3,lambda i,j: Ucols[j][i])
Kmat=Matrix(K,5,len(ker_basis),lambda i,j: ker_basis[j][i])

# w = [surface(8), params(3), P3(5)] = Amap*dp + Bmap*z
Amap=block_matrix(K,[[Spar],[identity_matrix(K,3)],[Upar]])
Bmap=block_matrix(K,[[zero_matrix(K,11,2)],[Kmat]])
Cmap=block_matrix(K,[[Spar],[identity_matrix(K,3)]])

print("E6P3ORD2|stage=hessians",flush=True)

# Hessian evaluated at base.  Rational differentiation is okay because all denominators are units.
def hessian_at(f,varnames):
    n=len(varnames)
    H=Matrix(K,n,n)
    for i,ni in enumerate(varnames):
        fi=f.derivative(d[ni])
        for j in range(i,n):
            v=eval_rf(fi.derivative(d[varnames[j]]))
            H[i,j]=v; H[j,i]=v
    return H

HF=[]
for i,f in enumerate(fiber):
    t0=time.time()
    HF.append(hessian_at(f,surf_names+param_names))
    print(f"E6P3ORD2|HF={i+1}/8|seconds={time.time()-t0:.2f}",flush=True)

HG=[]
for i,f in enumerate(p3):
    t0=time.time()
    HG.append(hessian_at(f,all_names))
    print(f"E6P3ORD2|HG={i+1}/8|seconds={time.time()-t0:.2f}",flush=True)

# Build quadratic obstruction forms using ONLY plain Python modular
# arithmetic. Sage/FFLAS/OpenBLAS SIGILLs on these tiny matrix products
# on this machine.

def arr(M):
    return [[int(M[i,j]) % MOD for j in range(M.ncols())]
            for i in range(M.nrows())]

def tr(A):
    return [list(row) for row in zip(*A)]

def mm_py(A,B):
    if not A:
        return []
    return [
        [
            sum(A[i][k]*B[k][j] for k in range(len(B))) % MOD
            for j in range(len(B[0]))
        ]
        for i in range(len(A))
    ]

def smul(c,A):
    c%=MOD
    return [[c*x % MOD for x in row] for row in A]

def sub_py(A,B):
    return [[(A[i][j]-B[i][j])%MOD
             for j in range(len(A[0]))]
            for i in range(len(A))]

def add_scaled(A,c,B):
    c%=MOD
    return [[(A[i][j]+c*B[i][j])%MOD
             for j in range(len(A[0]))]
            for i in range(len(A))]

def zero_py(m,n):
    return [[0]*n for _ in range(m)]

# Amap is 16x3 and Bmap 16x2.
Aa=arr(Amap)
Ba=arr(Bmap)
Ymap=[Aa[i]+Ba[i] for i in range(16)]

# Cmap is 11x3; append two zero columns.
Ca=arr(Cmap)
C5=[row+[0,0] for row in Ca]

HFa=[arr(H) for H in HF]
HGa=[arr(H) for H in HG]

inv2i=pow(2,-1,MOD)

# qF_i = 1/2 C5^T H_i C5
QF=[]
for i,H in enumerate(HFa):
    Q=mm_py(mm_py(tr(C5),H),C5)
    QF.append(smul(inv2i,Q))
    print(f"E6P3ORD2|QF={i+1}/8",flush=True)

# qG_i = 1/2 Ymap^T H_i Ymap
QG=[]
for i,H in enumerate(HGa):
    Q=mm_py(mm_py(tr(Ymap),H),Ymap)
    QG.append(smul(inv2i,Q))
    print(f"E6P3ORD2|QG={i+1}/8",flush=True)

# Effective second-order surface correction:
# qEff = qG - Gs*Js^-1*qF.
Tmat_py=mm_py(arr(Gs),arr(JsInv))

QEff=[]
for gi in range(8):
    Q=[row[:] for row in QG[gi]]
    for fi in range(8):
        c=Tmat_py[gi][fi]
        if c:
            Q=add_scaled(Q,-c,QF[fi])
    QEff.append(Q)

left_py=[[int(x)%MOD for x in lv] for lv in left_basis]

QObs=[]
for lv in left_py:
    Q=zero_py(5,5)
    for gi,c in enumerate(lv):
        if c:
            Q=add_scaled(Q,c,QEff[gi])
    QObs.append(Q)

nonzero=[
    i for i,Q in enumerate(QObs)
    if any(x for row in Q for x in row)
]

print(
    f"E6P3ORD2|obstruction_forms={len(QObs)}|nonzero={nonzero}",
    flush=True
)

for i,Q in enumerate(QObs):
    if i in nonzero:
        print(f"E6P3ORD2_Q|i={i}|matrix={Q}",flush=True)

def qeval(Q,y):
    total=0
    for i in range(5):
        if not y[i]:
            continue
        for j in range(5):
            if y[j]:
                total += Q[i][j]*y[i]*y[j]
    return total % MOD

survive_dirs=0
total_corrections=0
examples=[]
count=0
start=time.time()

for dr in range(MOD):
    for ds0 in range(MOD):
        for dx in range(MOD):
            dp=[dr,ds0,dx]
            found=0
            for z0 in range(MOD):
                for z1 in range(MOD):
                    y=dp+[z0,z1]
                    if all(qeval(Q,y)==0 for Q in QObs):
                        found+=1
            if found:
                survive_dirs+=1
                total_corrections+=found
                if len(examples)<args.max_examples:
                    examples.append((dr,ds0,dx,found))
            count+=1
            if args.progress and count%args.progress==0:
                print(f"E6P3ORD2|progress={count}/{MOD**3}|survive={survive_dirs}|corrections={total_corrections}|seconds={time.time()-start:.1f}",flush=True)

print(f"E6P3ORD2|done|directions={MOD**3}|survive={survive_dirs}|corrections={total_corrections}|seconds={time.time()-start:.1f}",flush=True)
for dr,ds0,dx,n in examples:
    print(f"E6P3ORD2_HIT|dr={dr}|ds0={ds0}|dx1={dx}|p3_corrections={n}",flush=True)
