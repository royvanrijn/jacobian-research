from sage.all import *
import argparse

ap=argparse.ArgumentParser(description="Scan first-order 31-adic parameter directions for P3 liftability.")
ap.add_argument("--progress",type=int,default=1000)
args=ap.parse_args()

K=GF(31)

# Surface/P1 solved variables
surf_names=["a1","a2","a4","lam","mu","s1","sl","sm"]
param_names=["r0","s0","x1"]

# P3 coeffs
p3_names=["u1","u2","v1","v2","v3"]

names=surf_names+param_names+p3_names
R=PolynomialRing(K,names,order="degrevlex")
d=R.gens_dict()

a1,a2,a4,lam,mu,s1,sl,sm=[d[n] for n in surf_names]
r0,s0,x1=[d[n] for n in param_names]
u1,u2,v1,v2,v3=[d[n] for n in p3_names]

# We build only the already-reduced P1/fiber system in one ring.
# Reconstruct a3,b4,b5,y2 symbolically via the same triangular formulas.
# Use auxiliary polynomial ring but collapse to R before differentiation.
RF=FractionField(R)
U=PolynomialRing(RF,["a3","b4","b5","y2"])
a3,b4,b5,y2=U.gens()
FU=FractionField(U)
Ut=PolynomialRing(FU,"t")
t=Ut.gen()

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
    out=FU(-e.subs({vv:0})/c)
    subs[vv]=out
    return out

solve_linear(Sec[7],b4)
solve_linear(Sec[6],b5)
solve_linear(Sec[5],y2)
solve_linear(Sec[4],a3)

def subpoly(P):
    return sum(settle(P[i])*t**i for i in range(P.degree()+1))

A=subpoly(A); B=subpoly(B)
Ap=A.derivative(t); Bp=B.derivative(t)
Delta=-16*(4*A**3+27*B**2)

fiber_fr=[
    Delta.derivative(t,2)(1),
    Delta.derivative(t,3)(1),
    A(L(lam))+3*L(sl)**2,
    B(L(lam))-2*L(sl)**3,
    Bp(L(lam))+L(sl)*Ap(L(lam)),
    A(L(mu))+3*L(sm)**2,
    B(L(mu))-2*L(sm)**3,
    Bp(L(mu))+L(sm)*Ap(L(mu)),
]

# Collapse Frac(U) -> Frac(R)
def collapse(e):
    e=settle(e)
    nu=U(e.numerator()); de=U(e.denominator())
    zero=(0,)*U.ngens()
    for mon in nu.dict():
        if tuple(mon)!=zero: raise RuntimeError("aux numerator survived")
    for mon in de.dict():
        if tuple(mon)!=zero: raise RuntimeError("aux denominator survived")
    n=RF(next((c for mon,c in nu.dict().items() if tuple(mon)==zero),0))
    q0=RF(next((c for mon,c in de.dict().items() if tuple(mon)==zero),0))
    return n/q0

fiber=[collapse(f) for f in fiber_fr]

# P3 equations after same surface reconstruction
X3=L(s0)+L(u1)*t+L(u2)*t**2
Y3=L(v1)*t+L(v2)*t**2+L(v3)*t**3+t**4
S3=Y3**2-X3**3-A*X3-B
p3=[collapse(S3[k]) for k in range(8)]

# Base point
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
    return K(num)/K(den)

# Linearized systems as ordinary Python integer arrays.
# Avoid Sage matrix_modn_dense_float / FFLAS/OpenBLAS: tiny matrices do
# not need optimized BLAS, and that backend SIGILLs on this machine.

MOD=31

def mi(x):
    return int(x) % MOD

def rank_mod(A):
    A=[[mi(x) for x in row] for row in A]
    if not A:
        return 0
    m=len(A); n=len(A[0])
    r=0
    for c in range(n):
        piv=next((i for i in range(r,m) if A[i][c]),None)
        if piv is None:
            continue
        A[r],A[piv]=A[piv],A[r]
        z=pow(A[r][c],-1,MOD)
        A[r]=[(x*z)%MOD for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                z=A[i][c]
                A[i]=[(A[i][j]-z*A[r][j])%MOD for j in range(n)]
        r+=1
        if r==m:
            break
    return r

def inv_mod(A):
    n=len(A)
    M=[
        [mi(x) for x in A[i]]
        + [1 if i==j else 0 for j in range(n)]
        for i in range(n)
    ]
    for c in range(n):
        piv=next((i for i in range(c,n) if M[i][c]),None)
        if piv is None:
            raise RuntimeError("singular matrix")
        M[c],M[piv]=M[piv],M[c]
        z=pow(M[c][c],-1,MOD)
        M[c]=[(x*z)%MOD for x in M[c]]
        for i in range(n):
            if i!=c and M[i][c]:
                z=M[i][c]
                M[i]=[(M[i][j]-z*M[c][j])%MOD for j in range(2*n)]
    return [row[n:] for row in M]

def mm(A,B):
    if not A:
        return []
    m=len(A); k=len(B); n=len(B[0])
    return [
        [
            sum(A[i][q]*B[q][j] for q in range(k)) % MOD
            for j in range(n)
        ]
        for i in range(m)
    ]

def msub(A,B):
    return [
        [(A[i][j]-B[i][j])%MOD for j in range(len(A[0]))]
        for i in range(len(A))
    ]

def mv(A,v):
    return [
        sum(A[i][j]*v[j] for j in range(len(v))) % MOD
        for i in range(len(A))
    ]

def transpose(A):
    return [list(x) for x in zip(*A)]

def nullspace_mod(A):
    A=[[mi(x) for x in row] for row in A]
    if not A:
        return []
    m=len(A); n=len(A[0])
    pivots=[]
    r=0
    for c in range(n):
        piv=next((i for i in range(r,m) if A[i][c]),None)
        if piv is None:
            continue
        A[r],A[piv]=A[piv],A[r]
        z=pow(A[r][c],-1,MOD)
        A[r]=[(x*z)%MOD for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                z=A[i][c]
                A[i]=[(A[i][j]-z*A[r][j])%MOD for j in range(n)]
        pivots.append(c)
        r+=1
        if r==m:
            break

    free=[c for c in range(n) if c not in pivots]
    basis=[]
    for f in free:
        v=[0]*n
        v[f]=1
        for row,pc in enumerate(pivots):
            v[pc]=(-A[row][f])%MOD
        basis.append(v)
    return basis

Js=[
    [mi(eval_rf(f.derivative(d[n]))) for n in surf_names]
    for f in fiber
]
Jp=[
    [mi(eval_rf(f.derivative(d[n]))) for n in param_names]
    for f in fiber
]

jsrank=rank_mod(Js)
if jsrank != 8:
    raise RuntimeError(f"surface Jacobian rank {jsrank} != 8")

Js_inv=inv_mod(Js)

G_s=[
    [mi(eval_rf(f.derivative(d[n]))) for n in surf_names]
    for f in p3
]
G_p=[
    [mi(eval_rf(f.derivative(d[n]))) for n in param_names]
    for f in p3
]
G_u=[
    [mi(eval_rf(f.derivative(d[n]))) for n in p3_names]
    for f in p3
]

# ds = -Js^-1 Jp dp
# P3 condition:
# G_u du = -(G_s ds + G_p dp)
#
# effective parameter matrix = G_p - G_s Js^-1 Jp
M_eff=msub(G_p, mm(mm(G_s,Js_inv),Jp))

gu_rank=rank_mod(G_u)

# left kernel of G_u == kernel(transpose(G_u))
left=nullspace_mod(transpose(G_u))

obs=[
    [
        sum(v[i]*M_eff[i][j] for i in range(len(v))) % MOD
        for j in range(3)
    ]
    for v in left
]

obs_rank=rank_mod(obs) if obs else 0

print(
    f"E6P3DIR|Js_rank={jsrank}|Gu_rank={gu_rank}"
    f"|Meff_rank={rank_mod(M_eff)}",
    flush=True,
)
print(
    f"E6P3DIR|left_dim={len(left)}|obs_rank={obs_rank}"
    f"|predicted_dim={3-obs_rank}",
    flush=True,
)
print(f"E6P3DIR|obs={obs}",flush=True)

# Enumerate to verify linear-algebra prediction.
for dr in K:
    for ds0 in K:
        for dx in K:
            dp=[int(dr),int(ds0),int(dx)]
            ok=all(
                sum(row[j]*dp[j] for j in range(3)) % MOD == 0
                for row in obs
            )
            if ok:
                hits+=1
                if len(examples)<20:
                    examples.append((int(dr),int(ds0),int(dx)))
            count+=1
            if args.progress and count%args.progress==0:
                print(f"E6P3DIR|progress={count}/{total}|hits={hits}",flush=True)

print(f"E6P3DIR|done|total={total}|hits={hits}",flush=True)
for e in examples:
    print(f"E6P3DIR_HIT|dr={e[0]}|ds0={e[1]}|dx1={e[2]}",flush=True)
