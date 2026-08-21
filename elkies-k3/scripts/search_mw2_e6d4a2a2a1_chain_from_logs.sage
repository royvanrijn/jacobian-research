from sage.all import *
from pathlib import Path
import argparse
import ast


ap=argparse.ArgumentParser(description="Validate the complete MW2 P1 chart and search P2.")
ap.add_argument("--p",type=int,default=23)
ap.add_argument("--dir",default="artifacts/local/elkies-k3/mw2-e6d4a2a2a1-p1")
ap.add_argument("--show-valid",action="store_true")
ap.add_argument("--progress-every",type=int,default=250)
ap.add_argument("--skip-p2",action="store_true")
ap.add_argument("--max-valid",type=int,default=0)
args=ap.parse_args()

p=args.p
K=GF(p)
Kt=PolynomialRing(K,"t")
t=Kt.gen()
root=Path(args.dir)
artifact=root/f"p{p}-p1.ms"
meta_path=artifact.with_suffix(".meta.txt")
all_names=["lam","mu","s","w","x1","x2","c1","c2","a2","a3","a4","a5"]
R=PolynomialRing(K,all_names,order="degrevlex")
d=R.gens_dict()
RF=FractionField(R)
derived=[]
inside=False
for line in meta_path.read_text().splitlines():
    if line=="DERIVED":
        inside=True
        continue
    if inside and " <- " in line:
        name,expression=line.split(" <- ",1)
        derived.append((name,RF(expression)))


def valuation_at(poly,point):
    factor=t-point
    value=0
    while poly and poly(point)==0:
        poly//=factor
        value+=1
    return value


def multiplicative_steps_polynomial(A,B,X,Y,fiber_point,node):
    P=PolynomialRing(K,("u","xx","yy"))
    u,xx,yy=P.gens()
    At=Kt(A(t+fiber_point)); Bt=Kt(B(t+fiber_point))
    Ap=sum(P(c)*u**i for i,c in enumerate(At.list()))
    Bp=sum(P(c)*u**i for i,c in enumerate(Bt.list()))
    surface=yy**2-(node+xx)**3-Ap*(node+xx)-Bp
    Xt=Kt(X(t+fiber_point)); Yt=Kt(Y(t+fiber_point))
    if Xt(0)!=node or Yt(0)!=0:
        return 0
    sx=Kt((Xt-node)//t); sy=Kt(Yt//t)
    surface=P(surface(u,u*xx,u*yy)//u**2)
    steps=1
    while True:
        cx,cy=K(sx(0)),K(sy(0))
        point={u:K(0),xx:cx,yy:cy}
        if any(surface.derivative(v).subs(point) for v in (u,xx,yy)):
            return steps
        sx=Kt((sx-cx)//t); sy=Kt((sy-cy)//t)
        surface=P(surface(u,cx+u*xx,cy+u*yy)//u**2)
        steps+=1
        if steps>5:
            return 99


def reconstruct(candidate):
    values={d[name]:K(value) for name,value in candidate.items()}
    try:
        for name,expression in derived:
            values[d[name]]=K(expression.subs(values))
    except (ZeroDivisionError,ValueError):
        return None
    if any(d[name] not in values for name in all_names):
        return None
    lam=values[d["lam"]]; mu=values[d["mu"]]; s=values[d["s"]]
    X=Kt([1,values[d["x1"]],values[d["x2"]]])
    c1=values[d["c1"]]; c2=values[d["c2"]]
    Y=(t-lam)*(t-mu)*(K(1)/(lam*mu)+c1*t+c2*t**2)
    A=sum(values[d[f"a{i}"]]*t**i for i in range(2,6))
    B=Y**2-X**3-A*X
    Delta=-16*(4*A**3+27*B**2)
    if A.degree()!=5 or B.degree()!=8 or Delta.degree()!=16:
        return None
    expected=[
        (K(0),6),(K(1),3),(lam,3),(mu,2),
    ]
    if any(valuation_at(Delta,point)!=order for point,order in expected):
        return None
    residual=Delta
    for point,order in expected:
        residual//=((t-point)**order)
    if residual.degree()!=2 or gcd(residual,residual.derivative()).degree():
        return None
    if any(residual(point)==0 for point,_ in expected):
        return None
    if X(1)==s and Y(1)==0:
        return None
    if multiplicative_steps_polynomial(A,B,X,Y,lam,X(lam))!=1:
        return None
    if multiplicative_steps_polynomial(A,B,X,Y,mu,X(mu))!=1:
        return None
    # IV*: ord_u(Abar,Bbar,Delta)=(3,4,8).  A polynomial section with
    # (ord_u xbar,ord_u ybar)=(2,2) meets a nonzero E6 component.
    if X.degree()>2 or Y.degree()!=4 or Y[4]==0 or B[8]!=Y[4]**2:
        return None
    return lam,mu,s,A,B,X,Y,Delta


def fast_square_roots(H,degree):
    if H.degree()!=2*degree:
        return []
    roots=[]
    for lead in K(H[2*degree]).sqrt(all=True):
        coeffs=[K(0)]*(degree+1); coeffs[degree]=lead
        for k in range(2*degree-1,degree-1,-1):
            j=k-degree
            known=K(0)
            for i in range(degree+1):
                h=k-i
                if 0<=h<=degree and i!=degree and h!=degree:
                    known+=coeffs[i]*coeffs[h]
            coeffs[j]=(H[k]-known)/(2*lead)
        M=Kt(coeffs)
        if M**2==H:
            roots.append(M)
    return roots


def local_poly(N,power,r,point,precision=7):
    S=PowerSeriesRing(K,"u",default_prec=precision)
    u=S.gen()
    numerator=sum(S(N[i])*(S(point)+u)**i for i in range(N.degree()+1))
    denominator=(S(point)+u-r)**power
    expansion=numerator/denominator
    return Kt([expansion[i] for i in range(precision)])


def rational_steps(A,B,N,M,r,point,node):
    localA=Kt(A(t+point)); localB=Kt(B(t+point))
    localX=local_poly(N,2,r,point); localY=local_poly(M,3,r,point)
    return multiplicative_steps_polynomial(localA,localB,localX,localY,K(0),node)


def find_p2(model,stats):
    lam,mu,s,A,B,P1X,P1Y,Delta=model
    z=polygen(K,"z")
    d4_roots=[
        K(c) for c,multiplicity in (z**3+A[2]*z+B[3]).roots()
        if 3*K(c)**2+A[2]!=0
    ]
    for r in K:
        if r in (K(0),K(1),lam,mu) or Delta(r)==0:
            continue
        q=t-r
        for c in d4_roots:
            n1=c*r**2
            for n2 in K:
                for n3 in K:
                    n4=s*(1-r)**2-n1-n2-n3
                    N=Kt([0,n1,n2,n3,n4])
                    if N(r)==0:
                        continue
                    H=N**3+A*N*q**4+B*q**6
                    for M in fast_square_roots(H,7):
                        stats["square_roots"]+=1
                        orientation="opposite" if M[7]==-P1Y[4] else "same"
                        stats[f"orientation_{orientation}"]+=1
                        if M(r)==0:
                            continue
                        if M[0]!=0 or M[1]!=0 or M(1)!=0:
                            continue
                        stats[f"node_{orientation}"]+=1
                        if rational_steps(A,B,N,M,r,K(1),s)!=1:
                            continue
                        # P2 is on the identity components at lam and mu.
                        if N(lam)==P1X(lam)*(lam-r)**2 and M(lam)==0:
                            continue
                        if N(mu)==P1X(mu)*(mu-r)**2 and M(mu)==0:
                            continue
                        pair=gcd(P1X*q**2-N,P1Y*q**3-M).degree()
                        stats[f"pair_{orientation}_{pair}"]+=1
                        if orientation!="opposite" or pair!=2:
                            continue
                        yield r,N,M


def candidates_from_logs():
    for log in sorted(root.glob("p1-seed*.scan.log")):
        meta=log.with_suffix("").with_suffix(".meta.txt")
        fixed={}
        for line in meta.read_text().splitlines():
            if line.startswith("values="):
                fixed=ast.literal_eval(line.split("=",1)[1])
        for line in log.read_text().splitlines():
            if not line.startswith("MW3A10SCAN_HIT|"):
                continue
            hit={}
            for item in line.split("|",1)[1].split(","):
                name,value=item.split("=",1); hit[name]=int(value)
            yield log.stem,{**fixed,**hit}


raw=valid=targets=0
seen=set()
from collections import defaultdict
stats=defaultdict(int)
for source,candidate in candidates_from_logs():
    raw+=1
    model=reconstruct(candidate)
    if model is None:
        continue
    key=tuple(tuple(poly.list()) if hasattr(poly,"list") else poly for poly in model[:-1])
    if key in seen:
        continue
    seen.add(key); valid+=1
    if args.show_valid or (args.progress_every and valid%args.progress_every==0):
        print(f"MW2CHAIN|valid_p1={valid}|source={source}",flush=True)
    for r,N,M in (() if args.skip_p2 else find_p2(model,stats)):
        targets+=1
        lam,mu,s,A,B,P1X,P1Y,Delta=model
        print(
            "MW2CHAIN_TARGET"
            +f"|lam={int(lam)}|mu={int(mu)}|s={int(s)}|pole={int(r)}"
            +"|A="+",".join(map(str,map(int,A.list())))
            +"|B="+",".join(map(str,map(int,B.list())))
            +"|P1X="+",".join(map(str,map(int,P1X.list())))
            +"|P1Y="+",".join(map(str,map(int,P1Y.list())))
            +"|N="+",".join(map(str,map(int,N.list())))
            +"|M="+",".join(map(str,map(int,M.list()))),flush=True,
        )
    if args.max_valid and valid>=args.max_valid:
        break

print(f"MW2CHAINSUMMARY|raw={raw}|valid_p1={valid}|targets={targets}",flush=True)
if not args.skip_p2:
    print(
        "MW2CHAINSTATS|"+"|".join(f"{key}={stats[key]}" for key in sorted(stats)),
        flush=True,
    )
