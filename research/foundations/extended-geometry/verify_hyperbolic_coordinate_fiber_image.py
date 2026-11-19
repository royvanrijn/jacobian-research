#!/usr/bin/env python3
"""Exact finite controls for the hyperbolic-coordinate polynomial fiber."""
if not __debug__:raise RuntimeError("Assertions must be enabled")
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import argparse,json,hashlib,platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PROOF=HERE/"HYPERBOLIC_COORDINATE_FIBER_IMAGE_WITNESSES.md"
LEDGER=HERE/"hyperbolic_coordinate_fiber_image_source_ledger.json"

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def trim(a):
    a=list(map(Q,a))
    while len(a)>1 and not a[-1]:a.pop()
    return a

def add(a,b,scale=1):
    return trim([(a[j]if j<len(a)else Q(0))+scale*(b[j]if j<len(b)else Q(0))
                 for j in range(max(len(a),len(b)))])

def multiply(a,b,order=None):
    length=len(a)+len(b)-1
    if order is not None:length=min(length,order+1)
    out=[Q(0)]*length
    for i,v in enumerate(a):
        for j,w in enumerate(b):
            if i+j<length:out[i+j]+=v*w
    return trim(out)

def power(a,n,order=None):
    out=[Q(1)]
    for _ in range(n):out=multiply(out,a,order)
    return out

def scale(a,v):return trim([v*c for c in a])

def pochhammer(a,r):
    out=Q(1)
    for j in range(r):out*=a+j
    return out

def kappa(m,a=1):
    if m%2:return Q(0)
    return sum((Q((-1)**j*comb(m,j),a*m+2*j+1)for j in range(m+1)),Q(0))

def mu(n,lam,r):
    return (-n*n*lam*lam)**r*pochhammer(Q(1,2*n),r)*pochhammer(1-Q(1,2*n),r)/(lam*pochhammer(Q(1,2),r)*factorial(r))

def badd(a,b,scale=1):
    out=dict(a)
    for p,v in b.items():out[p]=out.get(p,0)+scale*v
    return {p:v for p,v in out.items()if v}

def bscale(a,c):return {p:c*v for p,v in a.items()if c*v}

def bmul(a,b):
    out={}
    for(i,j),v in a.items():
        for(k,l),w in b.items():
            p=i+k,j+l;out[p]=out.get(p,0)+v*w
    return {p:v for p,v in out.items()if v}

def bpower(a,n):
    out={(0,0):1}
    for _ in range(n):out=bmul(out,a)
    return out

def beval(a,degree,xi1,xi2,z1,z2):
    return sum(v*xi1**i*xi2**(degree-i)*z1**j*z2**(degree-j)
               for(i,j),v in a.items())

def raw(C,D,G,M,a=1):
    terms=[(j,0,v)for j,v in enumerate(C)if v]+[(j,2,v)for j,v in enumerate(D)if v]
    P={(0,0):Q(1)};out=[]
    for m in range(1,M+1):
        temp={}
        for(i,t),v in P.items():
            for j,s,w in terms:
                if i+j<=M:temp[i+j,t+s]=temp.get((i+j,t+s),Q(0))+v*w
        P={key:v for key,v in temp.items()if v}
        pure=sum((v/Q(t+a*m+1)for(i,t),v in P.items()if i==m and not(t+a*m)%2),Q(0))
        marked=sum((v*coef(G,m-i)/Q(t+a*m+1)for(i,t),v in P.items()if i<=m and not(t+a*m)%2),Q(0))
        out.append((pure,marked))
    return out

R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}

def coef(a,j):return a[j]if j<len(a)else Q(0)

def evaluate(a,x):return sum((v*x**j for j,v in enumerate(a)),Q(0))
def divide(a,b):
    a=trim(a);b=trim(b);q=[Q(0)]*max(1,len(a)-len(b)+1)
    while a!=[0]and len(a)>=len(b):
        j=len(a)-len(b);v=a[-1]/b[-1];q[j]=v;a=add(a,[0]*j+scale(b,v),-1)
    return trim(q),a

def construct(N,a,ell,lam):
    assert N>=1 and a>=1 and a%2 and ell>=1 and ell%2 and lam
    U=[Q(1),Q(lam)];S=[Q(comb(2*N,j+1))*lam**j for j in range(2*N)]
    X0=scale(multiply(U,S),Q(1,2*N));C1=power(U,N+1)
    Qx=[Q(0)]
    for j in range(ell):Qx=add(Qx,power(U,2*N*j),(-1)**j)
    C=multiply(multiply(C1,power(X0,a-1)),power(Qx,a))
    D=scale(multiply(multiply(C1,power(X0,a+1)),power(Qx,a+2)),-1)
    d=(2*ell*(a+2)-1)*N+a+2;q=2*N*(ell-1)
    assert C[0]==1 and D[0]==-1 and Qx[0]==1
    assert len(Qx)-1==q and len(D)-1==d-a-1
    return C,D,Qx,X0,d,q

def original(C,D,N,a):
    e=len(D)-1;inner=bscale(bmul(W,bpower(R,e)),2)
    for j,v in enumerate(C):
        if j:inner=badd(inner,bmul(bpower(R,e+2-j),bpower(Z,j-1)),v)
    for j,v in enumerate(D):
        if j:inner=badd(inner,bmul(bpower(T,2),bmul(bpower(R,e-j),bpower(Z,j-1))),v)
    divided=bscale(inner,(2*N)**(a+1));F=bmul(bpower(T,a),divided)
    assert all(v.denominator==1 for v in F.values())
    return F,divided

def original_control(N,a,ell,lam,maximum):
    C,D,Qx,X0,d,q=construct(N,a,ell,Q(lam));F,divided=original(C,D,N,a)
    Gh={}
    for j,v in enumerate(Qx):Gh=badd(Gh,bmul(bpower(R,q-j),bpower(Z,j)),v)
    Gh=bmul(Z,Gh)
    radial=beval(F,d,1,1,1,-1)
    assert radial==(-1)**(N+1)*2**(a+2)*lam**((2*ell*(a+2)-1)*N+1)
    r=next(r for r in range(1,10)if evaluate(C,Q(r,2)))
    expected_T=(2*N)**(a+1)*(2*r)**(len(D)+1)*evaluate(C,Q(r,2))/Q(r*r)
    assert beval(divided,d-a,r,1,1,r)==expected_T and expected_T
    # Compare the ORIGINAL specialization with the literal divided-difference formula.
    specialized=[sum(v for(i,j),v in F.items()if i==d and d-j==k)for k in range(d+1)]
    YQ=scale(multiply(X0,Qx),2*N)
    square=power(YQ,2);assert square[0]==4*N*N
    expected=scale(multiply(multiply(multiply(power([1,lam],N+1),power(YQ,a-1)),Qx),square[1:]),-1)
    assert trim(specialized)==expected
    P={(0,0):Q(1)};rows=[]
    for m in range(1,maximum+1):
        P=bmul(P,F);md=d*m
        pure=sum(v*factorial(i)*factorial(md-i)for(i,j),v in P.items()if i==j)
        mixedP=bmul(P,Gh);gd=md+q+1
        mixed=sum(v*factorial(i)*factorial(gd-i)for(i,j),v in mixedP.items()if i==j)
        expected=Q(0)if m%2 else(2*N)**((a+1)*m)*kappa(m,a)*mu(N,Q(lam),m//2)
        assert pure==0 and Q(mixed,factorial(gd+1))==expected
        rows.append({'m':m,'support':len(P),'pure':0,'marked_factorial':gd+1,'normalized_marked':str(expected)})
    return {'N':N,'a':a,'ell':ell,'lambda':lam,'degree':d,'marker_degree':q+1,'phase_support':len(F),'radial_value':str(radial),'T_quotient_value':str(expected_T),'original_controls':rows}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path,help='Write exact finite-control receipt.')
    args=parser.parse_args()
    ledger=json.loads(LEDGER.read_text())
    pins=[PROOF,Path(__file__).resolve(),LEDGER]
    for item in ledger['retained_resources']:
        path=HERE/item['retained_path']
        assert digest(path)==item['sha256'],('Retained source changed',str(path))
        pins.append(path)
    cases=[];M=12
    for N in range(1,5):
        for a in[1,3,5]:
            for ell in[1,3,5]:
                for lam in[Q(1),Q(-2),Q(1,2)]:
                    C,D,Qx,X,d,q=construct(N,a,ell,lam);rows=raw(C,D,[0]+Qx,M,a)
                    for m,(pure,mixed)in enumerate(rows,1):
                        assert pure==0 and mixed==(Q(0)if m%2 else kappa(m,a)*mu(N,lam,m//2))
                    cases.append({'N':N,'a':a,'ell':ell,'lambda':str(lam),'degree':d,'marker_degree':q+1,'moments':M})
    originals=[original_control(*case)for case in[(1,1,1,1,4),(2,1,1,1,4),(1,1,3,1,4),(2,1,3,1,2),(2,1,3,-2,2),(2,1,5,1,2),(1,3,3,1,2),(2,3,3,1,1),(2,5,3,1,1)]]
    # Exact Laurent-divisibility and negative-exponent branch controls.
    fiber=[]
    for N in range(1,7):
        den=[1]+[0]*(2*N-1)+[1]
        for k in range(-18*N,18*N+1):
            num=[2]if k==0 else[1]+[0]*(abs(k)-1)+[1]
            quotient,rem=divide(num,den)
            allowed=bool(k)and abs(k)%(2*N)==0 and(abs(k)//(2*N))%2==1
            assert(rem==[0])==allowed
            if allowed:
                ell=abs(k)//(2*N)
                for a in[1,3,5]:
                    order=N+a+(0 if k>0 else a*k)
                    assert(k>0 or order<0 or(N,a,ell)==(1,1,1))
                fiber.append({'N':N,'k':k,'ell':ell,'positive':k>0})
    # M(U)=1+U^k is precisely compatible with the involution and unit identity.
    for k in[-12,-4,4,12]:
        for U in[Q(2),Q(3),Q(5,2)]:
            M1=1+U**k;M2=1+U**(-k)
            assert 1/M1+1/M2==1 and(1-M1)*(1-M2)==1
    C,D,Qx,X,d,q=construct(2,1,3,Q(1));badQ=Qx.copy();badQ[1]+=1
    Cbad=multiply(power([1,1],3),badQ);Dbad=scale(multiply(multiply(power([1,1],3),power(X,2)),power(badQ,3)),-1)
    assert any(pure for pure,mixed in raw(Cbad,Dbad,[1],M))
    # Missing Q in the fixed marker does not reproduce the old every-even coefficient.
    assert raw(C,D,[0,1],2)[1][1]!=kappa(2)*mu(2,Q(1),1)
    out={'passed':True,'scope':'Finite raw polynomial, exact Laurent-divisibility/unit and original four-variable controls. All-order purity, fixed-marker nonvanishing and exact fiber completeness use the separate written proof.',
         'raw_cases':cases,'original_cases':originals,'laurent_divisibility':fiber,
         'negative_controls':{'altered_Q_rejected':True,'omitted_marker_Q_rejected':True},
         'status':'passed','schema_version':1,'python':platform.python_version(),
         'sources':{str(path.relative_to(ROOT)):digest(path)for path in pins}}
    if args.write_receipt:args.write_receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'raw_cases':len(cases),'raw_moments':M*len(cases),'original_cases':len(originals),'original_powers':sum(len(c['original_controls'])for c in originals),'F37':originals[3]}))
if __name__=='__main__':main()
