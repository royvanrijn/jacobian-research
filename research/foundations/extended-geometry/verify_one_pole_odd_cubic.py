#!/usr/bin/env python3
"""Exact finite controls for the one-pole odd-cubic purity classification."""
if not __debug__:raise RuntimeError('Assertions must be enabled')
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import argparse,json,hashlib,platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PROOF=HERE/"ONE_POLE_ODD_CUBIC_CLASSIFICATION.md"
LEDGER=HERE/"one_pole_odd_cubic_source_ledger.json"

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def trim(a):
    a=list(map(Q,a))
    while len(a)>1 and not a[-1]:a.pop()
    return a

def add(a,b,c=1):return trim([(a[j]if j<len(a)else 0)+c*(b[j]if j<len(b)else 0)for j in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*v for v in a])
def mul(a,b,M=None):
    out=[Q(0)]*(len(a)+len(b)-1 if M is None else min(M+1,len(a)+len(b)-1))
    for i,v in enumerate(a):
        for j,w in enumerate(b):
            if i+j<len(out):out[i+j]+=v*w
    return trim(out)
def power(a,n,M=None):
    out=[Q(1)]
    for _ in range(n):out=mul(out,a,M)
    return out

def inverse(a,M):
    assert a[0]
    b=[1/Q(a[0])]
    for n in range(1,M+1):b.append(-sum((a[j]*b[n-j]for j in range(1,min(n,len(a)-1)+1)),Q(0))/a[0])
    return trim(b)

def compose(a,b,M=None):
    out=[Q(0)]
    for v in reversed(a):out=add(mul(out,b,M),[v])
    return out

def coef(a,j):return a[j]if j<len(a)else Q(0)
def evaluate(a,t):return sum((v*t**j for j,v in enumerate(a)),Q(0))
def divide_at1(a):
    assert sum(a)==0
    return trim([-sum(a[:j+1])for j in range(len(a)-1)])

def revert(a,M):
    assert a[0]==0 and a[1]
    b=[Q(0)]*(M+1)
    for j in range(1,M+1):b[j]=(int(j==1)-coef(compose(a,b,j),j))/a[1]
    assert compose(a,b,M)==[Q(0),Q(1)]
    return trim(b)

def construct(N,lam,d=None):
    assert N>=1 and lam
    if d is None:P=divide_at1([-1]+[0]*(2*N-1)+[1]);branch='inversion';degree=5*N+3
    else:
        assert N%d==0
        q=N//d
        twoW=[-1]+[0]*(d-1)+[2]
        P=divide_at1(add(power(twoW,q),[-1]));branch='fixed_zero';degree=3*N+3
    D_U=[Q(0)]*(N+3)+scale(power(P,2),-Q(1,4*N*N))
    U=[Q(1),Q(lam)]
    C=power(U,N+1);D=compose(D_U,U);Px=compose(P,U)
    assert C[0]==1 and D[0]==-1 and len(D)-1==degree-2
    assert evaluate(P,Q(1))==2*N
    return C,D,Px,degree

def kappa(m):
    return Q(0)if m%2 else sum((Q((-1)**j*comb(m,j),m+2*j+1)for j in range(m+1)),Q(0))

def raw(C,D,G,M):
    terms=[(j,0,v)for j,v in enumerate(C)if v]+[(j,2,v)for j,v in enumerate(D)if v]
    P={(0,0):Q(1)};out=[]
    for m in range(1,M+1):
        temp={}
        for(i,t),v in P.items():
            for j,s,w in terms:
                if i+j<=M:temp[i+j,t+s]=temp.get((i+j,t+s),Q(0))+v*w
        P={key:v for key,v in temp.items()if v}
        pure=sum((v/Q(t+m+1)for(i,t),v in P.items()if i==m and not(t+m)%2),Q(0))
        marked=sum((v*coef(G,m-i)/Q(t+m+1)for(i,t),v in P.items()if i<=m and not(t+m)%2),Q(0))
        out.append((pure,marked))
    return out

def badd(a,b,c=1):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,0)+c*v
    return {k:v for k,v in out.items()if v}
def bmul(a,b):
    out={}
    for(i,j),v in a.items():
        for(k,l),w in b.items():out[i+k,j+l]=out.get((i+k,j+l),0)+v*w
    return {k:v for k,v in out.items()if v}
def bpow(a,n):
    out={(0,0):Q(1)}
    for _ in range(n):out=bmul(out,a)
    return out

def beval(a,d,coords):
    u,v,w,z=coords
    return sum(c*u**i*v**(d-i)*w**j*z**(d-j)for(i,j),c in a.items())
R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}

def original(C,D,b):
    inner={k:2*v for k,v in bmul(W,bpow(R,b-2)).items()}
    for j,v in enumerate(C):
        if j:inner=badd(inner,bmul(bpow(R,b-j),bpow(Z,j-1)),v)
    for j,v in enumerate(D):
        if j:inner=badd(inner,bmul(bpow(T,2),bmul(bpow(R,b-2-j),bpow(Z,j-1))),v)
    return bmul(T,inner),inner

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path,help='Write exact finite-control receipt.')
    args=parser.parse_args()
    ledger=json.loads(LEDGER.read_text())
    pins=[PROOF,Path(__file__).resolve(),LEDGER]
    for resource in ledger['retained_resources']:
        source=HERE/resource['retained_path']
        assert digest(source)==resource['sha256'],('Retained source changed',str(source))
        pins.append(source)
    M=12;cases=[]
    for N in range(1,9):
        for d in [None]+[j for j in range(1,N+1)if N%j==0]:
            for lam in [Q(1),Q(-2),Q(1,2)]:
                C,D,P,b=construct(N,lam,d);U=[Q(1),lam];G=power(U,N)
                A=scale(mul(power(U,N),inverse(P,M),M),2*N);Cinv=inverse(C,M)
                rows=raw(C,D,G,M)
                residue=[]
                for m,(pure,mixed)in enumerate(rows,1):
                    common=mul(power(A,m+1,M),Cinv,M)
                    assert pure==kappa(m)*coef(common,m)==0
                    expect=kappa(m)*coef(mul(G,common,M),m)
                    assert mixed==expect
                    residue.append(str(mixed))
                if d is not None:
                    q=N//d
                    us=scale(add(power([1,-1],q),power([1,1],q),-1),Q(1,2*N)/lam)
                    assert all(not v for j,v in enumerate(us)if j%2==0)
                    den=add(power([1,-1],q-1),power([1,1],q-1))
                    Es=scale(inverse(mul([1,0,-1],den),M),2)
                    E=compose(Es,revert(us,M),M)
                    for m,(pure,mixed)in enumerate(rows,1):assert mixed==kappa(m)*coef(E,m)
                    assert evaluate(den,1)==evaluate(den,-1)==(2 if q==1 else 2**(q-1))
                    # Actual rational involution relations, at independent exact points.
                    for w in [Q(2),Q(3),Q(5,2)]:
                        sig=w/(2*w-1)
                        uw=lambda z:((2*z-1)**q-1)/(2*N*lam*z**q)
                        V=lambda z:(1-z**(-q))/(N*lam)
                        assert uw(sig)==-uw(w)and V(sig)==V(w)-2*uw(w)
                else:
                    for u in [Q(2),Q(3),Q(5,2)]:
                        uc=lambda z:(z**N-z**(-N))/(2*N*lam)
                        V=lambda z:(1-z**(-N))/(N*lam)
                        assert uc(1/u)==-uc(u)and V(1/u)==V(u)-2*uc(u)
                cases.append({'N':N,'d':d,'lambda':str(lam),'least_degree':b,'marker':'U^N','marked_controls':residue})
    originals=[]
    for N,d,lam,maxm in [(1,1,1,4),(2,1,1,6),(2,2,1,4),(3,1,1,4),(3,3,1,4),(4,2,1,4),(4,4,-2,4),(2,None,1,4)]:
        C,D,P,b=construct(N,Q(lam),d);F,divided=original(C,D,b);G=bpow(badd(R,Z,lam),N)
        assert beval(F,b,(1,1,1,-1))==8*D[-1]*(-1)**(b-3)and beval(F,b,(1,1,1,-1))
        r=next(j for j in range(1,8)if 2+lam*j)
        expected=Q((2*r)**(b-N-1)*(2*r+lam*r*r)**(N+1),r*r)
        assert beval(divided,b-1,(r,1,1,r))==expected and expected
        powerF={(0,0):Q(1)};rawrows=raw(C,D,power([Q(1),Q(lam)],N),maxm);controls=[]
        for m in range(1,maxm+1):
            powerF=bmul(powerF,F);degree=b*m
            pure=sum(v*factorial(i)*factorial(degree-i)for(i,j),v in powerF.items()if i==j)
            markedP=bmul(powerF,G);markedDegree=degree+N
            mixed=sum(v*factorial(i)*factorial(markedDegree-i)for(i,j),v in markedP.items()if i==j)
            assert pure==0 and Q(mixed,factorial(markedDegree+1))==rawrows[m-1][1]
            controls.append({'m':m,'support':len(powerF),'pure':0,'normalized_marked':str(rawrows[m-1][1])})
        originals.append({'N':N,'d':d,'lambda':lam,'degree':b,'phase_support':len(F),'controls':controls})
    # The old degree-9 point has zero x marker, but the uniform U^N marker is unsafe.
    C,D,P,b=construct(2,Q(1),1)
    assert D==scale(power([1,1],7),-1)
    assert all(mixed==0 for pure,mixed in raw(C,D,[0,1],M))
    assert raw(C,D,power([1,1],2),2)[1][1]
    altered=D.copy();altered[1]+=1
    assert any(pure for pure,mixed in raw(C,altered,[1],4))
    result={'passed':True,'scope':'Finite exact controls only. The written Galois/involution, norm, pole and endpoint arguments prove completeness and all orders.',
            'raw_cases':cases,'original_cases':originals,
            'negative_controls':{'altered_D_rejected':True,'uniform_x_marker_extrapolation_rejected':True},
            'status':'passed','schema_version':1,'python':platform.python_version(),
            'sources':{str(path.relative_to(ROOT)):digest(path)for path in pins}}
    if args.write_receipt:args.write_receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'passed':True,'raw_cases':len(cases),'raw_moments':M*len(cases),'original_cases':len(originals),'original_powers':sum(len(c['controls'])for c in originals)}))
if __name__=='__main__':main()
