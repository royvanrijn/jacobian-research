#!/usr/bin/env python3
"""Finite exact controls for the separate mixed-center structural proof."""
if not __debug__:raise RuntimeError('Assertions must be enabled')
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import argparse,json,hashlib,platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PROOF=HERE/"MIXED_CENTER_LOGARITHMIC_SEPARATION.md"
LEDGER=HERE/"mixed_center_logarithmic_separation_source_ledger.json"

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

def original(C,D,b):
    inner={k:2*v for k,v in bmul(W,bpow(R,b-2)).items()}
    for j,v in enumerate(C):
        if j:inner=badd(inner,bmul(bpow(R,b-j),bpow(Z,j-1)),v)
    for j,v in enumerate(D):
        if j:inner=badd(inner,bmul(bpow(T,2),bmul(bpow(R,b-2-j),bpow(Z,j-1))),v)
    return bmul(T,inner),inner

R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}

def divide(a,b):
    a=trim(a);b=trim(b);q=[Q(0)]*max(1,len(a)-len(b)+1)
    while a!=[0]and len(a)>=len(b):
        j=len(a)-len(b);v=a[-1]/b[-1];q[j]=v;a=add(a,[0]*j+scale(b,v),-1)
    return trim(q),a

def sqrt_series(lam,M):
    out=[Q(0)]*(M+1);v=Q(1);out[0]=1
    for r in range(1,M//2+1):v*=Q(3-2*r,2*r)*lam*lam;out[2*r]=v
    return trim(out)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path,help='Write exact bounded-control receipt.')
    args=parser.parse_args()
    ledger=json.loads(LEDGER.read_text())
    pins=[PROOF,Path(__file__).resolve(),LEDGER]
    for item in ledger['retained_resources']:
        path=HERE/item['retained_path']
        assert digest(path)==item['sha256'],('Retained source changed',str(path))
        pins.append(path)
    M=12;cases=[];originals=[]
    for branch in ['fixed_root','inversion']:
        for s in(range(1,5)if branch=='fixed_root'else range(4)):
            q=2*s if branch=='fixed_root'else 4*s
            Qpoly=scale(add([1],power([-1,2],q)),Q(1,2))if branch=='fixed_root'else[Q((-1)**(j//2))if j%2==0 else Q(0)for j in range(q+1)]
            assert Qpoly[0]==evaluate(Qpoly,1)==1 and len(Qpoly)-1==q
            if branch=='fixed_root':
                transformed=[Q(0)]
                for j,v in enumerate(Qpoly):transformed=add(transformed,[0]*j+scale(power([-1,2],q-j),v))
                assert transformed==Qpoly
            else:
                assert Qpoly==Qpoly[::-1]
                assert mul(Qpoly,[1,0,1])==[Q(1)]+[Q(0)]*(q+1)+[Q(1)]
            C_U=[0,0]+Qpoly
            D_U=scale([0]*4+power(Qpoly,3),-1)
            if branch=='inversion':D_U=scale(mul(D_U,[1,2,1]),Q(1,4))
            for lam in[Q(1),Q(-2),Q(1,2)]:
                U=[Q(1),lam];C=compose(C_U,U);D=compose(D_U,U)
                b=3*q+6 if branch=='fixed_root'else 3*q+8
                assert len(D)-1==b-2 and C[0]==1 and D[0]==-1
                A=U if branch=='fixed_root'else scale(mul(U,inverse(add(U,[1]),M),M),2)
                G=C if branch=='fixed_root'else U
                rows=raw(C,D,G,M);inverseC=inverse(C,M)
                if branch=='inversion':
                    Uv=add([0,lam],sqrt_series(lam,M))
                    marker_series=scale(mul(Uv,inverse(add([1],power(Uv,q+2,M)),M),M),2)
                records=[]
                for m,(pure,marked)in enumerate(rows,1):
                    common=mul(power(A,m+1,M),inverseC,M)
                    assert pure==kappa(m)*coef(common,m)==0
                    assert marked==kappa(m)*coef(mul(G,common,M),m)
                    expected=Q(0)if m%2 else(kappa(m)*(m+1)*lam**m if branch=='fixed_root'else kappa(m)*coef(marker_series,m))
                    assert marked==expected
                    if not m%2 and branch=='fixed_root':assert marked
                    records.append(str(marked))
                cases.append({'branch':branch,'s':s,'lambda':str(lam),'degree':b,'pure_controls':M,'marked_controls':records})
                if lam==1 and((branch=='fixed_root'and s<=2)or(branch=='inversion'and s<=1)):
                    maxm=4 if s<=1 else 2
                    F,divided=original(C,D,b);G_h={}
                    for j,v in enumerate(G):G_h=badd(G_h,bmul(bpow(R,len(G)-1-j),bpow(Z,j)),v)
                    P={(0,0):Q(1)};powers=[]
                    for m in range(1,maxm+1):
                        P=bmul(P,F);d=b*m
                        pure=sum(v*factorial(i)*factorial(d-i)for(i,j),v in P.items()if i==j)
                        mixed=bmul(P,G_h);md=d+len(G)-1
                        marked=sum(v*factorial(i)*factorial(md-i)for(i,j),v in mixed.items()if i==j)
                        assert pure==0 and Q(marked,factorial(md+1))==rows[m-1][1]
                        powers.append({'m':m,'support':len(P),'pure':0,'normalized_marked':str(rows[m-1][1])})
                    originals.append({'branch':branch,'s':s,'degree':b,'phase_support':len(F),'powers':powers})
    divisibility=[]
    for q in range(25):
        quotient,remainder=divide([1]+[0]*(q+1)+[1],[1,0,1])
        assert(remainder==[0])==(q%4==0)
        divisibility.append({'q':q,'divisible':remainder==[0]})
    for t in[Q(2),Q(3),Q(5,2)]:
        for a in[Q(1),Q(-1),Q(1,2)]:
            # L=Q(t),g=t^2: trace of f'/f equals logarithmic norm derivative.
            trace=1/(2*t*(t-a))+1/((-2*t)*(-t-a))
            assert trace==-1/(a*a-t*t)
    assert any(pure for pure,marked in raw([1,1],[-1],[1],4))
    out={'passed':True,'scope':'Exact bounded polynomial/divisibility/norm and raw moment/original contraction controls only; the separate proof supplies all-order separation and prime exclusion.',
         'cases':cases,'original_lifts':originals,'inversion_divisibility_controls':divisibility,
         'quadratic_norm_identity_controls':9,'negative_controls':{'linear_C_purity_rejected':True},
         'status':'passed','schema_version':1,'python':platform.python_version(),
         'sources':{str(path.relative_to(ROOT)):digest(path)for path in pins}}
    if args.write_receipt:args.write_receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'raw_cases':len(cases),'raw_moments':M*len(cases),'original_lifts':len(originals),'original_powers':sum(len(c['powers'])for c in originals)}))
if __name__=='__main__':main()
