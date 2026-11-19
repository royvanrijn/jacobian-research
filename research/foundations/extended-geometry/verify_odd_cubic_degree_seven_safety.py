#!/usr/bin/env python3
"""Exact finite algebraic controls for odd-cubic septic exclusion."""
if not __debug__:raise RuntimeError("Assertions must be enabled")
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import argparse,json,hashlib,platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PROOF=HERE/"ODD_CUBIC_DEGREE_SEVEN_SAFETY.md"
LEDGER=HERE/"odd_cubic_degree_seven_safety_source_ledger.json"

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

def kappa(m,a=1):
    if m%2:return Q(0)
    return sum((Q((-1)**j*comb(m,j),a*m+2*j+1)for j in range(m+1)),Q(0))

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

def coeff(a,j):return a[j]if j<len(a)else Q(0)
def evaluate(a,x):return sum((v*x**j for j,v in enumerate(a)),Q(0))
def divide(a,b):
    a=trim(a);b=trim(b);q=[Q(0)]*max(1,len(a)-len(b)+1)
    while a!=[0]and len(a)>=len(b):
        j=len(a)-len(b);v=a[-1]/b[-1];q[j]=v;a=add(a,[0]*j+scale(b,v),-1)
    return trim(q),a

def solve(A,b):
    n=len(b);M=[list(map(Q,row))+[Q(v)]for row,v in zip(A,b)]
    for j in range(n):
        k=next(k for k in range(j,n)if M[k][j]);M[j],M[k]=M[k],M[j]
        M[j]=[v/M[j][j]for v in M[j]]
        for k in range(n):
            if k!=j:
                v=M[k][j];M[k]=[x-v*y for x,y in zip(M[k],M[j])]
    return [row[-1]for row in M]

def partitions(n,cap=None):
    if n==0:yield();return
    for j in range(min(n,n if cap is None else cap),0,-1):
        for tail in partitions(n-j,j):yield(j,)+tail

def primitive_control(p):
    c=sum(p);h=len(p);C=[Q(1)]
    for r,m in enumerate(p,1):C=multiply(C,power([-r,1],m))
    C=scale(C,1/C[0]);columns=[];index=[]
    for r,m in enumerate(p,1):
        for j in range(1,m+1):
            q,rem=divide(C,power([-r,1],j));assert rem==[0]
            columns.append(q);index.append((r,j))
    a=solve([[coeff(q,k)for q in columns]for k in range(c)],[1]+[0]*(c-1))
    total=[Q(0)]
    for v,col in zip(a,columns):total=add(total,scale(col,v))
    assert total==[1]
    residues=[v for(r,j),v in zip(index,a)if j==1]
    if c>=2:assert sum(residues)==0
    B=[Q(1)]
    for r,m in enumerate(p,1):B=multiply(B,power([-r,1],m-1))
    A=[Q(0)]
    for(r,j),v in zip(index,a):
        if j>=2:
            q,rem=divide(B,power([-r,1],j-1));assert rem==[0]
            A=add(A,scale(q,-v/Q(j-1)))
    # Constant normalization R(0)=0 does not change its finite pole divisor.
    A=add(A,scale(B,-A[0]/B[0]))
    for r,m in enumerate(p,1):
        if m>=2:assert evaluate(A,Q(r))
    degree=max(len(A),len(B))-1
    assert degree==c-h
    if h==c:assert A==[0]
    if h==2:assert all(residues)and residues[0]==-residues[1]
    return {'partition':p,'C_coefficients':list(map(str,C)),'R_numerator':list(map(str,A)),'R_denominator':list(map(str,B)),'R_degree':degree,'residues':list(map(str,residues))}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path,help='Write exact finite algebraic-control receipt.')
    args=parser.parse_args()
    ledger=json.loads(LEDGER.read_text())
    pins=[PROOF,Path(__file__).resolve(),LEDGER]
    for item in ledger['retained_resources']:
        path=HERE/item['retained_path']
        assert digest(path)==item['sha256'],('Retained source changed',str(path))
        pins.append(path)
    patterns=[];primitives=[]
    for c in[1,3,5,7]:
        for p in partitions(c):
            h=len(p);D=c-h;primitives.append(primitive_control(p))
            for e in([5]if c<7 else[1,3,5]):
                assert max(c,e+2)==7 and c%2==e%2
                if h==c:reason='squarefree center excluded'
                elif h==1:
                    assert 7 not in[3*c,5*c-2];reason='one-pole degrees classified'
                elif h==2:reason='two-root logarithmic exclusion'
                else:
                    s=(3*c-e-2)//2;assert s>0
                    lower=s+1;upper=2*D
                    if lower>upper:reason='infinity zero degree exceeds pole sum'
                    else:
                        assert(c,e,h,D,s,lower,upper)==(7,5,3,4,7,8,8)
                        simple=sum(m==1 for m in p)
                        assert simple<=2 and simple<D
                        reason='maximal pole union forces simple-root capacity contradiction'
                patterns.append({'c':c,'e':e,'multiplicities':p,'h':h,'R_degree':D,'exclusion':reason})
    R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}
    Qpoly=badd(badd(bscale(bpower(R,3),4),bscale(bmul(bpower(R,2),Z),6)),badd(bscale(bmul(R,bpower(Z,2)),4),bpower(Z,3)))
    F=bmul(bmul(T,bpower(badd(R,Z),3)),badd(bscale(bmul(bpower(R,4),W),2),bmul(bpower(T,2),Qpoly),-1))
    assert beval(F,9,1,1,1,-1)==-8
    P={(0,0):Q(1)};controls=[]
    for m in range(1,5):
        P=bmul(P,F);d=9*m
        pure=sum(v*factorial(i)*factorial(d-i)for(i,j),v in P.items()if i==j)
        marked=sum(v*factorial(i+2)*factorial(d-i)for(i,j),v in P.items()if j==i+2)
        expected=Q(0)if m%2 else kappa(m)
        assert pure==0 and Q(marked,factorial(d+3))==expected
        controls.append({'m':m,'support':len(P),'pure':0,'normalized_Z2_marker':str(expected)})
    altered=badd(F,bpower(R,9))
    assert sum(v*factorial(i)*factorial(9-i)for(i,j),v in altered.items()if i==j)==factorial(10)
    out={'passed':True,'scope':'Exact finite degree/root-partition bookkeeping and rational partial-fraction/original-F9 controls. Field-involution and pole-divisor arguments in the separate proof establish universal septic exclusion; old all-order F9 supplies sharpness.',
         'degree_seven_patterns':patterns,'partial_fraction_controls':primitives,'retained_F9_original_controls':controls,
         'negative_controls':{'altered_F9_rejected':True},
         'status':'passed','schema_version':1,'python':platform.python_version(),
         'sources':{str(path.relative_to(ROOT)):digest(path)for path in pins}}
    if args.write_receipt:args.write_receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'degree_seven_patterns':len(patterns),'partial_fraction_controls':len(primitives),'original_F9_powers':len(controls)}))
if __name__=='__main__':main()
