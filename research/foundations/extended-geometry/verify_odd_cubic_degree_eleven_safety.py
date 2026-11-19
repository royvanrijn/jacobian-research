#!/usr/bin/env python3
"""Exact finite controls for odd-cubic degree-eleven exclusion."""
if not __debug__:raise RuntimeError("Assertions must be enabled")
from fractions import Fraction as Q
from math import comb,factorial,gcd
from pathlib import Path
import json,hashlib,argparse,sys
HERE=Path(__file__).resolve().parent

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

def local_filters(c,h,s,g,b,n):
    D=c-h
    if g<1 or D%g or n%g:return False
    # b<=g is necessary. b|g is deliberately NOT imposed: the cover need not be Galois.
    if b<1 or b>g or s%b or(c-1)%b:return False
    k=s//b
    if k%2==0 or n<max(s+1,g*(k+1))or n>2*D:return False
    if g==1 and(h%2 or n%2):return False
    return True

def poch(x,n):
    out=Q(1)
    for j in range(n):out*=x+j
    return out

def main(output):
    ledger_path=HERE/"odd_cubic_degree_eleven_safety_source_ledger.json"
    ledger=json.loads(ledger_path.read_text())
    for item in ledger["retained_resources"]:
        source=HERE/item["retained_path"]
        assert hashlib.sha256(source.read_bytes()).hexdigest()==item["sha256"],item["retained_path"]
    # A regression against the unjustified extra assumption b|g.
    assert local_filters(7,4,2,3,2,6)
    patterns=[];survivors=set();partial=[]
    for c in[1,3,5,7,9,11]:
        for p in partitions(c):
            h=len(p);D=c-h
            if c in[5,7]or(c in[9,11]and h<=4):partial.append(primitive_control(p))
            for e in([9]if c<11 else[1,3,5,7,9]):
                assert max(c,e+2)==11 and c%2==e%2
                if h==c:reason='squarefree excluded';cells=[]
                elif h==1:
                    assert 11 not in[3*c,5*c-2];reason='one-pole degree classification';cells=[]
                elif h==2:reason='two-root logarithmic exclusion';cells=[]
                else:
                    s=(3*c-e-2)//2
                    assert s>=D and 2*D<2*s+1
                    cells=[]
                    for g in range(1,D+1):
                        for b in range(1,g+1):
                            for n in range(s+1,2*D+1):
                                if not local_filters(c,h,s,g,b,n):continue
                                k=s//b;cell=(c,e,h,D,s,g,b,k,n);survivors.add(cell)
                                if g==1:
                                    if n==2*D:
                                        assert c>Q(3*h,2)
                                        why='rational-reflection maximal pole bound'
                                    else:
                                        assert c<Q(5*n,6)+Q(h,2)
                                        why='rational-reflection pair capacity'
                                elif c==5:
                                    assert(c,e,h,D,s,g,b,k,n)==(5,9,3,2,2,2,2,1,4)
                                    assert c<e and Q(c-e,2)==-2
                                    assert p in[(2,2,1),(3,1,1)]
                                    why='fixed infinity forces T=1 but its value is zero'
                                else:
                                    assert(c,e,h,D,s,g,b,k,n)==(11,3,3,8,14,2,2,7,16)
                                    assert sum(m==1 for m in p)<D
                                    why='disjoint simple-root pole capacity'
                                cells.append({'g':g,'b':b,'k':k,'degree_u':n,'last_exclusion':why})
                    reason='local ramification/parity/degree filters'if not cells else'final pole or fixed-infinity T-value exclusion'
                patterns.append({'c':c,'e':e,'multiplicities':p,'h':h,'R_degree':D,'exclusion':reason,'surviving_filter_cells':cells})
    expected={(5,9,3,2,2,2,2,1,4),(7,9,4,3,5,1,1,5,6),(11,3,3,8,14,2,2,7,16),(11,5,4,7,13,1,1,13,14),(11,9,4,7,11,1,1,11,12),(11,9,4,7,11,1,1,11,14)}
    assert survivors==expected
    capacities=0
    for m in range(1,12):
        for n in range(m,12):
            for k in range(n):
                if 2*k<=3*m:
                    assert Q(m+n)>=Q(5*k,3)+1;capacities+=1
    # Replay the RETAINED F13 upper bound in ORIGINAL variables.
    R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}
    Qh={}
    for j,v in enumerate([80,180,240,208,120,45,10,1]):Qh=badd(Qh,bmul(bpower(R,7-j),bpower(Z,j)),v)
    F=bmul(bmul(T,bpower(badd(R,Z),3)),badd(bscale(bmul(W,bpower(R,8)),32),bmul(bpower(T,2),Qh),-1))
    assert beval(F,13,1,1,1,-1)==-8
    P={(0,0):Q(1)};controls=[]
    for m in range(1,5):
        P=bmul(P,F);d=13*m
        pure=sum(v*factorial(i)*factorial(d-i)for(i,j),v in P.items()if i==j)
        marked=sum(v*factorial(i+1)*factorial(d-i)for(i,j),v in P.items()if j==i+1)
        if m%2:expect=Q(0)
        else:
            r=m//2;mu=(-4)**r*poch(Q(1,4),r)*poch(Q(3,4),r)/(poch(Q(1,2),r)*factorial(r))
            expect=16**m*kappa(m)*mu
        assert pure==0 and Q(marked,factorial(d+2))==expect
        controls.append({'m':m,'support':len(P),'pure':0,'normalized_Z_marker':str(expect)})
    altered=badd(F,bpower(R,13))
    assert sum(v*factorial(i)*factorial(13-i)for(i,j),v in altered.items()if i==j)==factorial(14)
    out={'passed':True,'scope':'Exact bounded integer/root-pattern and rational partial-fraction controls; original RETAINED F13 contractions support the existing upper bound. Stable-field, pullback-differential, zero-count, and fixed-infinity T-value arguments in the separate proof establish universal degree11 safety.',
         'degree_eleven_patterns':patterns,'six_filter_cells':sorted(survivors),'partial_fraction_controls':partial,
         'reflection_capacity_controls':capacities,'fixed_infinity_T_value_controls':{'c':5,'e':9,'T_growth':-2,'actual_value':0,'required_fixed_value':1},'retained_F13_original_controls':controls,
         'negative_controls':{'does_not_assume_b_divides_g':True,'altered_F13_rejected':True},
         'schema_version':1,'status':'passed','python':sys.version.split()[0],
         'sources':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in[HERE/'ODD_CUBIC_DEGREE_ELEVEN_SAFETY.md',Path(__file__),ledger_path]},
         'retained_inputs':ledger['retained_resources'],'independent_replay':False,'external_review':False,'formal_verification':False}
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'degree_eleven_patterns':len(patterns),'filter_cells':len(survivors),'partial_fraction_controls':len(partial),'reflection_capacity_controls':capacities,'original_F13_powers':len(controls)}))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-receipt',type=Path,required=True);args=parser.parse_args();main(args.write_receipt)
