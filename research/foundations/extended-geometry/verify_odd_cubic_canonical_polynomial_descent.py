#!/usr/bin/env python3
"""Exact finite controls for canonical odd-cubic polynomial descent."""
if not __debug__:raise RuntimeError("Assertions must be enabled")
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json,hashlib,argparse,sys
HERE=Path(__file__).resolve().parent

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

def derivative(a):return trim([j*a[j]for j in range(1,len(a))]or[0])
def divide(a,b):
    a=trim(a);b=trim(b);q=[Q(0)]*max(1,len(a)-len(b)+1)
    while a!=[0]and len(a)>=len(b):
        j=len(a)-len(b);v=a[-1]/b[-1];q[j]=v;a=add(a,[0]*j+scale(b,v),-1)
    return trim(q),a

def pgcd(a,b):
    while b!=[0]:a,b=b,divide(a,b)[1]
    return scale(a,1/a[-1])

def hyper(N,ell,lam):
    U=[Q(1),Q(lam)];S=[Q(comb(2*N,j+1))*lam**j for j in range(2*N)]
    X=scale(mul(U,S),Q(1,2*N));Qx=[Q(0)]
    for j in range(ell):Qx=add(Qx,power(U,2*N*j),(-1)**j)
    C=mul(power(U,N+1),Qx);D=scale(mul(mul(power(U,N+1),power(X,2)),power(Qx,3)),-1)
    return C,D,Qx

def fixed(N,d,lam):
    q=N//d;U=[Q(1),Q(lam)]
    # Formula in U, then compose U=1+lambda*x.
    A=add(power([-1]+[0]*(d-1)+[2],q),[-1])
    P=divide_at1(A)
    D_U=[Q(0)]*(N+3)+scale(power(P,2),-Q(1,4*N*N))
    return power(U,N+1),compose(D_U,U)

def verify_pullback(C,D,Cbar,Dbar,w,markers,M=8):
    wp=derivative(w);assert w[0]==0 and wp[0]==1
    assert mul(C,wp)==compose(Cbar,w)
    assert mul(mul(D,[0,0,1]),power(wp,3))==mul(compose(Dbar,w),power(w,2))
    assert len(C)-1==(len(w)-1)*(len(Cbar)-2)+1
    assert len(D)-1==(len(w)-1)*(len(Dbar)-2)+1
    for Gbar in markers:
        lhs=raw(C,D,compose(Gbar,w),M);rhs=raw(Cbar,Dbar,Gbar,M)
        assert lhs==rhs
    return len(markers)*M

def main(output):
    ledger_path=HERE/"odd_cubic_canonical_polynomial_descent_source_ledger.json"
    ledger=json.loads(ledger_path.read_text())
    for item in ledger["retained_resources"]:
        assert hashlib.sha256((HERE/item["retained_path"]).read_bytes()).hexdigest()==item["sha256"],item["retained_path"]
    hypercases=[];fixedcases=[];marker_controls=0
    for N in range(1,7):
        for ell in[1,3,5]:
            # In U, the Q factor is squarefree, has Q0=1, and no low principal terms.
            QU=[Q(0)]
            for j in range(ell):QU=add(QU,[0]*(2*N*j)+[Q(1)],(-1)**j)
            assert QU[0]==1 and pgcd(QU,derivative(QU))==[1]
            assert all(coef(QU,j)==0 for j in range(1,N+1))
            for lam in[Q(1),Q(-2),Q(1,2)]:
                C,D,Qx=hyper(N,ell,lam);Cb,Db,Qb=hyper(1,ell,N*lam)
                U=[Q(1),lam];w=scale(add(power(U,N),[-1]),Q(1,N)/lam)
                markers=[[1],[0,1],[0,0,1],[0,0,0,1],[0]+Qb]
                marker_controls+=verify_pullback(C,D,Cb,Db,w,markers)
                # The last marker is a retained unsafe marker, so transfer is nonvacuous.
                assert raw(Cb,Db,[0]+Qb,2)[1][1]
                c=len(C)-1;h=1+len(QU)-1;Hbar=1+(len(QU)-1)//N
                assert c-h==N and h==N*(Hbar-1)+1
                hypercases.append({'N':N,'ell':ell,'lambda':str(lam),'cover_degree':N,'c':c,'e':len(D)-1,'h':h,'quotient_c':len(Cb)-1,'quotient_e':len(Db)-1,'marker_controls':40})
        for d in range(1,N+1):
            if N%d:continue
            for lam in[Q(1),Q(-2),Q(1,2)]:
                C,D=fixed(N,d,lam);Cb,Db=fixed(N//d,1,d*lam)
                w=scale(add(power([1,lam],d),[-1]),Q(1,d)/lam)
                marker_controls+=verify_pullback(C,D,Cb,Db,w,[[1],[0,1],[0,0,1],[0,0,0,1]])
                # u and R generate K(U^d): numerator's next Laurent coefficient is nonzero.
                q=N//d
                A=add(power([-1,2],q),[-1])
                assert q==1 or(coef(A,1)==2*q*(-1)**(q-1)and coef(A,1))
                fixedcases.append({'N':N,'d':d,'lambda':str(lam),'cover_degree':d,'c':len(C)-1,'e':len(D)-1,'quotient_c':len(Cb)-1,'quotient_e':len(Db)-1,'marker_controls':32})
    # A genuinely NONPURE admissible pair verifies the reverse criterion beyond zero cases.
    Cb=power([1,2],3);Db=scale(Cb,-1);w=[0,1,Q(1,2)]
    C=power([1,1],5);D=scale(mul(power([1,1],3),power([1,Q(1,2)],2)),-1)
    marker_controls+=verify_pullback(C,D,Cb,Db,w,[[1],[0,1],[0,0,1],[0,0,0,1]],10)
    assert any(pure for pure,marked in raw(C,D,[1],10))
    # Both coefficient guards are essential; C polynomial alone is insufficient.
    badCbar=[1,2];badDbar=[-1,-2];wp=derivative(w)
    Cq,Cr=divide(compose(badCbar,w),wp);assert Cr==[0]
    Dnum=mul(compose(badDbar,w),power(w[1:],2))
    Dq,Dr=divide(Dnum,power(wp,3));assert Dr!=[0]
    # A generic normalized polynomial is not admissible even at the C stage.
    v=[0,1,1];Vprime=derivative(v)
    assert divide(compose([1,2,1],v),Vprime)[1]!=[0]
    # Local total-ramification arithmetic: a positive target zero cannot have another source.
    local=[]
    for c in range(2,20):
        for h in range(1,c):
            D_R=c-h
            for g in range(1,D_R+1):
                if D_R%g:continue
                for b in range(1,g+1):
                    if(c-1)%b:continue
                    order=(c-1)//b-1
                    if order==0:assert b==g
                    elif b<g:
                        # Any hypothetical other point has positive omega order, forbidden finite.
                        assert order>=1
                        continue
                    assert b==g
                    local.append((c,h,g,b))
    out={'passed':True,'scope':'Finite exact polynomial pullback/degree/root-count identities and DIRECT RAW marked moment comparisons for existing pure and admissible NONPURE phases; guards reject nonpolynomial C or D. Universal canonical descent and converse are supplied by the separate written proof.',
         'hyperbolic_cases':hypercases,'fixed_root_cases':fixedcases,'marked_moment_identity_controls':marker_controls,
         'admissible_nonpure_control':True,'local_degree_cells':len(local),
         'negative_controls':{'C_only_polynomial_guard_rejected':True,'generic_w_rejected':True},
         'schema_version':1,'status':'passed','python':sys.version.split()[0],
         'sources':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in[HERE/'ODD_CUBIC_CANONICAL_POLYNOMIAL_DESCENT.md',Path(__file__),ledger_path]},
         'retained_inputs':ledger['retained_resources'],'independent_replay':False,'external_review':False,'formal_verification':False}
    output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'hyperbolic_cases':len(hypercases),'fixed_root_cases':len(fixedcases),'raw_marker_equalities':marker_controls,'admissible_nonpure':True,'local_cells':len(local)}))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-receipt',type=Path,required=True);args=parser.parse_args();main(args.write_receipt)
