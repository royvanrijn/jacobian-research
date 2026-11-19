#!/usr/bin/env python3
"""Exact pole/divisor and original-contraction controls for the universal proof."""
if not __debug__:
    raise RuntimeError('Assertions must be enabled; -O is unsupported')
from fractions import Fraction as F
from pathlib import Path
from math import comb,factorial
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
ROOT=HERE

def trim(p):
    p=list(map(F,p))
    while p and not p[-1]:p.pop()
    return tuple(p)
def add(p,q,scale=1):
    out=[F(0)]*max(len(p),len(q))
    for j,v in enumerate(p):out[j]+=v
    for j,v in enumerate(q):out[j]+=scale*v
    return trim(out)
def mul(p,q):
    if not p or not q:return ()
    out=[F(0)]*(len(p)+len(q)-1)
    for i,v in enumerate(p):
        for j,w in enumerate(q):out[i+j]+=v*w
    return trim(out)
def power(p,n):
    out=(F(1),)
    for _ in range(n):out=mul(out,p)
    return out
def divmod_poly(p,q):
    assert q
    p=list(p);out=[F(0)]*max(0,len(p)-len(q)+1)
    while p and len(p)>=len(q):
        s=len(p)-len(q);c=p[-1]/q[-1];out[s]+=c
        for j,v in enumerate(q):p[s+j]-=c*v
        while p and not p[-1]:p.pop()
    return trim(out),trim(p)
def exactdiv(p,q):
    a,b=divmod_poly(p,q);assert not b;return a
def gcd_poly(p,q):
    while q:p,q=q,divmod_poly(p,q)[1]
    return tuple(v/p[-1] for v in p)if p else ()
def derivative(p):return trim([j*p[j] for j in range(1,len(p))])
def evaluate(p,x):
    out=F(0)
    for v in reversed(p):out=out*x+v
    return out
def shift(p,x):
    return trim([sum((p[j]*comb(j,k)*x**(j-k) for j in range(k,len(p))),F(0))for k in range(len(p))])
def fraction(n,d=(F(1),)):
    assert d
    if not n:return ((),(F(1),))
    g=gcd_poly(n,d);n=exactdiv(n,g);d=exactdiv(d,g);c=d[-1]
    return tuple(v/c for v in n),tuple(v/c for v in d)
def fadd(p,q,scale=1):return fraction(add(mul(p[0],q[1]),mul(q[0],p[1]),scale),mul(p[1],q[1]))
def fscale(p,c):return fraction(tuple(c*v for v in p[0]),p[1])
def fderivative(p):return fraction(add(mul(derivative(p[0]),p[1]),mul(p[0],derivative(p[1])),-1),power(p[1],2))
def degree(p):return max(len(p[0]),len(p[1]))-1 if p[0]else 0
def valuation(p,alpha):
    if not p:return 10**9
    count=0
    while evaluate(p,alpha)==0:p=exactdiv(p,(-alpha,F(1)));count+=1
    return count
def fpole(p,alpha):return max(valuation(p[1],alpha)-valuation(p[0],alpha),0)
def infinity_pole(p):return max(len(p[0])-len(p[1]),0)
def factored(roots):
    out=(F(1),)
    for alpha,m in roots.items():out=mul(out,power((F(1),-1/F(alpha)),m))
    return out

def rational_primitive(w,rootset):
    n,d=w;quotient,_=divmod_poly(n,d)
    integral=(F(0),)+tuple(v/F(j+1)for j,v in enumerate(quotient))
    r=fraction(trim(integral));logs=[]
    for alpha in rootset:
        m=valuation(d,alpha)
        if not m:continue
        rest=exactdiv(d,power((-alpha,F(1)),m));ns=shift(n,alpha);ds=shift(rest,alpha)
        series=[]
        for j in range(m):
            value=(ns[j]if j<len(ns)else F(0))-sum((ds[k]*series[j-k]for k in range(1,min(j,len(ds)-1)+1)),F(0))
            series.append(value/ds[0])
        for k in range(1,m+1):
            coefficient=series[m-k]
            if k==1:
                if coefficient:logs.append((alpha,coefficient))
            elif coefficient:
                r=fadd(r,fraction((-coefficient/F(k-1),),power((-alpha,F(1)),k-1)))
    r=fadd(r,fraction((-evaluate(r[0],F(0))/evaluate(r[1],F(0)),)))
    reconstructed=fderivative(r)
    for alpha,c in logs:reconstructed=fadd(reconstructed,fraction((c,),(-alpha,F(1))))
    assert reconstructed==w
    return r,logs

def divisor_control(a,pr,jr,hr):
    assert not set(pr)&set(jr)
    P,J,H=map(factored,(pr,jr,hr));p,j,h=len(P)-1,len(J)-1,len(H)-1
    u=fraction(mul((F(0),F(1)),power(P,a)),mul(H,power(J,a+2)))
    omega=fraction(power(P,a-1),mul(H,power(J,a+1))) if a>=1 else fraction((F(1),),mul(mul(H,J),P))
    roots=sorted(set(pr)|set(jr)|set(hr));R,logs=rational_primitive(omega,roots)
    rminus=fadd(R,u,-2)
    for alpha in roots:
        pa,ja,ha=pr.get(alpha,0),jr.get(alpha,0),hr.get(alpha,0)
        assert fpole(u,alpha)==max(ha+(a+2)*ja-a*pa,0)
        assert fpole(R,alpha)==max(ha+(a+1)*ja-(a-1)*pa-1,0)
    q=1+a*p-h-(a+2)*j
    assert infinity_pole(u)==max(q,0)
    assert infinity_pole(R)==max(q+j-p,0)
    compensating=any(hr.get(alpha,0)>=a*pa+1 for alpha,pa in pr.items())
    trivial=(p==j==h==0)
    if not compensating and not trivial:
        assert all(fpole(R,alpha)<fpole(u,alpha)for alpha in roots if fpole(u,alpha))
        if infinity_pole(u):assert infinity_pole(R)<infinity_pole(u)
        assert degree(rminus)>degree(R)
    minimal_D=h+2*p==a+3 and compensating
    if minimal_D and j:
        assert len(pr)==1 and list(pr.values())==[1]
        assert degree(rminus)-degree(R)>=sum(jr.values())+len(jr)-1>0
    return {'a':a,'degrees':{'P':p,'J':j,'H':h,'R':degree(R),'u':degree(u),'R_minus_2u':degree(rminus)},
            'compensating_P_root':compensating,'logarithmic_residues':len(logs),'minimal_D_nonconstant_J':bool(minimal_D and j)}

def raw_angular(a,m,q,C,D):
    if a*m%2:return F(0)
    index=m-q
    if index<0:return F(0)
    out=F(0)
    for j in range(m+1):
        poly=mul(power(C,m-j),power(D,j))
        coefficient=poly[index]if index<len(poly)else F(0)
        out+=F(comb(m,j),a*m+2*j+1)*coefficient
    return out

def mmul(p,q):
    out={}
    for m,v in p.items():
        for n,w in q.items():
            k=tuple(a+b for a,b in zip(m,n));out[k]=out.get(k,0)+v*w
    return {m:v for m,v in out.items()if v}
def madd(p,q,scale=1):
    out=dict(p)
    for m,v in q.items():out[m]=out.get(m,0)+scale*v
    return {m:v for m,v in out.items()if v}
def mpower(p,n):
    out={(0,0,0,0):1}
    for _ in range(n):out=mmul(out,p)
    return out
def contraction(p):
    return sum((v*factorial(i)*factorial(j)for(i,j,k,l),v in p.items()if i==k and j==l),0)
def original_controls():
    A={(1,0,1,0):1};B={(0,1,0,1):1}
    R=madd(A,B);T=madd(A,B,-1);Z={(1,0,0,1):1};W={(0,1,1,0):2}
    cubic=madd(mmul(W,mpower(R,2)),mmul(mpower(T,2),madd(R,Z,F(1,2))),-1)
    # Scale to the literal retained B=2WR²-T²(2R+Z).
    cubic={m:2*v for m,v in cubic.items()}
    out=[]
    for a in (0,1,2,3,5):
        phase=mmul(mmul(mpower(T,a),mpower(madd(R,Z),a+1)),cubic)
        d=2*a+4
        for m in (1,2):
            phasepower=mpower(phase,m);pure=contraction(phasepower);marked=contraction(mmul(Z,phasepower))
            kappa=sum((F((-1)**j*comb(m,j),a*m+2*j+1)for j in range(m+1)),F(0))if (a*m)%2==0 else F(0)
            assert pure==0 and marked==factorial(d*m+2)*kappa
            out.append({'a':a,'m':m,'phase_terms':len(phase),'power_terms':len(phasepower),'normalized_marked':str(F(marked,factorial(d*m+2)))})
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    rows=[]
    for a in (0,1,2,3,4,5,7):
        for h in range(a+3):rows.append(divisor_control(a,{-1:1},{},{-1:h}if h else{}))
        for pr,jr,hr in [({},{},{}),({},{-2:1},{}),({-1:1},{-2:1},{}),
             ({-1:1},{-2:1},{-3:1}),({-1:2},{-2:1},{-1:1}),
             ({-1:1,-2:1},{-3:1},{-1:1})]:rows.append(divisor_control(a,pr,jr,hr))
        for jr in ({-2:1},{-2:2},{-2:1,-3:1}):rows.append(divisor_control(a,{-1:1},jr,{-1:a+1}))
    angular=[]
    for a in (0,1,2,3,4,5,7,9):
        C=power((F(1),F(1)),a+1);D=tuple(-v for v in power((F(1),F(1)),a+3))
        for m in range(1,13):
            pure=raw_angular(a,m,0,C,D);marked=raw_angular(a,m,1,C,D)
            kappa=sum((F((-1)**j*comb(m,j),a*m+2*j+1)for j in range(m+1)),F(0))if (a*m)%2==0 else F(0)
            assert pure==0 and marked==kappa
            if (a*m)%2==0:assert kappa>0
            angular.append({'a':a,'m':m,'normalized_marked':str(marked)})
    # Reconstruct the UNIVERSAL even zero-core determinant numerator inZ[a].
    a1=(F(1),F(1));a3=(F(3),F(1))
    b1=(F(1),F(2));b3=(F(3),F(2));b5=(F(5),F(2))
    numerator=add(add(mul(power(a1,2),mul(b3,b5)),
                       mul(mul(a1,a3),mul(b1,b5)),-2),
                   mul(power(a3,2),mul(b1,b3)))
    assert numerator==(F(12),F(24),F(8))
    # A directly altered contraction expectation must be rejected by the
    # same exact equality used above, not by an approximate tolerance.
    originals=original_controls();assert any(F(x['normalized_marked']) for x in originals)
    try:assert F(originals[0]['normalized_marked'])==0
    except AssertionError:pass
    else:raise AssertionError('Corrupt marked moment accepted')
    ledger_path=HERE/'all_start_height_sharp_unsafe_degree_source_ledger.json'
    paths=[HERE/'ALL_START_HEIGHT_SHARP_UNSAFE_DEGREE.md',Path(__file__).resolve(),ledger_path]
    ledger=json.loads(ledger_path.read_text())
    for resource in ledger['retained_resources']:
        path=ROOT/resource['retained_path'];assert hashlib.sha256(path.read_bytes()).hexdigest()==resource['sha256'];paths.append(path)
    receipt={'passed':True,'scope':'Exact rational partial fractions, full pole-divisor degree comparisons, equality-case J-root gain budgets, direct raw angular contractions, and literal original four-variable controls of the retained witness. These finite controls do not replace the all-order proof or import a finite-prefix purity conclusion.',
             'divisor_controls':rows,'raw_angular_controls':angular,'original_four_variable_controls':originals,
             'even_zero_core_numerator_coefficients':[str(x)for x in numerator],
             'corrupted_marked_expectation_rejected':True,
             'pins':[{'path':str(p.relative_to(HERE)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in paths]}
    receipt.update({'schema_version':1,'status':'passed','independent_replay':False,'external_review':False,'formal_verification':False})
    output=args.output or HERE/'all_start_height_sharp_unsafe_degree_verification.json'
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'passed':True,'divisor_controls':len(rows),'raw_angular_controls':len(angular),'original_powers':len(originals),'output':str(output)}))
if __name__=='__main__':main()
