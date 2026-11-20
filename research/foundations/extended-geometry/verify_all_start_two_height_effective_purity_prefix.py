#!/usr/bin/env python3
"""Exact rational and integer controls for the uniform two-height prefix."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path

if not __debug__:
    raise SystemExit("Refusing -O: exact assertions must remain enabled.")


def trim(p):
    p=tuple(Q(c) for c in p)
    while p and p[-1]==0:p=p[:-1]
    return p


def add(p,q):
    return trim((p[i] if i<len(p) else 0)+(q[i] if i<len(q) else 0)
                for i in range(max(len(p),len(q))))


def scale(p,c):return trim(c*v for v in p)
def sub(p,q):return add(p,scale(q,-1))
def degree(p):return len(p)-1


def mul(p,q):
    if not p or not q:return ()
    out=[Q(0)]*(len(p)+len(q)-1)
    for i,pi in enumerate(p):
        for j,qj in enumerate(q):out[i+j]+=pi*qj
    return trim(out)


def power(p,n):
    assert n>=0
    out=(Q(1),)
    while n:
        if n%2:out=mul(out,p)
        p=mul(p,p);n//=2
    return out


def deriv(p):return trim(i*p[i] for i in range(1,len(p)))


def divmodp(p,q):
    assert q
    rem=list(p);out=[Q(0)]*max(0,len(p)-len(q)+1)
    while rem and len(rem)>=len(q):
        shift=len(rem)-len(q);c=rem[-1]/q[-1];out[shift]=c
        for i,qi in enumerate(q):rem[shift+i]-=c*qi
        rem=list(trim(rem))
    return trim(out),trim(rem)


def exactdiv(p,q):
    out,rem=divmodp(p,q);assert not rem
    return out


def gcdp(p,q):
    while q:p,q=q,divmodp(p,q)[1]
    return scale(p,1/p[-1]) if p else ()


class RF:
    def __init__(self,num=(Q(1),),den=(Q(1),)):
        num,den=trim(num),trim(den);assert den
        g=gcdp(num,den);num,den=exactdiv(num,g),exactdiv(den,g)
        self.num,self.den=scale(num,1/den[-1]),scale(den,1/den[-1])
    def __mul__(self,o):return RF(mul(self.num,o.num),mul(self.den,o.den))
    def __truediv__(self,o):return RF(mul(self.num,o.den),mul(self.den,o.num))
    def __add__(self,o):return RF(add(mul(self.num,o.den),mul(o.num,self.den)),mul(self.den,o.den))
    def __sub__(self,o):return self+RF(scale(o.num,-1),o.den)
    def __pow__(self,n):
        if n<0:return RF(power(self.den,-n),power(self.num,-n))
        return RF(power(self.num,n),power(self.den,n))
    def derivative(self):return RF(sub(mul(deriv(self.num),self.den),mul(self.num,deriv(self.den))),power(self.den,2))
    def map_degree(self):return max(degree(self.num),degree(self.den))
    def __eq__(self,o):return self.num==o.num and self.den==o.den


def moment(C,D,a,m):
    if a*m%2:return Q(0)
    ans=Q(0)
    for j in range(m+1):
        p=mul(power(C,m-j),power(D,j))
        if m<len(p):ans+=Q(comb(m,j),a*m+2*j+1)*p[m]
    return ans


def B(d):return 3*((d+5)**2//4)*(((d+3)**2//4)-1)-4


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name(
        "all_start_two_height_effective_purity_prefix_verification.json"))
    parser.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[3])
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    ledger_path=here/'all_start_two_height_effective_purity_prefix_source_ledger.json'
    ledger=json.loads(ledger_path.read_text())
    pins=[]
    for entry in ledger['retained_resources']:
        source=(args.repo_root/entry['path']) if entry['path'].startswith('research/') else (here/entry['path'])
        got=hashlib.sha256(source.read_bytes()).hexdigest()
        assert got==entry['sha256'],str(source)
        pins.append({'path':entry['path'],'sha256':got})
    xx=(Q(0),Q(1));one=(Q(1),);U=(Q(1),Q(1))
    fixtures=[(one,scale(one,-1)),(U,scale(one,-1)),
              (one,(Q(-1),Q(1))),(U,scale(U,-1)),
              (power(U,2),scale(power(U,4),-1)),
              ((Q(1),Q(2),Q(3)),(Q(-1),Q(1),Q(2)))]
    odd=even=raw=bad=0
    rows=[]
    for a in range(9):
        active=fixtures+[(power(U,a+1),scale(power(U,a+3),-1))]
        for C,D in active:
            c,e=degree(C),degree(D)
            d=max(c+a-1,e+a+1);assert d>=1
            if a%2:
                q=(a-1)//2
                f=RF(power(xx,2))*(RF(scale(D,-1))**a)/(RF(C)**(a+2))
                W=(RF(scale(D,-1))**q)/(RF(C)**(q+1))
                H=RF((Q(2),))*W/f.derivative()
                E=sub(add(scale(mul(C,D),2),scale(mul(xx,mul(C,deriv(D))),a)),
                      scale(mul(xx,mul(deriv(C),D)),a+2))
                Hformula=RF(scale(power(C,q+3),-2),mul(xx,mul(power(scale(D,-1),q),E)))
                assert H==Hformula
                n,h=f.map_degree(),H.map_degree()
                assert n<=max(2+a*e,(a+2)*c)
                assert h<=max(Q(a+5,2)*c,1+c+Q(a+1,2)*e)
                fixed=(n-1)*(n+4*h)-4
                L=d-a+1
                assert n<=(a+2)*L and h<=Q(a+5,2)*L
                assert fixed<=3*(a+4)*L*((a+2)*L-1)-4<=B(d)
                odd+=1
                status="odd_resultant_bound"
            else:
                u=RF(xx)*(RF(scale(D,-1))/RF(C))**(a//2)/RF(C)
                W2=(RF(scale(D,-1))/RF(C))**(a-1)/(RF(C)**2)
                residual=W2-u.derivative()**2
                n,w=u.map_degree(),W2.map_degree()
                fixed=w+4*n
                assert fixed<=B(d)
                if a:
                    L=d-a+1
                    assert n<=Q(a+2,2)*L and w<=(a+1)*L
                    assert fixed<=(3*a+5)*L
                else:assert fixed<=6*d+4
                status="pure" if not residual.num else "nonpure"
                even+=1
            first=None
            for m in range(1,9):
                value=moment(C,D,a,m);raw+=1
                if value and first is None:first=m
            if first is not None:
                assert first<=fixed<=B(d);bad+=1
            if not a%2:
                assert (status=="pure")== (first is None)
            rows.append({"a":a,"c":c,"e":e,"d":d,"status":status,
                         "actual_degree_bound":fixed,"uniform_prefix":B(d),
                         "first_bad_within_eight":first})
    inequalities=0
    for d in range(1,101):
        Nd,Ad=(d+3)**2//4,(d+5)**2//4
        assert Nd>=d+3 and Ad>=2*d+7
        assert B(d)>=6*d*d+33*d+38>=6*d+4
        assert B(d)>=Q((d+2)**2,4)+Q((d+3)**2,2)
        for a in range(d):
            L=d-a+1
            assert (a+2)*L<=Nd and (a+4)*L<=Ad
            if a%2:
                assert 3*(a+4)*L*((a+2)*L-1)-4<=B(d)
            inequalities+=1
    # Negative control: the even square criterion must detect this fixture.
    assert moment(U,scale(one,-1),0,1)==1
    proof=Path(__file__).with_name("ALL_START_TWO_HEIGHT_EFFECTIVE_PURITY_PREFIX.md")
    receipt={"status":"PASS","scope":"Exact finite rational-formula, raw moment and integer-bound controls; the all-order prefix and reduced-ideal results are the written proofs.",
             "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
             "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             "odd_H_formula_controls":odd,"even_square_criterion_controls":even,
             "direct_raw_moment_controls":raw,"nonpure_early_moment_controls":bad,
             "uniform_integer_degree_controls":inequalities,"fixtures":rows,
             "ordinary_moment_ideal_equality_claimed":False,
             "external_or_formal_review_claimed":False}
    assert receipt['proof_sha256']==ledger['proof_sha256']
    assert receipt['checker_sha256']==ledger['verification']['checker_sha256']
    receipt.update({'schema_version':1,'retained_source_pins':pins,
        'source_ledger_sha256':hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
        'independent_replay':False,'formal_verification':False,'external_review':False})
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k!="fixtures"},sort_keys=True))


if __name__=="__main__":main()
