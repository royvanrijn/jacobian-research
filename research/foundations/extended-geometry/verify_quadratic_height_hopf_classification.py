#!/usr/bin/env python3
"""Exact finite controls for QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md.

The separate note supplies the unrestricted classification and all-order
proofs. Raw interval expansions and original four-variable contractions
below independently control the displayed formulas and normalizations.
"""
if not __debug__:
    raise RuntimeError('This certificate requires Python assertions enabled.')

from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
from dataclasses import dataclass
from types import SimpleNamespace
import argparse,json,hashlib,platform

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def trim(a):
    a=list(map(Q,a))
    while len(a)>1 and not a[-1]:a.pop()
    return a or[Q(0)]
def coeff(a,n):return a[n]if 0<=n<len(a)else Q(0)
def padd(a,b,scale=1):
    return trim([coeff(a,j)+scale*coeff(b,j)for j in range(max(len(a),len(b)))])
def pmul(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):out[i+j]+=v*w
    return trim(out)
def ppow(a,n):
    out=[Q(1)]
    for _ in range(n):out=pmul(out,a)
    return out
def pdiv(a,b):
    a,b=trim(a),trim(b);assert b!=[0]
    out=[Q(0)]*max(1,len(a)-len(b)+1)
    while a!=[0]and len(a)>=len(b):
        s=len(a)-len(b);v=a[-1]/b[-1];out[s]=v
        a=padd(a,[Q(0)]*s+[v*c for c in b],-1)
    return trim(out),a
h=SimpleNamespace(trim=trim,add=padd,multiply=pmul,power=ppow,polynomial_division=pdiv)

@dataclass
class Balanced:
    degree:int
    # (i,j) means xi1^i xi2^(degree-i) z1^j z2^(degree-j).
    terms:dict
    def __post_init__(self):
        assert self.degree>=0
        assert all(0<=i<=self.degree and 0<=j<=self.degree
                   and type(v)is int for(i,j),v in self.terms.items())
        self.terms={p:v for p,v in self.terms.items()if v}
def badd(*polys):
    d=polys[0].degree;assert all(p.degree==d for p in polys)
    out={}
    for p in polys:
        for key,v in p.terms.items():out[key]=out.get(key,0)+v
    return Balanced(d,out)
def bscale(p,v):return Balanced(p.degree,{key:v*c for key,c in p.terms.items()})
def bmul(a,b):
    out={}
    for(i,j),v in a.terms.items():
        for(k,l),w in b.terms.items():
            key=i+k,j+l;out[key]=out.get(key,0)+v*w
    return Balanced(a.degree+b.degree,out)
def bpow(p,n):
    out=Balanced(0,{(0,0):1})
    for _ in range(n):out=bmul(out,p)
    return out
def contraction(p,marker_power=0):
    return sum(v*factorial(i+marker_power)*factorial(p.degree-i)
               for(i,j),v in p.terms.items()if j==i+marker_power)
def evaluate(p,a,b,c,d):
    return sum(v*a**i*b**(p.degree-i)*c**j*d**(p.degree-j)
               for(i,j),v in p.terms.items())
o=SimpleNamespace(Balanced=Balanced,add=badd,scale=bscale,multiply=bmul,power=bpow,
                  contraction=contraction,evaluate=evaluate,
                  R=Balanced(1,{(1,1):1,(0,0):1}),T=Balanced(1,{(1,1):1,(0,0):-1}),
                  Z=Balanced(1,{(1,0):1}),W=Balanced(1,{(0,1):2}))

def scale(a,v):return h.trim([v*c for c in a])
def shift(a):return [Q(0)]+list(a)
def pochhammer(a,m):
    z=Q(1)
    for j in range(m):z*=a+j
    return z
def kappa(m):return Q(4**m*factorial(m)**2,factorial(2*m+1))

def construct(n,c,lam,shear):
    c,lam=Q(c),Q(lam);U=[Q(1),lam]
    S=[Q(comb(n,j+1))*lam**j for j in range(n)]
    D=scale(h.multiply(S,h.power(U,n+2)),-1/(c*n*n*lam))
    K=2/(c*n*lam)
    B=scale(h.multiply(h.power(U,n+1),h.add([1],scale(S,Q(1,n)),-1)),K)
    B=h.add(B,shift(h.multiply(D,shear)))
    numerator=h.add(h.power(B,2),scale(h.power(U,2*n+2),K*K),-1)
    C,remainder=h.polynomial_division(numerator,scale(D,4))
    assert remainder==[0] and B[0]==0 and C[0]+D[0]==0
    A_num=h.power(U,n);A_den=scale(S,c*lam)
    left=h.multiply(h.add(scale(h.multiply(C,D),4),h.power(B,2),-1),A_den)
    right=scale(h.multiply(D,A_num),4)
    assert left==right
    return C,B,D

def raw_moments(C,B,D,maximum):
    terms=[(i,t,v)for t,P in enumerate([C,B,D])for i,v in enumerate(P)if v]
    current={(0,0):Q(1)};out=[]
    for m in range(1,maximum+1):
        nxt={}
        for (i,t),v in current.items():
            for j,u,w in terms:
                if i+j>maximum:continue
                key=i+j,t+u;nxt[key]=nxt.get(key,Q(0))+v*w
        current={k:v for k,v in nxt.items()if v}
        pure=sum((v/Q(t+1)for(i,t),v in current.items()if i==m and t%2==0),Q(0))
        mixed=sum((v/Q(t+1)for(i,t),v in current.items()if i==m-1 and t%2==0),Q(0))
        out.append((pure,mixed))
    return out

def original_degree_six():
    R,Z,W,T=o.R,o.Z,o.W,o.T
    quartic=o.add(o.scale(o.power(R,4),9),o.scale(o.multiply(Z,o.power(R,3)),32),
                  o.scale(o.multiply(o.power(Z,2),o.power(R,2)),56),
                  o.scale(o.multiply(o.power(Z,3),R),48),o.scale(o.power(Z,4),16))
    F=o.add(o.scale(o.multiply(W,o.power(R,5)),2),o.scale(o.power(R,6),3),
            o.scale(o.multiply(o.multiply(T,o.power(o.add(R,o.scale(Z,2)),2)),o.power(R,3)),2),
            o.scale(o.multiply(o.power(T,2),quartic),-1))
    assert F.degree==6 and o.evaluate(F,1,1,1,-1)==-64
    specialize=[sum(v for(i,j),v in F.terms.items()if 6-j==k)for k in range(7)]
    target=scale(h.multiply(h.multiply([0,1],[1,8,11]),[-1,-4,-3,4]),-4)
    assert specialize==target and specialize[0]==0 and specialize[1]==4
    P=o.Balanced(0,{(0,0):1});rows=[]
    for m in range(1,9):
        P=o.multiply(P,F)
        assert o.contraction(P)==0
        expected=Q(2**(2*m-1),2*m+1)
        assert Q(o.contraction(P,1),factorial(6*m+2))==expected
        rows.append({'m':m,'support':len(P.terms),'pure':0,'normalized_marked':str(expected)})
    return {'phase_degree':6,'radial_zero_value':-64,'simple_specialization_zero_derivative':4,
            'specialization_coefficients':specialize,'original_four_variable_controls':rows}

def singular_core_controls():
    # Substitute c=-d/3 and reduce b²=-4d²/15 in actual interval moments.
    data=[]
    for m in[1,2,3]:
        terms={(0,0,0):Q(1)}
        for _ in range(m):
            nxt={}
            for (b,d,t),v in terms.items():
                for B,D,T,w in[(1,0,1,Q(1)),(0,1,2,Q(1)),(0,1,0,-Q(1,3))]:
                    key=b+B,d+D,t+T;nxt[key]=nxt.get(key,Q(0))+v*w
            terms=nxt
        moment={}
        for(b,d,t),v in terms.items():
            if t%2:continue
            key=b,d;moment[key]=moment.get(key,Q(0))+v/Q(t+1)
        reduced={}
        for(b,d),v in moment.items():
            key=b%2,d+2*(b//2);reduced[key]=reduced.get(key,Q(0))+v*(-Q(4,15))**(b//2)
        reduced={k:v for k,v in reduced.items()if v}
        assert reduced==({(0,3):-Q(256,4725)}if m==3 else{})
        data.append({'m':m,'after_first_two_equations':{str(k):str(v)for k,v in reduced.items()}})
    return data

def safe_branch_controls():
    records=[]
    for a,B in[(Q(1),[0]),(Q(2),[0,1]),(Q(-3),[0,1,2]),(Q(1,2),[0,2,-1,0,3])]:
        C=h.add([a],scale(h.power(B,2),-1/(4*a)));D=[-a]
        raw=raw_moments(C,B,D,14)
        assert all(pure==0 for pure,mixed in raw)
        assert raw[0][1]==Q(2,3)*a
        assert all(mixed==0 for pure,mixed in raw[1:])
        b=len(h.trim(B))-1
        assert max(1,len(C)-2,b,len(D))==(2*b-1 if b else 1)
        records.append({'branch':'constant_center','A':str(a),'B':B,'moment_controls':len(raw)})
    for C,B,D in[([0,0,1],[0,0,2],[0,0,-3]),
                  ([0,0,1,2],[0,0,0,-2],[0,0,3,0,1])]:
        assert all(pure==mixed==0 for pure,mixed in raw_moments(C,B,D,14))
        records.append({'branch':'strictly_positive_weight','C':C,'B':B,'D':D,'moment_controls':14})
    return records

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path)
    args=parser.parse_args()
    records=[]
    for n in range(1,5):
        for shear in [[0],[2],[0,1]]:
            c,lam=(Q(1,4),Q(2))if n%2==0 else(Q(-2),Q(3,2))
            C,B,D=construct(n,c,lam,shear)
            raw=raw_moments(C,B,D,12)
            for m,(pure,mixed)in enumerate(raw,1):
                assert pure==0
                assert mixed==kappa(m)*pochhammer(Q(1,n),m)/(lam*c**m*factorial(m))
            b=len(B)-1 if B!=[0]else -1
            degree=max(1,len(C)-2,b,len(D))
            assert degree==(2*n+2 if b<=2*n+2 else 2*b-2*n-2) and degree%2==0
            records.append({'n':n,'c':str(c),'lambda':str(lam),'shear':shear,'C':list(map(str,C)),
                            'B':list(map(str,B)),'D':list(map(str,D)),'degree':degree,'moment_controls':len(raw)})
    C=[1,3];B=[0,2,8,8];D=[-1,-9,-32,-56,-48,-16]
    for m,(pure,mixed)in enumerate(raw_moments(C,B,D,12),1):
        assert pure==0 and mixed==Q(2**(2*m-1),2*m+1)
    alteredD=D.copy();alteredD[1]+=1
    assert raw_moments(C,B,alteredD,2)[0][0]!=0
    alteredB=B.copy();alteredB[1]+=1
    assert raw_moments(C,alteredB,D,2)[1][0]!=0
    pins=[HERE/'QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md',Path(__file__),
          HERE/'quadratic_height_hopf_classification_source_ledger.json']
    ledger=json.loads(pins[-1].read_text())
    for resource in ledger['retained_resources']:
        p=HERE/resource['retained_path']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==resource['sha256']
        pins.append(p)
    out={'status':'passed','python':platform.python_version(),
         'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in pins},
         'formal_family_cases':records,'degree_six':original_degree_six(),
         'zero_center_core':singular_core_controls(),'safe_branch_controls':safe_branch_controls(),
         'negative_controls':{'altered_D_rejected':True,'altered_B_rejected':True},
         'scope':'Bounded exact original interval and original four-variable contractions, parameter/divisibility/degree controls and singular-core reductions. The unrestricted classification, parity, all-order formula and safety use the separate written proof.'}
    if args.write_receipt:args.write_receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'pass':True,'family_cases':len(records),'family_moment_controls':12*len(records),'original_F6_controls':8,'degree_six_support':out['degree_six']['original_four_variable_controls'][-1]['support']}))

if __name__=='__main__':main()
