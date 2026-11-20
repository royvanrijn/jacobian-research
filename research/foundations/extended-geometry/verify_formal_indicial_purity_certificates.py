#!/usr/bin/env python3
"""Exact singular-point indicial and comparison-homogenization controls."""
if not __debug__:raise RuntimeError('Refusing -O: exact assertions are required')
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent

def trim(p):
    p=list(map(Q,p))
    while p and not p[-1]:p.pop()
    return tuple(p)
def add(p,q,scale=1):return trim([(p[i]if i<len(p)else 0)+scale*(q[i]if i<len(q)else 0)for i in range(max(len(p),len(q)))])
def mul(p,q):
    if not p or not q:return ()
    out=[Q(0)]*(len(p)+len(q)-1)
    for i,v in enumerate(p):
        for j,w in enumerate(q):out[i+j]+=v*w
    return trim(out)
def deriv(p):return trim([i*p[i]for i in range(1,len(p))])
def scale(p,c):return trim([c*v for v in p])
def power(p,n):
    out=(Q(1),)
    for _ in range(n):out=mul(out,p)
    return out
def evaluate(p,x):
    ans=Q(0)
    for v in reversed(p):ans=ans*x+v
    return ans
def fall(k):
    out=(Q(1),)
    for j in range(k):out=mul(out,(-j,1))
    return out

def indicial(M):
    orders={j:next(i for i,v in enumerate(p)if v)for j,p in enumerate(M)if p}
    assert orders
    shift=min(v-j for j,v in orders.items());I=()
    for j,v in orders.items():
        if v-j==shift:I=add(I,scale(fall(j),M[j][v]))
    assert I
    if len(I)==1:return shift,I,[],0
    radius=1+max(abs(v/I[-1])for v in I[:-1])
    roots=[n for n in range(radius.numerator//radius.denominator+1)if evaluate(I,n)==0]
    return shift,I,roots,max([0]+roots)

def infinity_indicial(M):
    degrees={j:len(p)-1 for j,p in enumerate(M)if p}
    shift=max(v-j for j,v in degrees.items());I=()
    for j,v in degrees.items():
        if v-j==shift:I=add(I,scale(fall(j),M[j][-1]))
    assert I
    if len(I)==1:return shift,I,[]
    radius=1+max(abs(v/I[-1])for v in I[:-1])
    roots=[n for n in range(radius.numerator//radius.denominator+1)if evaluate(I,n)==0]
    return shift,I,roots

def compose_multiply_z(L):
    out=[]
    for j in range(len(L)):
        out.append(add(mul((Q(0),Q(1)),L[j]),scale(L[j+1],j+1)if j+1<len(L)else()))
    return out

def compose_first(q,L):
    qp=deriv(q);out=[]
    for j in range(len(L)+1):
        pj=L[j]if j<len(L)else();prev=L[j-1]if j else()
        out.append(add(add(mul(q,prev),mul(q,deriv(pj))),mul(qp,pj),-1))
    while out and not out[-1]:out.pop()
    return out

def homogenize(L):
    q=scale(L[0],-1)
    return compose_first(q,L)if q else L

def apply(M,p):
    out=();current=p
    for coefficient in M:
        out=add(out,mul(coefficient,current));current=deriv(current)
    return out

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output','--write-receipt',dest='output',type=Path);args=parser.parse_args()
    delayed=[]
    for n in (1,2,7,31):
        M=[(-Q(n),),(Q(0),Q(1))]
        _,I,roots,B=indicial(M)
        assert roots==[n]and B==n and not apply(M,(Q(0),)*n+(Q(1),))
        # Multiplication by a power of z cannot change the cutoff.
        shifted=[(Q(0),)*5+p if p else()for p in M]
        assert indicial(shifted)[1:]==(I,roots,B)
        delayed.append({'n':n,'cutoff':B,'indicial':[str(v)for v in I]})
    irregular=[[(Q(-1),), (Q(0),Q(0),Q(1))],
               [(Q(1),),(Q(0),Q(1)),(Q(0),Q(0),Q(0),Q(1))]]
    for M in irregular:
        assert indicial(M)[2]==[]and indicial(M)[3]==0
        assert all(apply(M,(Q(0),)*n+(Q(1),))for n in range(21))
    # Actual normalized sphere period p=t, G=atanh(z)/z.
    Lt=[(Q(0),Q(-2)),(Q(2),Q(0),Q(-4)),(Q(0),Q(1),Q(0),Q(-1))]
    Mt=homogenize(Lt)
    # The product in the proof differs by scalar2 because q=2z.
    shift,It,rt,Bt=indicial(Mt)
    expected=scale(mul(mul((Q(0),Q(1)),(Q(-2),Q(1))),(Q(1),Q(1))),2)
    assert It==expected and rt==[0,2]and Bt==2
    Gt=tuple(Q(1,n+1)if n%2==0 else Q(0)for n in range(121))
    assert all(v==0 for v in apply(Lt,Gt)[:115])
    assert Gt[1]==0 and Gt[2]==Q(1,3)
    # Verify noncommutative composition on independent polynomial probes.
    for n in range(12):
        probe=(Q(1),)*n+(Q(2),)
        q=scale(Lt[0],-1)
        assert apply(Mt,probe)==add(mul(q,deriv(apply(Lt,probe))),mul(deriv(q),apply(Lt,probe)),-1)
    # Verified companion 3-height operator: indicial computation only here.
    c0=scale((Q(7),Q(-84),Q(4)),-12)
    c1=scale((Q(-147),Q(-14),Q(-3612),Q(232)),-6)
    c2=scale((Q(0),Q(2121),Q(-94),Q(-3052),Q(232)),-12)
    c3=scale((Q(0),Q(0),Q(5103),Q(-378),Q(-1764),Q(152)),-8)
    c4=scale(mul((Q(0),Q(0),Q(0),Q(1)),mul((Q(-21),Q(2)),(Q(-27),Q(0),Q(4)))),-16)
    L3=[c0,c1,c2,c3,c4];M3=homogenize(L3);shift3,I3,r3,B3=indicial(M3)
    target=(Q(1),)
    for factor in ((Q(0),Q(1)),(Q(-1),Q(1)),(Q(-1),Q(2)),(Q(-7),Q(6)),(Q(1),Q(6))):target=mul(target,factor)
    assert I3==scale(target,-10584)and shift3==-2 and r3==[0,1]and B3==1
    mu2=2*(Q(1,3)-Q(1,7));assert mu2==Q(8,21)
    # Actual fixed-marker polynomial vs nonpolynomial series with mean0.
    Ly=[(Q(-1),),(Q(0),Q(1))]
    assert infinity_indicial(Ly)[2]==[1]
    assert not apply(Ly,(Q(0),Q(1,3)))
    LA=compose_first((Q(0),Q(1)),compose_multiply_z(Lt))
    _,Iinf,rootsinf=infinity_indicial(LA)
    expected_inf=scale(mul(mul((Q(1),Q(1)),(Q(2),Q(1))),(Q(3),Q(1))),-1)
    assert Iinf==expected_inf and rootsinf==[]
    assert indicial(LA)[2]==[1]
    At=tuple(Q(1,n+2)if n%2==1 else Q(0)for n in range(120))
    assert At[0]==0 and At[1]==Q(1,3)
    assert all(v==0 for v in apply(LA,At)[:110])
    # e^z is a synthetic nonpolynomial-solution control, not a sphere claim.
    assert infinity_indicial([(Q(-1),),(Q(1),)])[2]==[]
    # Cancellation-preserving min-degree Euler normalization, includingd0.
    marker_controls=0
    for n in range(1,4):
        for d in range(3):
            for m in range(6):
                radial_degree=d*m
                common=Q(factorial(radial_degree+n-1),factorial(n-1))
                # F=R^d, G=R-n: sphere markerpieces1 and-n aggregate.
                actual=Q(factorial(radial_degree+n),factorial(n-1))-n*common
                assert actual==common*(radial_degree+n-n)
                if d==0:assert actual==0
                marker_controls+=1
                # Unbalanced G=z1, beta=e1 balances it toxi1*z1.
                coefficient=Q(factorial(radial_degree+n),factorial(n))
                balanced_factor=Q(factorial(radial_degree+n),factorial(n-1))
                assert coefficient==balanced_factor*Q(1,n)
                marker_controls+=1
    # Bad coefficient-specific shortcut would miss delayed z31.
    assert delayed[-1]['cutoff']==31>len([(Q(-31),),(Q(0),Q(1))])
    pins={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in[HERE/'FORMAL_INDICIAL_PURITY_CERTIFICATES.md',Path(__file__)]}
    out={'passed':True,'scope':'Exact operator algebra and finite controls supplement the unrestricted formal proof. Companion3height annihilator validity must be verified separately by its rational telescoping certificate.',
         'delayed_operator_controls':delayed,'irregular_root_free_controls':2,
         'sphere_t':{'cutoff':Bt,'indicial':[str(v)for v in It],'moment2':str(Gt[2])},
         'three_height':{'J_cutoff':B3,'sphere_cutoff':2*B3,'indicial':[str(v)for v in I3],'sphere_moment2':str(mu2)},
         'noncommutative_composition_controls':12,'general_marker_normalization_controls':marker_controls,'fixed_marker_controls':{'pure_y_marker_x':'z/3 polynomial','t_marker_t':'nonpolynomial, mean0, firstcoefficient1/3','infinity_indicial':[str(v)for v in Iinf]},'packet_sha256':pins}
    path=args.output or HERE/'formal_indicial_purity_certificates_verification.json'
    path.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
