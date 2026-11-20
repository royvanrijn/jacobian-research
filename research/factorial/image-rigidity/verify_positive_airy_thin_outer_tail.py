#!/usr/bin/env python3
"""Exact algebra/rational guards for thin positive Airy outer tails.

The written proof supplies all contour and integral arguments. No numerical
Airy sampling and no original-period remainder certification occurs here.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import argparse, json

if not __debug__:
    raise SystemExit('Unoptimized execution is required')
HERE=Path(__file__).resolve().parent

def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def scale(a,c):return {k:v*F(c) for k,v in a.items() if v*c}
def mul(a,b):
    out={}
    for k,v in a.items():
        for j,w in b.items():
            ij=tuple(x+y for x,y in zip(k,j));out[ij]=out.get(ij,F(0))+v*w
    return {k:v for k,v in out.items() if v}
def var(j,n=4):return {tuple(int(i==j) for i in range(n)):F(1)}
def d(a,i):
    out={}
    for k,v in a.items():
        if k[i]:
            kk=list(k);kk[i]-=1;out[tuple(kk)]=v*k[i]
    return out

def identities():
    p,q,x,y=[var(i) for i in range(4)]
    xx,yy=mul(x,x),mul(y,y)
    L=add(add(mul(p,add(yy,scale(xx,-1))),scale(mul(mul(q,x),y),2)),add(mul(x,yy),scale(mul(x,xx),-F(1,3))))
    target=add(scale(mul(p,yy),F(2,3)),scale(mul(x,yy),F(5,9)))
    residual=add(mul(add(p,scale(x,F(1,3))),add(scale(yy,F(1,3)),scale(xx,-1))),mul(x,add(scale(yy,F(1,3)),scale(mul(q,y),2))))
    assert L==add(target,residual)
    # Right-circle convexity: the second p derivative is exactly8p.
    rr,z=var(1),var(2)
    g=mul(p,add(scale(mul(p,p),F(4,3)),add(scale(mul(rr,rr),-1),scale(z,-1))))
    assert d(d(g,0),0)==scale(p,8)
    return 2

def guards():
    rho,rc=F(3,8),F(41,100)
    rows=[]
    def ck(name,a,b,relation='<'):
        a,b=F(a),F(b);assert (a<b if relation=='<' else a==b),(name,a,b)
        rows.append({'name':name,'left':str(a),'right':str(b),'relation':relation})
    ck('sqrt3 upper',3,F(9,5)**2)
    ck('sqrt3 lower',F(17,10)**2,3)
    ck('q modulus',F(1,25000),F(7,1000)**2)
    ck('s modulus',F(1,8),F(177,500)**2)
    ck('vertical cannot exit',F(1,8)+48*F(1,25000),rho**2)
    ck('large p y bound',F(5,6)*(rho-F(177,500)),F(7,400),'=')
    ck('large p loss',F(1,20000),F(2,3)*F(1,4)*F(7,400)**2)
    ck('small p s bound',F(1,16)+F(1,25000),F(251,1000)**2)
    ck('small p y bound',F(5,6)*(rho-F(251,1000)),F(31,300),'=')
    ck('small p x bound',(F(31,300)-F(21,500))*F(5,9),F(23,675),'=')
    ck('small p loss',F(1,5000),F(5,9)*F(23,675)*F(31,300)**2)
    ck('arc sqrt bound',4,F(6,5)**2*3)
    ck('whole Gaussian constant',(F(6,5))**2*12*4,81)
    ck('first weighted Gaussian coefficient',F(432,25),18)
    ck('large norm zero weight coefficient',F(16,3)*9,200)
    ck('large norm first weight coefficient',48+96,200)
    ck('compact cubic loss',F(1,16),F(5,18)*F(1,4))
    ck('compact bend start',12,24)
    ck('compact outer cubic dominates linear',16,24**2)
    ck('compact zero weight integral',F(6,5)*(48+2),60,'=')
    ck('compact first weight integral',F(6,5)*(48*(2+F(6,5)*24)+2*(2+F(6,5)*25)),1900)
    ck('compact normalized bound',1900*225*36000,2*10**10)
    ck('requested slack',2*10**10,10**15)
    ck('positive cut y',(rc-F(177,500))*F(5,6),F(7,150),'=')
    ck('positive cut imaginary sign',F(7,1000),F(7,150))
    ck('cos half angle lower',(F(24,25))**2,(2+F(17,10))/4)
    ck('sin half angle upper',(2-F(17,10))/4,F(1,3)**2)
    ck('sqrt small q',F(1,1000),F(4,125)**2)
    ck('negative right X large q',F(5,12)*(rc-F(177,500)),F(7,300),'=')
    ck('negative right re large q',F(7,300)*F(24,25)-rc/F(50),F(71,5000),'=')
    ck('negative right re small q',F(71,5000),F(5,12)*(rc-F(4,125))*F(24,25)-rc/3)
    ck('negative right imaginary sign',(F(5,6)*F(24,25)-F(1,3))*rc,F(7,15)*rc,'=')
    ck('negative incoming X',F(19,50)**2,rc**2-F(9,512))
    ck('negative incoming real part',F(8,25),F(19,50)*F(24,25)-F(531,4000)/3)
    ck('complex circular perturbation',2*rc*F(2,10**5),F(41,2500000),'=')
    ck('circular connector remaining margin',F(1,30000),F(1,20000)-F(41,2500000))
    ck('left-circle derivative',4*F(8,25)**2-rc**2-F(1,8),F(233,2000),'=')
    ck('left-circle remaining margin',F(1,10000),F(1,7000)-F(41,2500000))
    return rows

def negatives():
    # A false stronger large-p outer margin is actually above its proof bound.
    assert F(2,3)*F(1,4)*F(7,400)**2<F(1,19000)
    # The old unconstrained real-part rectangle cannot use the thin-q bound.
    p,q=F(1,100),F(3,100)
    assert q*q>F(1,25000) and 2*p*q>F(2,10**5)
    # Two small endpoint margins cannot absorb an arbitrary imaginary part.
    assert F(1,20000)-2*F(41,100)*F(1,1000)<0
    return 3

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'positive_airy_thin_outer_tail_verification.json');args=ap.parse_args()
    out={'schema':1,'scope':'thin canonical positive-Airy outer tails and conditional right-circle connector geometry only',
         'proof_sha256':sha256((HERE/'POSITIVE_AIRY_THIN_OUTER_TAIL_AND_CUT_GEOMETRY.md').read_bytes()).hexdigest(),
         'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'universal_identities':identities(),'rational_guards':guards(),'negative_controls':negatives(),
         'formal_verification':False,'external_review':False,'independent_replay':False}
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','identities':out['universal_identities'],'rational_guards':len(out['rational_guards']),'negative_controls':out['negative_controls']}))
if __name__=='__main__':main()
