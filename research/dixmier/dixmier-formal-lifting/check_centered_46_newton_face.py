#!/usr/bin/env python3
"""Exact necessary Newton-face controls for the centered order-(4,6) rectangle."""
from fractions import Fraction as Q
from pathlib import Path
import argparse, hashlib, json, platform, tempfile

HERE = Path(__file__).resolve().parent
NVAR = 9  # x,xi,A,B,C,D,E,K,t
ZERO = (0,) * NVAR
ONE = {ZERO: Q(1)}

def var(index):
    return {tuple(int(i == index) for i in range(NVAR)): Q(1)}
def clean(poly):
    return {e:c for e,c in poly.items() if c}
def add(*polys):
    out={}
    for poly in polys:
        for e,c in poly.items():out[e]=out.get(e,Q(0))+c
    return clean(out)
def scale(poly,c):return clean({e:a*c for e,a in poly.items()})
def mul(a,b):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            ef=tuple(i+j for i,j in zip(e,f));out[ef]=out.get(ef,Q(0))+c*d
    return clean(out)
def power(poly,n):
    out=ONE
    for _ in range(n):out=mul(out,poly)
    return out
def deriv(poly,index):
    out={}
    for e,c in poly.items():
        if e[index]:
            f=list(e);f[index]-=1;out[tuple(f)]=c*e[index]
    return out
def poisson(q,p):return add(mul(deriv(q,1),deriv(p,0)),scale(mul(deriv(q,0),deriv(p,1)),-1))
def mono(i,j):
    e=[0]*NVAR;e[0]=i;e[1]=j;return {tuple(e):Q(1)}
def xycoef(poly,i,j):return {tuple(0 if k<2 else n for k,n in enumerate(e)):c for e,c in poly.items() if e[:2]==(i,j)}
def substitute(poly,index,value):
    out={}
    for e,c in poly.items():
        rest=list(e);rest[index]=0
        out=add(out,mul({tuple(rest):c},power(value,e[index])))
    return out
def face(support,r):
    weight=max(Q(i)-r*j for i,j in support)
    return weight,{(i,j) for i,j in support if Q(i)-r*j==weight}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path);args=parser.parse_args()
    x,y,A,B,C,D,E,K,t=[var(i) for i in range(NVAR)]
    ps={(8,4)}|{(i,j) for i in range(8) for j in range(4)}
    qs={(12,6)}|{(i,j) for i in range(12) for j in range(6)}
    assert face(ps,Q(1,4))==(Q(7),{(8,4),(7,0)})
    assert face(qs,Q(1,4))==(Q(11),{(11,0)})
    pf=add(mono(8,4),mul(A,mono(7,0)))
    assert xycoef(poisson(mono(11,0),pf),18,3)==scale(ONE,-44)
    qs.remove((11,0));assert face(qs,Q(1,4))==(Q(43,4),{(11,1)})
    assert xycoef(poisson(mono(11,1),pf),18,4)==scale(ONE,-36)
    qs.remove((11,1));assert face(qs,Q(1,4))==(Q(21,2),{(12,6),(11,2)})
    qf=add(mono(12,6),mul(B,mono(11,2)))
    assert poisson(qf,pf)==add(mul(scale(add(scale(A,3),scale(B,-2)),14),mono(18,5)),mul(scale(mul(A,B),14),mono(17,1)))
    ps.remove((7,0));qs.remove((11,2))
    assert face(ps,Q(1,3))==(Q(20,3),{(8,4),(7,1)})
    assert face(qs,Q(1,3))==(Q(10),{(12,6),(11,3),(10,0)})
    pf=add(mono(8,4),mul(A,mono(7,1)))
    qf=add(mono(12,6),mul(B,mono(11,3)),mul(C,mono(10,0)))
    assert poisson(qf,pf)==add(mul(scale(add(scale(A,3),scale(B,-2)),10),mono(18,6)),mul(scale(add(mul(A,B),scale(C,-4)),10),mono(17,3)),mul(scale(mul(A,C),-10),mono(16,0)))
    ps.remove((7,1));qs.remove((11,3));qs.remove((10,0))
    assert face(ps,Q(1,2))==(Q(6),{(8,4),(7,2),(6,0)})
    assert face(qs,Q(1,2))==(Q(19,2),{(10,1)})
    pf=add(mono(8,4),mul(A,mono(7,2)),mul(B,mono(6,0)))
    assert xycoef(poisson(mono(10,1),pf),17,4)==scale(ONE,-32)
    qs.remove((10,1));assert face(qs,Q(1,2))==(Q(9),{(12,6),(11,4),(10,2),(9,0)})
    qf=add(mono(12,6),mul(C,mono(11,4)),mul(D,mono(10,2)),mul(E,mono(9,0)))
    u=mono(1,2)
    r3=add(scale(A,3),scale(C,-2))
    r2=add(mul(A,C),scale(B,6),scale(D,-4))
    r1=add(scale(mul(B,C),4),scale(mul(A,D),-1),scale(E,-6))
    r0=add(scale(mul(B,D),2),scale(mul(A,E),-3))
    expected=scale(mul(mono(15,1),add(mul(r3,power(u,3)),mul(r2,power(u,2)),mul(r1,u),r0)),6)
    assert poisson(qf,pf)==expected
    reduced=substitute(substitute(r0,5,add(scale(power(A,2),Q(3,8)),scale(B,Q(3,2)))),6,add(scale(mul(A,B),Q(3,4)),scale(power(A,3),Q(-1,16))))
    assert reduced==scale(power(add(B,scale(power(A,2),Q(-1,4))),2),3)
    pp=mul(mono(6,0),power(add(u,t),2));qq=mul(mono(9,0),power(add(u,t),3))
    assert poisson(qq,pp)=={}
    # The next half-weight spaces after the zero-Weyl-subprincipal gate.
    ps.remove((7,3));qs.remove((11,5))
    ptop=face(ps,Q(1,2))[0];qtop=face(qs,Q(1,2))[0]
    assert face({e for e in ps if Q(e[0])-Q(e[1],2)<ptop},Q(1,2))==(Q(11,2),{(6,1)})
    assert face({e for e in qs if Q(e[0])-Q(e[1],2)<qtop},Q(1,2))==(Q(17,2),{(10,3),(9,1)})
    s=mul(mono(3,0),add(u,t));pn=mul(K,mono(6,1));qn=add(mul(B,mono(10,3)),mul(C,mono(9,1)))
    R=add(scale(mul(s,pn),3),scale(qn,-2));M=add(scale(K,3),scale(B,-2));NN=add(scale(mul(K,t),3),scale(C,-2))
    expected=add(mul(scale(M,8),mono(13,4)),mul(add(scale(NN,14),scale(mul(M,t),-9)),mono(12,2)),mul(scale(mul(NN,t),-3),mono(11,0)))
    assert poisson(s,R)==expected
    assert add(poisson(power(s,3),pn),poisson(qn,power(s,2)))==mul(s,expected)
    assert poisson(s,add(scale(mul(s,pn),3),scale(scale(mul(s,pn),Q(3,2)),-2)))=={}
    assert Q(6)+Q(17,2)-Q(1,2)==14
    assert Q(6)+Q(9)-3*Q(1,2)==Q(27,2)<14
    def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
    receipt={'status':'PASS','claim':'Necessary polynomial Newton-face and next-half-face controls for centered order4/6 rectangle; method control only','environment':{'python':platform.python_version(),'arithmetic':'Independent Fraction sparse polynomial identities; complete rectangle monomial enumeration'},'initial_rectangles':[[8,4],[12,6]],'vanishings_P':[[7,0],[7,1]],'vanishings_Q':[[11,0],[11,1],[11,2],[11,3],[10,0],[10,1]],'weights':[[-1,4],[-1,3],[-1,2]],'critical_face':'p=x6(xxi²+t)²,q=x9(xxi²+t)³','field_condition':'(B-A²/4)²=0 implies B=A²/4 over fields','next_half_face':'underzeroWeylsubprincipal,q1=(3/2)s*p1; K remainsfree','hashes':{'CENTERED_46_RECTANGLE_NEWTON_FACE_GATE.md':sha(HERE/'CENTERED_46_RECTANGLE_NEWTON_FACE_GATE.md'),'check_centered_46_newton_face.py':sha(Path(__file__))},'scope':'Necessary conditions and method control, no exclusion of order4/6 case, new pair, operator square/cube, general normalization or DC1 theorem. No claim to go beyond published degree bounds.','formal_verification':False,'external_review':False}
    out=args.output_dir or Path(tempfile.mkdtemp(prefix='centered-46-newton-'));out.mkdir(parents=True,exist_ok=True);dest=out/'centered_46_newton_face_controls.json';dest.write_text(json.dumps(receipt,indent=2)+'\n')
    print('PASS complete rectangle face enumeration, exact coefficient gates and next-half-face identities');print('Receipt:',dest)

if __name__=='__main__':main()
