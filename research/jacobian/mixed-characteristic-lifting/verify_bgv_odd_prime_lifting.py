#!/usr/bin/env python3
"""Exact sparse controls for BGV odd-prime first obstruction and one Witt tower.

No third-party packages. The all-odd-prime assertion is proved in the note;
finite prime controls and finite Witt digits do not establish algebraization.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path

ZERO = {}
ONE = {(0, 0): 1}
X = {(1, 0): 1}
Y = {(0, 1): 1}

def clean(a): return {m: c for m, c in a.items() if c}
def add(a, b):
    out = dict(a)
    for m, c in b.items(): out[m] = out.get(m, 0) + c
    return clean(out)
def scale(a, c): return clean({m: v*c for m, v in a.items()})
def sub(a, b): return add(a, scale(b, -1))
def mul(a, b):
    out = {}
    for (i,j), c in a.items():
        for (k,l), d in b.items():
            m=(i+k,j+l);out[m]=out.get(m,0)+c*d
    return clean(out)
def power(a, n):
    out=ONE
    while n:
        if n & 1: out=mul(out,a)
        n//=2
        if n: a=mul(a,a)
    return out

def diff(a, axis):
    out={}
    for m,c in a.items():
        if m[axis]:
            n=list(m);n[axis]-=1;out[tuple(n)]=c*m[axis]
    return clean(out)
def jac(a,b): return sub(mul(diff(a,0),diff(b,1)),mul(diff(a,1),diff(b,0)))
def mod(a,p): return clean({m:c%p for m,c in a.items()})
def quotient(a,n):
    assert all(c%n==0 for c in a.values())
    return clean({m:c//n for m,c in a.items()})
def degree(a): return max((sum(m) for m in a),default=-1)
def serial(a): return [[i,j,c] for (i,j),c in sorted(a.items())]
def cartier(a,p):
    return clean({(i//p,j//p):c%p for (i,j),c in a.items()
                  if i%p==p-1 and j%p==p-1})
def family(p):
    z=add(X,mul(power(X,p),Y))
    P=add(z,power(mul(X,z),p))
    Q=add(Y,mul(power(X,p*(p-1)),power(z,p*(p-1)+1)))
    return z,P,Q

def primitive(h,p):
    """u dx+v dy with d(eta)=h dx dy, if Cartier(h)=0."""
    assert not cartier(h,p)
    u={};v={}
    for (i,j),c in mod(h,p).items():
        if (i+1)%p:
            v=add(v,{(i+1,j): c*pow(i+1,-1,p)})
        else:
            assert (j+1)%p
            u=add(u,{(i,j+1): -c*pow(j+1,-1,p)})
    u,v=mod(u,p),mod(v,p)
    assert mod(sub(diff(v,0),diff(u,1)),p)==mod(h,p)
    return u,v

def correction(P,Q,H,p):
    u,v=primitive(scale(H,-1),p)
    # eta=A dQ-B dP; determinant of this linear system is Jac(P,Q)=1.
    A=mod(add(scale(mul(diff(P,1),u),-1),mul(diff(P,0),v)),p)
    B=mod(add(scale(mul(diff(Q,1),u),-1),mul(diff(Q,0),v)),p)
    linear=add(jac(A,Q),jac(P,B))
    assert mod(add(H,linear),p)==ZERO
    return A,B

def mulmod(a,b,M): return mod(mul(a,b),M)
def powmod(a,n,M):
    out=ONE
    while n:
        if n&1:out=mulmod(out,a,M)
        n//=2
        if n:a=mulmod(a,a,M)
    return out

def radial_substitute(a,lam,M):
    out={};cache={0:ONE}
    for (i,j),c in a.items():
        d=i+j
        if d not in cache:cache[d]=powmod(lam,d,M)
        term={(u+i,v+j):c*t for (u,v),t in cache[d].items()}
        out=mod(add(out,term),M)
    return out

def euler_primitive(J,p,M):
    assert J.get((0,0),0)==1
    A={}
    for (i,j),c in J.items():
        d=i+j
        assert (d+2)%p
        A[i,j]=(2*c*pow(d+2,-1,M))%M
    A=clean(A)
    assert mod({m:(sum(m)+2)*c for m,c in A.items()},M)==mod(scale(J,2),M)
    assert mod(A,p)==ONE
    return A

def radial_tower(p,max_n):
    z,P,Q=family(p);J=jac(P,Q)
    assert {sum(m)%p for m in J} <= {0,p-1}
    assert all(i>0 for (i,j) in sub(J,ONE))
    lam=ONE;previous=None;rows=[]
    for n in range(2,max_n+1):
        M=p**n;A=euler_primitive(J,p,M)
        # One Hensel digit. Derivative of T^2 A(Tx,Ty)-1 is 2 mod p.
        E=mod(sub(mulmod(powmod(lam,2,M),radial_substitute(A,lam,M),M),ONE),M)
        digit=mod(scale(quotient(E,p**(n-1)),-pow(2,-1,p)),p)
        lam=mod(add(lam,scale(digit,p**(n-1))),M)
        assert mod(lam,p)==ONE
        assert all(i>0 for (i,j) in sub(lam,ONE))
        assert mod(mulmod(powmod(lam,2,M),radial_substitute(A,lam,M),M),M)==ONE
        FP=radial_substitute(P,lam,M);FQ=radial_substitute(Q,lam,M)
        assert mod(jac(FP,FQ),M)==ONE
        assert mod(FP,p)==mod(P,p) and mod(FQ,p)==mod(Q,p)
        assert {(i,j):c for (i,j),c in FP.items() if i==0}==ZERO
        assert {(i,j):c for (i,j),c in FQ.items() if i==0}==Y
        if previous:
            oldM,oldL,oldP,oldQ=previous
            assert mod(lam,oldM)==oldL
            assert mod(FP,oldM)==oldP and mod(FQ,oldM)==oldQ
        rows.append({'modulus':M,'lambda_degree':degree(lam),
                     'map_degrees':[degree(FP),degree(FQ)],
                     'map_terms':[len(FP),len(FQ)],
                     'branch_compatible':True,'jacobian_one':True,
                     'fixes_x_zero_line_pointwise':True,
                     'lambda_coefficients':serial(lam)})
        previous=M,lam,FP,FQ
    return rows

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write-receipt',type=Path)
    args=ap.parse_args()
    controls=[]
    for p in (2,3,5,7):
        z,P,Q=family(p)
        assert mod(jac(P,Q),p)==ONE
        H=mod(quotient(sub(jac(P,Q),ONE),p),p)
        # First-class simplification via the pth-power addition lemma.
        simple=add(mul(power(X,p-1),Y),
                   mul(power(X,p*p-1),power(z,p*p-p+1)))
        assert cartier(H,p)==cartier(simple,p)
        C=cartier(H,p)
        if p==2:
            assert C==add(ONE,add(power(X,3),mul(power(X,4),Y)))
        else:
            assert C==ZERO
            # Independent characteristic-zero degree controls; X=u, Y=P.
            zz=sub(Y,power(X,p))
            N=add(sub(power(zz,p+1),mul(X,power(zz,p-1))),mul(power(X,p*p),zz))
            assert max(i for i,j in N)==p*p
            assert {(i,j):c for (i,j),c in N.items() if i==p*p}=={(p*p,1):-p}
            assert {(i,j):c for (i,j),c in N.items() if i==0}=={(0,p+1):1}
            J=jac(P,Q);D=p**3+2*p*p+p-1
            assert degree(J)==D and D%2==1
            assert J[p**3+p*p+p-1,p*p]==p
        controls.append({'prime':p,'degrees':[degree(P),degree(Q)],
                         'jacobian_mod_p':1,'cartier_obstruction':serial(C),
                         'simplified_representative_checked':True})
    # Canonical digit convention: all coefficient representatives are 0,1,2;
    # integrate x first, then y only when the x exponent is 2 mod 3.
    p=3;z,P,Q=family(p);baseP,baseQ=P,Q;tower=[]
    for n in range(1,4):
        E=sub(jac(P,Q),ONE)
        H=mod(quotient(E,p**n),p);C=cartier(H,p)
        row={'input_modulus':p**n,'input_degrees':[degree(P),degree(Q)],
             'input_term_counts':[len(P),len(Q)],'error_term_count':len(H),
             'cartier_obstruction':serial(C)}
        if C:
            assert n==3
            row['conclusion']='This specific mod-27 representative has no mod-81 extension.'
            tower.append(row);break
        A,B=correction(baseP,baseQ,H,p)
        P=add(P,scale(A,p**n));Q=add(Q,scale(B,p**n))
        assert all(c%(p**(n+1))==0 for c in sub(jac(P,Q),ONE).values())
        row.update({'output_modulus':p**(n+1),'output_degrees':[degree(P),degree(Q)],
                    'correction_terms':[len(A),len(B)],
                    'digit_A':serial(A),'digit_B':serial(B),
                    'determinant_one_verified':True})
        tower.append(row)
    # Negative control: a monomial spanning top de Rham cohomology.
    assert cartier({(2,2):1},3)==ONE
    # Positive control: Artin-Schreier x-x^p, y has exact first error.
    AS=sub(X,power(X,3));H=mod(quotient(sub(jac(AS,Y),ONE),3),3)
    assert not cartier(H,3)
    A,B=correction(AS,Y,H,3)
    assert mod(jac(add(AS,scale(A,3)),add(Y,scale(B,3))),9)==ONE
    here=Path(__file__).resolve().parent
    source=here/'sources/bgv-2609.05746v1.pdf'
    report={'schema_version':1,'result':'PASS','kind':'internal_exact_replay',
            'python':platform.python_version(),'dependencies':'standard library only',
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'note_sha256':hashlib.sha256((here/'BGV_ODD_PRIME_LIFTING.md').read_bytes()).hexdigest(),
            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'prime_controls':controls,'canonical_p_three_tower':tower,
            'radial_p_three_tower':radial_tower(3,4),
            'negative_and_positive_controls':True,
            'scope':'Finite exact controls include canonical mod-9/mod-27 witnesses and a compatible radial tower modulo 9,27,81. The note proves compatible radial lifts for every odd prime and all Witt lengths, and proves their characteristic-zero limit is not polynomial. Canonical-gauge third-digit failure is for that recorded representative only. No polynomial characteristic-zero Keller counterexample or JC2 solution is claimed.',
            'formal_verification':False,'external_review':False,
            'independent_implementation':False}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.write_receipt:args.write_receipt.write_text(out)
    print(out,end='')

if __name__=='__main__': main()
