#!/usr/bin/env python3
"""All-order residue certificate and exact finite square-ordered quartic controls."""
from __future__ import annotations
import argparse
from functools import lru_cache
from pathlib import Path
import hashlib,json,platform,tempfile
import sympy as s
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    v,eta,g,S=s.symbols('V eta g s',nonzero=True)
    mu=s.symbols('mu')
    @lru_cache(None)
    def delta(f,k):
        if k==0:return f
        return s.factor(s.diff(delta(f,k-1),v)/v)
    def clean(op):return {i:s.factor(c) for i,c in op.items() if s.factor(c)!=0}
    def add(*ops):
        out={}
        for op in ops:
            for i,c in op.items():out[i]=out.get(i,0)+c
        return clean(out)
    def scale(op,c):return clean({i:c*a for i,a in op.items()})
    def mul(a,b,low=0):
        out={}
        for i,ca in a.items():
            for j,cb in b.items():
                for k in range(max(0,i+j-low+1)):
                    binom=s.binomial(i,k)
                    if binom==0:continue
                    n=i+j-k
                    out[n]=out.get(n,0)+ca*binom*eta**k*delta(cb,k)
        return clean(out)
    def comm(a,b):return add(mul(a,b),scale(mul(b,a),-1))
    def same(a,b):return add(a,scale(b,-1))=={}
    L={2:S*v**-3}
    P=add(mul(L,L),{2:g,0:mu*v})
    assert P==clean({4:S**2*v**-6,3:-6*eta*S**2*v**-8,
                     2:15*eta**2*S**2*v**-10+g,0:mu*v})
    B2=add(L,{0:g*v**3/(2*S)})
    # Independent coefficient recursion of the formal square root, only as far as needed.
    M=dict(B2)
    root_steps=[]
    for n in range(-1,-5,-1):
        target=n+2
        err=add(mul(M,M,target),scale(P,-1)).get(target,0)
        M[n]=s.factor(-err/(2*L[2]))
        assert add(mul(M,M,target),scale(P,-1)).get(target,0)==0
        root_steps.append({'root_order':n,'checked_square_order':target,'coefficient':str(M[n])})
    B6=mul(P,M)
    expected_B6={
        6:S**3*v**-9,
        5:-18*eta*S**3*v**-11,
        4:159*eta**2*S**3*v**-13+s.Rational(3,2)*g*S*v**-3,
        3:-780*eta**3*S**3*v**-15-s.Rational(9,2)*eta*g*S*v**-5,
        2:1800*eta**4*S**3*v**-17+s.Rational(21,2)*eta**2*g*S*v**-7+
          3*g**2*v**3/(8*S)+s.Rational(3,2)*mu*S*v**-2,
        1:s.Rational(3,2)*eta*mu*S*v**-4,
        0:18*eta**4*g*S*v**-11-3*eta**2*g**2*v**-1/(4*S)-
          2*eta**2*mu*S*v**-6-g**3*v**9/(16*S**3)+3*g*mu*v**4/(4*S)}
    assert same(B6,expected_B6)
    C2=comm(P,B2);C6=comm(P,B6)
    assert max(C2)==2 and max(C6)==2
    expected_C2={
        2:-3*eta**2*g*S*v**-7,
        1:12*eta**3*g*S*v**-9+3*eta*g**2*v/S-2*eta*mu*S*v**-4,
        0:36*eta**4*g*S*v**-11+3*eta**2*g**2*v**-1/(2*S)+eta**2*mu*S*v**-6}
    assert same(C2,expected_C2)
    expected_C6_order2=9*eta**2*g/(8*S*v**21)*(10208*eta**4*S**4-
            12*eta**2*g*S**2*v**10-11*g**2*v**20+4*mu*S**2*v**15)
    assert s.factor(C6[2]-expected_C6_order2)==0
    leading=s.expand(C6[2]).coeff(v,-21)
    assert s.factor(leading-11484*eta**6*g*S**3)==0
    assert s.expand(C2[2]).coeff(v,-21)==0
    coefficient_matrix=s.Matrix([[s.expand(C6[2]).coeff(v,-21),s.expand(C2[2]).coeff(v,-21)],
                                 [s.expand(C6[2]).coeff(v,-7),s.expand(C2[2]).coeff(v,-7)]])
    det=s.factor(coefficient_matrix.det())
    assert det==-34452*eta**8*g**2*S**4
    # Degree-eight polynomial subtraction retains a finite operator of order at most eight.
    P2=mul(P,P)
    assert max(P2)==8 and P2[8]==S**4*v**-12 and comm(P,P2)=={}
    assert [max(op) for op in [B2,P,B6,P2]]==[2,4,6,8]
    # Positive degeneration: when g=mu=0 the fractional powers are actual powers of L.
    zero=lambda op:clean({i:a.subs({g:0,mu:0}) for i,a in op.items()})
    assert same(zero(B2),L)
    assert same(zero(B6),mul(mul(L,L),L))
    assert zero(C2)=={} and zero(C6)=={}
    # Deliberate corruption produces forbidden higher commutator coefficients.
    bad=add(B6,{2:eta*g/v})
    bad_comm=comm(P,bad)
    assert max(bad_comm)>2 and not same(bad_comm,C6)
    # Leading rational logarithmic derivative has integral residue iff order is even.
    parity={m:bool((-s.Rational(3*m,2)).is_Integer) for m in range(1,11)}
    assert all(ok==(m%2==0) for m,ok in parity.items())
    # Universal polynomial certificate: both x and n remain indeterminates.
    x,n=s.symbols('x n')
    c=lambda y:y*(y-2)
    Qn=x*x/2+(21*n+23)*x/2+(245*n*n+511*n+246)/4
    rn=-(7*n+2)*(7*n+4)*(7*n+6)*(7*n+11)/4
    certificate=s.expand(c(x+7*(n+1))*c(x+7*(n+2))-rn-
              c(x+3)*Qn-c(x+7*(n+1))*Qn.subs(x,x+7))
    assert certificate==0
    corrupted_certificate=s.expand(c(x+7*(n+1))*c(x+7*(n+2))-rn-
              c(x+3)*(Qn+1)-c(x+7*(n+1))*(Qn.subs(x,x+7)+1))
    assert corrupted_certificate!=0
    k=s.symbols('k',integer=True,positive=True)
    recurrence_multiplier=-2*(7*k-6)*(7*k-5)*(7*k-4)*(14*k-3)
    assert s.factor(rn.subs(n,2*k-2)-recurrence_multiplier)==0
    assert s.expand(-(4*(1-14*k)-6)-(56*k+2))==0
    # Independent finite controls of the formal F equation and residue reduction.
    # They supplement the symbolic certificate, never substitute for induction.
    max_k=5
    f=[s.Rational(1,2)]
    def shifted_coefficient(order,shift):
        if order==0:return f[0]
        return sum(f[l]*s.binomial(order-1,l-1)*shift**(order-l)
                   for l in range(1,min(order,len(f)-1)+1))
    for order in range(1,4*max_k+2):
        lower_shift=sum(f[l]*s.binomial(order-1,l-1)*7**(order-l)
                        for l in range(1,order))
        previous_shift=shifted_coefficient(order-1,7)
        rhs=-2 if order==1 else 0
        f.append(s.factor((rhs-4*f[order-1]-(3*f[order-2] if order>=2 else 0)
                            -lower_shift+2*previous_shift)/2))
        actual=f[order]+4*f[order-1]+(3*f[order-2] if order>=2 else 0)+\
                      shifted_coefficient(order,7)-2*previous_shift
        assert s.factor(actual-rhs)==0
    assert f[1]==-s.Rational(3,2)
    def residue(poly,shift=0):
        poly=s.Poly(s.expand(poly),x)
        return s.factor(sum(coeff*shifted_coefficient(power[0]+1,shift)
                            for power,coeff in poly.terms()))
    def backwards(power,shift=0):
        return s.prod(c(x+shift-7*l) for l in range(power))
    Ds=[]
    Ddirect=[]
    predicted=-s.Rational(3,2)
    for index in range(max_k+1):
        pn=s.prod(c(x+7*l) for l in range(1,2*index+1))
        value=residue(pn)
        direct=sum(residue(backwards(2*index-j,3-7*j)*backwards(j),7*j)
                   for j in range(2*index+1))
        assert s.factor(value-direct)==0
        if index:
            predicted=s.factor(predicted*recurrence_multiplier.subs(k,index))
        assert value==predicted and value!=0
        Ds.append(str(value));Ddirect.append(str(direct))
    assert Ds==['-3/2','198','-7128000','2268414720000',
                 '-2920048606126080000','10552997261567530598400000']
    # The first two universal coefficients coincide with full independent Ore controls.
    assert s.factor(2*s.sympify(Ds[0])*eta**2*g*S-C2[2]*v**7)==0
    assert s.factor(58*s.sympify(Ds[1])*eta**6*g*S**3-leading)==0
    output=args.output or Path(tempfile.mkdtemp(prefix='square-quartic-'))/'square_ordered_quartic_controls.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    encode=lambda op:{str(i):str(c) for i,c in sorted(op.items(),reverse=True)}
    receipt={
        'scope':'All-order square-ordered quartic polynomial centralizer and no finite mate, with universal residue certificate and independent finite controls',
        'status':'PASS','environment':{'python':platform.python_version(),'sympy':s.__version__},
        'parameters':'eta,g,s nonzero; mu arbitrary; characteristic zero; all symbolic',
        'P':encode(P),'square_root_negative_coefficients':root_steps,
        'B2':encode(B2),'B6':encode(B6),'commutator_P_B2':encode(C2),'commutator_P_B6':encode(C6),
        'decisive_most_polar_A2_coefficient':str(leading),
        'A2_coefficient_independence_determinant':str(det),
        'rational_order_parity':parity,
        'centralizer':'k[P] in every finite differential order; normal-form and residue proof in canonical note',
        'universal_certificate':{'Q_n':str(Qn),'r_n':str(rn),'symbolic_remainder':str(certificate),
              'corrupted_identity_rejected':True,'D_recurrence_multiplier':str(recurrence_multiplier)},
        'formal_F_coefficients':[str(a) for a in f],
        'D_controls':{'k_range':[0,max_k],'functional_values':Ds,'direct_shifted_values':Ddirect},
        'global_bridge':'Every finite pole residue is polynomial in critical P; finite parameter-pole removal excludes every finite mate of the actual global square candidate',
        'positive_control':'g=mu=0 gives commuting B2=L and B6=L^3',
        'negative_control':'Corrupted B6 coefficient yields commutator order above two and is rejected',
        'hashes':{name:sha(HERE/name) for name in ['SQUARE_ORDERED_QUARTIC_FACE.md','check_square_ordered_quartic_face.py',
           'GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md','CRITICAL_LAURENT_FACE_OBSTRUCTION.md']},
        'coverage_boundary':'Universal polynomial identity checked with symbolic index; finite F/residue and full Ore controls independently support written all-order recurrence, normal-form, trace and parameter-pole arguments. No actual mate or bounded Long lift.',
        'independent_replay':False,'formal_verification':False,'external_review':False}
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print('Square-ordered quartic controls PASS')
    print('Verification output:',output)
if __name__=='__main__':main()
