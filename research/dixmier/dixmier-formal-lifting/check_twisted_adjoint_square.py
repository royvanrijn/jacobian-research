#!/usr/bin/env python3
"""Exact controls for the twisted-adjoint and nonsymmetric-square obstructions."""
from __future__ import annotations

import argparse
from math import comb
from pathlib import Path
import hashlib
import json
import tempfile

import sympy as sp

HERE=Path(__file__).resolve().parent
V=sp.symbols('V',nonzero=True)
a0,beta,mu,eta,sigma=sp.symbols('a0 beta mu eta sigma',nonzero=True)


def clean(op):
    return {key:sp.expand(c) for key,c in op.items() if sp.expand(c)!=0}


def add(left,right,factor=1):
    return clean({key:left.get(key,0)+factor*right.get(key,0)
                  for key in set(left)|set(right)})


def scale(op,c):
    return clean({key:c*value for key,value in op.items()})


def mul(left,right,minimum=0,Eta=eta):
    out={}
    for (i,e),a in left.items():
        for (j,f),b in right.items():
            falling=sp.Integer(1)
            for k in range(max(-1,i+j-minimum)+1):
                if i>=0 and k>i:
                    break
                key=(i+j-k,e+f-2*k)
                out[key]=out.get(key,0)+sp.binomial(i,k)*falling*Eta**k*a*b
                falling*=f-2*k
                if falling==0:
                    break
    return clean(out)


def power(op,n,Eta=eta):
    out={(0,0):sp.Integer(1)}
    for _ in range(n):
        out=mul(out,op,Eta=Eta)
    return out


def comm(left,right,Eta=eta):
    return add(mul(left,right,Eta=Eta),mul(right,left,Eta=Eta),-1)


def build(A0=a0,Beta=beta,Mu=mu,Eta=eta,Sigma=sigma):
    q={(1,0):sp.Integer(1),(2,-3):Beta}
    P=add(add(scale(mul(q,q,Eta=Eta),A0**2),{(0,1):Mu}),{(1,-2):Sigma})
    return q,P


def positive_half(q,P,k,A0=a0,Beta=beta,Eta=eta):
    R=scale(q,A0)
    for j in range(1,4*k+1):
        defect=add(mul(R,R,minimum=2-j,Eta=Eta),
                   {key:c for key,c in P.items() if key[0]>=2-j},-1)
        correction={(-j,e+3):-c/(2*A0*Beta)
                    for (i,e),c in defect.items() if i==2-j}
        R=add(R,correction)
    return mul(power(P,k,Eta=Eta),R,Eta=Eta)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()

    # Universal scalar coefficient calculation for an arbitrary quadratic q.
    D=lambda f:eta*sp.diff(f,V)/V
    aa,bb,cc=[sp.Function(name)(V) for name in ['aa','bb','cc']]
    twist=sp.cancel((D(aa)-bb)/aa)
    assert sp.simplify(2*D(aa)-2*twist*aa-bb-bb)==0
    extra=D(D(aa))-2*twist*D(aa)+(twist**2-D(twist))*aa-D(bb)+twist*bb
    assert sp.simplify(extra)==0

    q,P=build()
    Astar={(1,0):-sp.Integer(1),(0,-2):-3*eta,(0,3):-1/beta}
    def adjoint(op,star=Astar):
        out={}
        for (i,e),c in op.items():
            out=add(out,mul(power(star,i),{(0,e):c}))
        return out
    assert add(adjoint(q),q,-1)=={}
    assert add(adjoint(adjoint({(1,0):sp.Integer(1)})),{(1,0):1},-1)=={}
    assert comm({(0,1):1},Astar)=={(0,-1):eta}
    Pzero=build(Sigma=sp.Integer(0))[1]
    assert add(adjoint(Pzero),Pzero,-1)=={}
    broken=add(adjoint(P),P,-1)
    expected_break={(1,-2):-2*sigma,(0,-4):-sigma*eta,(0,1):-sigma/beta}
    assert add(broken,expected_break,-1)=={}
    skew=scale(add({(2,0):sp.Integer(1)},adjoint({(2,0):sp.Integer(1)}),-1),sp.Rational(1,2))
    assert max(i for i,_ in skew)==1
    assert add(adjoint(skew),skew)=={}

    B2=positive_half(q,P,0)
    B6=positive_half(q,P,1)
    Tmu={(2,-2):sp.Rational(3,2)*beta,(1,1):sp.Rational(3,2),
         (1,-4):sp.Rational(3,2)*beta*eta,
         (0,-1):sp.Rational(3,4)*eta,(0,-6):-2*beta*eta**2}
    Tsigma={(3,-5):sp.Rational(3,2)*beta,(2,-2):sp.Rational(3,2),
            (2,-7):-sp.Rational(21,4)*beta*eta,
            (1,-4):-sp.Rational(3,2)*eta,(1,-9):7*beta*eta**2,
            (0,-6):sp.Rational(9,8)*eta**2,
            (0,-11):-sp.Rational(21,2)*beta*eta**3}
    expected=add(add(add(scale(power(q,3),a0**3),scale(Tmu,a0*mu)),
                     scale(Tsigma,a0*sigma)),{(0,-1):3*sigma**2/(8*a0*beta)})
    assert add(B6,expected,-1)=={}
    C2,C6=comm(P,B2),comm(P,B6)
    assert C2[(2,-7)]==a0*beta*eta*sigma
    assert C6[(2,-21)]==-3828*a0**3*beta**3*eta**5*sigma
    assert set(e for i,e in C2 if i==2)=={-7}
    assert sp.expand(C6[(1,-3)].subs(sigma,0)+3*a0*beta*eta*mu**2)==0
    corrupted=dict(expected);corrupted[(0,-11)]+=a0*beta*eta**3*sigma
    assert add(corrupted,B6,-1)!={}
    assert comm(P,corrupted).get((2,-21),0)!=C6[(2,-21)]

    # A positive quartic Airy pair: a square leading coefficient permits odd order.
    Pairy={(4,0):sp.Integer(1),(0,2):sp.Rational(1,2)}
    assert comm(Pairy,{(1,0):-sp.Integer(1)})=={(0,0):eta}

    # Exact symbolic-index certificate; no induction inferred from samples.
    x,n=sp.symbols('x n')
    C=lambda y:y*(y-2)
    Q=x*x/2+(21*n+23)*x/2+(245*n*n+511*n+246)/4
    rn=-(7*n+2)*(7*n+4)*(7*n+6)*(7*n+11)/4
    certificate=(C(x+7*(n+1))*C(x+7*(n+2))-rn
                 -C(x+3)*Q-C(x+7*(n+1))*Q.subs(x,x+7))
    assert sp.expand(certificate)==0

    # The unique formal series for C(x+3)F(x)+C(x)F(x-7)=x.
    f={}
    for j in range(1,14):
        known=4*f.get(j-1,0)+3*f.get(j-2,0)
        known+=sum(f[r]*comb(j-1,r-1)*7**(j-r) for r in range(1,j))
        known-=2*sum(f[r]*comb(j-2,r-1)*7**(j-1-r) for r in range(1,j))
        f[j]=sp.Rational(1,2)*((1 if j==1 else 0)-known)
    assert f[1]==sp.Rational(1,2)
    moments=[]
    previous=None
    for k in range(4):
        p=sp.prod(C(x+7*j) for j in range(1,2*k+1))
        value=sp.expand(sum(coef*f[degree[0]+1]
                            for degree,coef in sp.Poly(p,x).terms()))
        if k:
            assert value==-2*(7*k-6)*(7*k-5)*(7*k-4)*(14*k-3)*previous
        previous=value
        moments.append(str(value))
    assert moments[1]=='-66'

    # Full tenth-order numerical hierarchy control, with sigma still symbolic.
    one=sp.Integer(1)
    q10,P10=build(one,one,one,one,sigma)
    B10=positive_half(q10,P10,2,one,one,one)
    C10=comm(P10,B10,Eta=one)
    predicted=114*sp.Integer(moments[2])*sigma
    assert C10[(2,-35)]==predicted
    assert min(e for i,e in C10 if i==2)==-35

    # The integer monomial-pole square family: symbolic ell and m.
    ell,m=sp.symbols('ell m',integer=True)
    H=sp.symbols('H')
    nell=ell+4
    qell={(1,0):sp.Integer(1),(2,-ell):beta}
    starell={(1,0):-sp.Integer(1),(0,-2):-ell*eta,(0,ell):-1/beta}
    assert add(adjoint(qell,starell),qell,-1)=={}
    Hop={(1,2):sp.Integer(1)}
    expected_H={(1,2):-sp.Integer(1),(0,0):-(ell+2)*eta,
                (0,ell+2):-1/beta}
    assert add(adjoint(Hop,starell),expected_H,-1)=={}
    def canonical_coefficient(op,weight):
        return sp.expand(sum(c*sp.prod(H-2*j*eta for j in range(i))
                             for (i,e),c in op.items()
                             if sp.simplify(e-2*i-weight)==0))
    Pell=add(scale(mul(qell,qell),a0**2),{(0,1):mu})
    Cpoly2=lambda z:z*(z-2*eta)
    Cn=Cpoly2(H-nell*eta)*Cpoly2(H)
    bottom=canonical_coefficient(Pell,-2*nell)
    assert sp.expand(bottom-a0**2*beta**2*Cn)==0
    reflection_degrees=[]
    for degree in range(4):
        T=mul({(0,m):sp.Integer(1)},power(Hop,degree))
        reflected=canonical_coefficient(adjoint(T,starell),m)
        assert sp.expand(reflected-(-H-(m+ell+2)*eta)**degree)==0
        reflection_degrees.append(degree)

    character_controls=[]
    t=sp.symbols('t')
    for N in [5,6,7,8,9,10,12,15]:
        character=(1+t**2)*(1+t**N)
        assert sp.gcd(character,sp.cyclotomic_poly(N,t))==1
        assert sp.rem(character,sp.cyclotomic_poly(2*N,t),t)==0
        assert sp.rem(character-2*(1+t**2),t**N-1,t)==0
        character_controls.append(N)

    nformal=sp.symbols('nformal',integer=True,positive=True)
    Cformal=Cpoly2(H-nformal*eta)*Cpoly2(H)
    boundary_controls=[]
    for degree in range(6):
        Dpoly=H**degree if degree==0 else H**degree+3*H**(degree-1)+1
        difference=(Cformal.subs(H,H+m*eta)*Dpoly
                    -Dpoly.subs(H,H-2*nformal*eta)*Cformal)
        coefficient=sp.Poly(sp.expand(difference),H).coeff_monomial(H**(degree+3))
        assert sp.expand(coefficient-eta*(4*m+2*nformal*degree))==0
        boundary=sp.Poly(sp.expand(difference.subs(m,2*nformal)),H)
        assert boundary.degree()==degree+3
        assert sp.expand(boundary.LC()-2*nformal*eta*(degree+4))==0
        boundary_controls.append(degree)

    names=['TWISTED_ADJOINT_SQUARE_OBSTRUCTION.md','check_twisted_adjoint_square.py',
           'SQUARE_ORDERED_QUARTIC_FACE.md','GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md',
           'CRITICAL_LAURENT_FACE_OBSTRUCTION.md']
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={'claim':'QAS1','status':'PASS',
             'universal_quadratic_adjoint_coefficient_identity':True,
             'rational_mixed_twist':str(-3*eta/V**2-V**3/beta),
             'Ore_relation_and_involution_controls':True,
             'mixed_square_selfadjoint':True,
             'sigma_perturbation_breaks_symmetry':True,
             'skew_part_odd_order_control':True,
             'complete_M6_positive_part_matches':True,
             'M2_sigma_A2_coefficient':str(C2[(2,-7)]),
             'M6_sigma_decisive_A2_coefficient':str(C6[(2,-21)]),
             'M6_sigma0_decisive_A1_coefficient':str(sp.factor(C6[(1,-3)].subs(sigma,0))),
             'corrupted_M6_rejected':True,
             'positive_square_leading_Airy_control':True,
             'symbolic_index_polynomial_certificate':True,
             'sigma_moments_D0_through_D3':moments,
             'M10_decisive_A2_coefficient':str(predicted),
             'integer_monomial_pole_family':{
                 'symbolic_ell_twisted_adjoint_and_bottom_coefficient':True,
                 'symbolic_m_adjoint_reflection_degrees':reflection_degrees,
                 'primitive_n_nonzero_character_examples':character_controls,
                 'primitive_2n_character_vanishes':True,
                 'symbolic_m_n_leading_difference_degrees':boundary_controls,
                 'm_equals_2n_boundary_degree_verified':True,
                 'mu_allowed_arbitrary':True,
                 'no_full_centralizer_claim':True,
             },
             'hashes':{name:sha(HERE/name) for name in names},
             'coverage_boundary':(
                 'Exact universal quadratic-adjoint identities, finite Ore '
                 'involution and symmetry controls, full-parameter M2/M6 '
                 'positive parts and commutator witnesses, symbolic-index '
                 'polynomial certificate, formal-series moments0..3, and '
                 'a full order10 control; symbolic-ell adjoint/bottom identities, '
                 'symbolic-m reflection and m=2n leading-degree controls for '
                 'the integer monomial-pole family. Finite scalar all-order exclusion '
                 'uses the written adjoint parity or recurrence/weight proof; '
                 'no quotient-division-ring theorem is certified.'),
             'independent_replay':False,'formal_verification':False,'external_review':False}
    output=args.output or Path(tempfile.mkdtemp(prefix='qas1-check-'))/'twisted_adjoint_square_controls.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'claim':'QAS1','status':'PASS','receipt':str(output)},indent=2))


if __name__=='__main__':
    main()
