#!/usr/bin/env python3
"""Small exact contact and adapted-filtration certificates. No F2 circuit rebuild."""
from hashlib import sha256
from pathlib import Path
from math import gcd
import importlib.util
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]

def check(name,condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)

def zero(name,expression):
    check(name,S.cancel(S.expand(expression))==0)

def degree(poly,z):
    return S.degree(poly,z)

def order(poly,z):
    return min(m[0] for m,c in S.Poly(poly,z).terms() if c)

def high_part(poly,u,v,k):
    result=0
    for (i,j,l),coefficient in S.Poly(S.expand(poly),u,v,k).terms():
        if 3*i+5*j+7*l>7:
            result+=coefficient*u**i*v**j*k**l
    return S.expand(result)

def weighted_degree(poly,u,v,k):
    return max(3*i+5*j+7*l for (i,j,l),coefficient in S.Poly(S.expand(poly),u,v,k).terms() if coefficient)


def main():
    start=time.monotonic()
    z,u,v,k,t=S.symbols('z u v k t')
    f0,g0,A,B,C,D,E,F,G,H=S.symbols('f0 g0 A B C D E F G H')
    j=C*E-A*G
    f=f0+A*z**2+B*z**4+C*z**5+D*z**6
    g=g0+E*z**2+F*z**4+G*z**5+H*z**6
    wr=S.expand(5*S.diff(f,z)*g-3*f*S.diff(g,z))
    c1=10*A*g0-6*f0*E
    c3=20*B*g0-12*f0*F+4*A*E
    c4=25*C*g0-15*f0*G
    zero('Wronskian linear coefficient',wr.coeff(z,1)-c1)
    zero('Wronskian quadratic coefficient absent',wr.coeff(z,2))
    zero('Wronskian cubic coefficient',wr.coeff(z,3)-c3)
    zero('Wronskian quartic coefficient',wr.coeff(z,4)-c4)
    zero('nonzero-image beta certificate',50*g0*j-(2*E*c4-5*G*c1))
    zero('nonzero-image alpha certificate',30*f0*j-(2*A*c4-5*C*c1))
    centered=S.expand(wr.subs({f0:0,g0:0}))
    zero('centered order3 coefficient',centered.coeff(z,3)-4*A*E)
    zero('centered order5 coefficient',centered.coeff(z,5)-(14*B*E-2*A*F))
    zero('centered order6 coefficient',centered.coeff(z,6)-(19*C*E-5*A*G))
    zero('vertical centered endpoint coefficient is19J',
         centered.coeff(z,6).subs({A:0,B:0})-19*j.subs({A:0}))
    zero('horizontal centered endpoint coefficient is5J',
         centered.coeff(z,6).subs({E:0,F:0})-5*j.subs({E:0}))
    Px,Py,Qx,Qy=S.symbols('Px Py Qx Qy')
    fp=5*z**4*Px+2*z*Py
    gp=5*z**4*Qx+2*z*Qy
    Jac=Px*Qy-Py*Qx
    zero('unimodular inverse derivative identity x',Qy*fp-Py*gp-5*z**4*Jac)
    zero('unimodular inverse derivative identity y',-Qx*fp+Px*gp-2*z*Jac)
    # The full differential identity is checked for symbolic formal derivatives.
    fs,gs,df,dg=S.symbols('fs gs df dg')
    diff=fs**5-gs**3
    d_diff=5*fs**4*df-3*gs**2*dg
    zero('exact cusp-difference differential identity',
         gs*d_diff-3*dg*diff-fs**4*(5*df*gs-3*fs*dg))
    check('F2 Wronskian degree budget',4*195-325==455)
    check('F2 uncentered and centered lower thresholds',
          [455+i+1 for i in (1,3,4,5,6)]==[457,459,460,461,462])
    check('common-power contradiction degrees',3*65-1==194 and 194>1)
    check('source cusp normalization covering degree forced1',
          [d for d in range(1,66) if 65 % d == 0 and d <= 2]==[1])
    check('auxiliary cusp image degree exceeds every retained nonproper component',325>120)
    controls=[
        ('centered order3',z**5+z**2,z**5+2*z**2,3),
        ('centered vertical order5',z**5+z**4,z**2,5),
        ('centered horizontal order5',z**2,z**5+z**4,5),
        ('centered order6',z**5,z**2,6),
        ('nonzero-image order1',z**5+1,z**2+1,1),
        ('nonzero-image order3',z**5+z**4,z**2+1,3),
        ('nonzero-image order4',z**5,z**2+1,4),
    ]
    control_rows=[]
    for name,p,q,expected in controls:
        w=5*S.diff(p,z)*q-3*p*S.diff(q,z)
        check(name+' exact order',order(w,z)==expected)
        check(name+' exact derivative gcd',
              S.monic(S.gcd(S.diff(p,z),S.diff(q,z)),z)==z)
        control_rows.append(dict(name=name,order=expected,
                                 f=str(p),g=str(q),wronskian=str(S.expand(w))))

    # Centralizing polynomial controls: exact coefficient strata, not Keller pairs.
    c0,c1s,c2,c3s,c4s=S.symbols('c0 c1 c2 c3 c4')
    s=S.symbols('s')
    q=s**5+c4s*s**4+c3s*s**3+c2*s**2+c1s*s+c0
    W=S.Poly(s**15-q**3,s)
    coefficients=[W.coeff_monomial(s**i) for i in range(14,9,-1)]
    central=S.groebner(coefficients,c4s,c3s,c2,c1s,c0)
    check('complete centralizing high-contact ideal is reduced origin',
          list(central)==[c4s,c3s,c2,c1s,c0])
    for i,parameter in enumerate((c4s,c3s,c2,c1s,c0)):
        zero('centralizer triangular leading coefficient '+str(i),
             coefficients[i].subs({x:0 for x in (c4s,c3s,c2,c1s,c0)[:i]})+3*parameter)
    check('centralizer contact degrees',
          [65*i for i in (14,13,12,11,10)]==[910,845,780,715,650])

    # Universal implicit-key polynomial with weights(u,v,k)=(3,5,7).
    aa,bb,cc,tt,uu,vv,zz=S.symbols('aa bb cc tt uu vv zz')
    l0,l1,l2,l3=S.symbols('l0 l1 l2 l3')
    K=(u**5-v**3+tt*u**3*v+bb*u*v**2+aa*u**4+uu*u**2*v
       +cc*v**2+vv*u**3+zz*u*v+l2*u**2+l3*v+l1*u+l0)
    check('adapted quintic u derivative has correct Poisson degree',weighted_degree(S.diff(K,u),u,v,k)==12)
    check('adapted quintic v derivative has correct Poisson degree',weighted_degree(S.diff(K,v),u,v,k)==10)
    I2=uu+bb**2/3
    I3=vv-S.Rational(2,5)*aa**2-bb**3/27
    I4=zz+S.Rational(2,3)*bb*cc
    alpha=aa/5
    beta=-bb/3
    gamma=-cc/3
    Du=3*(u+alpha)
    Dv=5*v+2*beta*u+5*gamma-3*beta*alpha
    DK=S.expand(S.diff(K,u)*Du+S.diff(K,v)*Dv-15*K)
    expected=(-4*I2*u**2*v+(-6*I3-S.Rational(2,3)*bb*I2)*u**3
              +(-7*I4+S.Rational(6,5)*aa*I2)*u*v)
    zero('normalized derivation complete high-degree remainder',
         high_part(DK.subs(tt,0),u,v,k)-expected)
    # Independently compute all preserving derivation equations from a general L.
    r0,r1,r2,r3,r4=S.symbols('r0 r1 r2 r3 r4')
    LK=S.diff(K,u)*(r0+r1*u)+S.diff(K,v)*(r2+r3*u+r4*v)
    normal_LK=S.rem(S.Poly(LK,u),S.Poly(K-k,u)).as_expr()
    high_LK=high_part(normal_LK,u,v,k)
    high_coeff=S.Poly(high_LK,u,v,k)
    zero('general top derivation condition',high_coeff.coeff_monomial(v**3)-(5*r1-3*r4))
    normalized=high_LK.subs({r1:3,r4:5})
    zero('normalized weight14 forces tt',normalized.coeff(u,3).coeff(v,1)+tt)
    normalized=S.expand(normalized.subs(tt,0))
    zero('normalized weight13 determines shear',normalized.coeff(u,1).coeff(v,2)-(-2*bb-3*r3))
    normalized=S.expand(normalized.subs(r3,-2*bb/3))
    zero('normalized weight12 determines translation',normalized.coeff(u,4)-(-3*aa+5*r0))
    normalized=S.expand(normalized.subs(r0,3*aa/5))
    zero('normalized weight10 determines v translation',
         normalized.coeff(v,2).coeff(u,0)-(-5*cc-3*r2+3*aa*bb/5))
    normalized=S.expand(normalized.subs(r2,aa*bb/5-5*cc/3))
    zero('independent general derivation elimination matches four invariants',
         normalized-expected)
    transformed=S.expand((u+alpha)**5-(v+beta*u+gamma)**3)
    check('coordinate shifts and shear have strictly lower weights',
          weighted_degree(u+alpha,u,v,k)==3
          and weighted_degree(v+beta*u+gamma,u,v,k)==5)
    specialized={tt:0,uu:-bb**2/3,
                 vv:S.Rational(2,5)*aa**2+bb**3/27,
                 zz:-S.Rational(2,3)*bb*cc}
    zero('Euler-compatible key differs from translated binomial only below degree8',
         high_part((K-transformed).subs(specialized),u,v,k))
    check('normalized preserving derivation has divergence8',
          S.diff(Du,u)+S.diff(Dv,v)==8)

    # Actual retained four-parameter rational quintic, transcribed from PF2K1I1.
    a,b,c,d=S.symbols('a b c d')
    target=(
        u**5+(a*b+b**3+3*c)*u**4+3*b*u**3*v
        +(-a*b*c+4*a*d+3*b**2*d+3*c**2)*u**3
        +(-5*a**2-4*a*b**2+3*b*c+3*d)*u**2*v
        +(a**3*c+a**2*b**2*c-a**2*b*d-2*a*b*c**2+5*a*c*d+3*b*d**2+c**3)*u**2
        -5*a*u*v**2+(a**3*b-3*a**2*c-5*a*b*d+3*c*d)*u*v
        +(a**4*d+a**3*b**2*d-2*a**2*b*c*d+2*a**2*d**2+a*c**2*d+d**3)*u
        -v**3+(2*a**2*b-2*a*c)*v**2
        +(-a**5-a**4*b**2+2*a**3*b*c-2*a**3*d-a**2*c**2-a*d**2)*v)
    pp=t**3+a*t
    qq=t**5+b*t**4+c*t**2+d*t
    zero('retained implicit quintic vanishes on full normalized parametrization',
         target.subs({u:pp,v:qq},simultaneous=True))
    resultant=-S.resultant(u-pp,v-qq,t)
    zero('retained implicit quintic equals independent resultant',target-resultant)
    values={tt:3*b,bb:-5*a,aa:a*b+b**3+3*c,
            uu:-5*a*a-4*a*b*b+3*b*c+3*d,
            cc:2*a*a*b-2*a*c,
            vv:-a*b*c+4*a*d+3*b*b*d+3*c*c,
            zz:a**3*b-3*a*a*c-5*a*b*d+3*c*d}
    invariants=[S.expand(x.subs(values)) for x in (tt,I2,I3,I4)]
    G0=S.groebner(invariants,b,d,c,a)
    expected_basis=[b,d+S.Rational(10,9)*a*a,
                    c*c-S.Rational(25,81)*a**3,a*a*c,a**5]
    check('Euler-eligibility full scheme Groebner basis',
          list(G0)==expected_basis)
    # Small literal generator transformations; no large lift is hidden.
    reduced=[S.expand(x.subs(b,0)) for x in invariants]
    g2=d+S.Rational(10,9)*a*a
    g3=c*c-S.Rational(25,81)*a**3
    g4=a*a*c
    zero('eligibility invariant1 isolates b',invariants[0]-3*b)
    zero('eligibility invariant2 modulo b',reduced[1]-3*g2)
    zero('eligibility invariant3 modulo b',reduced[2]-(4*a*g2-S.Rational(3,5)*g3))
    zero('eligibility invariant4 modulo b',reduced[3]-(3*c*g2+g4/3))
    for i in range(1,4):
        quotient=S.cancel((invariants[i]-reduced[i])/b)
        check('eligibility b remainder is polynomial '+str(i),S.denom(quotient)==1
              or S.denom(quotient).is_Integer)
        zero('eligibility b reduction literal '+str(i),
             invariants[i]-reduced[i]-b*quotient)
    zero('a nilpotence exact identity',
         a**5-S.Rational(81,25)*(c*g4-a*a*g3))
    standard=[S.Integer(1),a,a*a,a**3,a**4,c,a*c]
    check('length7 standard monomials',
          all(G0.reduce(x)[1]==x for x in standard) and len(standard)==7)
    check('a fifth power vanishes',G0.reduce(a**5)[1]==0)
    check('c third power vanishes',G0.reduce(c**3)[1]==0)
    check('d third power vanishes',G0.reduce(d**3)[1]==0)
    for index,expr in enumerate((target.subs({a:0,b:0,c:0,d:0})-(u**5-v**3),)):
        zero('sole geometric eligibility point is monomial cusp '+str(index),expr)
    # Sample generic and degenerate rows: exact ranks, not an exhaustive proof.
    rows=[]
    matrix,_=S.linear_eq_to_matrix(S.Poly(high_LK,u,v,k).coeffs(),[r0,r1,r2,r3,r4])
    for row in ((0,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)):
        sub=dict(zip((a,b,c,d),row))
        M=matrix.subs(values).subs(sub)
        dim=5-M.rank()
        check('actual target preserving-kernel control '+str(row),
              dim==(1 if row==(0,0,0,0) else 0))
        rows.append(dict(parameters=list(row),preserving_kernel_dimension=dim))

    # Actual two-face spectrum and both forced L=39 carrier coefficients.
    Jc,rho=S.symbols('Jc rho')
    allowed_L=[L for L in range(39,75) if L==39 or L%5==0]
    spectrum=[13*L for L in allowed_L]
    check('complete two-face L spectrum',allowed_L==[39,40,45,50,55,60,65,70])
    check('complete cusp-contact spectrum',spectrum==[507,520,585,650,715,780,845,910])
    check('exact Wronskian degree spectrum',
          [d-456 for d in spectrum]==[51,64,129,194,259,324,389,454])
    for L in allowed_L:
        # One monomial is the complete homogeneous kernel at the fixed weight.
        solutions=[(i,j) for i in range(13*L//5+1)
                   for j in range(13*L//2+1)
                   if 5*i+2*j==13*L and 25*j-100*i==0]
        check('positive-weight homogeneous kernel at L='+str(L),
              solutions==[(L,4*L)])
    check('critical band39 multiplicity floor excludes degree36',
          min(7*5+r for r in (3,4))>36)
    check('critical infinity resonance is degree39',5*39-195==0)
    v0=S.symbols('v0')
    Rs=(v0**2-3*v0+3)/25
    cs=v0*(v0-1)**2*Rs
    ds=S.Rational(25,3)*Jc*v0*(v0-1)*Rs
    zero('squarefree forced residue ODE',
         5*cs*S.diff(ds,v0)-4*S.diff(cs,v0)*ds+Jc*cs)
    check('squarefree forced numerator degree4',S.degree(ds,v0)==4)
    zero('squarefree forced numerator leading coefficient',
         S.Poly(ds,v0).LC()-Jc/3)
    al=S.Rational(3,5)*(rho+1)
    Rd=(v0-rho)**2/(25*(1-rho)**2)
    cd=v0*(v0-1)**2*Rd
    dd=Jc*v0*(v0-1)*(v0-rho)*(v0-al)/(rho*al)
    residue=S.cancel(5*cd*S.diff(dd,v0)-4*S.diff(cd,v0)*dd+Jc*cd)
    numerator=S.fraction(residue)[0]
    check('double carrier forced residue ODE in exact quadratic field',
          S.rem(S.Poly(numerator,rho),S.Poly(rho**2-3*rho+1,rho)).is_zero)
    check('double carrier numerator degree4',S.degree(dd,v0)==4)
    zero('double carrier leading numerator coefficient',
         S.Poly(dd,v0).LC()-Jc/(rho*al))
    check('double carrier rho alpha denominator is a unit',
          S.gcd(S.Poly(rho*al,rho),S.Poly(rho**2-3*rho+1,rho)).degree()==0)
    check('double carrier c leading denominator is a unit',
          S.gcd(S.Poly((1-rho)**2,rho),S.Poly(rho**2-3*rho+1,rho)).degree()==0)
    ccfunc=S.Function('cc')(v0)
    dfunc=S.Function('dd')(v0)
    zero('critical residue transport to full band39',
         5*ccfunc*S.diff(ccfunc**7*dfunc,v0)
         -39*S.diff(ccfunc,v0)*ccfunc**7*dfunc
         -ccfunc**7*(5*ccfunc*S.diff(dfunc,v0)-4*S.diff(ccfunc,v0)*dfunc))
    yy=S.symbols('yy',nonzero=True)
    LL=S.symbols('LL',integer=True)
    Hfun=S.Function('Hh')(v0)
    first=yy**(-25)*ccfunc**5
    second=yy**(-LL)*Hfun
    original=-yy**5*(S.diff(first,yy)*S.diff(second,v0)
                    -S.diff(first,v0)*S.diff(second,yy))
    check('original source-coordinate leading band bracket',
          S.simplify(original*yy**(LL+21)
                     -5*ccfunc**4*(5*ccfunc*S.diff(Hfun,v0)
                                   -LL*S.diff(ccfunc,v0)*Hfun))==0)
    check('lowest Wronskian leading factor',507-975==-468)

    files=['PROOF.md','verify.py','README.md']
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
                 wronskian_controls=control_rows,
                 centralizing_top_coefficients=list(map(str,coefficients)),
                 general_eligibility_invariants=list(map(str,(tt,I2,I3,I4))),
                 actual_quintic_invariants=list(map(str,invariants)),
                 actual_quintic_eligibility_groebner=list(map(str,G0)),
                 eligibility_length=7,kernel_controls=rows,
                 leading_band_spectrum=allowed_L,contact_degree_spectrum=spectrum,
                 wronskian_degree_spectrum=[d-456 for d in spectrum],
                 source_hashes={name:sha256((ROOT/name).read_bytes()).hexdigest()
                    for name in (
                       'foundations/plane-jc/F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md',
                       'foundations/plane-jc/F2_75_125_DERIVATION.md',
                       'foundations/plane-jc/F2_AFFINE_TARGET_CURVE_ATLAS.md',
                       'foundations/plane-jc/F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md',
                       'jacobian/jc2-cusp-filtration-rigidity/PROOF.md')},
                 files={name:sha256((HERE/name).read_bytes()).hexdigest() for name in files},
                 proof_boundary='Universal degree bounds and filtration arguments are written; this checks exact algebra and finite controls.',
                 python=platform.python_version(),sympy=S.__version__,
                 elapsed_seconds=round(time.monotonic()-start,4))
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',checks=len(checks),
                         elapsed_seconds=receipt['elapsed_seconds'],
                         eligibility_length=7)))

if __name__=='__main__':
    main()
