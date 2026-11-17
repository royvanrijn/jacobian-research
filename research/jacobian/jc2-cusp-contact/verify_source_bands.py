#!/usr/bin/env python3
"""Selected ACTUAL contact507 source equations; exact linear algebra, no full DAG."""
from pathlib import Path
from hashlib import sha256
import json
import platform
import time
import sympy as S
from sympy.polys.matrices import DomainMatrix
H=Path(__file__).resolve().parent
RROOT=H.parents[1]
checks=[]
def check(name,c):
    if not c: raise AssertionError(name)
    checks.append(name)
def zero(name,e):check(name,S.expand(e)==0)
v,q=S.symbols('v q')
def source_spaces(n):
    # P_(15-n) constraints after earlier common-root coefficients vanish.
    m=39-n
    amin=(15-n+4)//5; amax=(90-n)//6
    aord=(72-5*n+11)//12
    bk=v**max(0,amin-2)*(v-1)**max(0,aord-4)
    ks=S.symbols('k:'+str(amax-10-S.degree(bk,v)+1))
    k=bk*sum(x*v**j for j,x in enumerate(ks))
    wmin=(m+4)//5; wmax=(507+2*m)//15
    word=(5*m-15+11)//12
    bw=v**wmin*(v-1)**word
    ws=S.symbols('w:'+str(wmax-S.degree(bw,v)+1))
    K=bw*sum(x*v**j for j,x in enumerate(ws))
    return k,K,list(ws)+list(ks)
def linear(c,d,n,k,K):
    HH=c**7*d
    return S.Poly(S.expand(5*c*S.diff(K,v)-(39-n)*S.diff(c,v)*K
       +(5-n)*k*S.diff(HH,v)-39*S.diff(k,v)*HH+8*c**7*k),v)
def main():
    start=time.monotonic()
    CC=S.Function('C')(q,v);EE=S.Function('E')(q,v);J=S.symbols('J')
    pp=q**-15*CC**3;gg=q**-25*CC**5;ww=q**-39*EE
    bracket=q**5*(S.diff(gg,v)*S.diff(ww,q)-S.diff(gg,q)*S.diff(ww,v))
    pde=(5*CC-q*S.diff(CC,q))*S.diff(EE,v)-S.diff(CC,v)*(39*EE-q*S.diff(EE,q))+J*CC**8
    check('exact two-variable Keller PDE identity',
          S.simplify(bracket+5*J*pp**4-5*q**-60*CC**4*pde)==0)
    # Reproduce the retained primitive residue calculation supporting polynomial
    # root continuation; the universal induction is stated in SOURCE_BANDS.md.
    aa,bb,ss=S.symbols('aa bb ss')
    series=S.series((1+aa*ss+bb*ss**2)**S.Rational(5,3),ss,0,6).removeO()
    E4=aa**4-9*aa*aa*bb+27*bb*bb
    E5=7*aa**4-60*aa*aa*bb+135*bb*bb
    zero('double-prime fourth residue',series.coeff(ss,4)-S.Rational(5,243)*E4)
    zero('double-prime fifth residue',series.coeff(ss,5)+aa*E5/729)
    zero('double-prime residue resultant',S.resultant(E4,E5,bb)-1701*aa**8)
    check('all first-seven primitive residues precede target',
          all(5*n<36 for n in range(1,8)))
    c=v*(v-1)**2*(v*v-3*v+3)/25
    d=v*(v-1)*(v*v-3*v+3)/3
    zero('squarefree carrier ODE',5*c*S.diff(d,v)-4*S.diff(c,v)*d+c)
    rows=[]
    for n,expected in ((1,0),(2,0),(3,1)):
        k,K,unknown=source_spaces(n)
        expr=linear(c,d,n,k,K)
        M=S.linear_eq_to_matrix(expr.all_coeffs(),unknown)[0]
        null=M.nullspace()
        check('squarefree nullity n='+str(n),len(null)==expected)
        check('squarefree exact rank n='+str(n),M.rank()==len(unknown)-expected)
        if n==3:
            kk=v*(v-1)*(5*v*v-15*v+12)/5
            KK=5*c**7*v*(13*v*v-39*v+31)
            check('squarefree complete third-band generator',
                  linear(c,d,n,kk,KK).is_zero)
            vector=null[0]
            kv=S.expand(k.subs(dict(zip(unknown,vector))))
            Kv=S.expand(K.subs(dict(zip(unknown,vector))))
            ratio=S.cancel(kv/kk)
            check('squarefree third-band kernel has displayed generator',
                  not ratio.has(v) and ratio!=0 and S.expand(Kv-ratio*KK)==0)
        rows.append(dict(n=n,shape=list(M.shape),rank=len(unknown)-expected,
                         nullity=expected,k_support=str(k),K_support=str(K)))
    # Actual squarefree quadratic return at relative order6.
    k3=v*(v-1)*(5*v*v-15*v+12)/5
    K3=5*c**7*v*(13*v*v-39*v+31)
    k6,K6,uu=source_spaces(6)
    lin6=linear(c,d,6,k6,K6)
    force6=S.Poly(S.expand(2*k3*S.diff(K3,v)-36*S.diff(k3,v)*K3+28*c**6*k3*k3),v)
    degree6=max(lin6.degree(),force6.degree())
    M6=S.Matrix([[lin6.nth(j).coeff(u) for u in uu] for j in range(degree6+1)])
    b6=S.Matrix([force6.nth(j) for j in range(degree6+1)])
    D6=DomainMatrix.from_Matrix(M6).convert_to(S.QQ)
    A6=DomainMatrix.from_Matrix(M6.row_join(b6)).convert_to(S.QQ)
    check('squarefree order6 coefficient rank22',D6.rank()==22)
    check('squarefree order6 quadratic forcing is in the image',A6.rank()==22)
    k6p=-3*(2*v-3)
    K6p=c**6*v*(v-1)*(221*v**5-1326*v**4+3027*v**3-3090*v*v+1137*v+72)
    check('squarefree order6 explicit particular solution',
          S.Poly(S.expand(linear(c,d,6,k6p,K6p).as_expr()+force6.as_expr()),v).is_zero)
    check('squarefree order6 particular source and terminal spaces',
          S.Poly(S.cancel(k6p),v).degree()<=4 and
          S.Poly(S.cancel(K6p/(v**7*(v-1)**13)),v).degree()<=18)
    tt=S.symbols('tt')
    Cjet=c+tt*q**3*k3+tt*tt*q**6*k6p
    Ejet=c**7*d+tt*q**3*K3+tt*tt*q**6*K6p
    # Coefficient identities suffice; do not expand the full degree48 expression.
    check('squarefree selected nonzero-parameter jet through order7',
          linear(c,d,3,k3,K3).is_zero and
          S.Poly(S.expand(linear(c,d,6,k6p,K6p).as_expr()+force6.as_expr()),v).is_zero)
    # Double carrier: degree-one surviving k, not a claimed exclusion.
    r=(3+S.sqrt(5))/2
    dom=S.QQ.algebraic_field(S.sqrt(5))
    def canon(f):return S.Poly(S.expand(f),v,extension=S.sqrt(5))
    RR=(v-r)**2/(25*(1-r)**2)
    c=canon(v*(v-1)**2*RR).as_expr()
    al=3*(r+1)/5
    d=canon(v*(v-1)*(v-r)*(v-al)/(r*al)).as_expr()
    check('double carrier ODE',canon(5*c*S.diff(d,v)-4*S.diff(c,v)*d+c).is_zero)
    ka=(54*r-261)/142
    ea=(-6438+2397*r)/142
    eb=(9265-3332*r)/142
    ec=(15*r-45)/2
    k1=S.expand(v*(v-1)**2*(ka*v+1))
    e1=v*v*(v-1)**3*(v-r)*(ea*v*v+eb*v+ec)
    K1=canon(c**6*e1).as_expr()
    check('double explicit first-band solution',canon(linear(c,d,1,k1,K1).as_expr()).is_zero)
    k,K,unknown=source_spaces(1)
    eq=linear(c,d,1,k,K)
    M=S.linear_eq_to_matrix(eq.all_coeffs(),unknown)[0]
    DM=DomainMatrix.from_Matrix(M).convert_to(dom)
    check('double first-band exact nullity one',DM.rank()==len(unknown)-1)
    # Nonlinear order2: k1,K1 both multiplied by t; forcing is t^2 times b.
    k2,K2,unknown=source_spaces(2)
    lin=linear(c,d,2,k2,K2)
    forcing=canon(4*k1*S.diff(K1,v)-38*S.diff(k1,v)*K1+28*c**6*k1*k1)
    degree=max(lin.degree(),forcing.degree())
    M=S.Matrix([[S.simplify(lin.nth(j).coeff(u)) for u in unknown]
                for j in range(degree+1)])
    b=S.Matrix([forcing.nth(j) for j in range(degree+1)])
    DM=DomainMatrix.from_Matrix(M).convert_to(dom)
    DA=DomainMatrix.from_Matrix(M.row_join(b)).convert_to(dom)
    rank=DM.rank()
    check('double nonlinear second-band coefficient rank17',rank==17)
    check('double quadratic forcing is in the image',DA.rank()==rank)
    reduced,pivots=DA.rref()
    rr=reduced.to_Matrix()
    solution=[S.Integer(0)]*len(unknown)
    for row,col in enumerate(pivots):
        assert col<len(unknown)
        solution[col]=-rr[row,len(unknown)]
    witness=DomainMatrix.from_Matrix(M*S.Matrix(solution)+b).convert_to(dom)
    check('double quadratic forcing explicit particular solution',witness.is_zero_matrix)
    check('double order2 retains one homogeneous coordinate',len(unknown)-rank==1)
    # Symbolically verify the order2 forcing formula without a finite sample.
    x0,x1,x2,y0,y1,y2=S.symbols('x0 x1 x2 y0 y1 y2')
    # Treat coefficients as functions of v to retain all derivatives.
    X=[S.Function('x'+str(j))(v) for j in range(3)]
    Y=[S.Function('y'+str(j))(v) for j in range(3)]
    C=sum(X[j]*q**j for j in range(3)); E=sum(Y[j]*q**j for j in range(3))
    expr=(5*C-q*S.diff(C,q))*S.diff(E,v)-S.diff(C,v)*(39*E-q*S.diff(E,q))+C**8
    coeff=S.expand(expr).coeff(q,2)
    target=5*X[0]*S.diff(Y[2],v)-37*S.diff(X[0],v)*Y[2]+3*X[2]*S.diff(Y[0],v)-39*S.diff(X[2],v)*Y[0]+8*X[0]**7*X[2]+4*X[1]*S.diff(Y[1],v)-38*S.diff(X[1],v)*Y[1]+28*X[0]**6*X[1]**2
    zero('universal quadratic coefficient formula',coeff-target)
    receipt=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
       squarefree_source_rows=rows,
       squarefree_P14_zero=True,squarefree_P13_zero=True,
       squarefree_third_band_dimension=1,
       squarefree_quadratic_return=dict(order=6,shape=list(M6.shape),rank=22,augmented_rank=22,homogeneous_dimension=2,k_particular=str(k6p),K_particular=str(K6p)),
       double_first_band_dimension=1,
       double_second_band=dict(shape=list(M.shape),rank=rank,augmented_rank=rank,
          particular_solution=list(map(str,solution)),homogeneous_dimension=1),
       files={n:sha256((H/n).read_bytes()).hexdigest() for n in
              ('SOURCE_BANDS.md','verify_source_bands.py')},
       source_hashes={n:sha256((RROOT/n).read_bytes()).hexdigest() for n in
              ('foundations/plane-jc/F2_75_125_DERIVATION.md',
               'foundations/plane-jc/cas/classify_f2_75_125_layers.py',
               'foundations/plane-jc/F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md')},
       boundary='Actual necessary contact507 source equations, with prior polynomial root continuation. Nonzero P14/P13 is excluded on squarefree carrier only. Double order2 solvability is not a full Keller point.',
       python=platform.python_version(),sympy=S.__version__,
       elapsed_seconds=round(time.monotonic()-start,4))
    (H/'source-bands-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',checks=len(checks),
       squarefree_forced_zero_bands=[14,13],double_order2='SOLVABLE',
       elapsed_seconds=receipt['elapsed_seconds'])))
if __name__=='__main__':main()

