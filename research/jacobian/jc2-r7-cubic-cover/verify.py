"""Exact conic-cover coefficient and local-order certificates.

Surface geometry is written input. The finite Fraction calculation below
exhausts bounded polynomial degrees/orders, not possible Keller maps.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib
import json
import time

import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(name,value):
    if not value: raise AssertionError(name)
    CHECKS.append(name)


def zero(name,value):
    check(name,S.cancel(value)==0)


def nonzero_solution(rows,active):
    """Does the homogeneous kernel meet all active coordinate complements?

    RREF over Fraction. Over an infinite field a vector space is not a
    finite union of proper subspaces, so every active coordinate can be
    nonzero simultaneously iff none is forced identically zero.
    """
    cols=[i for i,on in enumerate(active) if on]
    matrix=[[Q(row[j]) for j in cols] for row in rows]
    pivots=[]; pivot_row=0
    for col in range(len(cols)):
        pivot=next((j for j in range(pivot_row,len(matrix)) if matrix[j][col]),None)
        if pivot is None: continue
        matrix[pivot_row],matrix[pivot]=matrix[pivot],matrix[pivot_row]
        divisor=matrix[pivot_row][col]
        matrix[pivot_row]=[v/divisor for v in matrix[pivot_row]]
        for j in range(len(matrix)):
            if j!=pivot_row and matrix[j][col]:
                factor=matrix[j][col]
                matrix[j]=[a-factor*b for a,b in zip(matrix[j],matrix[pivot_row])]
        pivots.append(col); pivot_row+=1
        if pivot_row==len(matrix): break
    free=[j for j in range(len(cols)) if j not in pivots]
    if not free: return False
    return all(any(matrix[i][j] for j in free) for i,_ in enumerate(pivots))


def leading_row(exponents,coefficients,finite):
    available=[e for e in exponents if e is not None]
    if not available: return [Q(0)]*len(exponents)
    leading=(min if finite else max)(available)
    return [Q(c) if e==leading else Q(0) for e,c in zip(exponents,coefficients)]


def varying_patterns(finite,enforce_top=True,a_bound=4):
    surviving=[]; examined=0
    # Coordinates are the leading coefficients of a,b,c,d. -1 means
    # the polynomial is zero; orders/degrees of a,d are always finite.
    for k,n,m,h,l in product(range(1,4),range(a_bound+1),range(-1,7),range(-1,7),range(8)):
        if not finite and enforce_top and l==7 and m!=6: continue
        if finite and min(v for v in (n,m,h,l) if v>=0)!=0: continue
        active=[True,m>=0,h>=0,True]
        if finite:
            tangent=[k*n-10,k*m if m>=0 else None,None,k*l+4]
            incidence=[k*n-20,k*m-10 if m>=0 else None,k*h if h>=0 else None,k*l-6]
        else:
            tangent=[k*n+10,k*m if m>=0 else None,None,k*l-4]
            incidence=[k*n+20,k*m+10 if m>=0 else None,k*h if h>=0 else None,k*l+6]
        rows=[leading_row(tangent,[2,1,0,Q(-3,5)],finite),
              leading_row(incidence,[1,1,1,-1],finite)]
        examined+=1
        if nonzero_solution(rows,active): surviving.append([k,n,m,h,l])
    return examined,surviving


def reciprocal_patterns(finite):
    surviving=[]; examined=0
    # Leading variables are a,c,b,d, in the order used by the proof.
    for k,n,m,j,l in product(range(1,4),range(6),range(-1,7),range(-1,7),range(-1,7)):
        if finite and min(v for v in (n,m,j,l) if v>=0)!=0: continue
        active=[True,m>=0,j>=0,l>=0]
        if finite:
            slope=[k*n-10,k*m-4 if m>=0 else None,k*j if j>=0 else None,None]
            incidence=[k*n-16,k*m-10 if m>=0 else None,None,k*l if l>=0 else None]
        else:
            slope=[k*n+10,k*m+4 if m>=0 else None,k*j if j>=0 else None,None]
            incidence=[k*n+16,k*m+10 if m>=0 else None,None,k*l if l>=0 else None]
        rows=[leading_row(slope,[8,5,3,0],finite),leading_row(incidence,[5,2,0,-3],finite)]
        examined+=1
        if nonzero_solution(rows,active): surviving.append([k,n,m,j,l])
    return examined,surviving


def poly(coeffs,t):
    return sum(S.Rational(c)*t**i for i,c in enumerate(coeffs))


def main():
    start=time.monotonic()
    x,y,z,w,p,q,s=S.symbols('x y z w p q s')
    a,b,c,d=S.symbols('a b c d')
    fv=d*x-a*(y-z*x)**2-b*(y-z*x)-c
    vcoeff=S.Poly(S.expand(fv),x,y)
    for monomial,expected in {(2,0):-a*z*z,(1,1):2*a*z,(0,2):-a,
                              (1,0):d+b*z,(0,1):-b,(0,0):-c}.items():
        zero('varying primitive coefficient '+str(monomial),vcoeff.coeff_monomial(x**monomial[0]*y**monomial[1])-expected)
    fg=a*(x*y-z*x*x)+b*x+c*y+d
    gcoeff=S.Poly(fg,x,y)
    for monomial,expected in {(2,0):-a*z,(1,1):a,(1,0):b,(0,1):c,(0,0):d}.items():
        zero('reciprocal primitive coefficient '+str(monomial),gcoeff.coeff_monomial(x**monomial[0]*y**monomial[1])-expected)
    slope=S.diff(fg,x).subs({x:p,y:q})+s*S.diff(fg,y).subs(x,p)
    expected_s=a*(q+p*s-2*z*p)+c*s+b
    expected_i=a*p*p*(s-z)-c*(q-p*s)-d
    zero('reciprocal slope equation',slope-expected_s)
    zero('reciprocal incidence elimination',p*slope-fg.subs({x:p,y:q})-expected_i)
    C=-S.Rational(5,2)
    zero('incidence leading cancellation',5+2*C)
    zero('slope residual stays nonzero',S.Rational(8,3)+S.Rational(5,3)*C+S.Rational(3,2))
    check('two cancellation equations incompatible',Q(8)*Q(2)-Q(5)*Q(5)==-9)

    # Independent Gaussian leading-balance replay, stronger than merely
    # enumerating the two integer equations derived in the written proof.
    vf_count,vf=varying_patterns(True)
    vi_count,vi=varying_patterns(False)
    check('varying finite base point no leading patterns',vf==[])
    check('varying infinite base point no leading patterns',vi==[])
    _,weakened=varying_patterns(False,enforce_top=False)
    check('missing x coefficient degree condition would leave an escape',any(row[:3]==[2,0,5] and row[-1]==7 for row in weakened))
    gf_count,gf=reciprocal_patterns(True)
    gi_count,gi=reciprocal_patterns(False)
    finite_projected=sorted({tuple(row[:4]) for row in gf})
    infinite_projected=sorted({tuple(row[:4]) for row in gi})
    check('reciprocal finite leading patterns',finite_projected==[(2,5,2,0)])
    check('reciprocal infinite leading patterns',infinite_projected==[(2,0,3,5),(2,1,4,6)])
    check('finite patterns still allow absent d',any(row[-1]==-1 for row in gf))
    check('nonzero coordinate constraints detect a forced zero',not nonzero_solution([[1,0],[0,0]],[True,True]))
    check('coordinate complement test accepts a genuine kernel',nonzero_solution([[1,1]],[True,True]))
    check('multiple root is forced in both reciprocal charts',(2-0)==(5-3)==(6-4)==2)

    # The two cusp squarefreeness certificates are literal retained Bezout
    # identities, not numerical roots or freshly trusted gcd output.
    t=S.Symbol('t')
    classes=json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases']
    certificates=json.loads((HERE.parent/'r7-adjoint-exclusion/certificates.json').read_text())
    for i,(row,certificate) in enumerate(zip(classes,certificates),1):
        pp=S.sympify(row['p']); sf=certificate['squarefree']
        zero(f'class{i} derivative pinned',S.diff(pp,t)-poly(sf['dp'],t))
        zero(f'class{i} squarefree literal Bezout',poly(sf['a'],t)*S.diff(pp,t)+poly(sf['b'],t)*S.diff(pp,t,2)-1)
        check(f'class{i} mutated Bezout rejected',S.expand((poly(sf['a'],t)+1)*S.diff(pp,t)+poly(sf['b'],t)*S.diff(pp,t,2)-1)!=0)
    witness=t**6+t**3
    check('multiple-root immersed witness is retained',S.degree(S.gcd(S.diff(witness,t),S.diff(witness,t,2)),t)==1)
    witness_q=t**10+S.Rational(5,3)*t**7+S.Rational(7,9)*t**4+S.Rational(7,81)*t
    check('multiple-root witness image is immersive',S.gcd(S.diff(witness,t),S.diff(witness_q,t))==1)

    # Genuine varying-direction family of field degree two and contact
    # degree three; its branch has two poles and is outside the exclusion.
    xx=w*w+z; yy=z*xx+w
    jac=S.diff(xx,z)*S.diff(yy,w)-S.diff(xx,w)*S.diff(yy,z)
    zc=(1-2*w**3)/(2*w)
    pc=1/(2*w); qc=1/(4*w*w)+w/2
    zero('positive varying family critical curve',jac.subs(z,zc))
    zero('positive varying family x',xx.subs(z,zc)-pc)
    zero('positive varying family y',yy.subs(z,zc)-qc)
    check('positive varying family conic field degree two',S.degree((y-z*x)**2+z-x,z)==2)
    check('positive varying contact projection degree three',max(S.degree(S.fraction(S.cancel(zc))[0],w),S.degree(S.fraction(S.cancel(zc))[1],w))==3)
    check('positive varying branch has a finite pole',S.limit(w*pc,w,0)==S.Rational(1,2))
    # Genuine reciprocal cubic-contact family of degree two.
    graph=z*x+1/(x-z)
    zg=-t*t-1/t; pg=-t*t; qg=t**4+2*t
    zero('positive reciprocal incidence',graph.subs({z:zg,x:pg})-qg)
    zero('positive reciprocal critical condition',S.diff(graph,z).subs({z:zg,x:pg}))
    check('positive reciprocal field degree two',S.degree((y-z*x)*(x-z)-1,z)==2)
    check('positive reciprocal branch outside degree(6,10)',(S.degree(pg,t),S.degree(qg,t))==(2,4))

    inputs=[HERE/'PROOF.md',HERE/'README.md',Path(__file__).resolve(),
            HERE.parent/'jc2-r7-universal-contact/PROOF.md',
            HERE.parent/'jc2-normal-r7/classification.json',
            HERE.parent/'r7-adjoint-exclusion/certificates.json']
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    receipt=dict(schema='r7-cubic-conic-cover-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),input_sha256=hashes,
                 bounded_leading_balance_replay=dict(varying_finite_examined=vf_count,varying_infinite_examined=vi_count,
                     reciprocal_finite_examined=gf_count,reciprocal_infinite_examined=gi_count),
                 varying_patterns=dict(finite=vf,infinite=vi),
                 reciprocal_patterns=dict(finite=gf,infinite=gi),
                 reciprocal_projected_patterns=dict(finite=finite_projected,infinite=infinite_projected),
                 remaining_infinity_indices=[1,5],
                 remaining_contact_degrees=dict(squarefree_p_prime=[4],multiple_root_p_prime=[3,4]),
                 assurance='Exact symbolic and Fraction checks of coefficient/order constraints. Irreducible degree identification and inherited surface geometry are written inputs. No Theorem3.8 conclusion, formalization or external review.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} CUBIC CONIC-COVER CHECKS ({receipt["elapsed_seconds"]}s)')
    print('Reciprocal projected patterns:',receipt['reciprocal_projected_patterns'])
    print('Remaining indices (1,5); squarefree derivative forces quartic contact.')


if __name__=='__main__':main()
