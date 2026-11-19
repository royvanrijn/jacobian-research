#!/usr/bin/env python3
"""Exact finite support for the all-order hyperbolic odd-height family."""
if not __debug__:raise RuntimeError('Assertions must be enabled')
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
import argparse,json,hashlib,platform

HERE=Path(__file__).resolve().parent
PROOF=HERE/"HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md"
LEDGER=HERE/"hyperbolic_odd_height_image_source_ledger.json"

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def trim(a):
    a=list(map(Q,a))
    while len(a)>1 and not a[-1]:a.pop()
    return a
def add(a,b,scale=1):
    return trim([(a[j]if j<len(a)else Q(0))+scale*(b[j]if j<len(b)else Q(0))
                 for j in range(max(len(a),len(b)))])
def multiply(a,b,order=None):
    length=len(a)+len(b)-1
    if order is not None:length=min(length,order+1)
    out=[Q(0)]*length
    for i,v in enumerate(a):
        for j,w in enumerate(b):
            if i+j<length:out[i+j]+=v*w
    return trim(out)
def power(a,n,order=None):
    out=[Q(1)]
    for _ in range(n):out=multiply(out,a,order)
    return out
def scale(a,v):return trim([v*c for c in a])
def pochhammer(a,r):
    out=Q(1)
    for j in range(r):out*=a+j
    return out
def kappa(m,a=1):
    if m%2:return Q(0)
    return sum((Q((-1)**j*comb(m,j),a*m+2*j+1)for j in range(m+1)),Q(0))

def construct_angular(n,lam,a=1):
    assert type(n)is int and n>=1 and lam and type(a)is int and a>=1 and a%2
    lam=Q(lam);U=[Q(1),lam]
    S=[Q(comb(2*n,j+1))*lam**j for j in range(2*n)]
    baseC=power(U,n+1)
    X=scale(multiply(U,S),Q(1,2*n))
    C=multiply(baseC,power(X,a-1))
    D=scale(multiply(baseC,power(X,a+1)),-1)
    assert C[0]==1 and D[0]==-1 and X[0]==1
    assert len(C)-1==(2*a-1)*n+1 and len(X)-1==2*n and len(D)-1==(2*a+3)*n+1
    return C,D,X

def raw_angular(C,D,maximum,a=1):
    terms=[(i,0,c)for i,c in enumerate(C)if c]+[(i,2,c)for i,c in enumerate(D)if c]
    P={(0,0):Q(1)};rows=[]
    for m in range(1,maximum+1):
        out={}
        for(i,t),v in P.items():
            for j,s,w in terms:
                if i+j>maximum:continue
                key=i+j,t+s;out[key]=out.get(key,Q(0))+v*w
        P={k:v for k,v in out.items()if v}
        pure=sum((v/Q(t+a*m+1)for(i,t),v in P.items()if i==m and(t+a*m)%2==0),Q(0))
        marked=sum((v/Q(t+a*m+1)for(i,t),v in P.items()if i==m-1 and(t+a*m)%2==0),Q(0))
        rows.append((pure,marked))
    return rows

def mu(n,lam,r):
    return (-n*n*lam*lam)**r*pochhammer(Q(1,2*n),r)*pochhammer(1-Q(1,2*n),r)/(lam*pochhammer(Q(1,2),r)*factorial(r))

def badd(a,b,scale=1):
    out=dict(a)
    for p,v in b.items():out[p]=out.get(p,0)+scale*v
    return {p:v for p,v in out.items()if v}
def bscale(a,c):return {p:c*v for p,v in a.items()if c*v}
def bmul(a,b):
    out={}
    for(i,j),v in a.items():
        for(k,l),w in b.items():
            p=i+k,j+l;out[p]=out.get(p,0)+v*w
    return {p:v for p,v in out.items()if v}
def bpower(a,n):
    out={(0,0):1}
    for _ in range(n):out=bmul(out,a)
    return out
def beval(a,degree,xi1,xi2,z1,z2):
    return sum(v*xi1**i*xi2**(degree-i)*z1**j*z2**(degree-j)
               for(i,j),v in a.items())

R={(1,1):1,(0,0):1};T={(1,1):1,(0,0):-1};Z={(1,0):1};W={(0,1):2}

def construct_original(n,lam,a=1):
    assert type(lam)is int and lam and n>=1 and a>=1 and a%2
    S=[comb(2*n,j+1)*lam**j for j in range(2*n)]
    Y=[(S[j]if j<len(S)else 0)+lam*(S[j-1]if j else 0)for j in range(2*n+1)]
    square=[0]*(4*n+1)
    for i,v in enumerate(Y):
        for j,w in enumerate(Y):square[i+j]+=v*w
    assert square[0]==4*n*n and square[-1]==lam**(4*n)
    quotient={}
    for j in range(1,4*n+1):
        term=bmul(bpower(R,4*n-j),bpower(Z,j-1))
        quotient=badd(quotient,term,square[j])
    bracket=badd(bscale(bmul(W,bpower(R,4*n)),8*n*n),bmul(bpower(T,2),quotient),-1)
    Y_h={}
    for j,v in enumerate(Y):Y_h=badd(Y_h,bmul(bpower(R,2*n-j),bpower(Z,j)),v)
    divided_T=bmul(bmul(bpower(badd(R,Z,lam),n+1),bpower(Y_h,a-1)),bracket)
    return bmul(bpower(T,a),divided_T),divided_T,Y,square,Y_h

def original_controls(n,lam,maximum,a=1):
    d=(2*a+3)*n+a+2;F,divided,Y,square,Y_h=construct_original(n,lam,a)
    radial=beval(F,d,1,1,1,-1)
    assert radial==2**(a+2)*(-1)**(n+1)*lam**((2*a+3)*n+1)
    r=next(r for r in range(1,4*n+4)if (2+lam*r)and(a==1 or beval(Y_h,2*n,r,1,1,r)))
    quotient=beval(divided,d-a,r,1,1,r)
    expected_T=8*n*n*2*(2*r)**(4*n)*((2*r)+lam*r*r)**(n+1)*beval(Y_h,2*n,r,1,1,r)**(a-1)
    assert quotient==expected_T and quotient
    P={(0,0):1};rows=[]
    for m in range(1,maximum+1):
        P=bmul(P,F);degree=d*m
        pure=sum(v*factorial(i)*factorial(degree-i)for(i,j),v in P.items()if i==j)
        marked=sum(v*factorial(i+1)*factorial(degree-i)for(i,j),v in P.items()if j==i+1)
        expected=Q(0)if m%2 else(2*n)**((a+1)*m)*kappa(m,a)*mu(n,Q(lam),m//2)
        assert pure==0 and Q(marked,factorial(degree+2))==expected
        if m==2:assert expected!=kappa(m,a)*mu(n,Q(lam),1)
        rows.append({'m':m,'support':len(P),'pure':0,'normalized_original_Z_marker':str(expected)})
    # Removing the final divided-difference coefficient destroys R-primitivity.
    removed=bmul(bpower(T,a+2),bmul(bpower(badd(R,Z,lam),n+1),bmul(bpower(Y_h,a-1),bpower(Z,4*n-1))))
    altered=badd(F,removed,square[-1])
    assert beval(altered,d,1,1,1,-1)==0
    # Exact original-coordinate specialization R=T=1,W=0,Z=z.
    specialized=[sum(v for(i,j),v in F.items()if i==d and d-j==k)for k in range(d+1)]
    expected_special=scale(multiply(multiply(power([1,lam],n+1),power(Y,a-1)),square[1:]),-1)
    assert trim(specialized)==trim(expected_special)
    # The general proof uses the unique simple negative-real root of Y(U)=2n.
    assert 1-(2*n+1)-2*n==-4*n and 2**(2*n+1)-6*n-2>=0
    return {'n':n,'a':a,'lambda':lam,'balanced_degree':d,'phase_support':len(F),'Y_coefficients':Y,
            'square_coefficients':square,'radial_zero_value':radial,
            'T_quotient_control_point':[r,1,1,r],'T_quotient_value':quotient,
            'nonpower_specialization_coefficients':list(map(str,trim(specialized))),
            'negative_real_root_bracket_polynomial_at1':-4*n,
            'negative_real_root_bracket_polynomial_at2':2**(2*n+1)-6*n-2,
            'original_four_variable_controls':rows,'unscaled_marker_formula_rejected':True,
            'omitted_top_lift_coefficient_rejected':True}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-receipt',type=Path,help='Write the exact bounded verification receipt.')
    args=parser.parse_args()
    ledger=json.loads(LEDGER.read_text())
    resources=[]
    for item in ledger['retained_resources']:
        source=HERE/item['retained_path']
        actual=digest(source)
        assert actual==item['sha256'],('Retained source changed',str(source),actual)
        resources.append({'path':item['retained_path'],'sha256':actual})
    cases=[];recurrence=[]
    for n in range(1,7):
        q=Q(n-1,n);a=Q(1)
        for r in range(1,13):
            a*= -Q((2*r-1)**2-q*q,(2*r)*(2*r-1))
            product=(-1)**r*pochhammer(Q(1,2*n),r)*pochhammer(1-Q(1,2*n),r)/(pochhammer(Q(1,2),r)*factorial(r))
            assert a==product and a
        recurrence.append({'n':n,'coefficient_controls':12,'final_coefficient':str(a)})
        for a in[1,3,5]:
            for lam in[Q(1),Q(-2),Q(1,2)]:
                C,D,X=construct_angular(n,lam,a)
                rows=raw_angular(C,D,14,a)
                for m,(pure,mixed)in enumerate(rows,1):
                    expected=Q(0)if m%2 else kappa(m,a)*mu(n,lam,m//2)
                    assert pure==0 and mixed==expected
                degree=(2*a+3)*n+a+2
                assert 2*degree-1==(2*n+1)*(2*a+3)
                cases.append({'n':n,'a':a,'lambda':str(lam),'degree':degree,'moment_controls':14,
                              'marked_even_controls':7})
    C,D,X=construct_angular(2,Q(1));bad=D.copy();bad[1]+=1
    assert any(pure for pure,mixed in raw_angular(C,bad,4))
    original=[original_controls(n,lam,m,a)for n,a,lam,m in[
        (1,1,1,4),(2,1,1,6),(2,1,-2,2),(3,1,1,2),(4,1,1,2),
        (1,3,1,2),(2,3,1,2),(3,3,1,2),(4,3,1,2),(4,5,1,2),(2,3,-2,2)]]
    family_degree_controls=[]
    for d in range(1,401,2):
        parameter_pairs=[(n,a)for n in range(1,d+1)for a in range(1,d+1,2)
                         if (2*a+3)*n+a+2==d]
        factor_pairs=[(u,v)for u in range(5,2*d,4)for v in range(5,2*d,4)if u*v==2*d-1]
        recovered=[((u-1)//2,(v-3)//2)for u,v in factor_pairs]
        assert sorted(parameter_pairs)==sorted(recovered)
        if parameter_pairs:family_degree_controls.append({'degree':d,'parameters':parameter_pairs,'factors':factor_pairs})
    out={'passed':True,'schema_version':1,'date':'2026-11-19',
         'proof':{'path':PROOF.name,'sha256':digest(PROOF)},
         'checker':{'path':Path(__file__).name,'sha256':digest(Path(__file__))},
         'source_ledger':{'path':LEDGER.name,'sha256':digest(LEDGER)},
         'retained_resources':resources,'python':platform.python_version(),
         'family_degree_equivalence_controls':family_degree_controls,
         'raw_angular_cases':cases,'coefficient_recurrence_controls':recurrence,
         'original_lifts':original,'negative_controls':{'altered_D_rejected':True},
         'scope':'Exact bounded raw interval, ODE/product-coefficient, lift degree/radial/T-multiplicity and ORIGINAL four-variable factorial controls. The separate written proof supplies all parameters and all moment orders; no finite prefix is treated as purity evidence.'}
    if args.write_receipt:
        args.write_receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':True,'raw_phase_cases':len(cases),'raw_moment_controls':14*len(cases),
                      'ODE_product_controls':12*len(recurrence),'original_cases':len(original),
                      'F13_support':original[1]['phase_support'],'F13_last_support':original[1]['original_four_variable_controls'][-1]['support']}))

if __name__=='__main__':main()
