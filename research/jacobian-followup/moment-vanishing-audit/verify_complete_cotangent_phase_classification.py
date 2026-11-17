#!/usr/bin/env python3
"""Exact support classifier from a written PSMT1 application;finite controls are not its proof."""
from pathlib import Path
import hashlib,json
import sympy as S
from verify_joint_safe_algebras import joint_moment

def parse_exact(H,s,t):
    H=S.sympify(H)
    if H.atoms(S.Float):raise ValueError('exact Q(i) coefficients required')
    try:p=S.Poly(H,s,t,domain=S.QQ_I)
    except Exception as error:raise ValueError('Q(i) coefficients with decided exact relations required') from error
    if p.nth(0,0)!=0:raise ValueError('zero constant coefficient required for this pure family')
    return p

def classify_phase(H,s,t):
    p=parse_exact(H,s,t);terms=[(a,c) for a,c in p.terms() if c!=0];charges=sorted({a[0]-a[1] for a,c in terms})
    if not charges:return dict(safe=True,pure_null=True,scope='this fixed cotangent algebra over Q(i)',charges=[],sufficient_linear_cutoff='0',bad_marker=None)
    safe=min(charges)>0 or max(charges)<0
    return dict(safe=safe,pure_null=True,scope='this fixed cotangent algebra over Q(i)',charges=charges,
        sufficient_linear_cutoff=str(S.Rational(2,min(abs(q) for q in charges))) if safe else None,
        bad_marker=None if safe else 'x*y',proof='strict source charge separation' if safe else 'quadratic marked one-radial functional and normalized prime separation',
        failure_subsequence=None if safe else 'k*p for every sufficiently large prime p;parameters not computed')

def marked_moment(H,s,t,m,max_total_degree=24):
    p=parse_exact(H,s,t)
    if not isinstance(m,int) or m<1:raise ValueError('positive integer moment index required')
    if p.is_zero:return S.S.Zero
    if p.total_degree()*m>max_total_degree:raise ValueError('configured exact coefficient limit;classifier remains separate')
    return S.expand(sum(c*S.factorial(2*a[0])/2**a[0] for a,c in S.Poly(S.expand(p.as_expr()**m),s,t).terms() if a[0]==a[1] and a[0]>0))

def safe_ideal(d):
    if not isinstance(d,int) or d<1 or d>4:raise ValueError('configured coefficient-degree gate')
    coeff={(i,j):S.Symbol('c_'+str(i)+'_'+str(j)) for i in range(d+1) for j in range(d+1-i) if i+j>0}
    plus=[v for (i,j),v in coeff.items() if i>j];minus=[v for (i,j),v in coeff.items() if i<j];neutral=[v for (i,j),v in coeff.items() if i==j]
    return coeff,neutral+[a*b for a in plus for b in minus],minus+neutral,plus+neutral

def vp(q,p):
    a,b=S.fraction(S.cancel(q));a=abs(int(a));b=abs(int(b));out=0
    while a and a%p==0:a//=p;out+=1
    while b%p==0:b//=p;out-=1
    return out

def main():
    here=Path(__file__).parent;s,t=S.symbols('s t');examples=[]
    for H,expected in [(0,True),(s,True),(t,True),(s*(1+s*t),True),(s**3+s**4*t,True),(s+t,False),(s*t,False),(s**3+t**3,False),((s*t)**2-6*s*t,False),((1+S.I)*s+(1-S.I)*t,False)]:
        q=classify_phase(H,s,t);assert q['safe']==expected
        examples.append(dict(H=str(H),classification=q))
    assert examples[4]['classification']['sufficient_linear_cutoff']=='2/3'
    delayed=(s*t)**2-6*s*t;assert marked_moment(delayed,s,t,1)==0 and marked_moment(delayed,s,t,2)==1656
    # Independent source Wick controls at small generator degree only.
    old=json.loads((here/'gaussian_keller_endomorphism_verification.json').read_text());x,y,z=S.symbols('x y z');pos=(x,y,z)
    inv=S.Matrix([[S.sympify(v) for v in row] for row in old['inverse_Jacobian']]);vs=inv[:,0];vt=-inv[:,1];cov=[1,-1,-S.Rational(1,2)];source=[]
    for H,m in [(s+t,2),(s*t,1),(s*s+t*t,2),(s*s*t,1)]:
        total=0
        for (i,j),c in S.Poly(S.expand(H**m),s,t).terms():
            total+=c*joint_moment(vs,vt,i,j,(1,1,0),(0,0,0),pos,cov)
        expected=marked_moment(H,s,t,m);assert S.expand(total-expected)==0
        source.append(dict(H=str(H),m=m,marked=str(total)))
    # A few rational normalized PS checks;the written product argument handles ALL odd primes.
    a=lambda n:S.factorial(2*n)/2**n;ps=[]
    for prime in [3,5]:
        minimum=100
        for n in [0,1]:
            for degree in range(n*prime,3*prime+1):assert vp(a(degree)/a(n*prime),prime)>=0
            for j in range(n+1,4):minimum=min(minimum,vp(a(j*prime)/a(n*prime),prime))
        assert minimum>=1;ps.append(dict(prime=prime,minimum_endpoint_valuation=minimum))
    ideals=[];aux=S.Symbol('_intersection')
    for d in [1,2]:
        coeff,J,Ipos,Ineg=safe_ideal(d);variables=list(coeff.values())
        elimination=S.groebner([aux*q for q in Ipos]+[(1-aux)*q for q in Ineg],aux,*variables)
        actual=[p.as_expr() for p in elimination.polys if not p.as_expr().has(aux)]
        expected=[p.as_expr() for p in S.groebner(J,*variables).polys]
        assert [p.as_expr() for p in S.groebner(actual,*variables).polys]==expected
        ideals.append(dict(degree=d,coefficient_count=len(coeff),safe_generators=[str(q) for q in expected]))
    # Count/charge nilpotent bound controls without expanding source phase powers.
    bounds=[]
    for d,L,delta in [(1,1,1),(1,2,1),(2,2,1),(3,2,2)]:
        m=2*delta+(d+1)*(L-1)+1;checked=0
        for n0 in range(m+1):
            for np in range(m-n0+1):
                nm=m-n0-np
                if n0+min(np,nm)>=L:continue
                assert max(np,nm)-d*min(np,nm)>2*delta;checked+=1
        bounds.append(dict(d=d,L=L,marker_degree=delta,test_index=m,nonzero_coefficient_count_patterns=checked,sufficient_slope=2+(d+1)*(L-1)))
    rejects=0
    for H in [1+s,S.Float('1.0')*s,S.Symbol('undecided')*s]:
        try:classify_phase(H,s,t)
        except ValueError:rejects+=1
    try:marked_moment(s+t,s,t,25)
    except ValueError:rejects+=1
    assert rejects==4
    src=here/'sources/long_torus_719ca23022821868320fc0a5b9824c051105ff7e.tex';raw=src.read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest();assert blob=='719ca23022821868320fc0a5b9824c051105ff7e'
    receipt=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sympy_version=S.__version__,
        exact_field_classifier_controls=examples,independent_source_marker_controls=source,balanced_first_zero_control=dict(H=str(delayed),first='0',second='1656'),
        finite_normalized_prime_controls=ps,ordinary_safe_ideal_controls=ideals,nilpotent_charge_bound_controls=bounds,guard_refusals=rejects,
        primary_git_blob_sha=blob,primary_bytes=len(raw),primary_sha256=hashlib.sha256(raw).hexdigest(),
        fixed_universal_bad_marker='x*y',field_common_slope='2,sharp for this entire fixed algebra',
        exact_specific_all_degree_field_classifier_implemented=True,guarded_marked_moment_extractor_implemented=True,explicit_safe_monomial_ideal_generator_implemented=True,
        supporting_checker_hashes={name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in ['verify_joint_safe_algebras.py','verify_weighted_cotangent_inverse.py','verify_gaussian_keller_endomorphisms.py','verify_polynomial_gaussian_symmetries.py','verify_actual_joint_factorial.py']},
        general_Gaussian_classifier_implemented=False,general_ordinary_ideal_or_Hom_solver_implemented=False,universal_theorem_verified_by_finite_controls=False,
        prime_nonvanishing_verified_by_finite_sweep=False,new_Gaussian_moment_preserving_map_constructed=False,minimum_GMC_dimension_or_degree_improved=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    (here/'complete_cotangent_phase_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['completed','checker_sha256','fixed_universal_bad_marker','field_common_slope','independent_source_marker_controls','balanced_first_zero_control']}))

if __name__=='__main__':main()
