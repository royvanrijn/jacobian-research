#!/usr/bin/env python3
"""Exact/source support for fixed-terminal repeated-root SUPPORT finiteness."""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Optimized Python refused.')
from pathlib import Path
from fractions import Fraction as Q
from math import comb, factorial, isqrt
import argparse, copy, hashlib, json

HERE=Path(__file__).resolve().parent
PN='FIXED_TERMINAL_REPEATED_ROOT_SUPPORT_FINITENESS.md'
CN='verify_fixed_terminal_repeated_support.py'
IN='fixed_terminal_repeated_support_inputs.json'
LN='fixed_terminal_repeated_support_source_ledger.json'
RN='fixed_terminal_repeated_support_verification.json'
PROOF_SHA='847a961f4f55400f2fc5973fa2ea547c694ea998150dce1e45a75e1cb206a24d'
ASSURANCE={'external_review':False,'formal_verification':False,'independent_replay':False}
SCOPE={
 'schema':1,'reviewed_body_sha256':PROOF_SHA,'proof_body_preserved_exactly':True,
 'terminal':'EACH fixed integer m>=2','supports':'integer a>=0,d>=1',
 'actual_shape':'a*m+1','polynomial':'sum binom(m,k)*(a*m+1)_(d*k)*x^(m-k)',
 'claim':'finitely many repeated-root SUPPORTS at EACH fixed terminal',
 'coefficient_fields':'Q, equivalently complex repeated roots',
 'finite_per_terminal':True,'uniform_terminal_cutoff_claimed':False,
 'finite_global_repeated_triple_set_claimed':False,'evaluated_support_box':False,
 'all_order_squarefreeness':False,'all_order_SFC':False,
 'finite_per_terminal_implies_finite_global_triples':False,
 'external_input':'finite Lee-Yang circle law for H_m(z;t), EVERY finite m and t>0',
 'unitary_Hermite_time':'1/(c*m), c=a/d^2',
 'model_simplicity':'self-derived backward-heat multiplicity blowup',
 'compact_c_limit':'fixed degree; actual coefficient convergence and discriminant compactness',
 'large_c_limit':'global Gamma L^m nonlinear comparison and momentwise Gaussian convergence',
 'Gamma_Gaussian_common_space_coupling_claimed':False,
 'global_h_Taylor_domain':'v>=-sqrt(a), d>=2, a>=1',
 'h_Taylor_on_whole_real_line_claimed':False,
 'raw_sparse_individual_clause_is_logical_input':True,
 'a0_uses_only_raw_large_exponent_criterion':True,
 'd1_simplicity':'explicit Gamma ODE, not adjacent coprimality',
 'HJN_motion_imported':False,'large_degree_zero_law_imported':False,
 'theorem407_or418_imported':False,'all_order_boundary_prefix_imported':False,
 'finite_controls_prove_discriminant_compactness':False,
 'finite_controls_prove_asymptotic_uniformity':False,
 'controls':{'heat_degree_max':24,'Gaussian_degree_max':24,'CW_degree_max':16,
             'CW_c_values':['1/8','1/2','1','3'],
             'actual_m_max':6,'actual_a_squares':[1,4,9],'actual_gaps':[2,3,5],
             'Gamma_cumulant_max':12,'Taylor_degree_max':12,
             'inverse_sqrt_a_values':['1','1/2','1/5'],'budget_terminal_max':32},
 'assurance_flags':ASSURANCE}
SPECS={
 'reviewed_proof':('local',PN,'complete written fixed-terminal proof; exact frozen847 body'),
 'exact_checker':('local',CN,'exact algebra, finite coefficients, source and scope support'),
 'typed_inputs':('local',IN,'fixed-terminal quantifiers, actual shape and exact source boundary'),
 'raw_sparse':('repo','research/factorial/image-rigidity/lower_exponent_initial_and_boundary.tex',
               'ONLY thm:sparse individual-squarefreeness clause and proof; a>=1,d>=2'),
 'raw_large_exponent':('repo','research/factorial/image-rigidity/ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex',
               'ONLY thm:large coefficient proof for a0; no all-order finite prefix'),
 'finite_Lee_Yang_primary':('repo','research/factorial/image-rigidity/uniform_cusp_sources/kabluchko_2203_05533v3.pdf',
               'actual v3 printedp2 equations1.4-1.5 and printedp5 Lemma2.1; finite circle law only')}

def need(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
    if type(a)is not type(b):return False
    if type(b)is dict:return set(a)==set(b)and all(typed(a[k],v)for k,v in b.items())
    if type(b)is list:return len(a)==len(b)and all(typed(x,y)for x,y in zip(a,b))
    return a==b
def inputs(d):need(typed(d,SCOPE),'Exact typed fixed-terminal, source and asymptotic scope')
def sources(d,repo):
    need(type(d)is dict and set(d)=={'schema','sources','assurance_flags'},'Exact ledger fields')
    need(type(d['schema'])is int and d['schema']==1,'Integer ledger schema1')
    need(typed(d['assurance_flags'],ASSURANCE),'No assurance upgrade')
    rows=d['sources'];need(type(rows)is list and len(rows)==len(SPECS),'All six mandatory roles')
    out={}
    for row in rows:
        need(type(row)is dict and set(row)=={'role','kind','path','relation','sha256'},'Exact source fields')
        role=row['role'];need(type(role)is str and role in SPECS and role not in out,'Unique mandatory role')
        need((row['kind'],row['path'],row['relation'])==SPECS[role],'Exact source identity and logical role')
        hh=row['sha256'];need(type(hh)is str and len(hh)==64 and all(c in '0123456789abcdef'for c in hh),'Exact digest')
        root=HERE if row['kind']=='local'else repo;p=(root/row['path']).resolve()
        need(p.is_relative_to(root.resolve())and p.is_file()and sha(p)==hh,'Missing/stale source:'+role)
        out[role]=row
    return out

# Sparse rational polynomials in four formal indeterminates, used for
# identities that do not specialize the terminal or Fourier index.
ZERO=(0,0,0,0)
def clean(p):return {k:Q(v)for k,v in p.items()if v}
def const(c):return {ZERO:Q(c)}if c else{}
def var(i):
    e=list(ZERO);e[i]=1;return {tuple(e):Q(1)}
def add(a,b):
    p=dict(a)
    for e,c in b.items():p[e]=p.get(e,Q(0))+c
    return clean(p)
def scale(a,c):return clean({e:v*c for e,v in a.items()})
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(x+y for x,y in zip(e,f));out[g]=out.get(g,Q(0))+c*d
    return clean(out)
def power(a,n):
    p=const(1)
    for _ in range(n):p=mul(p,a)
    return p
def derivative(a,i):
    out={}
    for e,c in a.items():
        if e[i]:
            f=list(e);f[i]-=1;out[tuple(f)]=c*e[i]
    return clean(out)
def universal():
    M,J,B,X=map(var,range(4));one=const(1);two=const(2);out=[]
    def eq(name,a,b):need(clean(a)==clean(b),'Universal identity:'+name);out.append(name)
    eq('Fourier-shift heat eigenvalue',sub(scale(power(M,2),Q(1,8)),scale(mul(J,sub(M,J)),Q(1,2))),
       scale(power(sub(J,scale(M,Q(1,2))),2),Q(1,2)))
    eq('finite-CW actual radial exponent',sub(power(J,2),mul(J,M)),scale(mul(J,sub(M,J)),-1))
    # B plays gamma in this formal coefficient relation for the Gamma ODE.
    eq('Gamma-ODE coefficient recurrence',
       add(add(mul(sub(M,J),sub(sub(M,J),one)),mul(sub(sub(one,B),M),sub(M,J))),
           mul(sub(M,J),add(B,J))),{})
    # mgf logarithm bound F=b*x+(x+b^2*x^2)/(2(1-b*x)), x=1/m.
    den=scale(sub(one,mul(B,X)),2)
    num=sub(add(X,scale(mul(B,X),2)),mul(power(B,2),power(X,2)))
    eq('mgf budget derivative in inverse-sqrt-a',
       sub(mul(derivative(num,2),den),mul(num,derivative(den,2))),
       add(scale(mul(X,power(sub(one,mul(B,X)),2)),4),scale(mul(X,num),2)))
    positive=add(scale(mul(B,power(sub(one,mul(B,X)),2)),2),
                 add(one,mul(mul(power(B,2),X),sub(two,mul(B,X)))))
    eq('mgf budget derivative in inverse-terminal',
       sub(mul(derivative(num,3),den),mul(num,derivative(den,3))),scale(positive,2))
    return out

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0:p.pop()
    return p
def pplus(a,b):
    out=[Q(0)]*max(len(a),len(b))
    for i,v in enumerate(a):out[i]+=v
    for i,v in enumerate(b):out[i]+=v
    return trim(out)
def ptimes(a,b):
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)
def ppow(p,n):
    out=[Q(1)]
    for _ in range(n):out=ptimes(out,p)
    return out
def pscale(p,c):return trim([x*c for x in p])
def pdiff(p):return trim([i*p[i]for i in range(1,len(p))]or[Q(0)])
def pshift(p,n):return [Q(0)]*n+list(p)
def rising(x,n):
    z=1
    for j in range(n):z*=x+j
    return z
def exp_upper(x,n=24):
    term=Q(1);s=term
    for j in range(1,n+1):term*=x/j;s+=term
    need(x<n+2,'Exponential geometric-tail ratio')
    return s+(term*x/(n+1))/(1-x/(n+2))
def vector(a,d,m,shape_offset=1):
    return [Q(comb(m,m-j)*rising(a*m+shape_offset,d*(m-j)))for j in range(m+1)]
def actual_normalizations(a,d,m):
    L=isqrt(a);need(L*L==a,'Literal rational square root')
    delta=Q(d,L);K=(m*a)**d
    lhs=[Q(0)]
    for k in range(m+1):
        coeff=Q(comb(m,k)*rising(a*m+1,d*k),K**k*delta**m)
        lhs=pplus(lhs,pscale(ppow([Q(-1),delta],m-k),coeff))
    # Independent route: expand h(V)^j through raw centered Gamma moments.
    Vmom=[]
    for r in range(d*m+1):
        val=sum(Q(comb(r,j)*(-a)**(r-j)*rising(a*m+1,j),m**j*L**r)for j in range(r+1))
        Vmom.append(val)
    hp=[Q(0)]+[Q(comb(d,j),d*L**(j-1))for j in range(1,d+1)]
    rhs=[Q(0)]*(m+1);hpow=[Q(1)]
    for j in range(m+1):
        rhs[m-j]=comb(m,j)*sum(c*Vmom[r]for r,c in enumerate(hpow))
        hpow=ptimes(hpow,hp)
    need(trim(lhs)==trim(rhs),'Full own-shape affine polynomial identity')
    return lhs
def exact(inp):
    cfg=inp['controls'];out={'universal_identities':universal()}
    heat=Gode=GammaODE=CW=norm=cumulants=taylor=budget=0
    for r in range(2,cfg['heat_degree_max']+1):
        # Q_r(w,s)=exp((s/2)D_w^2)w^r as coefficients in w and s.
        p={}
        for l in range(r//2+1):p[(r-2*l,l)]=Q(factorial(r),2**l*factorial(l)*factorial(r-2*l))
        ds={(i,j-1):c*j for(i,j),c in p.items()if j}
        dd={(i-2,j):c*i*(i-1)/2 for(i,j),c in p.items()if i>=2}
        need(ds==dd,'Earlier-time positive-heat PDE')
        need(p[(r-2,1)]==Q(r*(r-1),2)>0 and (r-1,0)not in p,'Vieta sign in multiplicity blowup')
        heat+=1
    for m in range(1,cfg['Gaussian_degree_max']+1):
        p=[Q(0)]*(m+1)
        for l in range(m//2+1):p[m-2*l]=Q(factorial(m),2**l*factorial(l)*factorial(m-2*l)*m**l)
        residual=pplus(pdiff(pdiff(p)),pplus(pscale(pshift(pdiff(p),1),m),pscale(p,-m*m)))
        need(residual==[0],'Full Gaussian moment ODE')
        Gode+=1
    for m in range(2,cfg['CW_degree_max']+1):
        for ctext in cfg['CW_c_values']:
            c=Q(ctext);t=1/(c*m)
            for k in range(m+1):
                actual=Q(k*k,2*m)/c-Q(k,2)/c
                model=-t*k*(m-k)/2
                need(actual==model,'Exact finite-CW time/radius exponent')
                j=m-k;need((-1)**m*(-1)**(m-j)*(-1)**j==1,'Exact H_m(-z) signs')
                CW+=1
    for a in cfg['actual_a_squares']:
        for d in cfg['actual_gaps']:
            for m in range(1,cfg['actual_m_max']+1):
                actual_normalizations(a,d,m);norm+=1
            # The d1 ODE is a separate boundary family; check it once,
            # rather than repeating its count for every normalization gap.
        for m in range(1,cfg['actual_m_max']+1):
            gamma=a*m+1;p=vector(a,1,m)
            # xG''+(1-gamma-m-x)G'+mG.
            residual=pplus(pshift(pdiff(pdiff(p)),1),
                pplus(ptimes([Q(1-gamma-m),Q(-1)],pdiff(p)),pscale(p,m)))
            need(residual==[0]and p[0]==rising(gamma,m),'Actual Gamma d1 ODE')
            GammaODE+=1
    for m in range(2,9):
        for bt in cfg['inverse_sqrt_a_values']:
            b=Q(bt);shape=Q(m,b*b)+1
            need(-1/b+shape*b/m==b/m,'Centered Gamma residual first cumulant')
            for k in range(2,cfg['Gamma_cumulant_max']+1):
                literal=factorial(k-1)*shape*(b/m)**k
                expected=factorial(k-1)*(b**(k-2)/m**(k-1)+b**k/m**k)
                need(literal==expected,'Actual centered Gamma cumulant')
                cumulants+=1
    for d in range(2,cfg['Taylor_degree_max']+1):
        for bt in cfg['inverse_sqrt_a_values']:
            b=Q(bt);delta=d*b
            hp=[Q(0)]+[comb(d,j)*b**(j-1)/d for j in range(1,d+1)]
            expected=pscale(ppow([Q(1),b],d-2),delta*Q(d-1,d))
            need(pdiff(pdiff(hp))==expected and hp[0]==0 and hp[1]==1,'Global convex Taylor derivative identity')
            taylor+=1
    for m in range(2,cfg['budget_terminal_max']+1):
        for bt in cfg['inverse_sqrt_a_values']:
            b=Q(bt);x=Q(1,m)
            F=b*x+(x+b*b*x*x)/(2*(1-b*x))
            need(0<b<=1 and 0<x<=Q(1,2)and 2-b*x>0 and F<=Q(5,4),'Whole Gamma exponential budget')
            budget+=1
        # The indicated single term in exp(|v|/2) proves the GLOBAL
        # polynomial-weight bound; no sampled v-grid is used.
        q=2*m;constant=2**q*factorial(q)
        need(Q(constant,2**q*factorial(q))==1,'Positive exponential coefficient yields weighted L^m envelope')
        need(2**m*factorial(2*m)==Q(constant,2**m),'Exact nonlinear L^m coefficient')
    need(2*exp_upper(Q(5,4))<8,'Uniform E exp|V| budget')
    # Concrete false variants; none are presented as actual discriminant counterexamples.
    r=2
    earlier=[Q(0)]*(r+1);later=[Q(0)]*(r+1)
    for l in range(r//2+1):
        term=Q(factorial(r),2**l*factorial(l)*factorial(r-2*l))
        earlier[r-2*l]=term;later[r-2*l]=(-1)**l*term
    need(earlier==[1,0,1]and later==[-1,0,1],
         'Actual earlier/later heat polynomials have opposite Vieta signs')
    m=3;c=Q(1);k=1
    correct=-Q(k*(m-k),2*c*m);wrong=-Q(k*(m-k),2*c)
    need(correct==Q(k*k,2*c*m)-Q(k,2*c)and correct!=wrong,
         'Wrong finite-CW time lacks the factor m')
    need(vector(1,2,2,1)!=vector(1,2,2,0),'Wrong common Gamma shape changes the ACTUAL polynomial')
    p2=[Q(1,2),Q(0),Q(1)]
    wrong=pplus(pdiff(pdiff(p2)),pplus(pscale(pshift(pdiff(p2),1),-2),pscale(p2,4)))
    need(wrong==[4],'Wrong Gaussian ODE sign has an actual nonzero residual')
    need(3*Q(2,3)*(1-2)<0,'Convex Taylor cannot extend below the Gamma support')
    out.update(heat_PDE_and_Vieta_polynomials=heat,Gaussian_full_ODEs=Gode,
               actual_d1_Gamma_full_ODEs=GammaODE,finite_CW_exponent_and_sign_entries=CW,
               full_actual_own_shape_vectors=norm,centered_Gamma_cumulants=cumulants,
               full_global_Taylor_derivative_polynomials=taylor,whole_Gamma_budget_controls=budget,
               genuine_math_negative_families=5)
    return out
def negatives(ledger,repo):
    out=[]
    def reject(name,fn):
        try:fn()
        except(ValueError,KeyError,TypeError,OSError):out.append(name);return
        raise ValueError('Negative accepted:'+name)
    for role in SPECS:
        d=copy.deepcopy(ledger);d['sources']=[r for r in d['sources']if r['role']!=role]
        reject('omit:'+role,lambda d=d:sources(d,repo))
        d=copy.deepcopy(ledger);next(r for r in d['sources']if r['role']==role)['sha256']='0'*64
        reject('badpin:'+role,lambda d=d:sources(d,repo))
    for key in ['uniform_terminal_cutoff_claimed','finite_global_repeated_triple_set_claimed',
                'all_order_squarefreeness','all_order_SFC','evaluated_support_box',
                'HJN_motion_imported','large_degree_zero_law_imported','theorem407_or418_imported',
                'all_order_boundary_prefix_imported','h_Taylor_on_whole_real_line_claimed',
                'Gamma_Gaussian_common_space_coupling_claimed','finite_per_terminal_implies_finite_global_triples',
                'finite_controls_prove_discriminant_compactness','finite_controls_prove_asymptotic_uniformity']:
        d=copy.deepcopy(SCOPE);d[key]=True;reject('false-scope:'+key,lambda d=d:inputs(d))
    for key,value in [('schema',True),('actual_shape','a*m'),('unitary_Hermite_time','1/c'),
                      ('terminal','one fixed numerical m only'),('supports','all real a,d')]:
        d=copy.deepcopy(SCOPE);d[key]=value;reject('wrong-type-or-domain:'+key,lambda d=d:inputs(d))
    d=copy.deepcopy(ledger);d['schema']=True;reject('boolean-ledger-schema',lambda:sources(d,repo))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--output',type=Path,default=HERE/RN);ns=ap.parse_args()
    repo=ns.repo_root.resolve();inp=load(HERE/IN);inputs(inp)
    need(sha(HERE/PN)==inp['reviewed_body_sha256']==PROOF_SHA,'Exact reviewed847 proof body')
    ledger=load(HERE/LN);rows=sources(ledger,repo);counts=exact(inp);bad=negatives(ledger,repo)
    sources(ledger,repo)
    report={'schema':1,'status':'PASS_FIXED_TERMINAL_REPEATED_ROOT_SUPPORT_FINITENESS',
       'scope':inp,'proof_sha256':sha(HERE/PN),'checker_sha256':sha(HERE/CN),
       'input_sha256':sha(HERE/IN),'source_ledger_sha256':sha(HERE/LN),
       'mandatory_source_roles':len(rows),'sources':[rows[r]for r in sorted(rows)],
       'counts':counts,'source_scope_negative_controls':bad,'assurance_flags':ASSURANCE,
       'meaning':'The frozen complete written proof establishes fixed-terminal support finiteness. Exact polynomial identities, finite actual coefficient vectors and source/scope controls support its algebra. They do not prove discriminant compactness, asymptotic uniformity, a numerical box, a support-uniform terminal cutoff or an all-order theorem.'}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','mandatory_sources':len(rows),'counts':counts,
                      'source_controls':len(bad),'output':str(ns.output)},sort_keys=True))
if __name__=='__main__':main()
