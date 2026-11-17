"""Original polynomial certificates for all-marker Gaussian safety.

This deliberately slow bounded enumerator is an independent check, not a
replacement for the agent's practical Singular integration implementation.
Uncapped classical holonomic existence implies termination. Capped failure
is UNKNOWN. A COMPLETE decision uses only original polynomial identities,
exact rational linear algebra, and finitely many original Wick moments.
"""
from pathlib import Path
import importlib.util
import sympy as sp

_spec=importlib.util.spec_from_file_location('ibp_side',Path(__file__).with_name('gaussian_ibp_annihilator_sidecheck_20261116.py'))
IBP=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(IBP)


def exponents(width,degree):
    if width==0:
        yield ()
    else:
        for a in range(degree+1):
            for tail in exponents(width-1,degree-a):yield (a,*tail)


def find_cover(phase_terms,covariance,max_degree=4,max_columns=300):
    n=len(covariance)
    C=sp.Matrix([[sp.Rational(c) for c in row] for row in covariance])
    if n<1 or C.shape!=(n,n) or C!=C.T or C.det()==0:raise ValueError('Invalid covariance')
    inv=C.inv();t=sp.Symbol('t');xs=sp.symbols('x0:'+str(n));variables=(t,*xs)
    P=sp.S.Zero
    for ex,c in phase_terms:
        if len(ex)!=n or any(not isinstance(e,int) or e<0 for e in ex):raise ValueError('Bad phase exponent')
        P+=sp.Rational(c)*sp.prod(x**a for x,a in zip(xs,ex))
    drift=[-sum(inv[i,j]*xs[j] for j in range(n))+t*sp.diff(P,xs[i]) for i in range(n)]
    theta=[sp.S.One]
    d=0
    while max_degree is None or d<=max_degree:
        while len(theta)<=d:
            B=theta[-1]
            theta.append(sp.expand(-n*B-sum(x*(sp.diff(B,x)+f*B) for x,f in zip(xs,drift))))
        remainder=[ex for ex in exponents(2*n+2,d) if sum(ex[n+1:2*n+1])>sum(ex[:n])]
        if max_columns is not None and d+1+len(remainder)>max_columns:
            return dict(outcome='UNKNOWN',reason='Euler-cover column budget exceeded',searched_degree=d)
        images=theta[:d+1]
        for ex in remainder:
            alpha,a,beta,b=ex[:n],ex[n],ex[n+1:2*n+1],ex[-1]
            B=P**b
            for i,power in enumerate(beta):
                for _ in range(power):B=sp.expand(sp.diff(B,xs[i])+drift[i]*B)
            images.append(sp.expand(t**a*sp.prod(x**p for x,p in zip(xs,alpha))*B))
        columns=[sp.Poly(B,*variables,domain=sp.QQ).as_dict() for B in images]
        keys=sorted(set().union(*(c.keys() for c in columns)))
        entries={(i,j):col[key] for j,col in enumerate(columns) for i,key in enumerate(keys) if key in col}
        matrix=sp.MutableSparseMatrix(len(keys),len(columns),entries)
        chosen=next((v for v in matrix.nullspace() if any(v[j] for j in range(d+1))),None)
        if chosen is not None:
            highest=max(j for j in range(d+1) if chosen[j])
            coeff=[sp.cancel(v/chosen[highest]) for v in chosen]
            relation=sp.expand(sum(c*B for c,B in zip(coeff,images)))
            if relation!=0:raise AssertionError('False original exponential Euler identity')
            z=sp.Symbol('z');bpoly=sp.Poly(sum(coeff[j]*z**j for j in range(d+1)),z,domain=sp.QQ)
            roots=sorted(int(a) for a in bpoly.ground_roots() if a.is_Integer and a>=0)
            if not roots or bpoly.eval(0)!=0:raise AssertionError('Gaussian unit requires b(0)=0')
            terms=[];slope=sp.S.Zero
            for j,ex in enumerate(remainder):
                c=coeff[d+1+j]
                if c:
                    drop=sum(ex[n+1:2*n+1])-sum(ex[:n]);gain=ex[n]-ex[-1]
                    slope=max(slope,sp.Rational(gain,drop))
                    terms.append(dict(powers=list(ex),coefficient=str(c),spatial_drop=drop,time_degree_gain=gain))
            return dict(outcome='COMPLETE',covariance=covariance,phase_terms=phase_terms,
                        b_coefficients=[str(bpoly.nth(j)) for j in range(bpoly.degree()+1)],
                        b_polynomial=str(bpoly.as_expr()),nonnegative_integer_roots=roots,
                        spatial_generator_cutoff=max(roots),remainder_terms=terms,
                        recurrence_degree_slope=str(slope),
                        original_exponential_identity_verified=True,
                        original_identity_residual=str(relation),search_degree=d)
        d+=1
    return dict(outcome='UNKNOWN',reason='Euler-cover degree budget exceeded',searched_degree=max_degree)


def decide(phase_terms,covariance,cover_degree_cap=4,ibp_degree_cap=5,max_low_markers=30):
    cover=find_cover(phase_terms,covariance,max_degree=cover_degree_cap)
    if cover['outcome']!='COMPLETE':return dict(outcome='UNKNOWN',cover=cover)
    n=len(covariance);L=cover['spatial_generator_cutoff']
    markers=list(exponents(n,L))
    if max_low_markers is not None and len(markers)>max_low_markers:
        return dict(outcome='UNKNOWN',reason='Low-marker budget exceeded',cover=cover)
    certificates=[];E=sp.S.Zero;C=sp.Rational(cover['recurrence_degree_slope'])
    for alpha in markers:
        cert=IBP.find(phase_terms,list(alpha),covariance,max_degree=ibp_degree_cap)
        certificates.append(cert)
        if cert['outcome']!='COMPLETE':
            return dict(outcome='UNKNOWN',cover=cover,low_marker_certificates=certificates)
        selected=cert['actual_period_test']
        if selected['outcome']=='NO_MATCH':
            return dict(outcome='UNSAFE',fixed_bad_marker_exponent=list(alpha),
                        cover=cover,low_marker_certificates=certificates,
                        reason='A certified original scalar period is not polynomial')
        polynomial=selected['selected_polynomial_images'][0]['coefficients']
        degree=max((j for j,c in enumerate(polynomial) if sp.Rational(c)),default=-1)
        E=max(E,degree-C*sum(alpha))
    return dict(outcome='SAFE',cover=cover,low_marker_certificates=certificates,
                affine_cutoff=dict(slope=str(C),intercept=str(E),
                                   statement='E[Q P^m]=0 whenever m > slope*deg(Q)+intercept'),
                all_marker_coverage='Proved by the original degree-lowering Euler/IBP recurrence',
                source_relation_completeness_assumed=False)
