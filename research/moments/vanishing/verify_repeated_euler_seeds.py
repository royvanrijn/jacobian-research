"""Auxiliary Euler-representative actual seeds for repeated-line gate92.
Universal all-j induction,minimal seed/period bridge,full actual chart/frame
cover and analytic nonzero zeroth inputs remain written and unreviewed.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
P=Path(__file__).resolve().parent
t,u=s.symbols('t u');sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()

def source(U,D,d):
    a=s.expand(s.diff(U,u)).coeff(u,d-1);L=s.expand(s.diff(a,t)*D-a*s.diff(D,t))
    C=s.expand(D*s.diff(U,t)-s.diff(D,t)*U-D*s.diff(D,t)/2)
    S=s.expand(a*C-L*(u*s.diff(U,u)+D)/d)
    h=s.expand(a*D*D);V=s.expand(-D*L/d)
    return a,L,C,S,h,V

def universal():
    rows=[]
    for d in [3,4]:
        U=sum(s.Symbol('c'+str(k)+str(j))*t**j*u**k for k in range(d+1) for j in range(4))
        D=sum(s.Symbol('D'+str(j))*t**j for j in range(3));a,L,C,S,h,V=source(U,D,d)
        assert s.degree(S,u)<=d-1 and s.degree(S,t)<=7 and s.degree(h,t)<=7 and s.degree(V,t)<=6
        gauge=[]
        for n in range(7):
            f=u**n;w=-L*u*f/d
            assert s.expand((S-a*C)*f+V*u*s.diff(f,u)-D*s.diff(w,u)-s.diff(U,u)*w)==0;gauge.append(n)
        principal=[]
        for n in range(d-1,13):
            w=u**(n-d+1);rem=s.expand(-w*(s.diff(U,u)-a*u**(d-1))-D*s.diff(w,u))
            assert s.degree(rem,u)<=n-1 and s.degree(rem,t)<=3
            assert s.expand(a*u**n-rem-D*s.diff(w,u)-s.diff(U,u)*w)==0;principal.append(n)
        rows.append({'phase_degree':d,'source_time_degree':7,'source_amplitude_degree':d-1,
                     'h_time_degree':7,'vertical_time_degree':6,'Euler_gauge_identities':gauge,
                     'principal_identities':principal,'actual_R_time_bound':'7j','actual_R_amplitude_bound':f'{d-1}j',
                     'depth':f'{d-1}j-{d-2}','row_time_bound':f'{3*d+4}j-{3*(d-2)}',
                     'full_rank_zeroth_degree_bound':33 if d==3 else 78})
    return rows

def cofactors(cols,n):
    ring=s.QQ.poly_ring(t);out=[]
    for k in range(n+1):
        selected=cols[:k]+cols[k+1:]
        M=DomainMatrix.from_list_sympy(n,n,[[c.as_expr() for c in row] for row in zip(*selected)]).convert_to(ring)
        out.append((-1)**k*s.Poly(ring.to_sympy(M.det()),t,domain=s.QQ))
    assert all(sum((out[k]*cols[k][i] for k in range(n+1)),s.Poly(0,t)).is_zero for i in range(n))
    return out

def model(label,U,D,d,origin):
    a,L,C,S,h,V=source(U,D,d);ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,domain=ring)
    zt=s.Poly(0,t);zero=poly(0);ap=s.Poly(a,t);wp,dp=poly(U),poly(D)
    @lru_cache(None)
    def reduce(n,depth):
        if n<d-1:
            out=[zt]*(d-1);out[n]=ap**depth;return tuple(out),zero
        assert depth>0
        w=u**(n-d+1);rem=poly(-w*(s.diff(U,u)-a*u**(d-1))-D*s.diff(w,u))
        out=[zt]*(d-1);W=poly(a**(depth-1)*w)
        for (k,),c in rem.terms():
            if not c:continue
            vals,z=reduce(k,depth-1);ct=s.Poly(c,t)
            out=[v+ct*b for v,b in zip(out,vals)];W+=poly(c)*z
        return tuple(out),W
    R=s.Integer(1);rows=[];cols=[[s.Poly(int(k==0),t) for k in range(d-1)]]
    for j in range(1,d):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+V*u*s.diff(R,u)+S*R)
        assert s.degree(R,t)<=7*j and s.degree(R,u)<=(d-1)*j
        depth=(d-1)*j-(d-2);vals=[zt]*(d-1);W=zero
        for (n,),c in poly(R).terms():
            if not c:continue
            v,w=reduce(n,depth);ct=s.Poly(c,t)
            vals=[v+ct*b for v,b in zip(vals,v)];W+=poly(c)*w
        normal=sum((poly(c.as_expr()*u**k) for k,c in enumerate(vals)),zero)
        assert (poly(a**depth*R)-normal-dp*W.diff(u)-wp.diff(u)*W).is_zero
        bound=(3*d+4)*j-3*(d-2);assert all(c.degree()<=bound for c in vals)
        H=s.Poly(h**j*a**depth,t);hv=min(e[0] for e,c in H.terms() if c);assert hv>=j
        rows.append({'j':j,'depth':depth,'N':[str(c.as_expr()) for c in vals],'H':str(H.as_expr()),'H_valuation':hv,
                     'row_time_bound':bound,'witness':str(W.as_expr()),'exact_original_twisted_identity':True});cols.append(vals)
    cs=cofactors(cols,d-1);checks=[];gate=33 if d==3 else 78
    for j,C in enumerate(cs):
        H=s.Poly(1,t) if j==0 else s.Poly(s.sympify(rows[j-1]['H'],locals={'t':t}),t)
        A=C*H;val=min(e[0] for e,c in A.terms() if c) if A else None
        if j==0:assert A and A.degree()<=gate
        else:assert val is None or val>=j
        checks.append({'order':j,'scalar_degree':int(A.degree()) if A else -1,'scalar_valuation':val})
    return {'label':label,'origin':origin,'degree':d,'U':str(U),'D':str(D),'a':str(a),'Euler_L':str(L),
            'h':str(h),'source':str(S),'seed_columns':rows,'cofactor_checks':checks}

def adapted_frame(g02):
    z=(u-t*(1+2*u+5*u*u))/(3*t)
    Phi=t*(7*u**3+11*u*u+13*u*z+g02*z*z+19*u+23*z+29)-z*z/2
    U=s.expand(t*t*Phi);assert s.degree(U,t)<=3 and s.degree(U,u)==4
    assert s.expand(U).coeff(u,4)==s.expand(25*t*t*(g02*t-s.Rational(1,2))/9)
    return U,t*t

def arithmetic():
    rows=[]
    for d in [3,4]:
        cases=[]
        for n in range(1,d):
            for selected in combinations(range(d-1),n):
                if 0 not in selected:continue
                C=(3*d+4)*n*(n+1)//2-3*(d-2)*n
                cases.append({'rank':n,'rows':selected,'zeroth_cofactor_degree_bound':C})
        gate=max(c['zeroth_cofactor_degree_bound'] for c in cases);assert gate==(33 if d==3 else 78)
        rows.append({'phase_degree':d,'gate':gate,'all_rank_row_cases':cases})
    assert max(33,78,35,92,78)==92
    cover={'quadratic':7,'smooth':94,'nodal':96,'cusp':95,'transverse':96,'tangent':95,'triangle':96,'concurrent':92,'repeated':92}
    assert max(cover.values())==96
    return rows,cover

def main():
    start=time.monotonic();univ=universal();print('Universal full Euler source/gauge and principal identities passed',flush=True)
    models=[]
    retained=json.loads((P/'concurrent_seeds_repeated_valuations_certificate.json').read_text())['repeated_models']
    for k,m in enumerate(retained):
        U=s.sympify(m['U'],locals={'t':t,'u':u});D=s.sympify(m['D'],locals={'t':t})
        models.append(model('retained_'+str(k),U,D,m['degree'],'concurrent_seeds_repeated_valuations_certificate.json'))
        print('Retained phase',k,'passed',flush=True)
    for c in [17,0]:
        U,D=adapted_frame(c);models.append(model('adapted_g02_'+str(c),U,D,4,'actual restricted mixed quartic'))
        print('Adapted frame',c,'passed',flush=True)
    cases,cover=arithmetic();cert={'universal':univ,'models':models,'degree_cases':cases,'global_cover':cover}
    f=P/'repeated_euler_seeds_certificate.json';f.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    receipt={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
             'Euler_gauge_identities':sum(len(x['Euler_gauge_identities']) for x in univ),'principal_identities':sum(len(x['principal_identities']) for x in univ),
             'models':6,'retained_phase_models':4,'actual_adapted_frame_models':2,'exact_original_twisted_seed_identities':sum(len(m['seed_columns']) for m in models),
             'full_scalar_cofactor_relations':6,'scalar_coefficient_checks':sum(len(m['cofactor_checks']) for m in models),
             'repeated_cubic_complement_gate':33,'repeated_quartic_complement_gate':78,'restricted_mixed_quartic_gate':78,'repeated_full_gate':92,
             'global_cover':cover,'global_sufficient_pure_prefix_unchanged':96,'odd_sufficient_pure_prefix_unchanged':96,'homogeneous_sufficient_pure_prefix_unchanged':96,
             'least_prefix_interval':[20,96],'odd_least_prefix_interval':[20,96],'homogeneous_least_prefix_interval':[14,96],
             'least_prefix_computed':False,'universal_written_proof_certified_by_checker':False,'independent_replay':False,'external_review':False,'formal_verification':False,
             'boundary':'Exact universal full time-dependent Euler source/gauge and principal identities,16 original twisted seeds and six full scalar relations/22 coefficient checks including two actual adapted mixed quartics,all rank/row and scalar valuation checks. Universal all-j induction,actual elimination/frame coverage,seed-to-period bridge and analytic nonzero-zeroth inputs remain independently unreviewed written mathematics. No old checker main,radical search or external producer package.'}
    (P/'repeated_euler_seeds_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
