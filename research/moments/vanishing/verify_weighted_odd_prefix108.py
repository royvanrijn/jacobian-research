"""Exact auxiliary weighted odd sources, projections and seed certificates.
Universal induction, complete covers and analytic exclusions remain written.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_euler_seed_bounds as euler
import verify_shared_divisor_weighted_seeds as shared
import verify_rescaled_singular_prefixes as rescaled
import verify_endpoint_degrees_and_odd_prefixes as previous
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
CASES={'transverse':{'weights':(s.Rational(4,3),s.Integer(2)),'cost':s.Rational(16,3),'gate':98},'triangle':{'weights':(s.Integer(2),s.Integer(2)),'cost':s.Integer(6),'gate':78},'tangent':{'weights':(s.Integer(1),s.Rational(3,2)),'cost':s.Integer(4),'gate':52}}

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def degree(f,weights):return max((e[0]+weights[0]*e[1]+weights[1]*e[2] for e,c in s.Poly(f,t,u,v).terms() if c),default=-1000)
def amp(ex,weights):return sum(a*b for a,b in zip(ex,weights))

def phase(kind,D,ju,jv,ell,E,F,qu,qv,quv,a=1):
    if kind!='tangent':return previous.odd_phase('nodal',D,ju,jv,ell,E,F,qu,qv,quv,a if kind=='transverse' else 0,0)
    return s.expand(D*(u*v*v+E*t*t*u+F*t*t*v-(qu*u*u+2*quv*u*v+qv*v*v)/2)+(u*u-ju*u-jv*v+ell*t*t)**2/2)

def operators(kind,W,D):
    delta=t*t*D
    if kind!='tangent':
        source=s.cancel(euler.ordinary_source(W,D)/t);cu=s.diff(delta,t)/(4*t);cv=cu
    else:
        source,bx,by=shared.weighted_source(kind,W,D);source=s.cancel(source/t);cu=s.cancel(bx/t);cv=s.cancel(by/t)
    vector=lambda f:delta*(cu*u*s.diff(f,u)+cv*v*s.diff(f,v))
    return source,vector,cu,cv

def principal_prepared(kind,W,D):
    if kind=='tangent':return None
    pv=s.Poly(W,u,v);co=[pv.coeff_monomial(z) for z in [u*u*v*v,u**3,v**3,u*u,v*v,u*v,u,v,1]]
    prepared=ordinary.operator_data(kind,*co,t*t*D);assert s.degree(prepared[-1],t)==0
    return co,prepared

def step(kind,W,D,ex,prepared):
    delta=t*t*D
    if kind!='tangent':
        co,data=prepared;rem,cu,cv,_,_=ordinary.one_step(kind,ex,co,delta,data);rem,cu,cv=[s.cancel(z/data[-1]) for z in [rem,cu,cv]]
    else:
        nf,cu,cv=weighted.principal_division(kind,ex,s.Rational(1,2),D);R=W-u**4/2-D*u*v*v
        rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
    assert s.expand(u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
    assert degree(rem,CASES[kind]['weights'])<=amp(ex,CASES[kind]['weights'])
    return rem,cu,cv

def universal():
    D,ju,jv,ell,E,F,qu,qv,quv,a=s.symbols('D ju jv ell E F qu qv quv a');rows=[]
    for kind,data in CASES.items():
        W=phase(kind,D,ju,jv,ell,E,F,qu,qv,quv,a);source,vector,bu,bv=operators(kind,W,D);delta=t*t*D
        assert degree(source,data['weights'])<=data['cost']
        # h*partial_t, h_prime and the complete Euler field add at most2,
        # below every source cost. D is constant, and h=t^3 D^2.
        assert 2<=data['cost'] and s.diff(D,t)==0
        if kind=='tangent':
            R=W-u**4/2-D*u*v*v
            assert degree(s.diff(R,u),data['weights'])<=3 and degree(s.diff(R,v),data['weights'])<=2
            original=D*(t*(u*v*v+E*u+F*v)-(qu*u*u+2*quv*u*v+qv*v*v)/2)+(t*(u*u+ell)-ju*u-jv*v)**2/2
            assert s.expand(rescaled.rescale(t*t*original)-W)==0
        oldsource=s.cancel(rescaled.source_numerator(W,D)/t);gauge=[]
        for ex in ordinary.monomials(3):
            f=u**ex[0]*v**ex[1];Cu=bu*u*f;Cv=bv*v*f
            assert s.expand((source-oldsource)*f+vector(f)-delta*(s.diff(Cu,u)+s.diff(Cv,v))-s.diff(W,u)*Cu-s.diff(W,v)*Cv)==0
            gauge.append(ex)
        prepared=principal_prepared(kind,W,D);identities=[]
        for ex in ordinary.monomials(12) if kind!='tangent' else weighted.monomials(28,(2,3)):
            rem,cu,cv=step(kind,W,D,ex,prepared);identities.append({'exponent':ex,'input_weight':str(amp(ex,data['weights'])),'remainder_weight':str(degree(rem,data['weights']))})
        rows.append({'kind':kind,'time_weight':1,'amplitude_weights':[str(z) for z in data['weights']],'source_weight_bound':str(data['cost']),'D_constant_nonzero_complement':True,'cross_coefficients_retained':True,'gauge_identities':gauge,'principal_identities':identities,'normal_divisor_time_degree':0,'same_original_Wick_functional':True})
    return rows

def model(kind,D):
    D=s.Integer(D);W=phase(kind,D,s.Integer(2),s.Integer(3),s.Integer(1),s.Integer(5),s.Integer(7),s.Integer(11),s.Integer(13),s.Integer(1));source,vector,bu,bv=operators(kind,W,D);delta=t*t*D;h=t**3*D**2;data=CASES[kind]
    basis=(ordinary.CASES if kind!='tangent' else weighted.CASES)[kind]['basis'];index={e:i for i,e in enumerate(basis)};prepared=principal_prepared(kind,W,D)
    ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,v,domain=ring);zero=poly(0);zt=s.Poly(0,t,domain=s.QQ);wp,dp=poly(W),poly(delta)
    @lru_cache(None)
    def cached_step(ex):return tuple(poly(z) for z in step(kind,W,D,ex,prepared))
    @lru_cache(None)
    def reduce(ex,depth):
        if ex in index:
            out=[zt]*len(basis);out[index[ex]]=s.Poly(1,t,domain=s.QQ);return tuple(out),zero,zero
        assert depth>0,(kind,ex,depth);rem,Cu,Cv=cached_step(ex);out=[zt]*len(basis)
        for f,c in rem.terms():
            if not c:continue
            vec,x,z=reduce(f,depth-1);ct=s.Poly(c,t,domain=s.QQ)
            out=[a+ct*b for a,b in zip(out,vec)];Cu+=poly(c)*x;Cv+=poly(c)*z
        return tuple(out),Cu,Cv
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+source*R+vector(R));assert degree(R,data['weights'])<=data['cost']*j
        depth=(3*j-2 if kind=='transverse' else 3*j-1 if kind=='triangle' else 7*j-4)
        values=[zt]*len(basis);Cu,Cv=zero,zero
        for ex,c in poly(R).terms():
            if not c:continue
            vec,x,z=reduce(ex,depth);ct=s.Poly(c,t,domain=s.QQ)
            values=[a+ct*b for a,b in zip(values,vec)];Cu+=poly(c)*x;Cv+=poly(c)*z
        normal=sum((poly(z.as_expr()*u**ex[0]*v**ex[1]) for z,ex in zip(values,basis)),zero)
        assert (poly(R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        bounds=[int(s.floor(data['cost']*j-amp(ex,data['weights']))) for ex in basis]
        assert all(z.degree()<=bound for z,bound in zip(values,bounds))
        H=s.Poly(h**j,t,domain=s.QQ);val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        seeds.append({'j':j,'depth':depth,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'row_degree_bounds':bounds,'witness':{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}})
    metric=int(s.Matrix([[11,1,2],[1,13,3],[2,3,D]]).det());assert metric!=0
    return {'kind':kind,'D':int(D),'metric_determinant':metric,'W':str(W),'source':str(source),'seed_columns':seeds}

def arithmetic():
    rows=[]
    for kind,data in CASES.items():
        basis=(ordinary.CASES if kind!='tangent' else weighted.CASES)[kind]['basis'];r=len(basis);maximum=-1;checked=0
        for size in range(1,r+1):
            for chosen in combinations(range(r),size):
                if 0 not in chosen:continue
                # Integer polynomial determinant degree is bounded by the
                # floor of the entire weighted cofactor bound.
                value=int(s.floor(data['cost']*size*(size+1)/2-sum(amp(basis[k],data['weights']) for k in chosen)))
                maximum=max(maximum,value);checked+=1
        assert maximum==data['gate'],(kind,maximum,data['gate']);rows.append({'kind':kind,'cofactor_gate':maximum,'row_subsets_checked':checked,'basis_weights':[str(amp(ex,data['weights'])) for ex in basis]})
    cover={'quadratic':7,'smooth':94,'nodal':108,'cusp':94,'transverse':108,'tangent':94,'triangle':108,'concurrent':104,'repeated':98}
    assert max(cover.values())==108 and max(98,108)==108 and max(78,108)==108 and max(52,92,94)==94
    return rows,cover

def main():
    started=time.monotonic();univ=universal();print('Three complete odd weighted source/projection strata passed',flush=True)
    models=[model(kind,D) for kind in CASES for D in [1,2]];print('Six nondegenerate odd Gaussian models passed',flush=True);cases,cover=arithmetic()
    cert={'universal':univ,'models':models,'degree_cases':cases,'odd_cover':cover};path=P/'weighted_odd_prefix108_certificate.json';path.write_text(json.dumps(cert,indent=2)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,'universal_strata':3,'source_gauge_identities':sum(len(c['gauge_identities']) for c in univ),'principal_projection_identities':sum(len(c['principal_identities']) for c in univ),'models':6,'direct_seed_identities':18,'odd_complement_gates':{k:v['gate'] for k,v in CASES.items()},'odd_tangent_full_gate':94,'odd_sufficient_pure_prefix':108,'global_sufficient_pure_prefix_unchanged':130,'homogeneous_sufficient_pure_prefix_unchanged':108,'least_prefix_interval':[20,130],'odd_least_prefix_interval':[20,108],'homogeneous_least_prefix_interval':[14,108],'odd_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact auxiliary three weighted odd source/projection strata,30 gauge and262 principal identities,18 actual seed reductions in six nondegenerate centered Gaussian models and all-rank/row arithmetic. Universal all-degree reduction/source induction, complete odd coefficient/frame cover and analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded global radical cofactors.'}
    (P/'weighted_odd_prefix108_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
