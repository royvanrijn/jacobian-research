#!/usr/bin/env python3
"""Exact divergence reductions/model certificates for ordinary-node prefixes."""
from pathlib import Path
from functools import lru_cache
import hashlib,json
import sympy as S

u,v,t=S.symbols('u v t')
c4,a3,d3,b2,f2,e11,ell,n0,k0,delta_symbol=S.symbols('c4 a3 d3 b2 f2 e11 ell n0 k0 delta')
PARAMETERS=[c4,a3,d3,b2,f2,e11,ell,n0,k0,delta_symbol]
PARAMETER_WEIGHTS=[4]*9+[3]
CASES={
 'nodal':{'basis':[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(0,3)],'rank':7,'max_source':7,'minimum_reducible':3,'steps':5,'B':66,'gate':1855},
 'transverse':{'basis':[(0,0),(1,0),(0,1),(2,0),(1,1),(0,2)],'rank':6,'max_source':6,'minimum_reducible':3,'steps':4,'B':54,'gate':1140},
 'triangle':{'basis':[(0,0),(1,0),(0,1),(1,1),(0,2)],'rank':5,'max_source':6,'minimum_reducible':2,'steps':5,'B':66,'gate':995}}


def monomials(limit):return [(i,j) for i in range(limit+1) for j in range(limit-i+1)]


def operator_data(kind,C,A,D,B,F,E,L,N,K,delta):
    V=C*u*u*v*v+A*u**3+D*v**3+B*u*u+F*v*v+E*u*v+L*u+N*v+K
    pairs=[(S.Integer(1),S.Integer(0)),(S.Integer(0),S.Integer(1)),(u,-v),
           (C*u*v,-C*v*v-S.Rational(3,2)*A*u-(B if kind!='nodal' else 0))]
    full=[S.expand(x*S.diff(V,u)+y*S.diff(V,v)+delta*(S.diff(x,u)+S.diff(y,v))) for x,y in pairs]
    p3=2*B*u*u-2*F*v*v if kind=='triangle' else 3*A*u**3-3*D*v**3
    principals=[2*C*u*v*v,2*C*u*u*v,p3,-3*C*D*v**4 if kind=='nodal' else -2*C*F*v**3]
    divisor=C*(B if kind=='triangle' else A)*(D if kind=='nodal' else F)
    return S.expand(V),pairs,full,principals,divisor


def principal_division(kind,exponent,C,A,D,B,F):
    i,j=exponent;zero=S.Integer(0);hs=[zero]*4
    if exponent in CASES[kind]['basis']:return u**i*v**j,hs
    if i>=1 and j>=2:
        hs[0]=u**(i-1)*v**(j-2)/(2*C);return zero,hs
    if kind!='triangle' and i>=2 and j>=1:
        hs[1]=u**(i-2)*v**(j-1)/(2*C);return zero,hs
    if kind=='nodal':
        if i==3 and j==0:hs[2]=1/(3*A);return D*v**3/A,hs
        if i>=4 and j==0:
            hs[2]=u**(i-3)/(3*A);hs[0]=D*u**(i-4)*v/(2*A*C);return zero,hs
        assert i==0 and j>=4
        hs[3]=-v**(j-4)/(3*C*D);return zero,hs
    if kind=='transverse':
        if i>=3 and j==0:hs[2]=u**(i-3)/(3*A);return zero,hs
        assert i==0 and j>=3
        hs[3]=-v**(j-3)/(2*C*F);return zero,hs
    if i>=2:
        hs[2]=u**(i-2)*v**j/(2*B)
        if i==2 and j==0:return F*v*v/B,hs
        if i==2 and j==1:hs[3]=-1/(2*B*C);return zero,hs
        assert i>=3 and j<=1
        hs[0]=F*u**(i-3)*v**j/(2*B*C);return zero,hs
    assert i==0 and j>=3
    hs[3]=-v**(j-3)/(2*C*F);return zero,hs


def one_step(kind,exponent,coefficients,delta,prepared=None):
    C,A,D,B,F,E,L,N,K=coefficients
    V,pairs,full,principals,divisor=prepared if prepared is not None else operator_data(kind,*coefficients,delta)
    nf,hs=principal_division(kind,exponent,C,A,D,B,F)
    cleared=[S.cancel(divisor*h) for h in hs]
    Cu=S.expand(sum(h*x for h,(x,y) in zip(cleared,pairs)))
    Cv=S.expand(sum(h*y for h,(x,y) in zip(cleared,pairs)))
    # Independent principal projection, with derivatives of the amplitude.
    remainder=S.expand(S.cancel(divisor*nf)-sum(h*(r-p) for h,r,p in zip(cleared,full,principals))
                       -delta*sum(x*S.diff(h,u)+y*S.diff(h,v) for h,(x,y) in zip(cleared,pairs)))
    return remainder,Cu,Cv,V,divisor


def parameter_degree(expr):
    poly=S.Poly(expr,*PARAMETERS)
    return max((sum(a*b for a,b in zip(exponent,PARAMETER_WEIGHTS)) for exponent,c in poly.terms() if c),default=-1)


def universal_checks():
    rows=[]
    for kind,data in CASES.items():
        coeffs=[c4,a3 if kind!='triangle' else S.Integer(0),d3 if kind=='nodal' else S.Integer(0),b2,f2,e11,ell,n0,k0]
        prepared=operator_data(kind,*coeffs,delta_symbol)
        V,pairs,full,principals,divisor=prepared
        assert all(S.Poly(r-p,u,v).total_degree()<S.Poly(p,u,v).total_degree() for r,p in zip(full,principals))
        count=0
        for exponent in monomials(data['max_source']):
            nf,hs=principal_division(kind,exponent,*coeffs[:5])
            assert S.expand(u**exponent[0]*v**exponent[1]-nf-sum(h*p for h,p in zip(hs,principals)))==0
            remainder,Cu,Cv,_,_=one_step(kind,exponent,coeffs,delta_symbol,prepared)
            defect=divisor*u**exponent[0]*v**exponent[1]-remainder-delta_symbol*(S.diff(Cu,u)+S.diff(Cv,v))-S.diff(V,u)*Cu-S.diff(V,v)*Cv
            assert S.expand(defect)==0
            assert all(e in data['basis'] or sum(e)<=sum(exponent)-1 for e,c in S.Poly(remainder,u,v).terms() if c)
            assert all(parameter_degree(c)<=12 for e,c in S.Poly(remainder,u,v).terms() if c)
            assert all(parameter_degree(c)<=8 for poly in [Cu,Cv] for e,c in S.Poly(poly,u,v).terms() if c)
            count+=1
        assert data['steps']==data['max_source']-data['minimum_reducible']+1
        assert data['B']==6+12*data['steps']
        assert data['gate']==data['B']*data['rank']*(data['rank']+1)//2+data['rank']
        rows.append({'kind':kind,'principal_and_divergence_identities':count,'cleared_remainder_t_degree':12,'cleared_cofactor_t_degree':8})
    return rows


def gaussian_phase(kind,singular):
    m=S.Integer(0) if singular else S.Integer(1)
    alpha=-S.Rational(1,2) if singular else S.Integer(1)
    D=m-2*alpha*t
    aa=S.Integer(1) if kind!='triangle' else S.Integer(0)
    dd=S.Integer(2) if kind=='nodal' else S.Integer(0)
    bb=S.Integer(0) if singular else S.Integer(3)
    cc=S.Integer(0) if singular and kind!='nodal' else S.Integer(5)
    G=aa*u**3+dd*v**3+bb*u*u+cc*v*v+7*u+11*v+2
    q=(11*u*u+13*v*v)/2+u*v+5*u+7*v+2
    L=t*(u*v+1)-2*u-3*v-1
    U=D*(t*G-q)+L*L/2
    V=S.expand(t*t*U.subs({u:u+3/t,v:v+2/t},simultaneous=True))
    pv=S.Poly(V,u,v)
    assert pv.coeff_monomial(u*u*v)==pv.coeff_monomial(u*v*v)==0
    terms=[u*u*v*v,u**3,v**3,u*u,v*v,u*v,u,v,S.Integer(1)]
    coefficients=[S.factor(pv.coeff_monomial(term)) for term in terms]
    divisor=coefficients[0]*(coefficients[3] if kind=='triangle' else coefficients[1])*(coefficients[2] if kind=='nodal' else coefficients[4])
    assert divisor!=0
    assert all(S.Poly(c,t).degree()<=4 for c in coefficients)
    metric=S.Matrix([[11,1,2],[1,13,3],[2,3,m]])
    assert metric.det()!=0
    return V,coefficients,D,t*t*D,int(metric.det()),str(divisor)


def model_connection(kind,singular):
    data=CASES[kind];ring=S.QQ.poly_ring(t)
    poly=lambda expr:S.Poly(expr,u,v,domain=ring)
    zero=poly(0);zero_t=S.Poly(0,t,domain=S.QQ)
    V,coeffs,D,delta,metric_det,_=gaussian_phase(kind,singular)
    prepared=operator_data(kind,*coeffs,delta);divisor=prepared[-1]
    At=S.Poly(S.expand(divisor),t,domain=S.QQ);Ap=poly(At.as_expr());dp=poly(delta)
    Vp=poly(V);Vu=Vp.diff(u);Vv=Vp.diff(v)
    basis=data['basis'];index={e:i for i,e in enumerate(basis)}
    @lru_cache(None)
    def step(exponent):
        remainder,Cu,Cv,_,_=one_step(kind,exponent,coeffs,delta,prepared)
        result=(poly(remainder),poly(Cu),poly(Cv))
        assert all(e in index or sum(e)<=sum(exponent)-1 for e,c in result[0].terms() if c)
        assert all(S.Poly(c,t).degree()<=12 for e,c in result[0].terms() if c)
        return result
    @lru_cache(None)
    def reduce_monomial(exponent,remaining):
        if exponent in index:
            values=[zero_t]*len(basis);values[index[exponent]]=At**remaining
            return tuple(values),zero,zero
        assert remaining>0,(kind,exponent)
        remainder,cu0,cv0=step(exponent)
        values=[zero_t]*len(basis);Cu=zero;Cv=zero
        for e,c in remainder.terms():
            if not c:continue
            vector,cu,cv=reduce_monomial(e,remaining-1)
            scalar_t=S.Poly(c,t,domain=S.QQ);scalar=poly(c)
            values=[x+scalar_t*y for x,y in zip(values,vector)]
            Cu+=scalar*cu;Cv+=scalar*cv
        Cu+=Ap**(remaining-1)*cu0;Cv+=Ap**(remaining-1)*cv0
        return tuple(values),Cu,Cv
    source=S.expand(S.diff(V,t)*delta-V*S.diff(delta,t)-t**4*D*S.diff(D,t)/2)
    assert S.Poly(source,t).degree()<=6
    sourcep=poly(source);H=S.expand(delta*delta*divisor**data['steps'])
    columns=[];certificates=[]
    for exponent in basis:
        f=sourcep*poly(u**exponent[0]*v**exponent[1])
        values=[zero_t]*len(basis);Cu=zero;Cv=zero
        for e,c in f.terms():
            if not c:continue
            vector,cu,cv=reduce_monomial(e,data['steps'])
            scalar_t=S.Poly(c,t,domain=S.QQ);scalar=poly(c)
            values=[x+scalar_t*y for x,y in zip(values,vector)]
            Cu+=scalar*cu;Cv+=scalar*cv
        normal=sum((poly(value.as_expr()*u**e[0]*v**e[1]) for e,value in zip(basis,values)),zero)
        defect=Ap**data['steps']*f-normal-dp*(Cu.diff(u)+Cv.diff(v))-Vu*Cu-Vv*Cv
        assert defect.is_zero
        assert all(value.degree()<=data['B'] for value in values)
        columns.append([str(value.as_expr()) for value in values])
        certificates.append({'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())})
    assert S.Poly(H,t).degree()<=data['B']
    return {'kind':kind,'singular_at_zero':singular,'metric_determinant':metric_det,
            'D':str(D),'delta':str(delta),'V':str(V),'divisor':str(divisor),'basis':[str(u**e[0]*v**e[1]) for e in basis],
            'source_numerator':str(source),'H':str(H),'N_columns':columns,'integration_by_parts_certificates':certificates,
            'connection_identity':'A^r*source*basis_j-N_j=delta*(partial_u Cu+partial_v Cv)+V_u*Cu+V_v*Cv',
            'H_degree':int(S.Poly(H,t).degree()),'N_max_degree':int(max(S.Poly(x,t).degree() for col in columns for x in col)),
            'columns_certified':len(basis),'B':data['B'],'rank':data['rank'],'prefix':data['gate']}


def main():
    p=Path(__file__).resolve().parent
    identities=universal_checks();print('Universal principal/divergence checks passed',flush=True)
    models=[]
    for kind in CASES:
        for singular in [False,True]:
            models.append(model_connection(kind,singular));print('Certified',kind,'singular' if singular else 'regular',flush=True)
    path=p/'ordinary_node_connection_models.json';path.write_text(json.dumps({'completed':True,'models':models},indent=2)+'\n')
    checkpoint=json.loads((p/'ordinary_node_phase_shift_checkpoint.json').read_text());assert checkpoint['completed']
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'model_certificate_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'normalization_checkpoint_sha256':hashlib.sha256((p/'ordinary_node_phase_shift_checkpoint.json').read_bytes()).hexdigest(),
            'universal_principal_divergence_checks':identities,
            'class_bounds':{kind:{key:row[key] for key in ['rank','max_source','minimum_reducible','steps','B','gate']} for kind,row in CASES.items()},
            'connection_models':[{key:row[key] for key in ['kind','singular_at_zero','metric_determinant','H_degree','N_max_degree','columns_certified','B','rank','prefix']} for row in models],
            'ordinary_isotropic_frame_exclusion_gate':888,
            'proposed_global_sufficient_pure_prefix':max([row['gate'] for row in CASES.values()]+[935,720,104,99,188,888,7]),
            'global_least_prefix_lower_bound':10,'global_least_prefix_computed':False,
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Universal principal/divergence identities and t-degree/ordinary-degree bounds; six exact Gaussian-normal-form rational model connections with all36 column integration-by-parts certificates, including D(0)=0 and vanishing quadratic coefficients. Full rank/normal-form coverage, all-order cofactor/seed argument, analytic nonzero zeroth coefficients and global1855 conclusion require the written proof; no high-order moment expansion, optimality, independent review or formal verification.'}
    (p/'ordinary_node_finite_prefixes_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','proposed_global_sufficient_pure_prefix','global_least_prefix_computed','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
