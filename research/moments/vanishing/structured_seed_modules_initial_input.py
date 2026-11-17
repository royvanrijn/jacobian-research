"""Exact auxiliary conic-module, mixed-parent Euler and tangent shear checks.

All analytic exclusions and full coefficient/rank coverage remain written inputs.
"""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import hashlib,json,time
import sympy as s
import verify_balanced_logarithmic_seed_bounds as balanced
import verify_euler_seed_bounds as euler
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_rescaled_singular_prefixes as rescaled
import verify_shared_divisor_weighted_seeds as shared
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
y=balanced.y

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def laurent(f):return balanced.laurent(f)
def total(f,*gens):return -1000 if f==0 else int(s.Poly(f,*gens).total_degree())

def coordinate_module_support():
    # Each coordinate is const*y+(c0+c1*t)+const*K/y.
    parts={(0,1,0),(0,0,0),(1,0,0),(0,-1,1)}
    support={(0,0,0)};rows=[]
    for m in range(4):
        if m:support={(a+d,b+e,c+f) for a,b,c in support for d,e,f in parts}
        for ell in range(3):
            for n,k,q in support:
                power=n+3-m
                assert q>=max(0,-k)
                if k<0:assert power+ell*(q+k)<=3+k
                else:assert power+ell*q<=3-k
            rows.append({'degree_G_monomial':m,'degree_K':ell,'support_terms_checked':len(support)})
    gaussian={(a+d,b+e,c+f) for a,b,c in parts for d,e,f in parts}
    for ell in range(3):
        for n,k,q in gaussian:
            assert q>=max(0,-k)
            assert n+ell*(q+min(k,0))<=(2+k if k<0 else 2-k)
    return rows

def conic_case(ell):
    wt,wy,base,highcost,lowcost=[(1,1,3,0,6),(2,1,8,3,3),(1,1,5,0,0)][ell]
    bj=lambda j:base*j+highcost*(3*j-2)+lowcost*max(0,2*j-3)
    rows=[]
    for size in range(1,7):
        for selected in combinations(range(-3,3),size):
            if 0 not in selected:continue
            gate=(sum(bj(j) for j in range(1,size+1))-wy*sum(selected))//wt
            rows.append({'rank':size,'rows':selected,'gate':gate})
    gate=max(x['gate'] for x in rows);assert gate==[216,199,108][ell]
    return {'degree_K':ell,'weight_t':wt,'weight_y':wy,'seed_weight_base':base,'high_cost':highcost,'low_cost':lowcost,'seed_degree_bounds':[bj(j) for j in range(1,7)],'gate':gate,'cofactor_cases':rows}

def conic_model(label,W,K):
    field=s.QQ.frac_field(t);raw=laurent(W);ell=int(s.Poly(K,t).degree());case=conic_case(ell)
    assert max(raw)==3 and min(raw)==-3
    assert not raw[3].has(t) and not s.cancel(raw[-3]/K**3).has(t)
    wt,wy=case['weight_t'],case['weight_y'];basis=list(range(-3,3));valid=set(basis)
    def polynomial(expr):
        expr=s.cancel(expr);assert not s.denom(expr).has(t)
        return s.Poly(expr,t,domain=s.QQ)
    def module(poly):
        for k,c in poly.items():
            co=polynomial(field.to_sympy(c))
            if k<0:co.exquo(s.Poly(K**(-k),t,domain=s.QQ))
    for k,c in raw.items():
        if k<0:assert polynomial(c/K**(-k)).degree()<=3+k
        else:assert polynomial(c).degree()<=3-k
    B=t*s.diff(K,t)-s.Rational(2,3)*K;h=t**3*K
    old=t*s.diff(W,t)-2*W-t*t;source=s.expand(K*old+B*y*s.diff(W,y))
    assert min(laurent(source))>=-2
    deriv={k:field.convert(k*c) for k,c in raw.items() if k}
    parse=lambda f:{k:field.convert(c) for k,c in laurent(f).items()}
    def add(out,k,c):
        out[k]=out.get(k,field.zero)+c
        if not out[k]:del out[k]
    def exact(poly):
        out={}
        for n,c in poly.items():
            add(out,n,field.convert(t*t*n)*c)
            for k,a in deriv.items():add(out,n+k,a*c)
        return out
    def reduce(poly):
        module(poly);out=poly.copy();witness={};low=high=0
        while any(k not in valid for k in out):
            n=max(out) if max(out)>2 else min(out)
            shift=n-3 if n>2 else n+3;lead=deriv[3] if n>2 else deriv[-3]
            factor=out[n]/lead;polynomial(field.to_sympy(factor))
            if n>2:high+=1
            else:low+=1
            add(witness,shift,factor)
            for k,c in exact({shift:factor}).items():add(out,k,-c)
            module(out)
        defect=poly.copy()
        for k,c in out.items():add(defect,k,-c)
        for k,c in exact(witness).items():add(defect,k,-c)
        assert not defect
        return [polynomial(field.to_sympy(out.get(k,field.zero))) for k in basis],witness,[high,low]
    # Unit connection and direct seeds have the module property; arbitrary
    # unscaled negative monomial basis columns need not have that property.
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        next=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+t*t*B*y*s.diff(R,y)+source*R)
        old_next=h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+K*old*R
        assert s.expand(next-old_next-B*(t*t*y*s.diff(R,y)+y*s.diff(W,y)*R))==0
        R=next;values,witness,counts=reduce(parse(R))
        assert counts[0]<=3*j-2 and counts[1]<=max(0,2*j-3)
        bnd=case['seed_degree_bounds'][j-1]
        assert all(wt*z.degree()+wy*k<=bnd for k,z in zip(basis,values) if not z.is_zero)
        H=s.Poly(s.expand(h**j),t,domain=s.QQ);assert wt*H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        seeds.append({'j':j,'degree_bound_weighted':bnd,'H':str(H.as_expr()),'H_valuation':val,'N':[str(z.as_expr()) for z in values],'reduction_counts':counts,'witness':{str(k):str(field.to_sympy(c)) for k,c in witness.items()}})
    return {'label':label,'degree_K':ell,'K':str(K),'W':str(W),'basis':basis,'polynomial_module_checked_after_every_cancellation':True,'seed_columns':seeds}

def parent_universal():
    z,T,a=s.symbols('z T a');g={(i,j):s.Symbol('g'+str(i)+str(j)) for i in range(4) for j in range(4-i)}
    F=z/t-a-z*z;W=s.expand(t*t*(t*sum(co*z**i*F**j for (i,j),co in g.items())-F*F/2));rows=[]
    for d in [5,6]:
        phase=W.subs(g[0,3],0) if d==5 else W
        lead=s.Poly(phase,z).coeff_monomial(z**d)
        assert not s.cancel(lead/t**3).has(t)
        old=t*s.diff(phase,t)-2*phase-t*t
        source=s.expand(old-t*t/d-z*s.diff(phase,z)/d)
        assert s.degree(source,z)<=d-1
        assert all(s.Poly(c,t).degree()<=3 for e,c in s.Poly(source,z).terms() if c)
        checks=[]
        for i in range(5):
            f=z**i;diff=(source-old)*f-t*t*z*s.diff(f,z)/d
            exact=-(t*t*s.diff(z*f,z)+s.diff(phase,z)*z*f)/d
            assert s.expand(diff-exact)==0;checks.append(i)
        rows.append({'d':d,'leading_coefficient_over_t_cubed':str(s.cancel(lead/t**3)),'connection_identity_checks':checks,'source_amplitude_degree_bound':d-1})
    return rows

def parent_model(label,W):
    ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,y,domain=ring)
    wp=poly(W);d=wp.degree();assert d in [5,6]
    derivative=wp.diff(y);A=s.Poly(derivative.LC(),t,domain=s.QQ);assert A.degree()==3
    assert not s.cancel(A.as_expr()/t**3).has(t)
    source=s.expand(t*s.diff(W,t)-2*W-t*t-t*t/d-y*s.diff(W,y)/d)
    assert s.degree(source,y)<=d-1
    ap=poly(A.as_expr());delta=poly(t*t);zero=poly(0);zt=s.Poly(0,t,domain=s.QQ)
    @lru_cache(None)
    def reduce(n,depth):
        if n<d-1:
            values=[zt]*(d-1);values[n]=A**depth;return tuple(values),zero
        assert depth>0;cu=poly(y**(n-d+1));rem=ap*poly(y**n)-delta*cu.diff(y)-derivative*cu
        assert rem.degree()<n;values=[zt]*(d-1);witness=ap**(depth-1)*cu
        for (k,),c in rem.terms():
            if not c:continue
            vec,w=reduce(k,depth-1);ct=s.Poly(c,t,domain=s.QQ)
            values=[a+ct*b for a,b in zip(values,vec)];witness+=poly(c)*w
        return tuple(values),witness
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        old=t*s.diff(W,t)-2*W-t*t;old_next=t**3*s.diff(R,t)-3*(j-1)*t*t*R+old*R
        next=s.expand(t**3*s.diff(R,t)-3*(j-1)*t*t*R-t*t*y*s.diff(R,y)/d+source*R)
        assert s.expand(next-old_next+(t*t*s.diff(y*R,y)+s.diff(W,y)*y*R)/d)==0
        R=next;assert s.degree(R,y)<=(d-1)*j
        depth=(d-1)*(j-1)+1;values=[zt]*(d-1);witness=zero
        for (n,),c in poly(R).terms():
            if not c:continue
            vec,w=reduce(n,depth);ct=s.Poly(c,t,domain=s.QQ)
            values=[a+ct*b for a,b in zip(values,vec)];witness+=poly(c)*w
        normal=sum((poly(z.as_expr()*y**k) for k,z in enumerate(values)),zero)
        assert (ap**depth*poly(R)-normal-delta*witness.diff(y)-derivative*witness).is_zero
        bnd=3*j+3*depth;assert all(z.degree()<=bnd for z in values)
        H=s.Poly(s.expand(t**(3*j)*A.as_expr()**depth),t,domain=s.QQ);assert H.degree()<=bnd
        seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'witness':str(witness.as_expr())})
    return {'label':label,'d':d,'W':str(W),'source':str(source),'seed_columns':seeds}

def tangent_phase(singular,zero_shear):
    D=t if singular else 1-2*t;m=0 if singular else 1;jv=0 if zero_shear else 3;ev=0 if zero_shear else 3
    ell=2*u+ev*v+1;J=2*u+jv*v+1
    g=3*u*u+5*u*v+7*v*v+11*u+13*v+17;q=(11*u*u+13*v*v)/2+u*v+5*u+7*v+2
    original=D*(t*(u*v*v+g)-q)+(t*(u*u+ell)-J)**2/2
    W=rescaled.rescale(t*t*original);metric=s.Matrix([[11,1,2],[1,13,jv],[2,jv,m]]).det();assert metric!=0
    return W,D,int(metric)

def tangent_source(W,D):
    delta=t*t*D;old,bx,by=shared.weighted_source('tangent',W,D)
    C=s.Poly(old,u,v).coeff_monomial(u*u*v);assert s.Poly(C,t).degree()<=3
    source=s.expand(D*old-C*u*s.diff(W,v)/2)
    assert total(source,t,u,v)<=7
    assert all(weighted.weight(e,(2,3))<=6 for e,c in s.Poly(source,u,v).terms() if c)
    return source,bx,by,C,old

def tangent_universal():
    d0,d1=s.symbols('d0 d1');D=d0+d1*t;delta=t*t*D
    W=u**4/2+D*u*v*v+sum(s.Symbol('a'+str(i)+'_'+str(j)+'_'+str(k))*t**k*u**i*v**j for i,j in weighted.monomials(7,(2,3)) for k in range(5-i-j))
    source,bx,by,C,old=tangent_source(W,D);checks=[]
    for i,j in ordinary.monomials(4):
        f=u**i*v**j;diff=(source-D*old)*f-delta*C*u*s.diff(f,v)/2;witness=-C*u*f/2
        assert s.expand(diff-delta*s.diff(witness,v)-s.diff(W,v)*witness)==0;checks.append([i,j])
    return {'source_joint_degree':7,'source_weight_bound':6,'shear_coefficient_C_degree_bound':3,'connection_identity_checks':checks,'new_seed_denominator':'delta^2*D'}

def tangent_model(singular,zero_shear):
    kind="tangent"
    W,D,metric=tangent_phase(singular,zero_shear);delta=t*t*D;source,bx,by,C,old=tangent_source(W,D);h=delta**2*D
    data=weighted.CASES[kind];basis=data['basis'];index={e:i for i,e in enumerate(basis)}
    k=1 if kind=='cusp' else 2;A=D**k
    ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,v,domain=ring);zero=poly(0);zt=s.Poly(0,t,domain=s.QQ)
    wp,dp,ap=map(poly,[W,delta,A]);At=s.Poly(A,t,domain=s.QQ);wu,wv=wp.diff(u),wp.diff(v)
    top=u**4/2+D*(v**3 if kind=='cusp' else u*v*v);remphase=poly(s.expand(W-top))
    @lru_cache(None)
    def step(e):
        nf,cu,cv=weighted.principal_division(kind,e,s.Rational(1,2),D)
        nf,cu,cv=[poly(s.cancel(A*z)) for z in [nf,cu,cv]]
        rem=nf-cu*remphase.diff(u)-cv*remphase.diff(v)-dp*(cu.diff(u)+cv.diff(v))
        assert (ap*poly(u**e[0]*v**e[1])-rem-dp*(cu.diff(u)+cv.diff(v))-wu*cu-wv*cv).is_zero
        assert all(f in index or weighted.weight(f,data['weights'])<=weighted.weight(e,data['weights'])-data['gap'] for f,c in rem.terms() if c)
        assert all(s.Poly(c,t).degree()+sum(f)<=sum(e)+k for f,c in rem.terms() if c)
        return rem,cu,cv
    @lru_cache(None)
    def reduce(e,depth):
        if e in index:
            out=[zt]*len(basis);out[index[e]]=At**depth;return tuple(out),zero,zero
        assert depth>0,(kind,e)
        rem,cu,cv=step(e);out=[zt]*len(basis);Cu=ap**(depth-1)*cu;Cv=ap**(depth-1)*cv
        for f,c in rem.terms():
            if not c:continue
            vec,x,y=reduce(f,depth-1);ct=s.Poly(c,t,domain=s.QQ);scalar=poly(c)
            out=[a+ct*b for a,b in zip(out,vec)];Cu+=scalar*x;Cv+=scalar*y
        return tuple(out),Cu,Cv
    def reduce_input(input,depth):
        out=[zt]*len(basis);Cu=zero;Cv=zero
        for e,c in poly(input).terms():
            if not c:continue
            vec,x,y=reduce(e,depth);ct=s.Poly(c,t,domain=s.QQ);scalar=poly(c)
            out=[a+ct*b for a,b in zip(out,vec)];Cu+=scalar*x;Cv+=scalar*y
        normal=sum((poly(z.as_expr()*u**e[0]*v**e[1]) for z,e in zip(out,basis)),zero)
        assert (ap**depth*poly(input)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wu*Cu-wv*Cv).is_zero
        return out,{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}
    columns=[]
    for e in basis:
        f=u**e[0]*v**e[1];input=source*f+delta*D*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))-delta*C*u*s.diff(f,v)/2
        maximum=6+weighted.weight(e,data['weights']);depth=(maximum-data['minimum_reducible'])//data['gap']+1
        values,witness=reduce_input(input,depth)
        assert all(z.degree()<=7+sum(e)+k*depth-sum(f) for z,f in zip(values,basis))
        columns.append({'depth':depth,'N':[str(z.as_expr()) for z in values],'witness':witness})
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        old_next=h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+D*old*R+delta*D*(bx*u*s.diff(R,u)+by*v*s.diff(R,v))
        next=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+source*R+delta*D*(bx*u*s.diff(R,u)+by*v*s.diff(R,v))-delta*C*u*s.diff(R,v)/2)
        Cu,Cv=s.Integer(0),-C*u*R/2
        assert s.expand(next-old_next-delta*(s.diff(Cu,u)+s.diff(Cv,v))-s.diff(W,u)*Cu-s.diff(W,v)*Cv)==0
        R=next;assert total(R,t,u,v)<=7*j
        assert all(weighted.weight(e,data['weights'])<=6*j for e,c in s.Poly(R,u,v).terms() if c)
        depth=6*j-4;bnd=7*j+k*depth
        values,witness=reduce_input(R,depth)
        assert all(z.degree()<=bnd-sum(e) for z,e in zip(values,basis))
        H=s.Poly(s.expand(h**j*A**depth),t,domain=s.QQ);assert H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=4*j
        seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'reduction_witness':witness,'gauge_witness':{'Cu':str(Cu),'Cv':str(Cv)}})
    return {'kind':kind,'singular_at_zero':singular,'zero_shear':zero_shear,'metric_determinant':metric,'W':str(W),'D':str(D),'A0':str(A),'basis':basis,'source':str(source),'C':str(C),'h':str(h),'beta_x':str(bx),'beta_y':str(by),'connection_columns':columns,'seed_columns':seeds}


def main():
    started=time.monotonic();support=coordinate_module_support();parent_checks=parent_universal();shear=tangent_universal();print('Universal conic support, parent Euler and tangent shear identities passed',flush=True)
    old=json.loads((P/'euler_seed_bounds_certificate.json').read_text())['conic_models'];K0=(t+1)*(5*t-1)/4-t*t
    phases=[(old[0]['label'],s.sympify(old[0]['W'],locals={'t':t,'y':y}),K0),(old[1]['label'],s.sympify(old[1]['W'],locals={'t':t,'y':y}),-2*t-7*t*t),(old[2]['label'],s.sympify(old[2]['W'],locals={'t':t,'y':y}),-2*t*t)]
    for label,K in [('constant',s.Integer(2)),('linear_regular',1+t),('linear_zero',t),('quadratic_regular',1+t+t*t),('double_zero',t*t),('repeated_root',(1-t)**2)]:
        W=1+t**3+sum((1+t**(3-k))*y**k+K**k*(2+t**(3-k))*y**(-k) for k in range(1,4))
        phases.append((label,s.expand(W),K))
    cm=[]
    for label,W,K in phases:cm.append(conic_model(label,W,K));print('Certified polynomial conic module',label,flush=True)
    z,T=s.symbols('z T');pm=[]
    for label,G in [('degree6',T**3),('degree5',z*T*T),('dense6',1+z+T+z*z+z*T+T*T+z**3+z*z*T+z*T*T+T**3)]:
        F=y/t-2-y*y;W=s.expand(t*t*(t*G.subs({z:y,T:F},simultaneous=True)-F*F/2));pm.append(parent_model(label,W));print('Certified parent Euler model',label,flush=True)
    tm=[]
    for singular in [False,True]:
        for zero_shear in [False,True]:tm.append(tangent_model(singular,zero_shear));print('Certified tangent shear',singular,zero_shear,flush=True)
    conic_cases=[conic_case(ell) for ell in range(3)]
    for d,expected in [(5,114),(6,210)]:assert sum(3*j+3*((d-1)*(j-1)+1) for j in range(1,d))==expected
    assert 19*15-8*5-6==239
    cover={'quadratic':7,'smooth':94,'nodal':227,'cusp':210,'transverse':220,'tangent':239,'triangle':216,'concurrent':104,'repeated':99};assert max(cover.values())==239
    retained=['euler_seed_bounds_certificate.json','euler_seed_bounds_verification.json','shared_divisor_weighted_seeds_verification.json','joint_degree_seed_bounds_verification.json','mixed_period_infinity_verification.json']
    cert={'coordinate_module_support':support,'conic_cases':conic_cases,'conic_models':cm,'parent_universal':parent_checks,'parent_models':pm,'tangent_universal':shear,'tangent_models':tm,'global_cover':cover,'retained_inputs':{n:sha(P/n) for n in retained}}
    f=P/'structured_seed_modules_certificate.json';f.write_text(json.dumps(cert,indent=2)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,
            'coordinate_support_cases':len(support),'conic_K_degree_cases':3,'conic_models':len(cm),'conic_seed_columns':sum(len(x['seed_columns']) for x in cm),'parent_symbolic_connection_checks':sum(len(x['connection_identity_checks']) for x in parent_checks),'parent_models':len(pm),'parent_seed_columns':sum(len(x['seed_columns']) for x in pm),'tangent_symbolic_connection_checks':len(shear['connection_identity_checks']),'tangent_models':len(tm),'tangent_connection_columns':sum(len(x['connection_columns']) for x in tm),'tangent_seed_columns':sum(len(x['seed_columns']) for x in tm),
            'conic_gates_by_K_degree':[216,199,108],'conic_full_gate':216,'mixed_parent_gate':210,'full_W_linear_gate':216,'tangent_complement_gate':239,'tangent_full_gate':239,'cusp_full_gate':210,
            'global_cover':cover,'global_sufficient_pure_prefix':239,'least_prefix_interval':[20,239],'odd_least_prefix_interval':[20,238],'all_homogeneous_sufficient_pure_prefix':216,'all_homogeneous_least_prefix_interval':[14,216],'binary_search_rounds_max':8,
            'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,
            'boundary':'Exact coordinate support/divisibility, per-cancellation polynomial Laurent module checks, parent Euler and tangent shear identities plus model seeds. Full coefficient/rank cover, all-degree reduction induction and retained analytic nonzero zeroth coefficients remain written and independently unreviewed. No expanded global radical cofactors or leastness computation.'}
    (P/'structured_seed_modules_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
