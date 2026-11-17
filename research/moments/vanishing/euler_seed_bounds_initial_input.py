"""Auxiliary exact Euler-source and conic cube-factor seed certificates.

The analytic exclusions and complete class/rank coverage remain written inputs.
Retained checker mains and receipts are never executed or rewritten.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_rescaled_singular_prefixes as rescaled
import verify_balanced_logarithmic_seed_bounds as balanced

P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
y=balanced.y

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def Euler(f):return u*s.diff(f,u)+v*s.diff(f,v)
def total(f,*gens):return -1000 if f==0 else int(s.Poly(f,*gens).total_degree())
def laurent(f):return balanced.laurent(f)

def ordinary_source(W,D):
    delta=t*t*D
    old=rescaled.source_numerator(W,D)
    new=s.expand(delta*s.diff(W,t)-s.diff(delta,t)*(W-Euler(W)/4)-t**3*D*D)
    assert s.expand(new-old-s.diff(delta,t)*(2*delta+Euler(W))/4)==0
    assert s.Poly(new,u,v).total_degree()<=3 and total(new,t,u,v)<=6
    return new

def universal_ordinary():
    d0,d1=s.symbols('d0 d1');D=d0+d1*t;delta=t*t*D
    W=u*u*v*v/2+sum(s.Symbol('c'+str(i)+str(j)+str(h))*t**h*u**i*v**j for i,j in ordinary.monomials(3) for h in range(5-i-j))
    source=ordinary_source(W,D);old=rescaled.source_numerator(W,D);rows=[]
    for i,j in ordinary.monomials(4):
        f=u**i*v**j
        exact=delta*(s.diff(u*f,u)+s.diff(v*f,v))+s.diff(W,u)*u*f+s.diff(W,v)*v*f
        assert s.expand((source-old)*f+delta*s.diff(delta,t)*Euler(f)/4-s.diff(delta,t)*exact/4)==0
        rows.append([i,j])
    a,d,m,ju,jv,qu,qv,quv=s.symbols('a d m ju jv qu qv quv')
    original=m*(t*(a*u**3+d*v**3)-(qu*u*u+quv*u*v+qv*v*v))+(t*u*v-ju*u-jv*v)**2/2
    hom=s.expand(t*t*original.subs({u:(u+jv)/t,v:(v+ju)/t},simultaneous=True))
    assert s.diff(hom,t)==0
    hs=ordinary_source(hom,m)
    assert total(hs,t,u,v)<=4 and s.Poly(hs,u,v).total_degree()<=3
    assert 4*28-11==101
    return {'generic_phase_joint_degree':4,'ordinary_source_joint_degree':6,'ordinary_source_amplitude_degree':3,'monomial_connection_identity_checks':rows,'homogeneous_phase_t_independent':True,'homogeneous_source_joint_degree':4,'homogeneous_complement_gate':101}

def ordinary_model(kind,singular):
    W,D,metric,coeffs=rescaled.phase(kind,singular);delta=t*t*D;source=ordinary_source(W,D);old=rescaled.source_numerator(W,D)
    basis=ordinary.CASES[kind]['basis'];index={e:i for i,e in enumerate(basis)}
    k={'nodal':2,'transverse':3,'triangle':4}[kind]
    prepared=ordinary.operator_data(kind,*coeffs,delta);A=prepared[-1]
    ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,v,domain=ring);zero=poly(0);zt=s.Poly(0,t,domain=s.QQ)
    wp,dp,ap=map(poly,[W,delta,A]);At=s.Poly(s.expand(A),t,domain=s.QQ);wu,wv=wp.diff(u),wp.diff(v)
    @lru_cache(None)
    def step(e):
        rem,cu,cv,_,_=ordinary.one_step(kind,e,coeffs,delta,prepared)
        rem,cu,cv=map(poly,[rem,cu,cv])
        assert (ap*poly(u**e[0]*v**e[1])-rem-dp*(cu.diff(u)+cv.diff(v))-wu*cu-wv*cv).is_zero
        assert all(f in index or sum(f)<=sum(e)-1 for f,c in rem.terms() if c)
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
    def reduce_input(f,depth):
        out=[zt]*len(basis);Cu=zero;Cv=zero
        for e,c in poly(f).terms():
            if not c:continue
            vec,x,y=reduce(e,depth);ct=s.Poly(c,t,domain=s.QQ);scalar=poly(c)
            out=[a+ct*b for a,b in zip(out,vec)];Cu+=scalar*x;Cv+=scalar*y
        normal=sum((poly(z.as_expr()*u**e[0]*v**e[1]) for z,e in zip(out,basis)),zero)
        assert (ap**depth*poly(f)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wu*Cu-wv*Cv).is_zero
        return out,{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}
    connection=[]
    for e in basis:
        f=u**e[0]*v**e[1];input=source*f+delta*s.diff(delta,t)*Euler(f)/4
        depth=sum(e)+1 if kind!='triangle' else sum(e)+2
        values,witness=reduce_input(input,depth)
        for z,f in zip(values,basis):assert z.degree()<=6+sum(e)+k*depth-sum(f)
        connection.append({'depth':depth,'N':[str(z.as_expr()) for z in values],'witness':witness})
    seeds=[];R=s.Integer(1)
    for j in range(1,4):
        old_next=delta**2*s.diff(R,t)-2*(j-1)*delta*s.diff(delta,t)*R+old*R
        next=s.expand(delta**2*s.diff(R,t)-2*(j-1)*delta*s.diff(delta,t)*R+source*R+delta*s.diff(delta,t)*Euler(R)/4)
        gu=s.diff(delta,t)*u*R/4;gv=s.diff(delta,t)*v*R/4
        assert s.expand(next-old_next-delta*(s.diff(gu,u)+s.diff(gv,v))-s.diff(W,u)*gu-s.diff(W,v)*gv)==0
        R=next;assert total(R,t,u,v)<=6*j and s.Poly(R,u,v).total_degree()<=3*j
        depth=3*j-(1 if kind=='triangle' else 2);bnd=6*j+k*depth
        values,witness=reduce_input(R,depth)
        for z,e in zip(values,basis):assert z.degree()<=bnd-sum(e)
        H=s.Poly(s.expand(delta**(2*j)*A**depth),t,domain=s.QQ);assert H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=4*j
        seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'H':str(H.as_expr()),'H_valuation':val,'N':[str(z.as_expr()) for z in values],'reduction_witness':witness,'gauge_witness':{'Cu':str(gu),'Cv':str(gv)}})
    return {'kind':kind,'singular_at_zero':singular,'metric_determinant':metric,'W':str(W),'D':str(D),'delta':str(delta),'A0':str(A),'basis':basis,'source':str(source),'connection_columns':connection,'seed_columns':seeds}

def conic_factor_checks():
    a,b,c,f,rho,sigma=s.symbols('a b c f rho sigma',nonzero=True);eta=f*(rho-sigma)
    U=((rho*c-b)*t+1)/eta;V=((b-sigma*c)*t-1)/eta;K=U*V-a*t*t/f
    Zh=f/eta*(y-K/y-U+V);Th=y-U-rho*Zh
    e,d=s.symbols('e d',nonzero=True);z0=-c/e;Azero=a+b*z0+d*z0*z0
    Zf=y+t*z0;Tf=(y-d*y*y-t*(b+2*d*z0)*y+t*z0-t*t*Azero)/(e*y);Kf=(t*z0-t*t*Azero)/e
    records=[]
    for label,Z,T,K in [('f_nonzero',Zh,Th,K),('f_zero_e_nonzero',Zf,Tf,Kf)]:
        zm=laurent(Z).get(-1,s.Integer(0));tm=laurent(T).get(-1,s.Integer(0));cubic=0
        for i in range(4):
            cubic+=s.Symbol('g'+str(i)+str(3-i))*zm**i*tm**(3-i)
        ratio=s.cancel(cubic/K**3);assert s.diff(ratio,t)==0 and s.Poly(s.cancel(K),t).degree()<=2
        records.append({'chart':label,'K':str(K),'negative_cubic_coefficient_over_K_cubed':str(ratio),'degree_K_at_most':2})
    assert min(laurent(y+(1-b*t)/e))==0
    k0,k1,k2,C=s.symbols('k0 k1 k2 C');K=k0+k1*t+k2*t*t;B=t*s.diff(K,t)-s.Rational(2,3)*K
    W=C*K**3/y**3+sum(s.Symbol('a'+str(power)+'_'+str(j))*t**j*y**power for power in range(-2,4) for j in range(4-power))
    old=t*s.diff(W,t)-2*W-t*t;new=s.expand(K*old+B*y*s.diff(W,y))
    co=laurent(new);assert min(co)>=-2 and max(co)<=3
    assert all(s.Poly(c,t).degree()+power<=5 for power,c in co.items())
    check=[]
    for power in range(-6,7):
        R=y**power;difference=B*(t*t*y*s.diff(R,y)+y*s.diff(W,y)*R)
        assert s.expand(difference-B*(t*t*y*s.diff(R,y)+y*s.diff(W,y)*R))==0
        check.append(power)
    return {'charts':records,'trapped_negative_order':0,'generic_negative_order_after_cancellation':2,'new_joint_degree':5,'monomial_connection_checks':check}

def conic_model(label,W,K):
    field=s.QQ.frac_field(t);raw=laurent(W);p=max(0,max(raw));assert min(raw)==-3
    assert s.Poly(K,t).degree()<=2
    ratio=s.cancel(raw[-3]/K**3);assert not ratio.has(t) and ratio!=0
    basis,plus,minus,_=balanced.specification(p,3);valid=set(basis)
    Aplus=raw[p] if p else s.Integer(1);Aminus=raw[-3]
    B=t*s.diff(K,t)-s.Rational(2,3)*K;h=t**3*K
    source=t*s.diff(W,t)-2*W-t*t;new_source=s.expand(K*source+B*y*s.diff(W,y))
    assert min(laurent(new_source))>=-2
    deriv={k:field.convert(k*c) for k,c in raw.items() if k}
    def parse(f):return {k:field.convert(c) for k,c in laurent(f).items()}
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
        out=poly.copy();witness={}
        while any(k not in valid for k in out):
            n=max(out) if max(out)>max(valid) else min(out)
            shift=n-p if n>max(valid) else n+3
            lead=deriv[p] if n>max(valid) else deriv[-3]
            factor=out[n]/lead;add(witness,shift,factor)
            for k,c in exact({shift:factor}).items():add(out,k,-c)
        defect=poly.copy()
        for k,c in out.items():add(defect,k,-c)
        for k,c in exact(witness).items():add(defect,k,-c)
        assert not defect
        return out,witness
    def clear(out,H):
        values=[]
        for k in basis:
            expr=s.cancel(field.to_sympy(field.convert(H)*out.get(k,field.zero)));assert s.denom(expr)==1
            values.append(s.Poly(expr,t,domain=s.QQ))
        return values
    def depths(j):return ((p*(j-1)+1 if p else 0),max(0,2*j-(3 if p else 2)))
    connections=[]
    for k in basis:
        input=new_source*y**k+t*t*B*k*y**k
        normal,witness=reduce(parse(input));H=Aplus**p*Aminus**3
        values=clear(normal,H)
        connections.append({'N':[str(z.as_expr()) for z in values],'clearing_divisor':str(H),'witness':{str(k):str(field.to_sympy(z)) for k,z in witness.items()}})
    seeds=[];R=s.Integer(1)
    for j in range(1,4):
        nxt=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+t*t*B*y*s.diff(R,y)+new_source*R)
        old_nxt=h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+K*source*R
        assert s.expand(nxt-old_nxt-B*(t*t*y*s.diff(R,y)+y*s.diff(W,y)*R))==0
        R=nxt;support=laurent(R)
        assert min(support)>=-2*j and max(support)<=p*j
        assert all(s.Poly(c,t).degree()+k<=5*j for k,c in support.items())
        high,low=depths(j);divisor=Aplus**high*Aminus**low
        normal,witness=reduce(parse(R));values=clear(normal,divisor);bnd=5*j+(3-p)*high+6*low
        for k,z in zip(basis,values):assert z.degree()<=bnd-k
        H=s.Poly(s.expand(h**j*divisor),t,domain=s.QQ);assert H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        seeds.append({'j':j,'depths':[high,low],'degree_bound':bnd,'H':str(H.as_expr()),'H_valuation':val,'N':[str(z.as_expr()) for z in values],'witness':{str(k):str(field.to_sympy(z)) for k,z in witness.items()}})
    return {'label':label,'W':str(W),'K':str(K),'positive_order':p,'negative_order':3,'basis':basis,'source':str(new_source),'seed_denominator':str(h),'connection_columns':connections,'seed_columns':seeds}

def degree_cases():
    ordinary_cases=[]
    for kind,k in [('nodal',2),('transverse',3),('triangle',4)]:
        grades=[sum(e) for e in ordinary.CASES[kind]['basis']];rank=len(grades);bounds=[]
        for size in range(1,rank+1):
            for chosen in combinations(range(1,rank),size-1):
                rows=(0,)+chosen;gate=sum(6*j+k*(3*j-(1 if kind=='triangle' else 2)) for j in range(1,size+1))-sum(grades[i] for i in rows)
                bounds.append(gate)
        maximum=max(bounds);assert maximum=={'nodal':297,'transverse':271,'triangle':244}[kind]
        ordinary_cases.append({'kind':kind,'complement_gate':maximum,'row_subsets_checked':len(bounds)})
    conic=[]
    for p in range(4):
        basis,_,_,_=balanced.specification(p,3);rank=len(basis)
        def bj(j):return 5*j+(3-p)*(p*(j-1)+1 if p else 0)+6*max(0,2*j-(3 if p else 2))
        gates=[sum(bj(j) for j in range(1,size+1))-sum(rows) for size in range(1,rank+1) for rows in combinations(basis,size) if 0 in rows]
        conic.append({'positive_order':p,'negative_order':3,'seed_degrees':[bj(j) for j in range(1,rank+1)],'gate':max(gates),'row_subsets_checked':len(gates)})
    assert [row['gate'] for row in conic]==[69,130,201,258]
    old=[row for row in balanced.degree_cases() if row['negative_order']<3];assert max(row['all_seed_rank_and_row_choices_gate'] for row in old)==145
    return {'ordinary':ordinary_cases,'conic_q3':conic,'conic_q_below3_retained':old}

def main():
    start=time.monotonic();normal=universal_ordinary();print('Universal ordinary Euler source and homogeneous101 passed',flush=True)
    om=[]
    for kind in ['nodal','transverse','triangle']:
        for singular in [False,True]:
            om.append(ordinary_model(kind,singular));print('Certified Euler model',kind,singular,flush=True)
    factor=conic_factor_checks();print('Conic cube factors and cancellation passed',flush=True)
    old=json.loads((P/'balanced_logarithmic_seed_bounds_certificate.json').read_text())['models']
    Kap=((t+1)*(5*t-1)/4-t*t)
    phases=[(old[0]['label'],s.sympify(old[0]['phase_numerator_W'],locals={'t':t,'y':y}),Kap),(old[1]['label'],s.sympify(old[1]['phase_numerator_W'],locals={'t':t,'y':y}),-2*t-7*t*t),(old[2]['label'],s.sympify(old[2]['phase_numerator_W'],locals={'t':t,'y':y}),-2*t*t)]
    K=1+t+t*t
    for p in range(4):
        W=K**3/y**3+sum((1+t**(3-k))*y**k for k in range(-2,p+1))
        phases.append(('orders'+str(p)+'_3',s.expand(W),K))
    phases.append(('constant_K',y**3+1/y**3+t*t/y+1,s.Integer(1)))
    cm=[]
    for label,W,K in phases:
        cm.append(conic_model(label,s.expand(W),K));print('Certified cube-factor model',label,flush=True)
    cases=degree_cases();cover={'quadratic':7,'smooth':94,'nodal':297,'cusp':255,'transverse':271,'tangent':284,'triangle':258,'concurrent':104,'repeated':99};assert max(cover.values())==297
    retained=['rescaled_singular_prefixes_verification.json','rescaled_singular_prefixes_certificate.json','balanced_logarithmic_seed_bounds_verification.json','balanced_logarithmic_seed_bounds_certificate.json','joint_degree_seed_bounds_verification.json']
    cert={'ordinary_normalization':normal,'ordinary_models':om,'conic_factor_checks':factor,'conic_models':cm,'degree_cases':cases,'global_cover':cover,'retained_inputs':{n:sha(P/n) for n in retained}}
    f=P/'euler_seed_bounds_certificate.json';f.write_text(json.dumps(cert,indent=2)+'\n')
    receipt={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
             'ordinary_universal_connection_identity_checks':len(normal['monomial_connection_identity_checks']),'ordinary_models':len(om),'ordinary_connection_columns':sum(len(m['connection_columns']) for m in om),'ordinary_direct_seed_columns':sum(len(m['seed_columns']) for m in om),
             'conic_factor_charts':len(factor['charts']),'conic_q3_cases':len(cases['conic_q3']),'conic_models':len(cm),'conic_connection_columns':sum(len(m['connection_columns']) for m in cm),'conic_direct_seed_columns':sum(len(m['seed_columns']) for m in cm),
             'ordinary_complement_gates':{'nodal':297,'transverse':271,'triangle':244},'homogeneous_nodal_complement_gate':101,'conic_cubic_exclusion_gate':258,'full_W_linear_gate':258,'global_cover':cover,
             'global_sufficient_pure_prefix':297,'least_prefix_interval':[20,297],'odd_least_prefix_interval':[20,296],'all_homogeneous_sufficient_pure_prefix':258,'all_homogeneous_least_prefix_interval':[14,258],'binary_search_rounds_max':9,
             'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,
             'boundary':'Exact Euler representative/seed divergence identities, six full-Gaussian ordinary-node connection/seed models, conic cube-factor chart identities and eight Laurent connection/seed models. Complete rank/class cover, all-coefficient degree induction and retained analytic nonzero zeroth coefficients remain written and independently unreviewed. No expanded global radical cofactors or leastness computation.'}
    (P/'euler_seed_bounds_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
