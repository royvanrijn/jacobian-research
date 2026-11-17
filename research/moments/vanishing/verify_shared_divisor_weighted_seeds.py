"""Exact auxiliary shared-D and weighted-Euler seed certificates.

The full rank/class and analytic seed exclusions remain retained written inputs.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_rescaled_singular_prefixes as rescaled
import verify_euler_seed_bounds as euler

P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def total(f,*gens):return -1000 if f==0 else int(s.Poly(f,*gens).total_degree())

def universal_shared_divisions():
    D,Le,Lf,J,L,N,K,delta=s.symbols('D Le Lf J L N K delta')
    a,d=s.symbols('a d',nonzero=True);C=s.Rational(1,2)
    pars=[D,Le,Lf,J,L,N,K,delta];weights=[1,1,1,2,3,3,4,3]
    def degree(f):
        return max((e[0]+e[1]+sum(x*y for x,y in zip(e[2:],weights)) for e,c in s.Poly(f,u,v,*pars).terms() if c),default=-1000)
    out=[]
    for kind,k in [('nodal',1),('transverse',2),('triangle',3)]:
        co=[C,a*D if kind!='triangle' else 0,d*D if kind=='nodal' else 0,D*Le,D*Lf,J,L,N,K]
        prepared=ordinary.operator_data(kind,*co,delta);oldA=prepared[-1]
        A=D if kind=='nodal' else D*Lf if kind=='transverse' else D*Le*Lf
        rows=[]
        for exponent in ordinary.monomials(8):
            rem,cu,cv,W,_=ordinary.one_step(kind,exponent,co,delta,prepared)
            rem,cu,cv=[s.cancel(f*A/oldA) for f in [rem,cu,cv]]
            # Polynomial in D and the axis-linear forms, allowing constant coefficient field C(a,d).
            for f in [rem,cu,cv]:s.Poly(f,u,v,*pars)
            assert s.expand(A*u**exponent[0]*v**exponent[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
            assert degree(rem)<=sum(exponent)+k
            assert degree(cu)<=sum(exponent)+k-3 and degree(cv)<=sum(exponent)+k-3
            assert all(e in ordinary.CASES[kind]['basis'] or sum(e)<=sum(exponent)-1 for e,c in s.Poly(rem,u,v).terms() if c)
            rows.append({'exponent':exponent,'joint_degree':degree(rem),'bound':sum(exponent)+k})
        out.append({'kind':kind,'A0':str(A),'joint_cost':k,'identities':rows})
    checkpoint=json.loads((P/'ordinary_node_phase_shift_checkpoint.json').read_text())
    W=rescaled.rescale(s.sympify(checkpoint['V'],locals={'t':t,'u':u,'v':v}));den=s.sympify(checkpoint['D'],locals={'t':t})
    pv=s.Poly(W,u,v);axis={}
    for label,mon,deg in [('u_cubic',u**3,0),('v_cubic',v**3,0),('u_quadratic',u*u,1),('v_quadratic',v*v,1)]:
        quotient=s.cancel(pv.coeff_monomial(mon)/den);assert s.Poly(quotient,t).degree()<=deg
        axis[label]={'coefficient_over_D':str(quotient),'degree_bound':deg}
    return {'divisions':out,'universal_axis_factor_checks':axis}

def shared_model_replay():
    old=json.loads((P/'euler_seed_bounds_certificate.json').read_text());out=[]
    local={'t':t,'u':u,'v':v};read=lambda f:s.sympify(f,locals=local)
    for model in old['ordinary_models']:
        kind=model['kind'];k={'nodal':1,'transverse':2,'triangle':3}[kind]
        W,D,delta,oldA,source=[read(model[n]) for n in ['W','D','delta','A0','source']]
        pv=s.Poly(W,u,v);Le=s.cancel(pv.coeff_monomial(u*u)/D);Lf=s.cancel(pv.coeff_monomial(v*v)/D)
        assert s.Poly(Le,t).degree()<=1 and s.Poly(Lf,t).degree()<=1
        A=D if kind=='nodal' else D*Lf if kind=='transverse' else D*Le*Lf
        assert A!=0 and s.Poly(A,t).degree()<=k
        ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,v,domain=ring)
        wp,dp,ap=map(poly,[W,delta,A]);wu,wv=wp.diff(u),wp.diff(v)
        def exquo(expr,divisor):
            return s.Poly(expr,u,v,t,domain=s.QQ).exquo(s.Poly(divisor,u,v,t,domain=s.QQ)).as_expr()
        def replay(input,row,depth):
            scale=s.cancel((oldA/A)**depth);assert not s.denom(scale).has(t)
            values=[s.Poly(exquo(read(f),scale),t,domain=s.QQ) for f in row['N']]
            witness=row.get('reduction_witness',row.get('witness'))
            Cu,Cv=[poly(exquo(read(witness[n]),scale)) for n in ['Cu','Cv']]
            normal=sum((poly(z.as_expr()*u**e[0]*v**e[1]) for z,e in zip(values,model['basis'])),poly(0))
            assert (ap**depth*poly(input)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wu*Cu-wv*Cv).is_zero
            return values,{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}
        columns=[]
        for exponent,row in zip(model['basis'],model['connection_columns']):
            f=u**exponent[0]*v**exponent[1];input=source*f+delta*s.diff(delta,t)*euler.Euler(f)/4
            values,witness=replay(input,row,row['depth'])
            assert all(z.degree()<=6+sum(exponent)+k*row['depth']-sum(e) for z,e in zip(values,model['basis']))
            columns.append({'depth':row['depth'],'N':[str(z.as_expr()) for z in values],'witness':witness})
        R=s.Integer(1);seeds=[]
        for row in model['seed_columns']:
            j=row['j'];R=s.expand(delta**2*s.diff(R,t)-2*(j-1)*delta*s.diff(delta,t)*R+delta*s.diff(delta,t)*euler.Euler(R)/4+source*R)
            depth=row['depth'];values,witness=replay(R,row,depth);bnd=6*j+k*depth
            assert all(z.degree()<=bnd-sum(e) for z,e in zip(values,model['basis']))
            H=s.Poly(s.expand(delta**(2*j)*A**depth),t,domain=s.QQ);assert H.degree()<=bnd
            val=min(e[0] for e,c in H.terms() if c);assert val>=4*j
            seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'witness':witness})
        out.append({'kind':kind,'singular_at_zero':model['singular_at_zero'],'metric_determinant':model['metric_determinant'],'W':str(W),'D':str(D),'A0':str(A),'old_A0_over_new_A0':str(s.cancel(oldA/A)),'basis':model['basis'],'connection_columns':columns,'seed_columns':seeds})
    return out

def weighted_source(kind,W,D):
    delta=t*t*D;old=rescaled.source_numerator(W,D)
    bx=s.diff(delta,t)/4
    by=2*t*D/3 if kind=='cusp' else t*D-s.diff(delta,t)/8
    new=s.expand(old+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v))
    assert total(new,t,u,v)<=6
    data=weighted.CASES[kind]
    assert all(weighted.weight(e,data['weights'])<=data['lower_weight'] for e,c in s.Poly(new,u,v).terms() if c)
    return new,bx,by

def universal_weighted_sources():
    d0,d1=s.symbols('d0 d1');D=d0+d1*t;delta=t*t*D;rows=[]
    for kind in ['cusp','tangent']:
        data=weighted.CASES[kind];W=u**4/2+D*(v**3 if kind=='cusp' else u*v*v)
        W+=sum(s.Symbol('a'+str(i)+'_'+str(j)+'_'+str(h))*t**h*u**i*v**j for i,j in weighted.monomials(data['lower_weight'],data['weights']) for h in range(5-i-j))
        source,bx,by=weighted_source(kind,W,D);old=rescaled.source_numerator(W,D);checks=[]
        for i,j in ordinary.monomials(4):
            f=u**i*v**j;Cu=bx*u*f;Cv=by*v*f
            diff=(source-old)*f+delta*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))
            assert s.expand(diff-delta*(s.diff(Cu,u)+s.diff(Cv,v))-s.diff(W,u)*Cu-s.diff(W,v)*Cv)==0
            checks.append([i,j])
        rows.append({'kind':kind,'beta_x':str(bx),'beta_y':str(by),'source_joint_degree':6,'source_weight_at_most':data['lower_weight'],'connection_identity_checks':checks})
    return rows

def weighted_model(kind,singular):
    W,D,metric,_=rescaled.phase(kind,singular);delta=t*t*D;source,bx,by=weighted_source(kind,W,D);old=rescaled.source_numerator(W,D)
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
        f=u**e[0]*v**e[1];input=source*f+delta*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))
        maximum=data['lower_weight']+weighted.weight(e,data['weights']);depth=(maximum-data['minimum_reducible'])//data['gap']+1
        values,witness=reduce_input(input,depth)
        assert all(z.degree()<=6+sum(e)+k*depth-sum(f) for z,f in zip(values,basis))
        columns.append({'depth':depth,'N':[str(z.as_expr()) for z in values],'witness':witness})
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        old_next=delta**2*s.diff(R,t)-2*(j-1)*delta*s.diff(delta,t)*R+old*R
        next=s.expand(delta**2*s.diff(R,t)-2*(j-1)*delta*s.diff(delta,t)*R+source*R+delta*(bx*u*s.diff(R,u)+by*v*s.diff(R,v)))
        Cu,Cv=bx*u*R,by*v*R
        assert s.expand(next-old_next-delta*(s.diff(Cu,u)+s.diff(Cv,v))-s.diff(W,u)*Cu-s.diff(W,v)*Cv)==0
        R=next;assert total(R,t,u,v)<=6*j
        assert all(weighted.weight(e,data['weights'])<=data['lower_weight']*j for e,c in s.Poly(R,u,v).terms() if c)
        depth=(5*j-3 if kind=='cusp' else 7*j-4);bnd=6*j+k*depth
        values,witness=reduce_input(R,depth)
        assert all(z.degree()<=bnd-sum(e) for z,e in zip(values,basis))
        H=s.Poly(s.expand(delta**(2*j)*A**depth),t,domain=s.QQ);assert H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=4*j
        seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'reduction_witness':witness,'gauge_witness':{'Cu':str(Cu),'Cv':str(Cv)}})
    return {'kind':kind,'singular_at_zero':singular,'metric_determinant':metric,'W':str(W),'D':str(D),'A0':str(A),'basis':basis,'source':str(source),'beta_x':str(bx),'beta_y':str(by),'connection_columns':columns,'seed_columns':seeds}

def degree_cases():
    rows=[]
    specifications=[('nodal',ordinary.CASES['nodal']['basis'],9,2,227),('transverse',ordinary.CASES['transverse']['basis'],12,4,220),('triangle',ordinary.CASES['triangle']['basis'],15,3,204),('cusp',weighted.CASES['cusp']['basis'],11,3,204),('tangent',weighted.CASES['tangent']['basis'],20,8,254)]
    for kind,basis,slope,offset,expected in specifications:
        grades=[sum(e) for e in basis];rank=len(basis);cases=[]
        for size in range(1,rank+1):
            for chosen in combinations(range(1,rank),size-1):
                selected=(0,)+chosen;gate=sum(slope*j-offset for j in range(1,size+1))-sum(grades[i] for i in selected)
                cases.append({'rank':size,'rows':selected,'gate':gate})
        assert max(x['gate'] for x in cases)==expected
        rows.append({'kind':kind,'slope':slope,'offset':offset,'gate':expected,'cases':cases})
    return rows

def main():
    started=time.monotonic();shared=universal_shared_divisions();print('135 symbolic shared-D divergence identities and four axis factors passed',flush=True)
    om=shared_model_replay();print('Six tightened ordinary connection/seed models replayed',flush=True)
    sources=universal_weighted_sources();print('30 symbolic weighted Euler connection identities passed',flush=True)
    wm=[]
    for kind in ['cusp','tangent']:
        for singular in [False,True]:
            wm.append(weighted_model(kind,singular));print('Certified weighted Euler model',kind,singular,flush=True)
    cases=degree_cases();cover={'quadratic':7,'smooth':94,'nodal':258,'cusp':255,'transverse':258,'tangent':255,'triangle':258,'concurrent':104,'repeated':99};assert max(cover.values())==258
    retained=['ordinary_node_phase_shift_checkpoint.json','euler_seed_bounds_certificate.json','euler_seed_bounds_verification.json','joint_degree_seed_bounds_verification.json','rescaled_singular_prefixes_certificate.json']
    cert={'shared_divisions':shared,'ordinary_models':om,'weighted_sources':sources,'weighted_models':wm,'degree_cases':cases,'global_cover':cover,'retained_inputs':{n:sha(P/n) for n in retained}}
    f=P/'shared_divisor_weighted_seeds_certificate.json';f.write_text(json.dumps(cert,indent=2)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,
            'shared_symbolic_divergence_identities':135,'axis_factor_checks':4,'ordinary_models_replayed':len(om),'ordinary_connection_columns_replayed':sum(len(m['connection_columns']) for m in om),'ordinary_seed_columns_replayed':sum(len(m['seed_columns']) for m in om),
            'weighted_symbolic_connection_identities':sum(len(x['connection_identity_checks']) for x in sources),'weighted_models':len(wm),'weighted_connection_columns':sum(len(m['connection_columns']) for m in wm),'weighted_seed_columns':sum(len(m['seed_columns']) for m in wm),
            'ordinary_complement_gates':{'nodal':227,'transverse':220,'triangle':204},'ordinary_full_gates':{'nodal':258,'transverse':258,'triangle':258},'cusp_complement_gate':204,'cusp_full_gate':255,'tangent_complement_gate':254,'tangent_full_gate':255,
            'global_cover':cover,'global_sufficient_pure_prefix':258,'least_prefix_interval':[20,258],'odd_least_prefix_interval':[20,258],'all_homogeneous_sufficient_pure_prefix':258,'all_homogeneous_least_prefix_interval':[14,258],'binary_search_rounds_max':8,
            'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,
            'boundary':'Exact135 shared-D symbolic divergence identities, universal axis factorization, six tighter ordinary connection/seed model replays,30 weighted Euler source identities and four full-Gaussian weighted connection/seed models. Complete class/rank coverage, all-coefficient degree induction and retained analytic nonzero zeroth coefficients remain written and independently unreviewed. No expanded global radical cofactors or numerical leastness computation.'}
    (P/'shared_divisor_weighted_seeds_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
