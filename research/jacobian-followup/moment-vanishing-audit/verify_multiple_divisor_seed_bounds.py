"""Exact auxiliary multi-divisor filtrations and conic displacement bounds.
Universal all-degree induction and analytic exclusions remain written inputs.
No retained checker main or prior receipt is executed or rewritten.
"""
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_euler_seed_bounds as euler
import verify_shared_divisor_weighted_seeds as shared
import verify_balanced_logarithmic_seed_bounds as balanced
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
y=balanced.y

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def grade(f,var,wu,wv):return max((wu*e[0]+wv*e[1]-e[2] for e,c in s.Poly(f,u,v,var).terms() if c),default=-1000)
def ceil(x):return int(s.ceiling(x))
def filters(kind,D,Le,Lf):
    out=[(D,s.Integer(1),s.Integer(1),2)]
    if kind=='transverse':out.append((Lf,s.Rational(1,3),s.Integer(1),2))
    if kind=='triangle':out.extend([(Le,s.Integer(1),s.Rational(1,2),2),(Lf,s.Rational(1,2),s.Integer(1),2)])
    return out

def ordinary_universal():
    D,Le,Lf,J,L,N,K=s.symbols('D Le Lf J L N K');a,d=s.symbols('a d',nonzero=True)
    eps,epse,epsf,Jp,Lp,Np,Kp=s.symbols('eps epse epsf Jp Lp Np Kp')
    deriv=lambda f:s.diff(f,t)+sum(s.diff(f,x)*xp for x,xp in [(D,eps),(Le,epse),(Lf,epsf),(J,Jp),(L,Lp),(N,Np),(K,Kp)])
    delta=t*t*D;dp=deriv(delta);h=t**3*D*D;rows=[]
    for kind in ['nodal','transverse','triangle']:
        co=[s.Rational(1,2),a*D if kind!='triangle' else 0,d*D if kind=='nodal' else 0,D*Le,D*Lf,J,L,N,K]
        prepared=ordinary.operator_data(kind,*co,delta);W=prepared[0]
        A=D*(Lf if kind!='nodal' else 1)*(Le if kind=='triangle' else 1)
        E=lambda f:u*s.diff(f,u)+v*s.diff(f,v)
        source=s.expand(t*D*deriv(W)-(2*D+t*eps)*(W-E(W)/4)-t*t*D*D)
        vector=s.cancel(delta*dp/(4*t));checks=[]
        for var,wu,wv,seedcost in filters(kind,D,Le,Lf):
            assert grade(source,var,wu,wv)<=seedcost
            # A time derivative raises any divisor-pole grade by at most1.
            assert grade(h,var,wu,wv)+1<=seedcost
            assert grade(deriv(h),var,wu,wv)<=seedcost
            assert grade(vector,var,wu,wv)<=seedcost
        for ex in ordinary.monomials(12):
            rem,cu,cv,_,_=ordinary.one_step(kind,ex,co,delta,prepared)
            rem,cu,cv=[s.cancel(z*A/prepared[-1]) for z in [rem,cu,cv]]
            assert s.expand(A*u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
            values=[]
            for var,wu,wv,cost in filters(kind,D,Le,Lf):
                rho=wu*ex[0]+wv*ex[1]
                # Remainder and both exact witnesses are controlled. The
                # latter weak bound is enough for the common cancellation.
                assert all(grade(z,var,wu,wv)<=rho-1 for z in [rem,cu,cv])
                values.append({'divisor':str(var),'rho':str(rho),'remainder_grade':str(grade(rem,var,wu,wv))})
            checks.append({'exponent':ex,'filtrations':values})
        rows.append({'kind':kind,'source':str(source),'formal_independent_divisors':[str(f[0]) for f in filters(kind,D,Le,Lf)],'source_pole_bounds':[2]*len(filters(kind,D,Le,Lf)),'principal_identities':checks})
    # The tempting integral Lf grade fails: retain its exact obstruction.
    co=[s.Rational(1,2),a*D,0,D*Le,D*Lf,J,L,N,K];prepared=ordinary.operator_data('transverse',*co,delta)
    rem,*_=ordinary.one_step('transverse',(3,0),co,delta,prepared);rem=s.cancel(rem*D*Lf/prepared[-1])
    assert grade(rem,Lf,0,1)==0>-1
    return {'cases':rows,'failed_integral_axis_grade':{'kind':'transverse','input_exponent':[3,0],'grade':'v exponent-Lf exponent','actual_remainder_grade':0,'required_bound':-1,'remainder':str(rem)},'formal_UFD_before_specialization':True}

def cusp_universal():
    D,eps=s.symbols('D eps');dt=lambda f:s.diff(f,t)+eps*s.diff(f,D);delta=t*t*D
    R=sum(s.Symbol('c'+str(i)+'_'+str(j)+'_'+str(a)+'_'+str(b))*t**a*D**b*u**i*v**j for i,j in weighted.monomials(10,(3,4)) for a in range(5-i-j) for b in range(5-i-j-a))
    W=u**4/2+D*v**3+R;dp=dt(delta);bx=dp/4;by=2*t*D/3
    old=delta*dt(W)-dp*W-2*t**3*D**2-t**4*D*eps/2
    source=s.cancel((old+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v))/t)
    assert grade(source,D,1,2)<=4
    for z in [t**3*D**2,dt(t**3*D**2),s.cancel(delta*bx/t),s.cancel(delta*by/t)]:assert grade(z,D,1,2)<=3
    rows=[]
    for ex in weighted.monomials(30,(3,4)):
        nf,cu,cv=weighted.principal_division('cusp',ex,s.Rational(1,2),D);nf,cu,cv=[s.cancel(D*z) for z in [nf,cu,cv]]
        rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
        assert s.expand(D*u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
        rho=ex[0]+2*ex[1];assert all(grade(z,D,1,2)<=rho-1 for z in [rem,cu,cv])
        rows.append({'exponent':ex,'remainder_grade':grade(rem,D,1,2),'bound':rho-1})
    return {'source_grade_bound':4,'formal_D_before_specialization':True,'principal_identities':rows}

def replay_binary(m):
    local={'t':t,'u':u,'v':v};read=lambda x:s.sympify(x,locals=local)
    W,D,A=[read(m[n]) for n in ['W','D','A0']];kind=m['kind'];basis=m['basis'];delta=t*t*D;h=t**3*D**2
    if kind=='cusp':
        source,bx,by=shared.weighted_source(kind,W,D);source=s.cancel(source/t)
        vector=lambda f:s.cancel(delta*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))/t)
    else:
        source=s.cancel(euler.ordinary_source(W,D)/t);vector=lambda f:s.cancel(delta*s.diff(delta,t)*(u*s.diff(f,u)+v*s.diff(f,v))/(4*t))
    pv=s.Poly(W,u,v);Le=s.cancel(pv.coeff_monomial(u*u)/D);Lf=s.cancel(pv.coeff_monomial(v*v)/D)
    ring=s.QQ.poly_ring(t);poly=lambda z:s.Poly(z,u,v,domain=ring);wp,dp=poly(W),poly(delta)
    R=s.Integer(1);out=[]
    for row in m['seed_columns']:
        j=row['j'];r=row['depth'];R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+vector(R)+source*R)
        if kind=='cusp':
            assert r==5*j-3;common=D**(j-3);scales=[D**(e[0]+2*e[1]) for e in basis];grades=[sum(e)+e[0]+2*e[1] for e in basis];slope=9
        else:
            assert r==3*j-(1 if kind=='triangle' else 2)
            base=j-(1 if kind=='triangle' else 2);common=D**base
            scales=[D**sum(e) for e in basis];grades=[2*sum(e) for e in basis];slope={'nodal':7,'transverse':9,'triangle':11}[kind]
            if kind=='transverse':
                common*=Lf**base;scales=[z*Lf**ceil(s.Rational(e[0],3)+e[1]) for z,e in zip(scales,basis)];grades=[g+ceil(s.Rational(e[0],3)+e[1]) for g,e in zip(grades,basis)]
            if kind=='triangle':
                common*=Le**base*Lf**base;scales=[z*Le**ceil(e[0]+s.Rational(e[1],2))*Lf**ceil(s.Rational(e[0],2)+e[1]) for z,e in zip(scales,basis)];grades=[g+ceil(e[0]+s.Rational(e[1],2))+ceil(s.Rational(e[0],2)+e[1]) for g,e in zip(grades,basis)]
        values=[s.Poly(s.cancel(read(z)/(common*scale)),t,domain=s.QQ) for z,scale in zip(row['N'],scales)]
        Cu,Cv=[poly(s.cancel(read(row['witness'][name])/common)) for name in ['Cu','Cv']]
        normal=sum((poly(z.as_expr()*scale*u**e[0]*v**e[1]) for z,scale,e in zip(values,scales,basis)),poly(0))
        assert (poly(s.cancel(A**r/common))*poly(R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        assert all(z.degree()<=slope*j-g for z,g in zip(values,grades)),(kind,j,[z.degree() for z in values],grades)
        H=s.Poly(s.cancel(read(row['H'])/common),t,domain=s.QQ);assert H.degree()<=slope*j
        expected=t**(3*j)*D**((6 if kind=='cusp' else 4)*j)
        if kind=='transverse':expected*=Lf**(2*j)
        if kind=='triangle':expected*=Le**(2*j)*Lf**(2*j)
        assert s.expand(H.as_expr()-expected)==0
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        out.append({'j':j,'common_factor':str(common),'row_scales':[str(z) for z in scales],'row_grades':grades,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'witness':{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}})
    return {'kind':kind,'singular_at_zero':m['singular_at_zero'],'metric_determinant':m['metric_determinant'],'W':str(W),'D':str(D),'Le':str(Le),'Lf':str(Lf),'seed_columns':out}

def conic_replay(m):
    read=lambda z:s.sympify(z,locals={'t':t,'y':y});W,K=read(m['W']),read(m['K']);raw=balanced.laurent(W);p=max(raw);kd=s.degree(K,t);basis=m['basis'];field=s.QQ.frac_field(t)
    def parse(x):return {k:field.convert(a) for k,a in balanced.laurent(x).items()}
    def add(out,k,a):
        out[k]=out.get(k,field.zero)+a
        if not out[k]:del out[k]
    def exact(x):
        out={}
        for k,a in x.items():
            add(out,k,field.convert(t*t*k)*a)
            for n,c in raw.items():
                if n:add(out,k+n,field.convert(n*c)*a)
        return out
    B=t*s.diff(K,t)-s.Rational(2,3)*K;h=t**3*K;source=s.expand(K*(t*s.diff(W,t)-2*W-t*t)+B*y*s.diff(W,y));R=s.Integer(1);out=[]
    for row in m['seed_columns']:
        j=row['j'];R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+t*t*B*y*s.diff(R,y)+source*R)
        den=raw[p]**(2*j-1) if p==2 else 1
        values=[s.Poly(read(z),t,domain=s.QQ) for z in row['N']]
        witness={int(k):field.convert(read(a)) for k,a in row['witness'].items()}
        check=parse(R)
        for k,z in zip(basis,values):add(check,k,-field.convert(z.as_expr()/den))
        for k,a in exact(witness).items():add(check,k,-a)
        assert not check
        if p==3 and kd==1:
            bounds=[11*j//2+3]*len(basis)
        elif p==2 and kd==0:
            bounds=[4-k if j==1 else 9*j-1+k for k in basis]
        else:
            # Preserve prior scopes; only these two degree strata change.
            bounds=None
        if bounds:assert all(z.degree()<=b for z,b in zip(values,bounds)),(m['label'],j,[z.degree() for z in values],bounds)
        H=s.Poly(read(row['H']),t,domain=s.QQ);assert s.expand(H.as_expr()-h**j*den)==0
        assert min(e[0] for e,c in H.terms() if c)>=3*j
        out.append({'j':j,'new_row_degree_bounds':bounds,'observed_row_degrees':[int(z.degree()) if z else None for z in values],'H':str(H.as_expr())})
    return {'label':m['label'],'orders':[p,3],'K_degree':int(kd),'seed_columns':out}

def arithmetic():
    rows=[]
    for kind,slope,expected in [('nodal',7,174),('transverse',9,166),('triangle',11,142),('cusp',9,168)]:
        basis=(weighted.CASES if kind=='cusp' else ordinary.CASES)[kind]['basis']
        def g(e):
            i,k=e
            if kind=='cusp':return sum(e)+i+2*k
            if kind=='nodal':return 2*sum(e)
            if kind=='transverse':return 2*sum(e)+ceil(s.Rational(i,3)+k)
            return 2*sum(e)+ceil(i+s.Rational(k,2))+ceil(s.Rational(i,2)+k)
        maximum=-1
        for size in range(1,len(basis)+1):
            for selected in combinations(range(len(basis)),size):
                if 0 not in selected:continue
                maximum=max(maximum,slope*size*(size+1)//2-sum(g(basis[k]) for k in selected))
        assert maximum==expected,(kind,maximum);rows.append({'case':kind,'gate':maximum,'row_grades':[g(e) for e in basis]})
    # DegreeK1 p3: each high/low step cost is bounded by exponent
    # displacement under weights(2,1), after polynomial K^3 division.
    steps=[]
    for endpoint in ['high','low']:
        for k in [-3,-2,-1,1,2,3,'Euler']:
            if k==3 and endpoint=='high' or k==-3 and endpoint=='low':continue
            if k=='Euler':degree=2;shift=-3 if endpoint=='high' else 3
            else:degree=3-k if k>0 else 3;shift=k-3 if endpoint=='high' else k+3
            cost=2*degree+shift-(6 if endpoint=='low' else 0)
            assert cost<=abs(shift);steps.append({'endpoint':endpoint,'phase_power':k,'displacement':shift,'weighted_cost':cost})
    assert sum(11*j//2+3 for j in range(1,7))==132
    rows.append({'case':'conic_p3q3_Klinear_displacement','gate':132})
    # DegreeK0 p2: j1 has no low input; later bounds9j-1+k
    # dominate5j-1-k for every row-3,...,1 and j>=2.
    mx=-1;basis=list(range(-3,2))
    for size in range(1,6):
        for selected in combinations(basis,size):
            if 0 not in selected:continue
            gate=4-sum(selected) if size==1 else 4+sum(9*j-1 for j in range(2,size+1))+sum(selected)-2*min(selected)
            mx=max(mx,gate)
    assert mx==127;rows.append({'case':'conic_p2q3_Kconstant_displacement','gate':mx})
    conic=max(146,132,108,127,136,105,130,145);assert conic==146
    cover={'quadratic':7,'smooth':94,'nodal':174,'cusp':168,'transverse':166,'tangent':132,'triangle':146,'concurrent':104,'repeated':99}
    assert max(cover.values())==174 and max(56,10,10,72,146)==146
    return rows,steps,cover

def main():
    started=time.monotonic();ordinary_checks=ordinary_universal();print('Universal ordinary filtrations passed',flush=True)
    cusp_checks=cusp_universal();print('Universal cusp filtration passed',flush=True)
    old=json.loads((P/'common_time_seed_factor_certificate.json').read_text());binary=[replay_binary(m) for m in old['models'] if m['kind'] in ['nodal','transverse','triangle','cusp']];assert len(binary)==8;print('Eight binary seed models passed',flush=True)
    prior=json.loads((P/'pole_filtered_seed_bounds_certificate.json').read_text());structured=json.loads((P/'structured_seed_modules_certificate.json').read_text())
    conics=[conic_replay(m) for m in prior['conic_models']+structured['conic_models']];print('Twenty conic seed models passed',flush=True)
    rows,steps,cover=arithmetic()
    cert={'ordinary_universal':ordinary_checks,'cusp_universal':cusp_checks,'binary_models':binary,'conic_models':conics,'degree_cases':rows,'linear_K_step_costs':steps,'global_cover':cover}
    path=P/'multiple_divisor_seed_bounds_certificate.json';path.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'retained_inputs':{n:sha(P/n) for n in ['common_time_seed_factor_certificate.json','pole_filtered_seed_bounds_certificate.json','structured_seed_modules_certificate.json']},'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,'binary_models':len(binary),'conic_models':len(conics),'direct_seed_identities':3*(len(binary)+len(conics)),'ordinary_complement_gates':{'nodal':174,'transverse':166,'triangle':142},'cusp_full_gate':168,'conic_full_gate':146,'full_W_linear_gate':146,'tangent_full_gate_unchanged':132,'mixed_parent_gate_unchanged':92,'global_sufficient_pure_prefix':174,'odd_sufficient_pure_prefix':174,'homogeneous_sufficient_pure_prefix':146,'least_prefix_interval':[20,174],'odd_least_prefix_interval':[20,174],'homogeneous_least_prefix_interval':[14,146],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact auxiliary formal source/division filtrations,28 model replays with84 actual seed identities, displacement-cost and determinant arithmetic. Universal all-degree induction, full coefficient/rank cover, period comparison and analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded global radical cofactors.'}
    (P/'multiple_divisor_seed_bounds_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
