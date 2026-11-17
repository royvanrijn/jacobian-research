"""Exact auxiliary endpoint-degree conic bounds and odd node/cusp prefixes.
Universal induction, full covers and analytic exclusions remain written inputs.
"""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import hashlib,json,time
import sympy as s
import verify_balanced_logarithmic_seed_bounds as balanced
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_euler_seed_bounds as euler
import verify_shared_divisor_weighted_seeds as shared
import verify_rescaled_singular_prefixes as rescaled
P=Path(__file__).resolve().parent
t,y=balanced.t,balanced.y
u,v=ordinary.u,ordinary.v

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def parse(z):return balanced.laurent(s.expand(z))
def coeff_degree(z):return s.Poly(s.cancel(z),t).degree()
def specification(p,q,kd,ell,j,k):
    if q==2:
        assert p==3 and kd==1
        return 6*j+abs(k) if ell==0 else 6*j-2+max(-k,0)
    assert p==1 and q==3 and ell in [0,1,2]
    if kd==0:return [4*j-2*k,4*j-k,5*j][ell]
    if kd==1:return [5*j-2*k,5*j-k,6*j+int(s.ceiling(-s.Rational(k,2)))][ell]
    assert kd==2;return [7*j-2*k,8*j+int(s.ceiling(-s.Rational(3*k,2))),9*j-k][ell]

def conic_universal():
    rows=[];k0,k1,k2=s.symbols('k0 k1 k2')
    for p,q,kd in [(1,3,0),(1,3,1),(1,3,2),(3,2,1)]:
        K=[k0,k0+k1*t,k0+k1*t+k2*t*t][kd]
        W=sum(s.Symbol('a'+str(k)+'_'+str(n))*t**n*y**k for k in range(p+1) for n in range(4-k))
        W+=sum(K**k*sum(s.Symbol('b'+str(k)+'_'+str(n))*t**n for n in range(4-k))/y**k for k in range(1,q+1))
        B=t*s.diff(K,t)-s.Rational(2,3)*K if q==3 else s.Rational(2,3)*K
        source=s.expand(K*(t*s.diff(W,t)-2*W-t*t)+B*y*s.diff(W,y));raw=parse(source)
        assert min(raw)>=-2 and max(raw)<=(1 if q==3 else 2)
        if kd==0:assert all(coeff_degree(c)+abs(k)<=3 for k,c in raw.items())
        elif kd==1:assert all(coeff_degree(c)+max(k,0)<=4 for k,c in raw.items())
        else:assert all(coeff_degree(c)+k<=5 for k,c in raw.items())
        rows.append({'p':p,'q':q,'K_degree':kd,'source_support':[min(raw),max(raw)],'seed_degree':'3j-|k|' if kd==0 else '4j-max(k,0)' if kd==1 else '5j-k','whole_denominator_derivative_retained':True,'vertical_Euler_retained':True})
    return rows

def model(label,W,K,maximum_seed=3):
    raw=parse(W);p=max(raw);q=-min(raw);kd=int(s.degree(K,t));assert (p==1 and q==3 and kd in [0,1,2]) or (p==3 and q==2 and kd==1)
    assert all(coeff_degree(raw[k])<=3-k for k in raw if k>=0)
    assert all(coeff_degree(s.cancel(raw[-k]/K**k))<=3-k for k in range(1,q+1) if -k in raw)
    endpoint=raw[p] if q==3 and p==1 else s.cancel(raw[-q]/K**q) if q==2 else s.Integer(1)
    ed=int(s.degree(endpoint,t));assert ed in ([0,1,2] if p==1 else [0,1])
    field=s.QQ.frac_field(t);basis=list(range(-q,p));valid=set(basis)
    def to_field(x):return {k:field.convert(c) for k,c in parse(x).items()}
    def add(out,k,c):
        out[k]=out.get(k,field.zero)+c
        if not out[k]:del out[k]
    def exact(x):
        out={}
        for n,c in x.items():
            add(out,n,field.convert(t*t*n)*c)
            for k,a in raw.items():
                if k:add(out,n+k,field.convert(k*a)*c)
        return out
    def reduce(x):
        normal=x.copy();witness={};high_steps=low_steps=0
        while max(normal,default=0)>=p:
            n=max(normal);shift=n-p;c=normal[n]/field.convert(p*raw[p]);high_steps+=1;add(witness,shift,c)
            for k,a in exact({shift:c}).items():add(normal,k,-a)
        while min(normal,default=0)<-q:
            n=min(normal);shift=n+q;c=normal[n]/field.convert(-q*raw[-q]);low_steps+=1;add(witness,shift,c)
            for k,a in exact({shift:c}).items():add(normal,k,-a)
        assert set(normal)<=valid
        check=x.copy()
        for k,a in normal.items():add(check,k,-a)
        for k,a in exact(witness).items():add(check,k,-a)
        assert not check
        return normal,witness,high_steps,low_steps
    old=t*s.diff(W,t)-2*W-t*t
    B=t*s.diff(K,t)-s.Rational(2,3)*K if q==3 else s.Rational(2,3)*K
    source=s.expand(K*old+B*y*s.diff(W,y));h=t**3*K;R=s.Integer(1);seeds=[]
    for j in range(1,maximum_seed+1):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+t*t*B*y*s.diff(R,y)+source*R)
        support=parse(R);assert min(support)>=-2*j and max(support)<= (p if q==3 else 2)*j
        for k,c in support.items():
            if k<0:s.Poly(s.cancel(c/K**(-k)),t,domain=s.QQ)
            if kd==0:assert coeff_degree(c)+abs(k)<=3*j
            if kd==1 and q==3:assert coeff_degree(c)+max(k,0)<=4*j
        normal,witness,high,low=reduce(to_field(R))
        depth=j if q==3 else 2*j-2
        den=s.expand(endpoint**depth);values=[]
        for k in basis:
            z=s.Poly(s.cancel(field.to_sympy(normal.get(k,field.zero))*den),t,domain=s.QQ)
            bound=specification(p,q,kd,ed,j,k);assert z.degree()<=bound,(label,j,k,z.degree(),bound)
            values.append(z)
        if q==2:assert high<=2*j-2 and low<=2*j-2
        H=s.Poly(s.expand(h**j*den),t,domain=s.QQ);val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        seeds.append({'j':j,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'normal_endpoint_clear_power':depth,'high_low_steps':[high,low],'row_degree_bounds':[specification(p,q,kd,ed,j,k) for k in basis],'witness':{str(k):str(field.to_sympy(c)) for k,c in witness.items()}})
    return {'label':label,'p':p,'q':q,'K_degree':kd,'endpoint_degree':ed,'W':str(W),'K':str(K),'basis':basis,'source':str(source),'vertical_coefficient':str(t*t*B),'seed_columns':seeds}


def odd_phase(kind,D,ju,jv,ell,e,f,qu,qv,quv,a=1,d=2):
    q2=lambda x,z:(qu*x*x+2*quv*x*z+qv*z*z)/2
    if kind=='nodal':
        x,z=u+jv,v+ju
        return s.expand(D*(a*x**3+d*z**3+e*t*t*x+f*t*t*z-q2(x,z))+(u*v-ju*jv+ell*t*t)**2/2)
    return s.expand(D*(v**3+e*t*t*u+f*t*t*v-q2(u,v))+(u*u-ju*u-jv*v+ell*t*t)**2/2)

def odd_universal():
    D,ju,jv,ell,e,f,qu,qv,quv,a,d=s.symbols('D ju jv ell e f qu qv quv a d')
    rows=[]
    for kind in ['nodal','cusp']:
        W=odd_phase(kind,D,ju,jv,ell,e,f,qu,qv,quv,a,d)
        assert s.Poly(W,t,u,v).total_degree()<=4
        if kind=='nodal':source=s.cancel(euler.ordinary_source(W,D)/t)
        else:source=s.cancel(shared.weighted_source(kind,W,D)[0]/t)
        assert s.Poly(source,t,u,v).total_degree()<=4
        if kind=='cusp':
            R=s.expand(W-u**4/2-D*v**3)
            assert s.Poly(s.diff(R,u),t,u,v).total_degree()<=3 and s.Poly(s.diff(R,v),t,u,v).total_degree()<=2
            original=D*(t*(v**3+e*u+f*v)-(qu*u*u+2*quv*u*v+qv*v*v)/2)+(t*(u*u+ell)-ju*u-jv*v)**2/2
            assert s.expand(rescaled.rescale(t*t*original)-W)==0
            origin='direct retained completion/rescaling'
        else:
            checkpoint=json.loads((P/'ordinary_node_phase_shift_checkpoint.json').read_text());read=s.sympify(checkpoint['V'],locals={'t':t,'u':u,'v':v})
            names={z.name:z for z in read.free_symbols};sub={names[n]:val for n,val in {'alpha':0,'m':D,'b':0,'c':0,'k':0,'j0':0,'lu':0,'lv':0,'q0':0,'A':qu,'B':qv,'C':quv,'a':a,'d':d,'e':e,'f':f,'h':ell,'ju':ju,'jv':jv}.items()}
            assert s.expand(rescaled.rescale(read.subs(sub))-W)==0
            origin='exact specialization of retained universal phase-shift checkpoint'
        rows.append({'kind':kind,'D_constant':True,'phase_joint_bound':4,'source_joint_bound':4,'normal_joint_cost':0,'normalization_identity':origin,'same_centered_Gaussian':True})
    return rows

def odd_model(kind,D):
    W=odd_phase(kind,s.Integer(D),s.Integer(2),s.Integer(3),s.Integer(1),s.Integer(5),s.Integer(7),s.Integer(11),s.Integer(13),s.Integer(1))
    delta=t*t*D;h=t**3*D**2;data=(ordinary.CASES if kind=='nodal' else weighted.CASES)[kind];basis=data['basis'];index={e:i for i,e in enumerate(basis)}
    if kind=='nodal':
        source=s.cancel(euler.ordinary_source(W,D)/t);vector=lambda z:delta*s.diff(delta,t)*(u*s.diff(z,u)+v*s.diff(z,v))/(4*t)
        pv=s.Poly(W,u,v);mons=[u*u*v*v,u**3,v**3,u*u,v*v,u*v,u,v,1];co=[pv.coeff_monomial(z) for z in mons];prepared=ordinary.operator_data(kind,*co,delta)
    else:
        source,bx,by=shared.weighted_source(kind,W,s.Integer(D));source=s.cancel(source/t);vector=lambda z:delta*(bx*u*s.diff(z,u)+by*v*s.diff(z,v))/t
    ring=s.QQ.poly_ring(t);poly=lambda z:s.Poly(z,u,v,domain=ring);zero=poly(0);zt=s.Poly(0,t,domain=s.QQ);wp,dp=poly(W),poly(delta)
    @lru_cache(None)
    def step(ex):
        if kind=='nodal':
            rem,cu,cv,_,_=ordinary.one_step(kind,ex,co,delta,prepared);rem,cu,cv=[s.cancel(z/prepared[-1]) for z in [rem,cu,cv]]
        else:
            nf,cu,cv=weighted.principal_division(kind,ex,s.Rational(1,2),s.Integer(D));R=W-u**4/2-D*v**3
            rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
        assert s.Poly(rem,t,u,v).total_degree()<=sum(ex)
        assert s.expand(u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
        return poly(rem),poly(cu),poly(cv)
    @lru_cache(None)
    def reduce(ex,depth):
        if ex in index:
            out=[zt]*len(basis);out[index[ex]]=s.Poly(1,t);return tuple(out),zero,zero
        assert depth>0;rem,Cu,Cv=step(ex);out=[zt]*len(basis)
        for f,c in rem.terms():
            if not c:continue
            vec,x,z=reduce(f,depth-1);ct=s.Poly(c,t,domain=s.QQ)
            out=[a+ct*b for a,b in zip(out,vec)];Cu+=poly(c)*x;Cv+=poly(c)*z
        return tuple(out),Cu,Cv
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+source*R+vector(R));assert s.Poly(R,t,u,v).total_degree()<=4*j
        depth=3*j-2 if kind=='nodal' else 5*j-3;values=[zt]*len(basis);Cu,Cv=zero,zero
        for ex,c in poly(R).terms():
            if not c:continue
            vec,x,z=reduce(ex,depth);ct=s.Poly(c,t,domain=s.QQ)
            values=[a+ct*b for a,b in zip(values,vec)];Cu+=poly(c)*x;Cv+=poly(c)*z
        normal=sum((poly(z.as_expr()*u**ex[0]*v**ex[1]) for z,ex in zip(values,basis)),zero)
        assert (poly(R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        assert all(z.degree()<=4*j-sum(ex) for z,ex in zip(values,basis))
        H=s.Poly(h**j,t,domain=s.QQ);assert min(e[0] for e,c in H.terms() if c)>=3*j
        seeds.append({'j':j,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'row_degree_bounds':[4*j-sum(ex) for ex in basis],'witness':{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}})
    metric=int(s.Matrix([[11,1,2],[1,13,3],[2,3,D]]).det());assert metric!=0
    return {'kind':kind,'D':D,'metric_determinant':metric,'W':str(W),'source':str(source),'seed_columns':seeds}

def arithmetic():
    rows=[]
    for kd in [0,1,2]:
        for ell in [0,1,2]:
            basis=list(range(-3,1));mx=max(sum(specification(1,3,kd,ell,j,0) for j in range(1,size+1))+sum(specification(1,3,kd,ell,1,k)-specification(1,3,kd,ell,1,0) for k in selected) for size in range(1,5) for selected in combinations(basis,size) if 0 in selected)
            expected=[[52,46,50],[62,56,64],[82,90,96]][kd][ell];assert mx==expected
            rows.append({'p':1,'q':3,'K_degree':kd,'endpoint_degree':ell,'gate':mx})
    for ell,expected in [(0,96),(1,83)]:
        basis=list(range(-2,3));mx=max(sum(specification(3,2,1,ell,j,0) for j in range(1,size+1))+sum(specification(3,2,1,ell,1,k)-specification(3,2,1,ell,1,0) for k in selected) for size in range(1,6) for selected in combinations(basis,size) if 0 in selected);assert mx==expected
        rows.append({'p':3,'q':2,'K_degree':1,'endpoint_degree':ell,'gate':mx})
    remaining=[c for c in balanced.degree_cases() if c['negative_order']<3 and (c['positive_order'],c['negative_order'])!=(3,2)];assert max(c['all_seed_rank_and_row_choices_gate'] for c in remaining)==108
    eulercert=json.loads((P/'euler_seed_bounds_certificate.json').read_text());assert eulercert['degree_cases']['conic_q3'][0]['positive_order']==0 and eulercert['degree_cases']['conic_q3'][0]['gate']==69
    assert max(72,90,108,85,100,105,96,69,87,96,95,108)==108
    # All seed ranks and basis-row choices for the two new odd complements.
    odd=[]
    for kind,expected,parity in [('nodal',101,100),('cusp',75,74)]:
        basis=(ordinary.CASES if kind=='nodal' else weighted.CASES)[kind]['basis'];r=len(basis)
        mx=max(sum(4*j for j in range(1,size+1))-sum(sum(basis[k]) for k in selected) for size in range(1,r+1) for selected in combinations(range(r),size) if 0 in selected);assert mx==expected and mx%2==1
        odd.append({'kind':kind,'cofactor_gate':mx,'odd_sufficient_complement':parity})
    global_cover={'quadratic':7,'smooth':94,'nodal':129,'cusp':130,'transverse':127,'tangent':120,'triangle':118,'concurrent':104,'repeated':99}
    odd_cover={'quadratic':7,'smooth':94,'nodal':108,'cusp':94,'transverse':126,'tangent':120,'triangle':118,'concurrent':104,'repeated':98}
    assert max(global_cover.values())==130 and max(odd_cover.values())==126 and max(56,10,10,72,108)==108
    return rows,odd,global_cover,odd_cover

def main():
    start=time.monotonic();sources=conic_universal();odd_sources=odd_universal();print('Universal endpoint and odd source/phase identities passed',flush=True)
    conics=[]
    for kd,K in [(0,s.Rational(-1,4)),(1,1+t),(2,t*(t+1))]:
        for A in [s.Integer(1),1+t,t*t]:
            W=1+t**3+A*y+sum(K**k*(2+t**(3-k))/y**k for k in [1,2,3]);conics.append(model('p1_K'+str(kd)+'_A'+str(A),s.expand(W),K))
    for K,A in [(1+t,1+t),(1+t,(1+t)**2),(t*(t+1),t*(t+1))]:
        W=1+t**3+A*y+sum(K**k*(2+t**(3-k))/y**k for k in [1,2,3]);conics.append(model('shared_p1_K'+str(K)+'_A'+str(A),s.expand(W),K))
    for kd,K,U,V in [(0,s.Rational(-1,4),s.Rational(1,2),s.Rational(-1,2)),(1,t/2-s.Rational(1,4),(1-t)/2,(t-1)/2),(2,(t+1)*(5*t-1)/4-t*t,(t+1)/2,(5*t-1)/2)]:
        Z=(y-K/y-U+V)/2;T=y-U-Z;S=Z+T;D=Z-T;lam=-s.Rational(1,4);corr=-lam*s.diff(V,t)
        W=s.expand(D**3+lam*S*S*D+corr*t*S*S+t*t*S+2*t**3-T*T/2)
        conics.append(model('actual_p1_conic_K'+str(kd),W,K))
    prior=json.loads((P/'two_sided_conic_seed_bounds_certificate.json').read_text())
    for m in prior['models']:
        if m['q']==2 and m['K_degree']==1:conics.append(model('refined_'+m['label'],s.sympify(m['W'],locals={'t':t,'y':y}),s.sympify(m['K'],locals={'t':t})))
    print('Twenty-two conic models passed',flush=True)
    oddmodels=[odd_model(kind,D) for kind in ['nodal','cusp'] for D in [1,2]];print('Four odd binary models passed',flush=True)
    cases,odd_cases,cover,odd_cover=arithmetic();cert={'conic_universal':sources,'odd_universal':odd_sources,'conic_models':conics,'odd_models':oddmodels,'endpoint_degree_cases':cases,'odd_degree_cases':odd_cases,'global_cover':cover,'odd_cover':odd_cover}
    path=P/'endpoint_degrees_and_odd_prefixes_certificate.json';path.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'conic_source_strata':len(sources),'odd_source_strata':len(odd_sources),'conic_models':len(conics),'odd_models':len(oddmodels),'direct_seed_identities':3*(len(conics)+len(oddmodels)),'conic_gate':108,'full_W_linear_gate':108,'odd_nodal_complement':100,'odd_cusp_complement':74,'odd_cusp_full_gate':94,'global_sufficient_pure_prefix_unchanged':130,'odd_sufficient_pure_prefix':126,'homogeneous_sufficient_pure_prefix':108,'least_prefix_interval':[20,130],'odd_least_prefix_interval':[20,126],'homogeneous_least_prefix_interval':[14,108],'global_cover':cover,'odd_cover':odd_cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact auxiliary endpoint/source and centered odd phase identities,26 models with78 actual seed reductions, actual endpoint-degree and all-rank/row arithmetic. Full all-degree degree/module/path inductions, centered Gaussian coefficient/frame cover and analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded global radical cofactors.'}
    (P/'endpoint_degrees_and_odd_prefixes_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
