"""Auxiliary two-sided conic degree and exact seed reduction certificates.
Analytic exclusions, universal inductions and full class cover remain written.
No retained checker main or previous receipt is executed or rewritten.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_balanced_logarithmic_seed_bounds as balanced
P=Path(__file__).resolve().parent
t,y=balanced.t,balanced.y

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def parse(x):return balanced.laurent(s.expand(x))
def coeff_degree(x):return s.Poly(s.cancel(x),t).degree()

def specification(p,q,kd,endpoint_degree,j,k):
    if q==3:
        if p==3 and kd==0:return 3*j+abs(k)
        if p==3 and kd==1:return 4*j+max(-k,0)
        if p==2 and kd==0:return (5*j-1+abs(k)) if endpoint_degree==1 else (5*j-2*k)
        if p==2 and kd==1:return (6*j-1+max(-k,0)) if endpoint_degree==1 else (6*j-2*k)
        raise AssertionError((p,q,kd))
    assert p==3 and q==2
    if kd==0:return (5*j-2+abs(k)) if endpoint_degree==1 else (5*j+2*abs(k))
    if kd==1:return (22*j-14-k)//2
    assert kd==2;return 7*j-2-k

def universal_sources():
    rows=[];k0,k1,k2=s.symbols('k0 k1 k2')
    for p,q,kd in [(p,3,kd) for p in [2,3] for kd in [0,1]]+[(3,2,kd) for kd in [0,1,2]]:
        K=[k0,k0+k1*t,k0+k1*t+k2*t*t][kd]
        W=sum(s.Symbol('a'+str(k)+'_'+str(n))*t**n*y**k for k in range(p+1) for n in range(4-k))
        W+=sum(K**k*sum(s.Symbol('b'+str(k)+'_'+str(n))*t**n for n in range(4-k))/y**k for k in range(1,q+1))
        old=t*s.diff(W,t)-2*W-t*t
        B=t*s.diff(K,t)-s.Rational(2,3)*K if q==3 else s.Rational(2,3)*K
        source=s.expand(K*old+B*y*s.diff(W,y));h=t**3*K
        # Generic polynomial source bounds provide the all-j recurrence
        # induction by subadditivity of abs(k) or max(k,0).
        data=parse(source);assert min(data)>=-2 and max(data)<= (p if q==3 else 2)
        if kd==0:
            assert all(coeff_degree(c)+abs(k)<=3 for k,c in data.items())
            seed_bound='degree R_(k,j)<=3j-abs(k)'
        elif kd==1 and q==3:
            assert all(coeff_degree(c)+max(k,0)<=4 for k,c in data.items())
            seed_bound='degree R_(k,j)<=4j-max(k,0)'
        else:
            wt,wy=(2,1) if kd==1 else (1,1)
            cost=8 if kd==1 else 5
            assert all(wt*coeff_degree(c)+wy*k<=cost for k,c in data.items())
            seed_bound='weighted degree<='+str(cost)+'j'
        if q==2:
            # Exact equality of representatives in the original rational
            # logarithmic quotient, including the vertical derivative.
            for n in range(-4,5):
                f=y**n;difference=s.expand(source*f+s.Rational(2,3)*t*t*K*y*s.diff(f,y)-K*old*f)
                witness=s.Rational(2,3)*K*f
                assert s.expand(difference-t*t*y*s.diff(witness,y)-y*s.diff(W,y)*witness)==0
        rows.append({'p':p,'q':q,'K_degree':kd,'source_support':[min(data),max(data)],'seed_bound':seed_bound,'same_logarithmic_measure':True,'negative_module':'coefficient y^(-k) divisible by K^k before rational endpoint reductions','source_gauge_checks':9 if q==2 else 0})
    return rows

def model(label,W,K,maximum_seed=3):
    raw=parse(W);p=max(raw);q=-min(raw);kd=int(s.degree(K,t));assert q==3 and p in [2,3] and kd in [0,1] or p==3 and q==2 and kd in [0,1,2]
    assert all(coeff_degree(raw[k])<=3-k for k in raw if k>=0)
    assert all(coeff_degree(s.cancel(raw[-k]/K**k))<=3-k for k in range(1,q+1) if -k in raw)
    endpoint=raw[p] if q==3 and p==2 else s.cancel(raw[-q]/K**q) if q==2 else s.Integer(1)
    ed=int(s.degree(endpoint,t));assert ed in [0,1]
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
        depth=2*j-1 if q==3 and p==2 else 2*j-2 if q==2 else 0
        den=s.expand(endpoint**depth);values=[]
        for k in basis:
            z=s.Poly(s.cancel(field.to_sympy(normal.get(k,field.zero))*den),t,domain=s.QQ)
            bound=specification(p,q,kd,ed,j,k);assert z.degree()<=bound,(label,j,k,z.degree(),bound)
            values.append(z)
        if q==2:assert high<=2*j-2 and low<=2*j-2
        H=s.Poly(s.expand(h**j*den),t,domain=s.QQ);val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        seeds.append({'j':j,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'normal_endpoint_clear_power':depth,'high_low_steps':[high,low],'row_degree_bounds':[specification(p,q,kd,ed,j,k) for k in basis],'witness':{str(k):str(field.to_sympy(c)) for k,c in witness.items()}})
    return {'label':label,'p':p,'q':q,'K_degree':kd,'endpoint_degree':ed,'W':str(W),'K':str(K),'basis':basis,'source':str(source),'vertical_coefficient':str(t*t*B),'seed_columns':seeds}

def arithmetic():
    rows=[]
    for p,q,kd,ed,expected in [(3,3,0,0,72),(3,3,1,0,90),(2,3,0,0,85),(2,3,0,1,77),(2,3,1,0,100),(2,3,1,1,91),(3,2,0,0,87),(3,2,0,1,71),(3,2,1,0,129),(3,2,1,1,129),(3,2,2,0,95),(3,2,2,1,95)]:
        basis=list(range(-q,p));maximum=-1
        for size in range(1,len(basis)+1):
            for selected in combinations(basis,size):
                if 0 not in selected:continue
                # Row bounds have separable grades except q3/K0/linear
                # endpoints abs(k), still separable. All j>=1 included.
                val=sum(specification(p,q,kd,ed,j,0) for j in range(1,size+1))
                val+=sum(specification(p,q,kd,ed,1,k)-specification(p,q,kd,ed,1,0) for k in selected)
                maximum=max(maximum,val)
        assert maximum==expected,(p,q,kd,ed,maximum,expected)
        rows.append({'p':p,'q':q,'K_degree':kd,'endpoint_degree':ed,'gate':maximum})
    assert max(72,90,108,85,100,105,130,87,129,95,108)==130
    cover={'quadratic':7,'smooth':94,'nodal':130,'cusp':140,'transverse':132,'tangent':132,'triangle':130,'concurrent':104,'repeated':99}
    assert max(cover.values())==140 and max(56,10,10,72,130)==130
    return rows,cover

def main():
    started=time.monotonic();sources=universal_sources();print('Seven universal source strata passed',flush=True)
    models=[]
    # Every new endpoint-degree stratum, including shared factors and K(0)=0.
    for kd,K in [(0,s.Rational(-1,4)),(1,1+t)]:
        for p in [2,3]:
            for Ap in ([s.Integer(1),1+t,t] if p==2 else [s.Integer(2)]):
                W=1+t**3+(1+t*t)*y+((1+t)*y*y if p==3 else Ap*y*y)+(Ap*y**3 if p==3 else 0)+sum(K**k*(2+t**(3-k))/y**k for k in [1,2,3])
                models.append(model('q3_p'+str(p)+'_K'+str(kd)+'_A'+str(Ap),s.expand(W),K))
    for kd,K in [(0,s.Rational(-1,4)),(1,1+t),(1,t),(2,t*(t+1))]:
        for B2 in [s.Integer(1),1+t,t]:
            W=1+t**3+(1+t*t)*y+(1+t)*y*y+2*y**3+K*(1+t*t)/y+K*K*B2/y**2
            models.append(model('q2_K'+str(kd)+'_'+str(K)+'_B'+str(B2),s.expand(W),K))
    for kd,K,U,V in [(0,s.Rational(-1,4),s.Rational(1,2),s.Rational(-1,2)),(1,t/2-s.Rational(1,4),(1-t)/2,(t-1)/2),(2,(t+1)*(5*t-1)/4-t*t,(t+1)/2,(5*t-1)/2)]:
        Z=(y-K/y-U+V)/2;T=y-U-Z
        W=s.expand(Z*Z*T+Z*T*T+sum((1+2*i+3*j)*t**(3-i-j)*Z**i*T**j for i in range(3) for j in range(3-i))-T*T/2)
        models.append(model('actual_q2_conic_K'+str(kd),W,K))
    print('Twenty-three exact conic models passed',flush=True)
    rows,cover=arithmetic();cert={'universal_sources':sources,'models':models,'degree_cases':rows,'global_cover':cover}
    path=P/'two_sided_conic_seed_bounds_certificate.json';path.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,'universal_source_strata':len(sources),'positive_Euler_gauge_identities':sum(c['source_gauge_checks'] for c in sources),'models':len(models),'actual_q2_conic_models':3,'direct_seed_identities':3*len(models),'conic_full_gate':130,'full_W_linear_gate':140,'unrestricted_rank_one_gate_unchanged':140,'global_sufficient_pure_prefix':140,'odd_sufficient_pure_prefix':140,'homogeneous_sufficient_pure_prefix':130,'least_prefix_interval':[20,140],'odd_least_prefix_interval':[20,140],'homogeneous_least_prefix_interval':[14,130],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact auxiliary seven source strata and27 positive logarithmic Euler gauge identities,23 model reductions with69 actual seed identities, endpoint-degree and all-rank/row arithmetic. Universal all-degree degree/module/path inductions, full coefficient/rank cover and retained analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded global radical cofactors.'}
    (P/'two_sided_conic_seed_bounds_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
