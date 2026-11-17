"""Exact auxiliary checks for joint-degree seeds; retained analytic proofs are external inputs.

Never executes or rewrites the retained checkers' mains or receipts.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_squarefree_singular_finite_prefixes as weighted

P=Path(__file__).resolve().parent
t,u,v,y,T,Z=s.symbols('t u v y T Z')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def degree(f,*gens):
    return -1000 if f==0 else int(s.Poly(f,*gens).total_degree())

def weighted_identities():
    out=[]
    d0,d1=s.symbols('d0 d1');D=d0+d1*t;delta=t*t*D
    for kind,k in [('cusp',1),('tangent',2)]:
        data=weighted.CASES[kind];a=s.Rational(1,2);b=D;A=a*b**k
        top=a*u**4+b*(v**3 if kind=='cusp' else u*v*v)
        R=sum(s.Symbol('r_'+str(i)+'_'+str(j)+'_'+str(h))*t**h*u**i*v**j
              for i,j in weighted.monomials(data['lower_weight'],data['weights'])
              for h in range(5-i-j))
        phase=top+R;rows=[]
        for e in weighted.monomials(2*data['top_weight'],data['weights']):
            nf,cu,cv=weighted.principal_division(kind,e,a,b)
            nf,cu,cv=[s.cancel(A*f) for f in [nf,cu,cv]]
            remainder=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
            defect=A*u**e[0]*v**e[1]-remainder-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(phase,u)*cu-s.diff(phase,v)*cv
            assert s.expand(defect)==0
            assert degree(cu,t,u,v)<=sum(e)+k-3 and degree(cv,t,u,v)<=sum(e)+k-3
            assert degree(remainder,t,u,v)<=sum(e)+k
            assert all(f in data['basis'] or weighted.weight(f,data['weights'])<=weighted.weight(e,data['weights'])-data['gap'] for f,c in s.Poly(remainder,u,v).terms() if c)
            rows.append({'exponent':e,'remainder_total_degree':degree(remainder,t,u,v),'bound':sum(e)+k})
        out.append({'kind':kind,'principal_divisor_degree':k,'identities':rows})
    return out

def replay_weighted_models():
    old=json.loads((P/'rescaled_singular_prefixes_certificate.json').read_text());out=[]
    for m in old['models']:
        kind=m['kind']
        if kind not in ['cusp','tangent']:continue
        data=weighted.CASES[kind];k=1 if kind=='cusp' else 2
        local={'t':t,'u':u,'v':v};read=lambda f:s.sympify(f,locals=local)
        phase,D,A=map(read,[m['V'],m['D'],m['A0']]);delta=t*t*D
        source=read(m['source']);assert degree(phase,t,u,v)<=4 and degree(source,t,u,v)<=6
        assert s.expand(source-(s.diff(phase,t)*delta-phase*s.diff(delta,t)-2*t**3*D**2-t**4*D*s.diff(D,t)/2))==0
        ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,v,domain=ring)
        vp,dp,ap=map(poly,[phase,delta,A]);vu,vv=vp.diff(u),vp.diff(v)
        def verify(f,row,depth):
            values=[read(z) for z in row['N']];cu,cv=[poly(read(row['witness'][z])) for z in ['Cu','Cv']]
            normal=sum((poly(z*u**e[0]*v**e[1]) for z,e in zip(values,data['basis'])),poly(0))
            assert (ap**depth*poly(f)-normal-dp*(cu.diff(u)+cv.diff(v))-vu*cu-vv*cv).is_zero
            return values
        for e,row in zip(data['basis'],m['connection_columns']):
            values=verify(source*u**e[0]*v**e[1],row,m['connection_depth'])
            for f,value in zip(data['basis'],values):assert degree(value,t)<=6+sum(e)+k*m['connection_depth']-sum(f)
        R=s.Integer(1);seeds=[]
        for row in m['direct_seed_columns']:
            j=row['j'];R=s.expand(delta**2*s.diff(R,t)-2*(j-1)*delta*s.diff(delta,t)*R+source*R)
            depth=(6*j-3 if kind=='cusp' else 8*j-4);assert row['depth']==depth
            values=verify(R,row,depth);bnd=6*j+k*depth
            for e,value in zip(data['basis'],values):assert degree(value,t)<=bnd-sum(e)
            assert degree(read(row['H']),t)<=bnd
            seeds.append({'j':j,'degree_bound':bnd,'row_degrees':[degree(x,t) for x in values]})
        out.append({'kind':kind,'singular_at_zero':m['singular_at_zero'],'metric_determinant':m['metric_determinant'],'connection_columns_replayed':len(m['connection_columns']),'direct_seed_columns_replayed':seeds})
    return out

def universal_rank_one():
    a,b,c,f=s.symbols('a b c f');q=1-b*t;A=a+c*T+f*T*T
    G=sum(s.Symbol('g'+str(i)+str(j))*Z**i*T**j for i in range(4) for j in range(4-i))
    old=s.expand(sum(s.Symbol('g'+str(i)+str(j))*t**(i+1)*A**i*T**j/q**i for i in range(4) for j in range(4-i))-T*T/2)
    ahat=a*t*t+c*t*y+f*y*y
    W=s.expand(sum(s.Symbol('g'+str(i)+str(j))*t**(3-i-j)*q**(3-i)*ahat**i*y**j for i in range(4) for j in range(4-i))-q**3*y*y/2)
    delta=t*t*q**3
    assert s.cancel(W/delta-old.subs(T,y/t))==0
    source=s.expand(s.diff(W,t)*delta-W*s.diff(delta,t)-delta**2/t+b*delta**2/q)
    assert degree(W,t,y)<=6 and degree(delta,t)<=5 and degree(source,t,y)<=10
    pulled=(s.diff(old,t)+b/q).subs(T,y/t)
    exact=(delta+y*s.diff(W,y))/(t*delta)
    assert s.cancel(source/delta**2-pulled+exact)==0
    assert all(degree(co,t)<=6-e[0] for e,co in s.Poly(W,y).terms() if co)
    return {'phase_total_degree':6,'delta_degree':5,'source_total_degree':10,'coordinate':'y=t*T','prefactor':'1/(t*q)','source':'W_t*delta-W*delta_prime-delta^2/t+b*delta^2/q','gauge_identity':'new_source-pulled_old_source=-Dop(y)/(t*delta)','coefficient_degree':'deg_t W_k<=6-k'}

def one_variable_model(label,phase,delta,source,seed_denominator,total_bound,lead_bound,kind):
    lead_bound=int(lead_bound)
    ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,y,domain=ring)
    wp,dp,sp=map(poly,[phase,delta,source]);derivative=wp.diff(y);d=wp.degree();assert 2<=d<=6
    A=s.Poly(derivative.LC(),t,domain=s.QQ);ap=poly(A.as_expr());basis=list(range(d-1));zt=s.Poly(0,t,domain=s.QQ);zero=poly(0)
    assert A.degree()<=lead_bound
    @lru_cache(None)
    def step(n):
        power=n-d+1;cu=poly(y**power)
        remainder=ap*poly(y**n)-dp*cu.diff(y)-derivative*cu
        assert remainder.degree()<n
        if kind=='rank_one':assert degree(remainder.as_expr(),t,y)<=n+lead_bound
        else:assert all(degree(co,t)<=lead_bound for e,co in remainder.terms() if co)
        return remainder,cu
    @lru_cache(None)
    def reduce(n,depth):
        if n<d-1:
            values=[zt]*len(basis);values[n]=A**depth
            return tuple(values),zero
        assert depth>0,(label,n,depth)
        remainder,cu=step(n);values=[zt]*len(basis);witness=ap**(depth-1)*cu
        for (j,),co in remainder.terms():
            if not co:continue
            vector,w=reduce(j,depth-1);ct=s.Poly(co,t,domain=s.QQ)
            values=[x+ct*z for x,z in zip(values,vector)];witness+=poly(co)*w
        return tuple(values),witness
    def reduce_input(f,depth):
        values=[zt]*len(basis);witness=zero
        for (n,),co in poly(f).terms():
            if not co:continue
            vector,w=reduce(n,depth);ct=s.Poly(co,t,domain=s.QQ)
            values=[x+ct*z for x,z in zip(values,vector)];witness+=poly(co)*w
        normal=sum((poly(z.as_expr()*y**n) for n,z in enumerate(values)),zero)
        assert (ap**depth*poly(f)-normal-dp*witness.diff(y)-derivative*witness).is_zero
        return values,witness
    connections=[]
    for n in basis:
        values,w=reduce_input(source*y**n,d)
        for k,z in enumerate(values):assert z.degree()<=total_bound+(n if kind=='rank_one' else 0)+lead_bound*d-(k if kind=='rank_one' else 0)
        connections.append({'N':[str(z.as_expr()) for z in values],'witness':str(w.as_expr())})
    h=seed_denominator;R=s.Integer(1);seeds=[]
    for j in range(1,4):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+source*R)
        if kind=='rank_one':assert degree(R,t,y)<=total_bound*j
        else:assert all(degree(co,t)<=total_bound*j for e,co in s.Poly(R,y).terms() if co)
        assert s.Poly(R,y).degree()<=d*j
        depth=d*(j-1)+2;values,w=reduce_input(R,depth);bnd=total_bound*j+lead_bound*depth
        for n,z in enumerate(values):assert z.degree()<=bnd-(n if kind=='rank_one' else 0)
        H=s.Poly(s.expand(h**j*A.as_expr()**depth),t,domain=s.QQ);assert H.degree()<=bnd
        val=min(e[0] for e,co in H.terms() if co);assert val>=(4*j if kind=='rank_one' else 3*j)
        seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'witness':str(w.as_expr())})
    return {'label':label,'kind':kind,'d':d,'W':str(phase),'delta':str(delta),'source':str(source),'seed_denominator':str(h),'A':str(A.as_expr()),'connection_columns':connections,'seed_columns':seeds}

def polynomial_models():
    out=[]
    rank_models=[('dense_b1',1,1,2,1,{(i,j):s.Integer(2+3*i+5*j) for i in range(4) for j in range(4-i)}),
                 ('dense_b0',2,0,-1,1,{(i,j):s.Integer(1+i+j) for i in range(4) for j in range(4-i)}),
                 ('degree5',0,1,1,1,{(2,1):s.Integer(1)}),('degree4',1,1,1,1,{(2,0):s.Integer(1)}),
                 ('degree3',1,1,1,1,{(0,3):s.Integer(1)}),('flattened',1,1,0,1,{(1,0):-s.Rational(1,2),(0,2):-s.Rational(1,2)}),
                 ('safe_null',2,0,-1,1,{})]
    for label,a,b,c,f,G in rank_models:
        q=1-b*t;ah=a*t*t+c*t*y+f*y*y
        W=s.expand(sum(co*t**(3-i-j)*q**(3-i)*ah**i*y**j for (i,j),co in G.items())-q**3*y*y/2)
        delta=t*t*q**3;source=s.expand(s.diff(W,t)*delta-W*s.diff(delta,t)-delta**2/t+b*delta**2/q)
        d=s.degree(W,y);out.append(one_variable_model(label,W,delta,source,delta**2,10,6-d,'rank_one'))
    for label,G in [('mixed6',T**3),('mixed5',Z*T*T),('mixed4',T*T),('mixed_constant',0)]:
        F=y/t-2-y*y;Phi=s.expand(t*s.sympify(G).subs({Z:y,T:F},simultaneous=True)-F*F/2)
        W=s.expand(t*t*Phi);source=s.expand(t*s.diff(W,t)-2*W-t*t)
        assert all(degree(co,t)<=3 for e,co in s.Poly(W,y).terms() if co)
        assert all(degree(co,t)<=3 for e,co in s.Poly(source,y).terms() if co)
        out.append(one_variable_model(label,W,t*t,source,t**3,3,3,'mixed_parent'))
    return out

def degree_cases():
    out=[]
    for kind,rank,grades,slope,offset in [('cusp',6,[0,1,2,1,2,3],12,3),('tangent',5,[0,1,2,1,2],22,8)]:
        vals=[]
        for r in range(1,rank+1):
            for rows in combinations(range(1,rank),r-1):
                selected=(0,)+rows;bound=sum(slope*j-offset for j in range(1,r+1))-sum(grades[k] for k in selected)
                vals.append({'rank':r,'rows':selected,'bound':bound})
        gate=max(x['bound'] for x in vals);assert gate==({'cusp':225,'tangent':284}[kind])
        out.append({'kind':kind,'gate':gate,'cases':vals})
    rank_cases=[];mixed_cases=[]
    for d in range(2,7):
        k=6-d;r=d-1;bj=lambda j:10*j+k*(d*(j-1)+2)
        values=[sum(bj(j) for j in range(1,h+1))-sum(rows) for h in range(1,r+1) for rows0 in combinations(range(1,r),h-1) for rows in [(0,)+rows0]]
        gate=max(values);rank_cases.append({'d':d,'leading_degree_bound':k,'gate':gate,'row_subsets_checked':len(values)})
    assert [x['gate'] for x in rank_cases]==[18,50,93,132,140]
    for d in [4,5,6]:
        gate=sum(3*j+3*(d*(j-1)+2) for j in range(1,d));mixed_cases.append({'d':d,'gate':gate})
    assert [x['gate'] for x in mixed_cases]==[72,144,255]
    return {'weighted':out,'rank_one':rank_cases,'mixed_parent':mixed_cases}

def main():
    start=time.monotonic();universal=weighted_identities();print('Joint-degree weighted identities passed',flush=True)
    models=replay_weighted_models();print('Retained four weighted connection/seed models replayed',flush=True)
    rank_normal=universal_rank_one();print('Universal rank-one rescaling/source identity passed',flush=True)
    polynomial=polynomial_models();print('Polynomial model connections and first three seeds passed',flush=True)
    cases=degree_cases();cover={'quadratic':7,'smooth':94,'nodal':353,'cusp':255,'transverse':336,'tangent':284,'triangle':336,'concurrent':104,'repeated':99}
    assert max(cover.values())==353
    retained=['rescaled_singular_prefixes_certificate.json','rescaled_singular_prefixes_verification.json','rank_one_wt2_verification.json','mixed_period_infinity_verification.json','balanced_logarithmic_seed_bounds_verification.json']
    cert={'universal_weighted':universal,'weighted_models_replayed':models,'rank_one_normalization':rank_normal,'polynomial_models':polynomial,'degree_cases':cases,'global_cover':cover,'retained_receipts':{n:sha(P/n) for n in retained}}
    f=P/'joint_degree_seed_bounds_certificate.json';f.write_text(json.dumps(cert,indent=2)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
            'weighted_symbolic_identities':sum(len(x['identities']) for x in universal),'weighted_models_replayed':len(models),'weighted_connection_columns_replayed':sum(x['connection_columns_replayed'] for x in models),'weighted_direct_seed_columns_replayed':sum(len(x['direct_seed_columns_replayed']) for x in models),
            'polynomial_models':len(polynomial),'polynomial_connection_columns':sum(len(x['connection_columns']) for x in polynomial),'polynomial_direct_seed_columns':sum(len(x['seed_columns']) for x in polynomial),
            'rank_one_gate':140,'mixed_parent_gate':255,'cusp_complement_gate':225,'cusp_full_gate':255,'tangent_complement_gate':284,'tangent_full_gate':284,'full_W_linear_gate':336,
            'global_cover':cover,'global_sufficient_pure_prefix':353,'least_prefix_interval':[20,353],'odd_least_prefix_interval':[20,352],'all_homogeneous_sufficient_pure_prefix':336,'all_homogeneous_least_prefix_interval':[14,336],'binary_search_rounds_max':9,
            'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,
            'boundary':'Exact joint-degree divergence identities, retained weighted model replay, universal rank-one coordinate/source identities and eleven polynomial connection/seed models. Analytic exclusions, all-coefficient reduction induction, class cover and nonzero zeroth coefficients remain written and independently unreviewed. No expanded global radical cofactors or leastness computation.'}
    (P/'joint_degree_seed_bounds_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
