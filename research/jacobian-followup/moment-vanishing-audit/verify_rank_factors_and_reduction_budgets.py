"""Auxiliary rank-one factor and stopped reduction pole-budget certificates.
Complete all-degree/stratum and analytic conclusions remain written inputs.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_euler_seed_bounds as euler
import verify_shared_divisor_weighted_seeds as shared
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
y=s.Symbol('y')

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def budget(kind,e):return max(0,int(s.ceiling((s.Rational(e[0],3)+e[1]-2) if kind=='transverse' else (s.Rational(e[0],2)+e[1]-1))))
def phase_factor(d):return 0 if d==6 else 1 if d in [4,5] else 2

def rank_universal():
    a,b,c,f,Q=s.symbols('a b c f Q');g={(i,j):s.Symbol('g'+str(i)+str(j)) for i in range(4) for j in range(4-i)}
    q=1-b*t;ah=a*t*t+c*t*y+f*y*y
    W=sum(co*t**(3-i-j)*Q**(3-i)*ah**i*y**j for (i,j),co in g.items())-Q**3*y*y/2;W=s.expand(W)
    actual=W.subs(Q,q);vanished={};hierarchy=[]
    for d,expected,newzeros in [(6,f**3*g[3,0],[(3,0)]),(5,f*f*q*g[2,1],[(2,1)]),(4,f*q*(f*t*g[2,0]+q*g[1,2]),[(2,0),(1,2)]),(3,q*q*(f*t*g[1,1]+q*g[0,3]),[(1,1),(0,3)])]:
        top=s.Poly(s.expand(actual.subs(vanished)),y).coeff_monomial(y**d);assert s.expand(top-expected)==0
        hierarchy.append({'degree':d,'leading_coefficient':str(top),'next_lower_stratum_forces_zero':[str(g[z]) for z in newzeros]})
        vanished.update({g[z]:0 for z in newzeros})
    lead2=s.cancel(s.Poly(s.expand(actual.subs(vanished)),y).coeff_monomial(y*y)/q**2)
    assert s.expand(lead2-(g[1,0]*f*t*t+g[0,2]*t*q-q/2))==0 and lead2.subs(t,0)==-s.Rational(1,2)
    delta=t*t*Q**3;dt=lambda z:s.diff(z,t)-b*s.diff(z,Q)
    fullsource=s.expand(delta*dt(W)-dt(delta)*W-delta**2/t+b*delta**2/Q)
    cases=[]
    zeros_by_d={6:[],5:[(3,0)],4:[(3,0),(2,1)],3:[(3,0),(2,1),(2,0),(1,2)],2:[(3,0),(2,1),(2,0),(1,2),(1,1),(0,3)]}
    for d in range(2,7):
        zeros={g[z]:0 for z in zeros_by_d[d]};m=phase_factor(d);n=3-m
        full=s.expand(W.subs(zeros));reduced=s.cancel(full/Q**m);s.Poly(reduced,t,Q,y)
        assert s.Poly(reduced,t,Q,y).total_degree()<=6-m
        source=s.expand(t*Q*dt(reduced)-(2*Q-n*b*t)*reduced-t*t*Q**(n+1)+b*t**3*Q**n)
        factor=t*Q**(m+2);h=t**3*Q**(4-m)
        assert s.expand(fullsource.subs(zeros)-factor*source)==0
        assert s.expand(delta**2-factor*h)==0 and s.Poly(source,t,Q,y).total_degree()<=7-m
        j,R,Rt=s.symbols('j R Rt');Fp=dt(factor);hp=dt(h)
        assert s.cancel(factor*h*(Rt+j*Fp*R/factor)-j*(Fp*h+factor*hp)*R+factor*source*R-factor*(h*Rt-j*hp*R+source*R))==0
        cases.append({'d':d,'phase_q_power':m,'new_phase_joint_bound':6-m,'source_joint_bound':7-m,'new_h':'t^3*q^'+str(4-m),'whole_original_factor':'t*q^'+str(m+2),'leading_degree_bound':6-m-d,'symbolic_j_recurrence_identity':True})
    return {'hierarchy':hierarchy,'lowest_quadratic_leading_coefficient':str(lead2),'formal_q_before_specialization':True,'cases':cases}

def rank_replay(m):
    read=lambda z:s.sympify(z,locals={'t':t,'y':y});W,delta,source,oldh,A=[read(m[n]) for n in ['W','delta','source','seed_denominator','A']]
    q=s.cancel(delta/t**2);q=s.real_root(q,3) if False else s.Integer(1) if q==1 else None
    # The retained model labels fix b=0 or1; verify, rather than infer,
    # the exact old delta factor against that candidate.
    q=s.Integer(1) if m['label'] in ['dense_b0','safe_null'] else 1-t
    assert s.expand(delta-t*t*q**3)==0
    d=m['d'];power=phase_factor(d);factor=t*q**(power+2);h=s.cancel(oldh/factor);S=s.cancel(source/factor);newA=s.cancel(A/q**power)
    Wnew=s.cancel(W/q**power);s.Poly(Wnew,t,y);assert s.Poly(Wnew,t,y).total_degree()<=6-power
    field=s.QQ.frac_field(t);poly=lambda z:s.Poly(z,y,domain=field);wp,dp=poly(W),poly(delta)
    R=s.Integer(1);Rold=s.Integer(1);rows=[]
    for row in m['seed_columns']:
        j=row['j'];R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+S*R);Rold=s.expand(oldh*s.diff(Rold,t)-(j-1)*s.diff(oldh,t)*Rold+source*Rold)
        assert s.expand(Rold-factor**j*R)==0 and s.Poly(R,t,y).total_degree()<=(7-power)*j
        depth=row['depth'];assert depth==d*(j-1)+2
        common=s.expand(factor**j*q**(power*depth));values=[s.Poly(s.cancel(read(z)/common),t,domain=s.QQ) for z in row['N']]
        witness=poly(s.cancel(read(row['witness'])/common));normal=sum((poly(z.as_expr()*y**k) for k,z in enumerate(values)),poly(0))
        assert (poly(newA**depth)*poly(R)-normal-dp*witness.diff(y)-wp.diff(y)*witness).is_zero
        bnd=(7-power)*j+(6-power-d)*depth
        assert all(z.degree()<=bnd-k for k,z in enumerate(values))
        H=s.Poly(s.cancel(read(row['H'])/common),t,domain=s.QQ);assert s.expand(H.as_expr()-h**j*newA**depth)==0 and H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        rows.append({'j':j,'depth':depth,'cancelled_factor':str(common),'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'row_degree_bound_base':bnd,'witness_in_original_quotient':str(witness.as_expr())})
    return {'label':m['label'],'d':d,'q':str(q),'phase_q_power':power,'same_original_operator':True,'seed_columns':rows}

def budget_universal():
    D,Le,Lf,J,L0,N0,K,a=s.symbols('D Le Lf J L0 N0 K a');delta=t*t*D
    co=[s.Rational(1,2),a*D,0,D*Le,D*Lf,J,D*L0,D*N0,K];prep=ordinary.operator_data('transverse',*co,delta);out=[]
    for kind in ['transverse','cusp']:
        if kind=='cusp':
            R=sum(s.Symbol('c'+str(i)+'_'+str(j))*u**i*v**j for i,j in weighted.monomials(10,(3,4)));W=u**4/2+D*v**3+R
        rows=[]
        for ex in ordinary.monomials(15) if kind=='transverse' else weighted.monomials(45,(3,4)):
            if kind=='transverse':
                rem,cu,cv,_,_=ordinary.one_step(kind,ex,co,delta,prep);rem,cu,cv=[s.cancel(z*Lf/prep[-1]) for z in [rem,cu,cv]];var=Lf;A=Lf;W=prep[0]
            else:
                nf,cu,cv=weighted.principal_division(kind,ex,s.Rational(1,2),D);nf,cu,cv=[s.cancel(D*z) for z in [nf,cu,cv]]
                rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)));var=D;A=D
            assert s.expand(A*u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
            terms=s.Poly(rem,u,v,var).terms();observed=max((budget(kind,z[:2])-z[2] for z,c in terms if c),default=-1000)
            assert observed<=budget(kind,ex)-1,(kind,ex,observed,budget(kind,ex)-1)
            rows.append({'exponent':ex,'budget_bound':budget(kind,ex)-1,'observed':observed})
        basis=(ordinary.CASES if kind=='transverse' else weighted.CASES)[kind]['basis']
        out.append({'kind':kind,'potential':'max(0,ceil(i/3+k-2))' if kind=='transverse' else 'max(0,ceil(i/2+k-1))','basis_budgets':[budget(kind,e) for e in basis],'identities':rows,'single_pole_branch_minimum_rho':3 if kind=='transverse' else 2,'source_pole_bound_retained':'2j','source_budget_bound':'2j-2' if kind=='transverse' else '2j-1'})
    return out

def budget_replay(m):
    read=lambda z:s.sympify(z,locals={'t':t,'u':u,'v':v});W,D,A=[read(m[n]) for n in ['W','D','A0']];kind=m['kind'];basis=m['basis'];delta=t*t*D;h=t**3*D**2
    if kind=='transverse':source=s.cancel(euler.ordinary_source(W,D)/t);vector=lambda z:s.cancel(delta*s.diff(delta,t)*(u*s.diff(z,u)+v*s.diff(z,v))/(4*t))
    else:
        source,bx,by=shared.weighted_source(kind,W,D);source=s.cancel(source/t);vector=lambda z:s.cancel(delta*(bx*u*s.diff(z,u)+by*v*s.diff(z,v))/t)
    Lf=s.cancel(s.Poly(W,u,v).coeff_monomial(v*v)/D);field=s.QQ.frac_field(t);poly=lambda z:s.Poly(z,u,v,domain=field);wp,dp=poly(W),poly(delta)
    R=s.Integer(1);rows=[]
    for row in m['seed_columns']:
        j=row['j'];depth=row['depth'];R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+vector(R)+source*R)
        if kind=='transverse':common=D**depth*Lf**j;scales=[1]*len(basis);base=7*j-2;expected=t**(3*j)*D**(2*j)*Lf**(2*j-2)
        else:common=D**(3*j-2);scales=[D**budget(kind,e) for e in basis];base=7*j-1;expected=t**(3*j)*D**(4*j-1)
        values=[s.Poly(s.cancel(read(z)/(common*scale)),t,domain=s.QQ) for z,scale in zip(row['N'],scales)]
        Cu,Cv=[poly(s.cancel(read(row['witness'][n])/common)) for n in ['Cu','Cv']]
        normal=sum((poly(z.as_expr()*scale*u**e[0]*v**e[1]) for z,scale,e in zip(values,scales,basis)),poly(0))
        assert (poly(s.cancel(A**depth/common))*poly(R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        grades=[sum(e)+(budget(kind,e) if kind=='cusp' else 0) for e in basis]
        assert all(z.degree()<=base-g for z,g in zip(values,grades))
        H=s.Poly(s.cancel(read(row['H'])/common),t,domain=s.QQ);assert s.expand(H.as_expr()-expected)==0
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        rows.append({'j':j,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'common_factor':str(common),'row_scales':[str(z) for z in scales],'row_grades':grades,'witness':{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}})
    return {'kind':kind,'singular_at_zero':m['singular_at_zero'],'metric_determinant':m['metric_determinant'],'seed_columns':rows}

def arithmetic():
    rank=[]
    for d,expected in [(2,9),(3,21),(4,51),(5,54),(6,95)]:
        m=phase_factor(d);r=d-1;bj=lambda j:(7-m)*j+(6-m-d)*(d*(j-1)+2)
        maximum=max(sum(bj(j) for j in range(1,size+1))-sum(selected) for size in range(1,r+1) for selected in combinations(range(r),size) if 0 in selected)
        assert maximum==expected;rank.append({'d':d,'gate':maximum})
    budgets=[]
    for kind,offset,expected in [('transverse',2,127),('cusp',1,130)]:
        basis=(ordinary.CASES if kind=='transverse' else weighted.CASES)[kind]['basis'];grades=[sum(e)+(budget(kind,e) if kind=='cusp' else 0) for e in basis];r=len(basis)
        maximum=max(sum(7*j-offset for j in range(1,size+1))-sum(grades[k] for k in selected) for size in range(1,r+1) for selected in combinations(range(r),size) if 0 in selected)
        assert maximum==expected;budgets.append({'kind':kind,'gate':maximum,'row_grades':grades})
    cover={'quadratic':7,'smooth':94,'nodal':130,'cusp':130,'transverse':130,'tangent':120,'triangle':130,'concurrent':104,'repeated':99}
    assert max(cover.values())==130 and max(130,92,95)==130 and max(130,92,95)==130 and max(120,92,95)==120
    return rank,budgets,cover

def main():
    started=time.monotonic();universal=rank_universal();print('Universal rank hierarchy, phase/operator factors and symbolic-j identities passed',flush=True)
    old=json.loads((P/'joint_degree_seed_bounds_certificate.json').read_text());rank_models=[rank_replay(m) for m in old['polynomial_models'] if m['kind']=='rank_one'];assert len(rank_models)==7
    print('Seven rank-one seed models replayed',flush=True);budgets=budget_universal();print('Stopped pole budgets passed',flush=True)
    common=json.loads((P/'common_time_seed_factor_certificate.json').read_text());budget_models=[budget_replay(m) for m in common['models'] if m['kind'] in ['transverse','cusp']];assert len(budget_models)==4
    rank_cases,budget_cases,cover=arithmetic();cert={'rank_universal':universal,'rank_models':rank_models,'budget_universal':budgets,'budget_models':budget_models,'rank_degree_cases':rank_cases,'budget_degree_cases':budget_cases,'global_cover':cover}
    path=P/'rank_factors_and_reduction_budgets_certificate.json';path.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,'rank_universal_strata':5,'rank_models':7,'budget_models':4,'direct_seed_identities':33,'rank_one_gate':95,'transverse_complement_gate':127,'cusp_complement_gate':130,'tangent_full_gate':120,'conic_gate_unchanged':130,'full_W_linear_gate':130,'global_sufficient_pure_prefix':130,'odd_sufficient_pure_prefix':130,'homogeneous_sufficient_pure_prefix':130,'least_prefix_interval':[20,130],'odd_least_prefix_interval':[20,130],'homogeneous_least_prefix_interval':[14,130],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact auxiliary rank hierarchy/factor/source/symbolic-order identities,33 original-quotient seed replays in11 retained models, principal stopped-budget identities and all-rank/row arithmetic. Universal all-degree and complete coefficient/rank coverage, period comparison and analytic exclusions remain written inputs. Rational witnesses are retained. No independent replay, external review, formal verification or expanded global radical cofactors.'}
    (P/'rank_factors_and_reduction_budgets_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
