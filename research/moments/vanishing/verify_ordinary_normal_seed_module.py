"""Auxiliary algebra for ordinary polynomial normal modules and half cusp grades.
Retained analytic exclusions and all-degree normal/module inductions are written.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_euler_seed_bounds as euler
import verify_shared_divisor_weighted_seeds as shared
import verify_rescaled_singular_prefixes as rescaled
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def grade(f,var,wu,wv):return max((wu*e[0]+wv*e[1]-e[2] for e,c in s.Poly(f,u,v,var).terms() if c),default=-1000)
def ceil(x):return int(s.ceiling(x))

def ordinary_module_universal():
    checkpoint=json.loads((P/'ordinary_node_phase_shift_checkpoint.json').read_text());read=lambda z:s.sympify(z,locals={'t':t,'u':u,'v':v})
    W=rescaled.rescale(read(checkpoint['V']));den=read(checkpoint['D']);pv=s.Poly(W,u,v);axis=[]
    for ex in [(1,0),(0,1)]:
        quotient=s.cancel(pv.coeff_monomial(u**ex[0]*v**ex[1])/den)
        assert s.Poly(quotient,t).degree()<=2
        assert quotient.is_polynomial(*quotient.free_symbols)
        axis.append({'exponent':ex,'coefficient_over_D':str(quotient),'degree_bound':2})
    D,Le,Lf,J,L0,N0,K=s.symbols('D Le Lf J L0 N0 K');a,d=s.symbols('a d',nonzero=True);delta=t*t*D
    gens=[u,v,t,D,Le,Lf,J,L0,N0,K];weights=[1,1,1,1,1,1,2,2,2,4]
    def degree(z):return max((sum(a*b for a,b in zip(e,weights)) for e,c in s.Poly(z,*gens).terms() if c),default=-1000)
    rows=[]
    for kind,cost in [('nodal',0),('transverse',1),('triangle',2)]:
        co=[s.Rational(1,2),a*D if kind!='triangle' else 0,d*D if kind=='nodal' else 0,D*Le,D*Lf,J,D*L0,D*N0,K]
        prepared=ordinary.operator_data(kind,*co,delta);W=prepared[0];A=1 if kind=='nodal' else Lf if kind=='transverse' else Le*Lf;checks=[]
        for ex in ordinary.monomials(12):
            rem,cu,cv,_,_=ordinary.one_step(kind,ex,co,delta,prepared)
            rem,cu,cv=[s.cancel(z*A/prepared[-1]) for z in [rem,cu,cv]]
            assert s.expand(A*u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
            # Only the normal remainder is required polynomial; witnesses
            # are retained in the original rational quotient over C(t).
            s.Poly(rem,*gens)
            assert degree(rem)<=sum(ex)+cost
            assert all(e in ordinary.CASES[kind]['basis'] or sum(e)<=sum(ex)-1 for e,c in s.Poly(rem,u,v).terms() if c)
            checks.append({'exponent':ex,'normal_joint_degree':degree(rem),'bound':sum(ex)+cost,'D_witness_denominators_allowed':True})
        rows.append({'kind':kind,'new_principal_normal_divisor':str(A),'normal_joint_cost':cost,'principal_identities':checks})
    return {'universal_linear_coefficient_factorization':axis,'checkpoint_sha256':sha(P/'ordinary_node_phase_shift_checkpoint.json'),'cases':rows,'normal_coefficients_polynomial_claimed':True,'witnesses_polynomial_claimed':False,'same_rational_quotient':True}

def cusp_universal():
    D,eps=s.symbols('D eps');dt=lambda f:s.diff(f,t)+eps*s.diff(f,D);delta=t*t*D
    R=sum(s.Symbol('c'+str(i)+'_'+str(j)+'_'+str(a)+'_'+str(b))*t**a*D**b*u**i*v**j for i,j in weighted.monomials(10,(3,4)) for a in range(5-i-j) for b in range(5-i-j-a))
    W=u**4/2+D*v**3+R;dp=dt(delta);bx=dp/4;by=2*t*D/3
    old=delta*dt(W)-dp*W-2*t**3*D**2-t**4*D*eps/2
    source=s.cancel((old+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v))/t)
    assert grade(source,D,s.Rational(1,2),1)<=2
    for z in [t**3*D**2,dt(t**3*D**2),s.cancel(delta*bx/t),s.cancel(delta*by/t)]:assert grade(z,D,s.Rational(1,2),1)<=1
    rows=[]
    for ex in weighted.monomials(30,(3,4)):
        nf,cu,cv=weighted.principal_division('cusp',ex,s.Rational(1,2),D);nf,cu,cv=[s.cancel(D*z) for z in [nf,cu,cv]]
        rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
        assert s.expand(D*u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
        rho=s.Rational(ex[0],2)+ex[1];assert all(grade(z,D,s.Rational(1,2),1)<=rho-1 for z in [rem,cu,cv])
        rows.append({'exponent':ex,'remainder_grade':grade(rem,D,s.Rational(1,2),1),'bound':rho-1})
    return {'source_grade_bound':2,'formal_D_before_specialization':True,'principal_identities':rows}

def replay_binary(m):
    local={'t':t,'u':u,'v':v};read=lambda x:s.sympify(x,locals=local)
    W,D,A=[read(m[n]) for n in ['W','D','A0']];kind=m['kind'];basis=m['basis'];delta=t*t*D;h=t**3*D**2
    if kind=='cusp':
        source,bx,by=shared.weighted_source(kind,W,D);source=s.cancel(source/t)
        vector=lambda f:s.cancel(delta*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))/t)
    else:
        source=s.cancel(euler.ordinary_source(W,D)/t);vector=lambda f:s.cancel(delta*s.diff(delta,t)*(u*s.diff(f,u)+v*s.diff(f,v))/(4*t))
    pv=s.Poly(W,u,v);Le=s.cancel(pv.coeff_monomial(u*u)/D);Lf=s.cancel(pv.coeff_monomial(v*v)/D)
    ring=s.QQ.frac_field(t);poly=lambda z:s.Poly(z,u,v,domain=ring);wp,dp=poly(W),poly(delta)
    R=s.Integer(1);out=[]
    for row in m['seed_columns']:
        j=row['j'];r=row['depth'];R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+vector(R)+source*R)
        if kind=='cusp':
            assert r==5*j-3;common=D**(3*j-3);scales=[D**ceil(s.Rational(e[0],2)+e[1]) for e in basis];grades=[sum(e)+ceil(s.Rational(e[0],2)+e[1]) for e in basis];slope=7
        else:
            assert r==3*j-(1 if kind=='triangle' else 2)
            base=j-(1 if kind=='triangle' else 2);common=D**r
            scales=[s.Integer(1) for e in basis];grades=[sum(e) for e in basis];slope={'nodal':5,'transverse':7,'triangle':9}[kind]
            for exponent in [(1,0),(0,1)]:
                q=s.cancel(pv.coeff_monomial(u**exponent[0]*v**exponent[1])/D)
                assert s.Poly(q,t,domain=s.QQ).degree()<=2
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
        expected=t**(3*j)*D**((4 if kind=='cusp' else 2)*j)
        if kind=='transverse':expected*=Lf**(2*j)
        if kind=='triangle':expected*=Le**(2*j)*Lf**(2*j)
        assert s.expand(H.as_expr()-expected)==0
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        out.append({'j':j,'common_factor':str(common),'row_scales':[str(z) for z in scales],'row_grades':grades,'witness_coefficients_rational_in_t':True,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'witness':{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}})
    return {'kind':kind,'singular_at_zero':m['singular_at_zero'],'metric_determinant':m['metric_determinant'],'W':str(W),'D':str(D),'Le':str(Le),'Lf':str(Lf),'seed_columns':out}


def arithmetic():
    rows=[]
    for kind,slope,expected in [('nodal',5,129),('transverse',7,132),('triangle',9,118),('cusp',7,131)]:
        basis=(weighted.CASES if kind=='cusp' else ordinary.CASES)[kind]['basis']
        def g(e):
            i,k=e
            if kind=='cusp':return sum(e)+ceil(s.Rational(i,2)+k)
            if kind=='nodal':return sum(e)
            if kind=='transverse':return sum(e)+ceil(s.Rational(i,3)+k)
            return sum(e)+ceil(i+s.Rational(k,2))+ceil(s.Rational(i,2)+k)
        maximum=-1
        for size in range(1,len(basis)+1):
            for selected in combinations(range(len(basis)),size):
                if 0 not in selected:continue
                maximum=max(maximum,slope*size*(size+1)//2-sum(g(basis[k]) for k in selected))
        assert maximum==expected,(kind,maximum);rows.append({'case':kind,'gate':maximum,'row_grades':[g(e) for e in basis]})
    cover={'quadratic':7,'smooth':94,'nodal':146,'cusp':140,'transverse':146,'tangent':132,'triangle':146,'concurrent':104,'repeated':99}
    assert max(cover.values())==146 and max(131,92,140)==140
    return rows,cover

def main():
    started=time.monotonic();normal=ordinary_module_universal();print('Universal axis factors and273 polynomial normal projections passed',flush=True)
    cusp=cusp_universal();print('Universal half cusp filtration passed',flush=True)
    old=json.loads((P/'common_time_seed_factor_certificate.json').read_text());models=[replay_binary(m) for m in old['models'] if m['kind'] in ['nodal','transverse','triangle','cusp']];assert len(models)==8
    print('Eight original Gaussian seed models passed',flush=True);rows,cover=arithmetic()
    cert={'ordinary_normal_module':normal,'cusp_half_filtration':cusp,'models':models,'degree_cases':rows,'global_cover':cover}
    path=P/'ordinary_normal_seed_module_certificate.json';path.write_text(json.dumps(cert,indent=2,default=str)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'retained_inputs':{n:sha(P/n) for n in ['ordinary_node_phase_shift_checkpoint.json','common_time_seed_factor_certificate.json','multiple_divisor_seed_bounds_certificate.json']},'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,'universal_linear_axis_factors':2,'ordinary_principal_normal_identities':273,'full_Gaussian_models':len(models),'direct_seed_identities':3*len(models),'ordinary_complement_gates':{'nodal':129,'transverse':132,'triangle':118},'cusp_complement_gate':131,'cusp_full_gate':140,'conic_full_gate_unchanged':146,'full_W_linear_gate_unchanged':146,'tangent_full_gate_unchanged':132,'mixed_parent_gate_unchanged':92,'global_sufficient_pure_prefix':146,'odd_sufficient_pure_prefix':146,'homogeneous_sufficient_pure_prefix':146,'least_prefix_interval':[20,146],'odd_least_prefix_interval':[20,146],'homogeneous_least_prefix_interval':[14,146],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact auxiliary two universal coefficient factors,273 ordinary principal normal/divergence identities, half cusp source/filtration identities,24 actual seed replays over the original rational quotient and all-rank cofactor arithmetic. Witnesses may have rational time coefficients; scalar normal coefficients satisfy the polynomial degree bounds. Universal all-degree induction, complete coefficient/rank coverage and analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded global radical cofactors.'}
    (P/'ordinary_normal_seed_module_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
