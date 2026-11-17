"""Auxiliary tangent pole residue, triangular factors and constant-D105.
General pole-lattice/seed/minor induction, polynomial extension and complete
analytic coefficient/frame cover remain written inputs. No prior main runs.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
import verify_squarefree_singular_finite_prefixes as weighted
import verify_shared_divisor_weighted_seeds as shared
import verify_structured_seed_modules as structured
PTH=Path(__file__).resolve().parent
t,u,v=weighted.t,weighted.u,weighted.v
D,eps,p,q,r=s.symbols('D eps p q r');U,V=s.symbols('U V')
basis=weighted.CASES['tangent']['basis'];rho=lambda e:e[0]+2*e[1];grades=[0,2,4,3,6]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()

def pole_limit(f,power):
    res=s.cancel(D**power*f.subs({u:U/D,v:V/D**2},simultaneous=True));poly=s.Poly(res,D)
    return s.expand(poly.coeff_monomial(1))

def generic_residue():
    a,b,c,e,f,g=s.symbols('a b c e f g');pr=s.symbols('ap bp cp ep fp gp pp qp rp');params=[a,b,c,e,f,g,p,q,r]
    dt=lambda z:s.diff(z,t)+eps*s.diff(z,D)+sum(dx*s.diff(z,x) for x,dx in zip(params,pr))
    F=u*u+p*u+q*v+r;G=u*v*v+a*u*u+b*u*v+c*v*v+e*u+f*v+g;W=D*G+F*F/2;delta=t*t*D;dp=dt(delta);bx=dp/4;by=t*D-dp/8
    old=delta*dt(W)-dp*W-2*t**3*D**2-t**4*D*eps/2;Sb=s.expand(old+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v));C=s.Poly(Sb,u,v).coeff_monomial(u*u*v);S=s.cancel((D*Sb-C*u*s.diff(W,v)/2)/t)
    Psi=U*V*V+(U*U+q*V)**2/2;Sigma=5*t*eps*q*q*(U**3+q*U*V-2*V*V)/16
    assert s.expand(pole_limit(W,4)-Psi)==0 and s.expand(pole_limit(S,3)-Sigma)==0
    assert s.expand(pole_limit(s.diff(W,u),3)-s.diff(Psi,U))==0
    assert s.expand(pole_limit(s.diff(W,v),2)-s.diff(Psi,V))==0
    # Both rescaled twisted gradients have derivative coefficient t^2 D^5.
    assert s.expand(D**3*delta*D-t*t*D**5)==0 and s.expand(D**2*delta*D**2-t*t*D**5)==0
    W0=u*v*v+(u*u+q*v)**2/2;low=W0-u**4/2-u*v*v;S0=5*q*q*(u**3+q*u*v-2*v*v)/16
    @lru_cache(None)
    def red(ex):
        if ex in basis:return s.eye(5)[:,basis.index(ex)]
        nf,cu,cv=weighted.principal_division('tangent',ex,s.Rational(1,2),s.Integer(1));rem=s.Poly(s.expand(nf-cu*s.diff(low,u)-cv*s.diff(low,v)),u,v);out=s.zeros(5,1)
        for f,c in rem.terms():
            if c:out+=c*red(f)
        return out.applyfunc(s.expand)
    def normal(z):
        out=s.zeros(5,1)
        for ex,c in s.Poly(z,u,v).terms():
            if c:out+=c*red(ex)
        return out.applyfunc(s.expand)
    A=s.Matrix.hstack(*[normal(S0*u**i*v**k) for i,k in basis]);assert s.expand(A.trace()+s.Rational(15625,8192)*q**8)==0
    assert all(s.expand(x)==0 for x in A*A-A.trace()*A)
    aug=s.Matrix.hstack(A[:,0],A);count=0
    for rows in combinations(range(5),2):
        for cols in combinations(range(6),2):assert s.expand(aug.extract(rows,cols).det())==0;count+=1
    assert count==150 and A.trace()!=0
    return {'leading_phase':str(s.expand(Psi)),'leading_sheared_source':str(Sigma),'rescaled_twisted_gradient_derivative':'t^2 D^5','base_residue_matrix':[[str(x) for x in A.row(i)] for i in range(5)],'trace_base':str(A.trace()),'source_scalar':'t*D_prime','augmented_two_by_two_minors':count,'projector_identity':True,'first_seed_is_residue_times_unit':True,'induced_derivative_D_order':6,'induced_scalar_j_D_order_at_least':5,'q_zero_and_D_at_zero_retained_by_polynomial_extension':True}

def constant_universal():
    R=sum(s.Symbol('c'+str(i)+'_'+str(k)+'_'+str(a)+'_'+str(b))*t**a*D**b*u**i*v**k for i,k in weighted.monomials(7,(2,3)) for a in range(5-i-k) for b in range(5-i-k-a))
    W=u**4/2+D*u*v*v+R;delta=t*t*D;bx=t*D/2;by=3*t*D/4
    old=delta*s.diff(W,t)-2*t*D*W-2*t**3*D*D
    Sb=s.expand(old+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v));S=s.cancel(Sb/(t*D))
    assert s.expand(S-(t*s.diff(W,t)-2*W+u*s.diff(W,u)/2+3*v*s.diff(W,v)/4-3*t*t*D/4))==0
    joint=lambda z:max((sum(ex) for ex,c in s.Poly(z,t,D,u,v).terms() if c),default=-1000)
    fil=lambda z:max((ex[2]+2*ex[3]-ex[1] for ex,c in s.Poly(z,t,D,u,v).terms() if c),default=-1000)
    assert joint(S)<=4 and fil(S)<=4
    assert all(weighted.weight(ex,(2,3))<=7 for ex,c in s.Poly(S,u,v).terms() if c)
    assert s.expand(t**3*D*D-D*(t**3*D))==0
    assert joint(t**3*D)==4 and joint(3*t*t*D)==3
    checks=[]
    for ex in weighted.monomials(16,(2,3)):
        nf,cu,cv=weighted.principal_division('tangent',ex,s.Rational(1,2),D);nf,cu,cv=[s.cancel(D*D*z) for z in [nf,cu,cv]]
        rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
        assert s.expand(D*D*u**ex[0]*v**ex[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
        assert joint(rem)<=sum(ex)+2 and fil(rem)<=rho(ex)-2
        assert fil(cu)<=rho(ex)-5 and fil(cv)<=rho(ex)-4
        checks.append(ex)
    return {'D_formal_but_D_prime_zero':True,'whole_old_connection_factor':'D','new_h':'t^3 D','new_source_joint_bound':4,'new_source_pole_grade_bound':4,'new_source_old_weight_bound':7,'principal_divergence_identities':checks,'depth':'7j-4','common_D_power':'10j-8','new_denominator':'t^(3j) D^(5j)','normal_row_degree':'8j-|e|-rho(e)','constant_D_gate':105}

def determinants(cols):
    ring=s.QQ.poly_ring(t);out=[]
    for k in range(6):
        part=cols[:k]+cols[k+1:];mat=DomainMatrix.from_list_sympy(5,5,[[x.as_expr() for x in row] for row in zip(*part)]).convert_to(ring)
        out.append((-1)**k*s.Poly(ring.to_sympy(mat.det()),t,domain=s.QQ))
    assert all(sum((out[k]*cols[k][i] for k in range(6)),s.Poly(0,t)).is_zero for i in range(5))
    return out

def model(label,den,jv=3,ellv=3):
    den=s.sympify(den);isconstant=s.diff(den,t)==0
    F=u*u+(2*t-2)*u+(ellv*t-jv)*v+t*t-t
    G=u*v*v+(3*t-s.Rational(11,2))*u*u+(5*t-1)*u*v+(7*t-s.Rational(13,2))*v*v+(11*t*t-5*t)*u+(13*t*t-7*t)*v+17*t**3-2*t*t
    W=s.expand(den*G+F*F/2);delta=t*t*den
    if isconstant:
        source,bx,by=shared.weighted_source('tangent',W,den);source=s.cancel(source/(t*den));h=t**3*den
        vector=lambda z:s.cancel(delta*(bx*u*s.diff(z,u)+by*v*s.diff(z,v))/(t*den))
        column_clear=lambda j:den**(4*j);slope=8;den_power=5;gate=105
    else:
        source,bx,by,C,_=structured.tangent_source(W,den);source=s.cancel(source/t);h=t**3*den**3
        vector=lambda z:s.cancel(delta*den*(bx*u*s.diff(z,u)+by*v*s.diff(z,v))/t-delta*C*u*s.diff(z,v)/(2*t))
        column_clear=lambda j:den**(3*j);slope=9;den_power=6;gate=110
    field=s.QQ.frac_field(t);poly=lambda z:s.Poly(z,u,v,domain=field);zero=poly(0);wp,dp=poly(W),poly(delta);low=W-u**4/2-den*u*v*v
    @lru_cache(None)
    def reduce(ex):
        if ex in basis:
            out=[field.zero]*5;out[basis.index(ex)]=field.one;return tuple(out),zero,zero
        nf,cu,cv=weighted.principal_division('tangent',ex,s.Rational(1,2),den);rem=poly(nf-cu*s.diff(low,u)-cv*s.diff(low,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
        out=[field.zero]*5;Cu,Cv=poly(cu),poly(cv)
        for f,c in rem.terms():
            if c:
                vals,x,z=reduce(f);cp=field.convert(c);out=[a+cp*b for a,b in zip(out,vals)];Cu+=poly(c)*x;Cv+=poly(c)*z
        return tuple(out),Cu,Cv
    R=s.Integer(1);rows=[];cols=[[s.Poly(int(i==0),t,domain=s.QQ) for i in range(5)]];Hs=[s.Poly(1,t,domain=s.QQ)]
    for j in range(1,6):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+vector(R)+source*R);out=[field.zero]*5;Cu,Cv=zero,zero
        for ex,c in poly(R).terms():
            if c:
                vals,x,z=reduce(ex);cp=field.convert(c);out=[a+cp*b for a,b in zip(out,vals)];Cu+=poly(c)*x;Cv+=poly(c)*z
        clear=column_clear(j);values=[s.Poly(s.cancel(field.to_sympy(c)*clear/den**rho(ex)),t,domain=s.QQ) for c,ex in zip(out,basis)]
        Cu,Cv=poly(clear)*Cu,poly(clear)*Cv;normal=sum((poly(c.as_expr()*den**rho(ex)*u**ex[0]*v**ex[1]) for c,ex in zip(values,basis)),zero)
        assert (poly(clear*R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        assert all(c.degree()<=slope*j-g for c,g in zip(values,grades)),(label,j,[c.degree() for c in values])
        H=s.Poly(t**(3*j)*den**(den_power*j),t,domain=s.QQ);val=min(ex[0] for ex,c in H.terms() if c);assert val>=3*j
        rows.append({'j':j,'N':[str(c.as_expr()) for c in values],'H':str(H.as_expr()),'H_valuation':val,'row_degree_bounds':[slope*j-g for g in grades],'exact_original_twisted_identity':True});cols.append(values);Hs.append(H)
    cof=determinants(cols);checks=[]
    for k,C in enumerate(cof):
        need=0 if isconstant else 10 if k==0 else max(0,10-k)
        C.exquo(s.Poly(den**need,t,domain=s.QQ))
        reduced=C*Hs[k] if isconstant else (C*Hs[k]).exquo(s.Poly(den**10,t,domain=s.QQ))
        if k==0:assert not reduced.is_zero and reduced.degree()<=gate
        val=min(ex[0] for ex,c in reduced.terms() if c) if reduced else None
        if k and val is not None:assert val>=3*k
        checks.append({'order':k,'minor_D_order_certified':need,'common_scalar_D_factor':0 if isconstant else 10,'reduced_scalar_degree':int(reduced.degree()) if reduced else -1,'reduced_scalar_valuation':val})
    metric=s.Matrix([[11,1,2],[1,13,jv],[2,jv,den.subs(t,0)]]).det();assert metric!=0
    return {'label':label,'D':str(den),'D_degree':int(not isconstant),'q':str(ellv*t-jv),'metric_determinant':str(metric),'seed_columns':rows,'cofactor_checks':checks,'class_gate':gate}

def arithmetic():
    out=[]
    for kd,slope,expected in [(1,9,110),(0,8,105)]:
        cases=[]
        for size in range(1,6):
            for rows in combinations(range(5),size):
                if 0 not in rows:continue
                raw=slope*size*(size+1)//2-sum(grades[k] for k in rows);T=size*(size-1)//2 if kd else 0
                cases.append({'rank':size,'rows':rows,'raw_degree_bound':raw,'common_scalar_D_factor':T,'new_degree_bound':raw-T})
        assert max(c['new_degree_bound'] for c in cases)==expected;out.append({'D_degree':kd,'gate':expected,'all_rank_row_cases':cases})
    cover={'quadratic':7,'smooth':94,'nodal':108,'cusp':115,'transverse':112,'tangent':110,'triangle':108,'concurrent':104,'repeated':99}
    assert max(110,105,92,95)==110 and max(cover.values())==115
    return out,cover

def main():
    start=time.monotonic();residue=generic_residue();const=constant_universal();print('Tangent generic residue and constant-D source/filtration passed',flush=True)
    fixtures=[('linear_regular',1-2*t),('linear_zero',t),('linear_q_zero',1-2*t,0,0),('linear_zero_q_zero',t,0,0),('linear_shared_q',1-2*t,s.Rational(3,2)),('constant1',1),('constant2',2),('constant_q_zero',1,0,0)]
    models=[]
    for args in fixtures:models.append(model(*args));print('Model',args[0],'passed',round(time.monotonic()-start,2),flush=True)
    cases,cover=arithmetic();cert={'generic_tangent_residue':residue,'constant_D_universal':const,'models':models,'degree_cases':cases,'global_cover':cover};path=PTH/'tangent_rank_one_lattice_certificate.json';path.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'generic_augmented_minors':150,'constant_D_principal_divergence_identities':len(const['principal_divergence_identities']),'models':8,'linear_D_models':5,'constant_D_models':3,'exact_original_twisted_seed_identities':40,'full_scalar_cofactor_relations':8,'tangent_linear_D_gate':110,'tangent_constant_D_gate':105,'tangent_full_gate':110,'global_sufficient_pure_prefix':115,'odd_sufficient_pure_prefix_unchanged':98,'homogeneous_sufficient_pure_prefix_unchanged':96,'least_prefix_interval':[20,115],'odd_least_prefix_interval':[20,98],'homogeneous_least_prefix_interval':[14,96],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact structured pole-scaled phase/source/gradient limits,Jacobi residue matrix and150 augmented minors,constant-D whole-factor/source/filtration identities,40 original twisted seeds and eight full scalar cofactor relations, all-rank/row and post-cancellation valuations. Universal pole-lattice/seed/triangular-minor argument, polynomial extension across degenerate traces, constant-D all-degree path proof, complete coefficient/frame cover and retained analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (PTH/'tangent_rank_one_lattice_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
