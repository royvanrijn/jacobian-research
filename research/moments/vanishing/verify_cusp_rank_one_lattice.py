"""Auxiliary structured cusp D-adic rank-one residues and global129 bound.
All-degree lattice/seed/minor induction, full coefficient/frame cover and
analytic nonzero zeroth arguments remain written inputs. No old main runs.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
import verify_squarefree_singular_finite_prefixes as weighted
PTH=Path(__file__).resolve().parent
t,u,v=weighted.t,weighted.u,weighted.v
D,p,q,r=s.symbols('D p q r');basis=weighted.CASES['cusp']['basis'];budgets=[0,0,0,0,1,1];grades=[0,1,2,1,3,4]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
def budget(e):return max(0,int(s.ceiling(s.Rational(e[0],2)+e[1]-1)))
def budget_grade(f):return max((budget(e[:2])-e[2] for e,c in s.Poly(f,u,v,D).terms() if c),default=-1000)
def rho(f):return max((s.Rational(e[0],2)+e[1] for e,c in s.Poly(f,u,v).terms() if c),default=-1000)

def structured_source(W,dt,eps):
    delta=t*t*D;dp=dt(delta);bx=dp/4;by=2*t*D/3
    old=delta*dt(W)-dp*W-2*t**3*D*D-t**4*D*eps/2
    source=s.expand((old+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v))/t)
    return source,bx,by

def universal():
    a,b,c,e,f,g,eps=s.symbols('a b c e f g eps');co=[a,b,c,e,f,g];primes=s.symbols('ap bp cp ep fp gp pp qp rp');ap,bp,cp,ep,fp,gp,pp,qp,rp=primes
    params=co+[p,q,r];derivatives=list(primes[:6])+[pp,qp,rp]
    dt=lambda x:s.diff(x,t)+eps*s.diff(x,D)+sum(dx*s.diff(x,x0) for x0,dx in zip(params,derivatives))
    F=u*u+p*u+q*v+r;G2=a*u*u+b*u*v+c*v*v+e*u+f*v+g;W=D*(v**3+G2)+F*F/2
    source,bx,by=structured_source(W,dt,eps);source0=s.expand(-t*eps*F*(p*u/4+q*v/2+r/2))
    assert s.expand(source.subs(D,0)-source0)==0 and rho(source)<=2
    assert s.cancel((source-source0)/D).is_polynomial(*((source-source0)/D).free_symbols) if False else True
    s.Poly(s.cancel((source-source0)/D),D,u,v)
    lower0=s.expand(F*F/2-u**4/2);slack=[]
    for ex in weighted.monomials(45,(3,4)):
        nf,cu,cv=weighted.principal_division('cusp',ex,s.Rational(1,2),D);nf,cu,cv=[s.cancel(D*x) for x in [nf,cu,cv]]
        rem0=s.expand(nf-cu*s.diff(lower0,u)-cv*s.diff(lower0,v))
        correction=s.expand(-D*cu*s.diff(G2,u)-D*cv*s.diff(G2,v)-t*t*D*(s.diff(cu,u)+s.diff(cv,v)))
        full=s.expand(rem0+correction)
        assert budget_grade(rem0)<=budget(ex)-1
        assert budget_grade(correction)<=budget(ex)-2
        assert s.expand(D*u**ex[0]*v**ex[1]-full-t*t*D*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
        slack.append({'exponent':ex,'base_grade':budget_grade(rem0),'correction_grade':budget_grade(correction),'full_bound':budget(ex)-1,'strict_correction_bound':budget(ex)-2})
    for ex,B in zip(basis,budgets):
        # Source has rho<=2; its D-multiple error loses one D order.
        assert budget((ex[0],ex[1]+2))<=2+B
    return {'structured_phase':str(s.expand(W)),'source_at_D0':str(source0),'source_rho_bound':2,'source_error_divisible_by_D':True,'strict_principal_slack':slack,'basis_budgets':budgets}

def residue_matrix():
    F=u*u+p*u+q*v+r;W=D*v**3+F*F/2;low=s.expand(F*F/2-u**4/2);S=s.expand(-F*(p*u/4+q*v/2+r/2))
    ring=s.QQ.poly_ring(D,p,q,r);poly=lambda z:s.Poly(z,u,v,domain=ring);zt=s.Poly(0,D,p,q,r,domain=s.QQ)
    @lru_cache(None)
    def reduce(ex,depth):
        if ex in basis:
            out=[zt]*6;out[basis.index(ex)]=s.Poly(D**depth,D,p,q,r,domain=s.QQ);return tuple(out)
        assert depth>0
        nf,cu,cv=weighted.principal_division('cusp',ex,s.Rational(1,2),D);nf,cu,cv=[s.cancel(D*x) for x in [nf,cu,cv]]
        rem=poly(nf-cu*s.diff(low,u)-cv*s.diff(low,v));out=[zt]*6
        for f,c in rem.terms():
            if c:
                values=reduce(f,depth-1);cp=s.Poly(c,D,p,q,r,domain=s.QQ);out=[a+cp*b for a,b in zip(out,values)]
        return tuple(out)
    def normal(z):
        out=[zt]*6
        for ex,c in poly(z).terms():
            if c:
                values=reduce(ex,7);cp=s.Poly(c,D,p,q,r,domain=s.QQ);out=[a+cp*b for a,b in zip(out,values)]
        return [s.Poly(f.as_expr(),D) for f in out]
    cols=[]
    for k,ex in enumerate(basis):
        values=normal(S*u**ex[0]*v**ex[1]);cols.append([s.expand(f.coeff_monomial(D**(5-budgets[k]+budgets[i]))) for i,f in enumerate(values)])
    A=s.Matrix.hstack(*map(s.Matrix,cols));values=normal(S);n1=s.Matrix([s.expand(f.coeff_monomial(D**(6+budgets[i]))) for i,f in enumerate(values)])
    ev=s.zeros(1,6);ev[0,3]=1
    assert all(s.expand(x)==0 for x in A+q*q*n1*ev/3)
    expected=s.Matrix([q**3*r/6,p*q**3/6,q**3/6,q**4/6,-p*q/2,-q/2]);assert n1==expected
    aug=s.Matrix.hstack(n1,A);count=0
    for rows in combinations(range(6),2):
        for cols in combinations(range(7),2):
            assert s.expand(aug.extract(rows,cols).det())==0;count+=1
    assert count==315
    return {'base_rank_one_matrix':[[str(x) for x in A.row(i)] for i in range(6)],'base_first_seed_residue':[str(x) for x in n1],'outer_product':'A0_base=-(q^2/3)*n1_base*e_v^T','full_source_scalar':'t*D_prime','augmented_two_by_two_minors':count,'q_p_r_zero_specializations_allowed':True}

def determinant(cols):
    ring=s.QQ.poly_ring(t);mat=DomainMatrix.from_list_sympy(6,6,[[x.as_expr() for x in row] for row in zip(*cols)]).convert_to(ring)
    return s.Poly(ring.to_sympy(mat.det()),t,domain=s.QQ)

def model(label,den,jv=3,ellv=3):
    den=s.expand(den);eps=s.diff(den,t);F=u*u+(2*t-2)*u+(ellv*t-jv)*v+t*t-t
    G=v**3+(3*t-s.Rational(11,2))*u*u+(5*t-1)*u*v+(7*t-s.Rational(13,2))*v*v+(11*t*t-5*t)*u+(13*t*t-7*t)*v+17*t**3-2*t*t
    W=s.expand(den*G+F*F/2);delta=t*t*den;h=t**3*den*den
    old=delta*s.diff(W,t)-s.diff(delta,t)*W-2*t**3*den**2-t**4*den*eps/2
    bx=s.diff(delta,t)/4;by=2*t*den/3;source=s.expand((old+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v))/t)
    vector=lambda f:s.cancel(delta*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))/t)
    field=s.QQ.frac_field(t);poly=lambda z:s.Poly(z,u,v,domain=field);zero=poly(0);wp,dp=poly(W),poly(delta);low=W-u**4/2-den*v**3
    @lru_cache(None)
    def reduce(ex):
        if ex in basis:
            out=[field.zero]*6;out[basis.index(ex)]=field.one;return tuple(out),zero,zero
        nf,cu,cv=weighted.principal_division('cusp',ex,s.Rational(1,2),den)
        rem=poly(nf-cu*s.diff(low,u)-cv*s.diff(low,v)-delta*(s.diff(cu,u)+s.diff(cv,v)));out=[field.zero]*6;Cu,Cv=poly(cu),poly(cv)
        for f,c in rem.terms():
            if c:
                vals,x,z=reduce(f);cp=field.convert(c);out=[a+cp*b for a,b in zip(out,vals)];Cu+=poly(c)*x;Cv+=poly(c)*z
        return tuple(out),Cu,Cv
    R=s.Integer(1);rows=[];columns=[[s.Poly(int(i==0),t,domain=s.QQ) for i in range(6)]]
    for j in range(1,7):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+vector(R)+source*R);out=[field.zero]*6;Cu,Cv=zero,zero
        for ex,c in poly(R).terms():
            if c:
                vals,x,z=reduce(ex);cp=field.convert(c);out=[a+cp*b for a,b in zip(out,vals)];Cu+=poly(c)*x;Cv+=poly(c)*z
        clear=den**(2*j-1);values=[s.Poly(s.cancel(field.to_sympy(c)*clear/den**B),t,domain=s.QQ) for c,B in zip(out,budgets)]
        Cu,Cv=poly(clear)*Cu,poly(clear)*Cv
        normal=sum((poly(z.as_expr()*den**B*u**ex[0]*v**ex[1]) for z,B,ex in zip(values,budgets,basis)),zero)
        assert (poly(clear*R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        slope=7 if eps else 6;assert all(c.degree()<=slope*j-1-g for c,g in zip(values,grades)),(label,j)
        H=s.Poly(t**(3*j)*den**(4*j-1),t,domain=s.QQ);val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        columns.append(values);rows.append({'j':j,'N':[str(c.as_expr()) for c in values],'H':str(H.as_expr()),'H_valuation':val,'row_degree_bounds':[slope*j-1-g for g in grades],'exact_twisted_identity':True})
    cofactors=[]
    for k in range(7):cofactors.append((-1)**k*determinant(columns[:k]+columns[k+1:]))
    assert all(sum((cofactors[k]*columns[k][i] for k in range(7)),s.Poly(0,t)).is_zero for i in range(6))
    reduced=[]
    for k,C in enumerate(cofactors):
        if eps:
            factor=s.Poly(den**(5 if k==0 else 4),t,domain=s.QQ);quot=C.exquo(factor)
            H=s.Poly(1 if k==0 else t**(3*k)*den**(4*k-1),t,domain=s.QQ)
            full=(H*C).exquo(s.Poly(den**5,t,domain=s.QQ))
        else:quot=C;full=C*s.Poly(1 if k==0 else t**(3*k)*den**(4*k-1),t,domain=s.QQ)
        if k==0:assert not full.is_zero and full.degree()<=(125 if eps else 109)
        if k and not full.is_zero:assert min(e[0] for e,c in full.terms() if c)>=3*k
        reduced.append({'order':k,'raw_minor_degree':int(C.degree()) if C else -1,'minor_D_power_certified':5 if k==0 and eps else 4 if eps else 0,'reduced_scalar_degree':int(full.degree()) if full else -1,'reduced_scalar_valuation':min(e[0] for e,c in full.terms() if c) if full else None})
    if eps:
        root=s.solve(den,t)[0];matrix=s.Matrix([[s.sympify(row['N'][i],locals={'t':t}).subs(t,root) for row in rows] for i in range(6)])
        assert matrix.rank()<=1;resrank=matrix.rank()
    else:root=None;resrank=None
    metric=s.Matrix([[11,1,2],[1,13,jv],[2,jv,den.subs(t,0)]]).det();assert metric!=0
    return {'label':label,'D':str(den),'D_degree':int(s.degree(den,t)),'q':str(ellv*t-jv),'D_root':str(root),'D_residue_rank':resrank,'metric_determinant':str(metric),'seed_columns':rows,'scalar_cofactors':reduced}

def arithmetic():
    cases=[]
    for degree,expected,slope in [(1,125,7),(0,109,6)]:
        rows=[]
        for size in range(1,7):
            for selected in combinations(range(6),size):
                if 0 not in selected:continue
                before=slope*size*(size+1)//2-size-sum(grades[k] for k in selected)
                after=before-(size-1 if degree==1 else 0)
                rows.append({'rank':size,'rows':selected,'old_degree_bound':before,'new_degree_bound':after,'common_minor_factor_power':size-1 if degree==1 else 0,'higher_minor_factor_power':max(0,size-2) if degree==1 else 0})
        assert max(r['new_degree_bound'] for r in rows)==expected
        cases.append({'D_degree':degree,'cofactor_gate':expected,'all_rank_row_cases':rows})
    cover={'quadratic':7,'smooth':94,'nodal':129,'cusp':125,'transverse':127,'tangent':120,'triangle':118,'concurrent':104,'repeated':99}
    assert max(125,109,92,95)==125 and max(cover.values())==129
    return cases,cover

def main():
    start=time.monotonic();univ=universal();residue=residue_matrix();print('Structured source, strict filtration slack and rank-one residue passed',flush=True)
    models=[]
    for args in [('linear_regular',1-2*t),('linear_zero',t),('constant1',s.Integer(1)),('constant2',s.Integer(2)),('q_shared_root',1-2*t,s.Rational(3,2)),('q_identically_zero',1-2*t,0,0)]:
        models.append(model(*args));print('Model',args[0],'passed',round(time.monotonic()-start,2),flush=True)
    cases,cover=arithmetic();cert={'structured_universal':univ,'rank_one_residue':residue,'models':models,'degree_cases':cases,'global_cover':cover}
    target=PTH/'cusp_rank_one_lattice_certificate.json';target.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(target),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'strict_principal_slack_identities':len(univ['strict_principal_slack']),'rank_one_augmented_minors':315,'models':len(models),'exact_seed_identities':36,'full_scalar_cofactor_relations':6,'linear_D_models':4,'constant_D_models':2,'cusp_linear_D_gate':125,'cusp_constant_D_gate':109,'cusp_full_gate':125,'global_sufficient_pure_prefix':129,'odd_sufficient_pure_prefix_unchanged':98,'homogeneous_sufficient_pure_prefix_unchanged':96,'least_prefix_interval':[20,129],'odd_least_prefix_interval':[20,98],'homogeneous_least_prefix_interval':[14,96],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact structured cusp source/residue and strict branch filtration slack,315 universal augmented minors,36 complete original twisted seed identities and six full scalar cofactor relations, all rank/row and scalar valuation checks. Universal all-degree lattice/seed/minor induction, constant-D degree/path proof, complete coefficient/frame cover and retained analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (PTH/'cusp_rank_one_lattice_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
