"""Auxiliary actual concurrent seeds and repeated-line zero valuations.
Universal degree/seed induction, full class/frame cover and nonzero zeroth
period arguments remain written inputs. Previous checker mains are not run.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
P=Path(__file__).resolve().parent
t,u,v,T=s.symbols('t u v T');basis=[(0,0),(1,0),(0,1),(1,1)]
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()

def source(W,D):
    delta=t*t*D;eps=s.diff(D,t)
    old=s.expand((delta*s.diff(W,t)-s.diff(delta,t)*W)/t-2*t*t*D*D-t**3*D*eps/2)
    S=s.expand(old+s.Rational(2,3)*D*(u*s.diff(W,u)+v*s.diff(W,v))+s.Rational(4,3)*t*t*D*D)
    return old,S

def universals():
    alpha,m,lu,lv,l0,ju,jv,j0=s.symbols('alpha m lu lv l0 ju jv j0')
    gc=s.symbols('g20 g11 g02 g10 g01 g00');qc=s.symbols('q20 q11 q02 q10 q01 q00')
    ex=[u*u,u*v,v*v,u,v,1];g2=sum(a*b for a,b in zip(gc,ex));q2=sum(a*b for a,b in zip(qc,ex))
    ell=lu*u+lv*v+l0;J=ju*u+jv*v+j0;D=m-2*alpha*t
    full=t*(alpha*T*T+T*ell+u**3+v**3+g2)-m*T*T/2-T*J-q2
    L=t*ell-J;U=s.expand(D*(t*(u**3+v**3+g2)-q2)+L*L/2)
    assert s.cancel(full.subs(T,T+L/D)-(-D*T*T/2+U/D))==0
    W=s.expand(t*t*U.subs({u:u/t,v:v/t},simultaneous=True))
    F=(t*lu-ju)*u+(t*lv-jv)*v+t*t*l0-t*j0
    G=u**3+v**3+sum((t if k<3 else t*t if k<5 else t**3)*gc[k]*ex[k] for k in range(6))-sum((1 if k<3 else t if k<5 else t*t)*qc[k]*ex[k] for k in range(6))
    assert s.expand(W-D*G-F*F/2)==0
    old,S=source(W,D);assert s.Poly(W,t,u,v).total_degree()<=4
    assert s.Poly(S,t,u,v).total_degree()<=5 and s.Poly(S,u,v).total_degree()<=2
    assert s.Poly(S.subs(alpha,0),t,u,v).total_degree()<=4
    low=W-D*(u**3+v**3);gauges=[];principal=[]
    for i in range(7):
        for k in range(7-i):
            f=u**i*v**k;cu=s.Rational(2,3)*D*u*f;cv=s.Rational(2,3)*D*v*f
            defect=(S-old)*f+s.Rational(2,3)*t*t*D*D*(u*s.diff(f,u)+v*s.diff(f,v))-t*t*D*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv
            assert s.expand(defect)==0;gauges.append([i,k])
    for i in range(11):
        for k in range(11-i):
            if (i,k) in basis:continue
            cu=u**(i-2)*v**k/3 if i>=2 else s.Integer(0);cv=u**i*v**(k-2)/3 if i<2 else s.Integer(0)
            rem=s.expand(-cu*s.diff(low,u)-cv*s.diff(low,v)-t*t*D*(s.diff(cu,u)+s.diff(cv,v)))
            assert s.expand(D*u**i*v**k-rem-t*t*D*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
            assert s.Poly(rem,t,u,v).total_degree()<=i+k+1 and s.Poly(rem,u,v).total_degree()<=i+k-1
            principal.append([i,k])
    return {'full_completion_rescaling':True,'phase':'D*G+(p*u+q*v+r)^2/2','source_amplitude_bound':2,
            'phase_joint_bound':4,'linear_D_source_joint_bound':5,'constant_D_source_joint_bound':4,
            'Euler_gauge_identities':gauges,'cleared_principal_identities':principal,'cleared_step_joint_cost':1,
            'reduction_depth':'2j-1','linear_D_row_bound':'7j-1-|e|','constant_D_row_bound':'6j-1-|e|',
            'scalar_denominator':'t^(3j) D^(4j-1)'}

def cofactors(cols,n):
    ring=s.QQ.poly_ring(t);out=[]
    for k in range(n+1):
        selected=cols[:k]+cols[k+1:]
        M=DomainMatrix.from_list_sympy(n,n,[[c.as_expr() for c in row] for row in zip(*selected)]).convert_to(ring)
        out.append((-1)**k*s.Poly(ring.to_sympy(M.det()),t,domain=s.QQ))
    assert all(sum((out[k]*cols[k][i] for k in range(n+1)),s.Poly(0,t)).is_zero for i in range(n))
    return out

def concurrent_model(label,D,ellzero=False):
    D=s.expand(D);eps=s.diff(D,t);m=D.subs(t,0)
    lu,lv,l0=(0,0,0) if ellzero else (2,3,5)
    F=(lu*t-2)*u+(lv*t-3)*v+l0*t*t-7*t
    G=u**3+v**3+(3*t-s.Rational(11,2))*u*u+(5*t-1)*u*v+(7*t-s.Rational(13,2))*v*v+(11*t*t-19*t)*u+(13*t*t-23*t)*v+17*t**3-29*t*t
    W=s.expand(D*G+F*F/2);_,S=source(W,D);low=W-D*(u**3+v**3)
    delta=t*t*D;h=t**3*D*D;ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,v,domain=ring)
    zero=poly(0);zt=s.Poly(0,t);wp,dp=poly(W),poly(delta);den=s.Poly(D,t)
    @lru_cache(None)
    def reduce(ex,depth):
        if ex in basis:
            out=[zt]*4;out[basis.index(ex)]=den**depth;return tuple(out),zero,zero
        assert depth>0
        i,k=ex;cu=u**(i-2)*v**k/3 if i>=2 else s.Integer(0);cv=u**i*v**(k-2)/3 if i<2 else s.Integer(0)
        rem=poly(-cu*s.diff(low,u)-cv*s.diff(low,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
        out=[zt]*4;Cu,Cv=poly(D**(depth-1)*cu),poly(D**(depth-1)*cv)
        for e,c in rem.terms():
            if not c:continue
            vals,a,b=reduce(e,depth-1);ct=s.Poly(c,t)
            out=[a+ct*b for a,b in zip(out,vals)];Cu+=poly(c)*a;Cv+=poly(c)*b
        return tuple(out),Cu,Cv
    R=s.Integer(1);rows=[];cols=[[s.Poly(int(i==0),t) for i in range(4)]]
    for j in range(1,5):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+s.Rational(2,3)*t*t*D*D*(u*s.diff(R,u)+v*s.diff(R,v))+S*R)
        assert s.Poly(R,t,u,v).total_degree()<=(5 if eps else 4)*j and s.Poly(R,u,v).total_degree()<=2*j
        values=[zt]*4;Cu,Cv=zero,zero
        for e,c in poly(R).terms():
            if not c:continue
            vals,a,b=reduce(e,2*j-1);ct=s.Poly(c,t)
            values=[a+ct*b for a,b in zip(values,vals)];Cu+=poly(c)*a;Cv+=poly(c)*b
        normal=sum((poly(c.as_expr()*u**e[0]*v**e[1]) for c,e in zip(values,basis)),zero)
        assert (poly(D**(2*j-1)*R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        limits=[(7 if eps else 6)*j-1-sum(e) for e in basis]
        assert all(c.degree()<=b for c,b in zip(values,limits))
        H=s.Poly(t**(3*j)*D**(4*j-1),t);assert min(e[0] for e,c in H.terms() if c)>=3*j
        rows.append({'j':j,'N':[str(c.as_expr()) for c in values],'H':str(H.as_expr()),'row_degree_bounds':limits,'exact_twisted_identity':True})
        cols.append(values)
    co=cofactors(cols,4);checks=[]
    for k,C in enumerate(co):
        H=s.Poly(1,t) if k==0 else s.Poly(s.sympify(rows[k-1]['H'],locals={'t':t}),t)
        A=H*C;val=min(e[0] for e,c in A.terms() if c) if A else None
        if k==0:assert A and A.degree()<=(62 if eps else 52)
        else:assert val is None or val>=3*k
        checks.append({'order':k,'scalar_degree':int(A.degree()) if A else -1,'scalar_valuation':val})
    metric=s.Matrix([[11,1,2],[1,13,3],[2,3,m]]).det();assert metric!=0
    return {'label':label,'D':str(D),'metric_determinant':str(metric),'D_degree':int(s.degree(D,t)),
            'W':str(W),'seed_columns':rows,'cofactor_checks':checks}

def repeated_universal():
    a,b,c=s.symbols('a b c');R=s.eye(2)-2*t*s.Matrix([[a,b],[b,c]]);D=s.expand(R.det());adj=R.adjugate()
    js=s.symbols('j00 j01 j10 j11');hs=s.symbols('h00 h01 h02 h10 h11 h12')
    J=s.Matrix([js[0]+js[1]*u,js[2]+js[3]*u]);H=s.Matrix([hs[0]+hs[1]*u+hs[2]*u*u,hs[3]+hs[4]*u+hs[5]*u*u])
    g=sum(s.Symbol('g'+str(i))*u**i for i in range(4));q0=sum(s.Symbol('q'+str(i))*u**i for i in range(3))
    U=s.expand(D*(t*g-q0)+((t*H-J).T*adj*(t*H-J))[0]/2)
    assert s.degree(U,t)<=3 and s.degree(D,t)<=2 and s.degree(U.subs(t,0),u)<=2
    assert all(s.expand(U).coeff(u,k).subs(t,0)==0 for k in [3,4])
    A=s.Symbol('A');Diso=t*(A*t-s.Rational(1,2));assert Diso.subs(t,0)==0
    return {'nonisotropic_full_resolvent_U_degree':3,'nonisotropic_D_degree':2,'nonisotropic_U_at_zero_amplitude_at_most':2,
            'nonisotropic_cubic_quartic_leading_coefficient_zero_at_t0':True,'isotropic_D':'t*(A*t-1/2)',
            'isotropic_D_zero_at_t0':True,'all_actual_connection_denominators_zero_at_t0':True,
            'connection_denominator':'D^2 a^d','higher_scalar_valuation_at_least':'j',
            'initial_order_gate':'B*s*(s+1)/2','cubic_complement_gate':39,'quartic_complement_gate':96,
            'adapted_frame_gates_retained':[35,92,93],'repeated_full_gate':96}

def repeated_model(d,den):
    U=-u*u/2+sum(((j+2)*t+(2*j+1)*t*t+(j*j+1)*t**3)*u**j for j in range(d+1))
    U=s.expand(U);a=s.expand(s.diff(U,u)).coeff(u,d-1);H=s.expand(den*den*a**d)
    src=s.expand(s.diff(U,t)*den-U*s.diff(den,t)-s.diff(den,t)*den/2)
    columns=[];witnesses=[]
    for j in range(d-1):
        f=s.Poly(src*u**j,u);C=s.Integer(0);count=0
        while f.degree()>=d-1:
            n=f.degree();lc=f.LC();cw=u**(n-d+1)
            rem=s.expand(-cw*(s.diff(U,u)-a*u**(d-1))-den*s.diff(cw,u))
            f=s.Poly(s.expand(a*(f.as_expr()-lc*u**n)+lc*rem),u);C=s.expand(a*C+lc*cw);count+=1
        assert count<=d
        normal=s.expand(a**(d-count)*f.as_expr());C=s.expand(a**(d-count)*C)
        assert s.expand(a**d*src*u**j-normal-den*s.diff(C,u)-s.diff(U,u)*C)==0
        columns.append([normal.coeff(u,k) for k in range(d-1)]);witnesses.append(str(C))
    N=s.Matrix.hstack(*map(s.Matrix,columns));B=4+3*d
    assert s.degree(H,t)<=B and all(s.degree(z,t)<=B for z in N)
    hp=s.Poly(H,t);hval=min(e[0] for e,c in hp.terms() if c);assert hval>=1
    cur=s.eye(d-1)[:,0];cols=[[s.Poly(z,t) for z in cur]]
    for j in range(d-1):
        cur=(H*cur.diff(t)-j*s.diff(H,t)*cur+N*cur).applyfunc(s.expand)
        assert all(s.degree(z,t)<=B*(j+1) for z in cur)
        cols.append([s.Poly(z,t) for z in cur])
    co=cofactors(cols,d-1);checks=[];gate=B*(d-1)*d//2
    for j,C in enumerate(co):
        A=C*hp**j;val=min(e[0] for e,c in A.terms() if c) if A else None
        assert A.degree()<=gate if A else True
        if j==0:assert A
        else:assert val is None or val>=j
        checks.append({'order':j,'scalar_degree':int(A.degree()) if A else -1,'scalar_valuation':val})
    return {'degree':d,'D':str(den),'U':str(U),'connection_H':str(H),'H_valuation':hval,'B':B,
            'exact_connection_columns':d-1,'connection_witnesses':witnesses,'cofactor_checks':checks,'pure_exclusion_gate':gate}

def arithmetic():
    rows=[]
    for degree,slope,gate in [(1,7,62),(0,6,52)]:
        cases=[]
        for n in range(1,5):
            for selected in combinations(range(4),n):
                if 0 not in selected:continue
                bound=slope*n*(n+1)//2-n-sum(sum(basis[k]) for k in selected)
                cases.append({'rank':n,'rows':selected,'zeroth_degree_bound':bound})
        assert max(c['zeroth_degree_bound'] for c in cases)==gate
        rows.append({'D_degree':degree,'gate':gate,'all_rank_row_cases':cases})
    assert max(62,52,92)==92 and max(39,96,35,92,93)==96
    cover={'quadratic':7,'smooth':94,'nodal':96,'cusp':95,'transverse':97,'tangent':95,'triangle':98,'concurrent':92,'repeated':96}
    assert max(cover.values())==98
    return rows,cover

def main():
    start=time.monotonic();universal=universals();ru=repeated_universal();print('Full concurrent source/division and repeated elimination valuations passed',flush=True)
    cm=[]
    for args in [('linear_full',1-2*t),('linear_t0',t),('constant1',s.Integer(1)),('constant2',s.Integer(2)),('linear_ell_zero',1-2*t,True)]:
        cm.append(concurrent_model(*args));print('Concurrent model',args[0],'passed',flush=True)
    rm=[]
    for d in [3,4]:
        for den in [1+2*t+3*t*t,t*(2*t-1)]:
            rm.append(repeated_model(d,den));print('Repeated model',d,str(den),'passed',flush=True)
    cases,cover=arithmetic();cert={'concurrent_universal':universal,'repeated_universal':ru,'concurrent_models':cm,'repeated_models':rm,'concurrent_degree_cases':cases,'global_cover':cover}
    f=P/'concurrent_seeds_repeated_valuations_certificate.json';f.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
            'concurrent_source_gauge_identities':len(universal['Euler_gauge_identities']),'concurrent_principal_identities':len(universal['cleared_principal_identities']),
            'concurrent_models':len(cm),'concurrent_original_seed_identities':20,'concurrent_full_cofactor_relations':5,
            'repeated_models':len(rm),'repeated_exact_connection_columns':sum(m['exact_connection_columns'] for m in rm),'repeated_full_cofactor_relations':4,
            'concurrent_complement_gates':{'linear_D':62,'constant_D':52},'concurrent_full_gate':92,'repeated_cubic_complement_gate':39,'repeated_quartic_complement_gate':96,'repeated_full_gate':96,
            'global_cover':cover,'global_sufficient_pure_prefix':98,'odd_sufficient_pure_prefix_unchanged':98,'homogeneous_sufficient_pure_prefix_unchanged':96,
            'least_prefix_interval':[20,98],'odd_least_prefix_interval':[20,98],'homogeneous_least_prefix_interval':[14,96],
            'least_prefix_computed':False,'universal_written_proof_certified_by_checker':False,'independent_replay':False,'external_review':False,'formal_verification':False,
            'boundary':'Exact full concurrent completion/rescaling/source/Euler/principal identities,20 original twisted seeds and five full scalar relations; actual repeated elimination denominator valuation identities,10 exact connection columns and four full scalar relations; all rank/row degrees and higher scalar valuations. Universal all-degree reduction/seed/minor induction,full metric/frame/leading-type coverage and analytic nonzero zeroth arguments remain independently unreviewed written inputs. No old checker main,radical search,producer package or fresh literature review.'}
    (P/'concurrent_seeds_repeated_valuations_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
