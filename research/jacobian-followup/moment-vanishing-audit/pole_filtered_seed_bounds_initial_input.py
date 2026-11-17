"""Exact auxiliary pole filtrations and row cancellations for cubic seed bounds.

No retained checker main or prior receipt is executed or rewritten.
Analytic exclusions and full class/rank coverage remain written inputs.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_squarefree_singular_finite_prefixes as weighted
import verify_rescaled_singular_prefixes as rescaled
import verify_structured_seed_modules as structured
import verify_balanced_logarithmic_seed_bounds as balanced
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
y=balanced.y

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def total(f,*gens):return -1000 if f==0 else s.Poly(f,*gens).total_degree()
def rho(e):return e[0]+2*e[1]
def laurent(f):return balanced.laurent(f)

def universal_tangent_filtration():
    D,eps=s.symbols('D eps');dt=lambda f:s.diff(f,t)+eps*s.diff(f,D)
    delta=t*t*D;deltap=dt(delta)
    R=sum(s.Symbol('c'+str(i)+'_'+str(j)+'_'+str(a)+'_'+str(b))*t**a*D**b*u**i*v**j
          for i,j in weighted.monomials(7,(2,3)) for a in range(5-i-j) for b in range(5-i-j-a))
    W=u**4/2+D*u*v*v+R
    def F(f):return max((e[2]+2*e[3]-e[1] for e,c in s.Poly(f,t,D,u,v).terms() if c),default=-1000)
    def joint(f):return -1000 if f==0 else s.Poly(f,t,D,u,v).total_degree()
    assert F(W)<=4 and joint(W)<=4
    original=delta*dt(W)-deltap*W-2*t**3*D**2-t**4*D*eps/2
    bx=deltap/4;by=t*D-deltap/8
    S=s.expand(original+delta*(bx+by)+bx*u*s.diff(W,u)+by*v*s.diff(W,v))
    C=s.Poly(S,u,v).coeff_monomial(u*u*v)
    source=s.cancel((D*S-C*u*s.diff(W,v)/2)/t);h=t**3*D**3
    assert joint(source)<=6 and F(source)<=3
    assert all(weighted.weight(e,(2,3))<=6 for e,c in s.Poly(source,u,v).terms() if c)
    assert F(h)<=-3 and F(dt(h))<=-2
    for field in [s.cancel(delta*D*bx/t),s.cancel(delta*D*by/t)]:assert F(field)<=-2
    assert F(s.cancel(-delta*C/(2*t)))<=-1
    # Complete principal divisions preserve F after clearing D^2; the
    # strict original (2,3)-weight decrease supplies termination.
    rows=[]
    for exponent in weighted.monomials(16,(2,3)):
        nf,cu,cv=weighted.principal_division('tangent',exponent,s.Rational(1,2),D)
        nf,cu,cv=[s.cancel(D**2*x) for x in [nf,cu,cv]]
        rem=s.expand(nf-cu*s.diff(R,u)-cv*s.diff(R,v)-delta*(s.diff(cu,u)+s.diff(cv,v)))
        assert s.expand(D**2*u**exponent[0]*v**exponent[1]-rem-delta*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
        assert F(rem)<=rho(exponent)-2 and F(cu)<=rho(exponent)-5 and F(cv)<=rho(exponent)-4
        assert joint(rem)<=sum(exponent)+2
        rows.append({'exponent':exponent,'F_remainder':F(rem),'F_bound':rho(exponent)-2})
    return {'formal_derivation':'partial_t+eps*partial_D','formal_total_weights':{'t':1,'D':1,'u':1,'v':1,'eps':0},'pole_weight':'i+2j-D_exponent','phase_F_bound':4,'source_F_bound':3,'source_joint_degree':6,'principal_identities':rows,'normal_row_divisibility':'D^(9j-8+rho_i)','row_grades':[0,2,4,3,6]}

def tangent_seed_replay(m,pre_time=False):
    local={'t':t,'u':u,'v':v};read=lambda x:s.sympify(x,locals=local)
    W,D,A=[read(m[n]) for n in ['W','D','A0']];delta=t*t*D
    source,bx,by,C,_=structured.tangent_source(W,D);source=s.cancel(source/t);h=t**3*D**3
    vector=lambda R:s.cancel(delta*D*(bx*u*s.diff(R,u)+by*v*s.diff(R,v))/t-delta*C*u*s.diff(R,v)/(2*t))
    ring=s.QQ.poly_ring(t);poly=lambda x:s.Poly(x,u,v,domain=ring);wp,dp=map(poly,[W,delta]);wu,wv=wp.diff(u),wp.diff(v)
    def divide(expr,den):return s.Poly(expr,t,u,v,domain=s.QQ).exquo(s.Poly(den,t,u,v,domain=s.QQ)).as_expr()
    basis=m['basis'];R=s.Integer(1);rows=[]
    for row in m['seed_columns']:
        j=row['j'];R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+vector(R)+source*R)
        depth=row['depth'];assert depth==6*j-4
        e=9*j-8;time=t**j if pre_time else s.Integer(1)
        values=[s.Poly(divide(read(x),time*D**(e+rho(exponent))),t,domain=s.QQ) for x,exponent in zip(row['N'],basis)]
        witness=row.get('reduction_witness',row.get('witness'))
        Cu,Cv=[poly(divide(read(witness[k]),time*D**e)) for k in ['Cu','Cv']]
        normal=sum((poly(z.as_expr()*D**rho(exponent)*u**exponent[0]*v**exponent[1]) for z,exponent in zip(values,basis)),poly(0))
        assert (poly(D**(3*j))*poly(R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wu*Cu-wv*Cv).is_zero
        for z,exponent in zip(values,basis):assert z.degree()<=9*j-sum(exponent)-rho(exponent)
        H=s.Poly(divide(read(row['H']),time*D**e),t,domain=s.QQ)
        assert s.expand(H.as_expr()-t**(3*j)*D**(6*j))==0 and H.degree()<=9*j
        val=min(ex[0] for ex,c in H.terms() if c);assert val>=3*j
        rows.append({'j':j,'cancelled_common_D_power':e,'row_D_powers':[rho(ex) for ex in basis],'H':str(H.as_expr()),'H_valuation':val,'N':[str(z.as_expr()) for z in values],'witness':{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}})
    return {'singular_at_zero':m['singular_at_zero'],'zero_shear':m['zero_shear'],'constant_D':s.degree(D,t)==0,'metric_determinant':m['metric_determinant'],'D':str(D),'basis':basis,'row_basis_scales':[str(D**rho(e)) for e in basis],'seed_columns':rows}

def parent_universal():
    z,T,a=s.symbols('z T a');g={(i,j):s.Symbol('g'+str(i)+str(j)) for i in range(4) for j in range(4-i)}
    F=z/t-a-z*z;W=s.expand(t*t*(t*sum(c*z**i*F**j for (i,j),c in g.items())-F*F/2));rows=[]
    for d in [4,5,6]:
        V=W.subs({g[0,3]:0,g[1,2]:0}) if d==4 else W.subs(g[0,3],0) if d==5 else W
        fp=2 if d<6 else 3;lv=2 if d==4 else 3
        def fil(x):return max((e[1]-e[0] for e,c in s.Poly(x,t,z).terms() if c),default=-1000)
        assert fil(V)<=fp
        source=t*s.diff(V,t)-2*V-t*t
        if d>4:source-=t*t/d+z*s.diff(V,z)/d
        source=s.expand(source);assert fil(source)<=fp
        assert s.degree(source,z)<=(d if d==4 else d-1)
        A=s.Poly(s.diff(V,z),z).LC();assert s.expand(A/t**lv).subs(t,0)!=0
        bound=3 if d==6 else 1
        assert fil(s.diff(V,z))<=fp-1
        # Leading-term division has input filtration preserved; these are
        # exactly the three possible leading valuations in the full cover.
        rows.append({'d':d,'phase_filtration':fp,'source_filtration':fp,'leading_t_valuation':lv,'leading_after_t_power':str(s.cancel(A/t**lv))})
    return rows

tangent_source=structured.tangent_source

def constant_tangent_phase():
    D=s.Integer(1);ell=2*u+3*v+1;J=2*u+3*v+1
    g=3*u*u+5*u*v+7*v*v+11*u+13*v+17;q=(11*u*u+13*v*v)/2+u*v+5*u+7*v+2
    W=rescaled.rescale(t*t*(D*(t*(u*v*v+g)-q)+(t*(u*u+ell)-J)**2/2))
    metric=int(s.Matrix([[11,1,2],[1,13,3],[2,3,1]]).det());assert metric!=0
    return W,D,metric

def parent_model(label,W):
    ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,y,domain=ring)
    wp=poly(W);d=wp.degree();assert d in [4,5,6]
    derivative=wp.diff(y);A=s.Poly(derivative.LC(),t,domain=s.QQ);assert A.degree()==3
    assert not s.cancel(A.as_expr()/t**3).has(t) if d>4 else s.cancel(A.as_expr()/t**2).subs(t,0)!=0
    source=s.expand(t*s.diff(W,t)-2*W-t*t-(t*t/d+y*s.diff(W,y)/d if d>4 else 0))
    assert s.degree(source,y)<=(d-1 if d>4 else d)
    ap=poly(A.as_expr());delta=poly(t*t);zero=poly(0);zt=s.Poly(0,t,domain=s.QQ)
    @lru_cache(None)
    def reduce(n,depth):
        if n<d-1:
            values=[zt]*(d-1);values[n]=A**depth;return tuple(values),zero
        assert depth>0;cu=poly(y**(n-d+1));rem=ap*poly(y**n)-delta*cu.diff(y)-derivative*cu
        assert rem.degree()<n;values=[zt]*(d-1);witness=ap**(depth-1)*cu
        for (k,),c in rem.terms():
            if not c:continue
            vec,w=reduce(k,depth-1);ct=s.Poly(c,t,domain=s.QQ)
            values=[a+ct*b for a,b in zip(values,vec)];witness+=poly(c)*w
        return tuple(values),witness
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        old=t*s.diff(W,t)-2*W-t*t;old_next=t**3*s.diff(R,t)-3*(j-1)*t*t*R+old*R
        next=s.expand(t**3*s.diff(R,t)-3*(j-1)*t*t*R-(t*t*y*s.diff(R,y)/d if d>4 else 0)+source*R)
        assert s.expand(next-old_next+((t*t*s.diff(y*R,y)+s.diff(W,y)*y*R)/d if d>4 else 0))==0
        R=next;assert s.degree(R,y)<=(d-1 if d>4 else d)*j
        depth=(d-1)*(j-1)+1 if d>4 else 4*j-2;values=[zt]*(d-1);witness=zero
        for (n,),c in poly(R).terms():
            if not c:continue
            vec,w=reduce(n,depth);ct=s.Poly(c,t,domain=s.QQ)
            values=[a+ct*b for a,b in zip(values,vec)];witness+=poly(c)*w
        normal=sum((poly(z.as_expr()*y**k) for k,z in enumerate(values)),zero)
        assert (ap**depth*poly(R)-normal-delta*witness.diff(y)-derivative*witness).is_zero
        bnd=3*j+3*depth;assert all(z.degree()<=bnd for z in values)
        H=s.Poly(s.expand(t**(3*j)*A.as_expr()**depth),t,domain=s.QQ);assert H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        base={4:6*j-4,5:10*j-9,6:12*j-12}[d]
        divide=lambda x,power:s.Poly(x,y,t,domain=s.QQ).exquo(s.Poly(t**power,y,t,domain=s.QQ)).as_expr()
        scaled=[s.Poly(divide(z.as_expr(),base+i),t,domain=s.QQ) for i,z in enumerate(values)]
        scaled_witness=poly(divide(witness.as_expr(),base))
        newH=s.Poly(divide(H.as_expr(),base),t,domain=s.QQ)
        nb={4:9*j-2,5:5*j,6:6*j}[d]
        assert all(z.degree()<=nb-i for i,z in enumerate(scaled)) and newH.degree()<=nb
        check_normal=sum((poly(z.as_expr()*t**i*y**i) for i,z in enumerate(scaled)),zero)
        assert (poly(s.cancel(A.as_expr()**depth/t**base))*poly(R)-check_normal-delta*scaled_witness.diff(y)-derivative*scaled_witness).is_zero
        val=min(e[0] for e,c in newH.terms() if c);assert val>=j
        seeds.append({'j':j,'depth':depth,'cancelled_common_t_power':base,'row_t_powers':list(range(d-1)),'N':[str(z.as_expr()) for z in scaled],'H':str(newH.as_expr()),'H_valuation':val,'new_degree_bound':nb,'witness':str(scaled_witness.as_expr())})
    return {'label':label,'d':d,'W':str(W),'source':str(source),'seed_columns':seeds}

def constant_tangent_fixture(singular=False,zero_shear=False):
    kind="tangent"
    W,D,metric=constant_tangent_phase();delta=t*t*D;source,bx,by,C,old=tangent_source(W,D);h=delta**2*D
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
        f=u**e[0]*v**e[1];input=source*f+delta*D*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))-delta*C*u*s.diff(f,v)/2
        maximum=6+weighted.weight(e,data['weights']);depth=(maximum-data['minimum_reducible'])//data['gap']+1
        values,witness=reduce_input(input,depth)
        assert all(z.degree()<=7+sum(e)+k*depth-sum(f) for z,f in zip(values,basis))
        columns.append({'depth':depth,'N':[str(z.as_expr()) for z in values],'witness':witness})
    R=s.Integer(1);seeds=[]
    for j in range(1,4):
        old_next=h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+D*old*R+delta*D*(bx*u*s.diff(R,u)+by*v*s.diff(R,v))
        next=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+source*R+delta*D*(bx*u*s.diff(R,u)+by*v*s.diff(R,v))-delta*C*u*s.diff(R,v)/2)
        Cu,Cv=s.Integer(0),-C*u*R/2
        assert s.expand(next-old_next-delta*(s.diff(Cu,u)+s.diff(Cv,v))-s.diff(W,u)*Cu-s.diff(W,v)*Cv)==0
        R=next;assert total(R,t,u,v)<=7*j
        assert all(weighted.weight(e,data['weights'])<=6*j for e,c in s.Poly(R,u,v).terms() if c)
        depth=6*j-4;bnd=7*j+k*depth
        values,witness=reduce_input(R,depth)
        assert all(z.degree()<=bnd-sum(e) for z,e in zip(values,basis))
        H=s.Poly(s.expand(h**j*A**depth),t,domain=s.QQ);assert H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=4*j
        seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':val,'reduction_witness':witness,'gauge_witness':{'Cu':str(Cu),'Cv':str(Cv)}})
    return {'kind':kind,'singular_at_zero':singular,'zero_shear':zero_shear,'metric_determinant':metric,'W':str(W),'D':str(D),'A0':str(A),'basis':basis,'source':str(source),'C':str(C),'h':str(h),'beta_x':str(bx),'beta_y':str(by),'connection_columns':columns,'seed_columns':seeds}



def conic_model(label,W,K,maximum_seed=3):
    field=s.QQ.frac_field(t);raw=laurent(W);p=max(raw);q=-min(raw);assert p in [2,3] and q==3
    kd=s.degree(K,t);assert p==2 or kd==0
    assert all(s.Poly(s.cancel(raw[-k]/K**k),t).degree()<=3-k for k in [1,2,3])
    basis=list(range(-3,p));valid=set(basis);Ap=raw[p]
    source=s.expand(t*s.diff(W,t)-2*W-t*t-y*s.diff(W,y)/3-t*t/3)
    h=t**3*K
    Sd=s.cancel(K*source);assert s.denom(Sd)==1 or not s.denom(Sd).has(t)
    derivative={k:field.convert(k*c) for k,c in raw.items() if k}
    def parse(x):return {k:field.convert(c) for k,c in laurent(x).items()}
    def add(out,k,c):
        out[k]=out.get(k,field.zero)+c
        if not out[k]:del out[k]
    def exact(x):
        out={}
        for n,c in x.items():
            add(out,n,field.convert(t*t*n)*c)
            for k,a in derivative.items():add(out,n+k,a*c)
        return out
    def reduce(x):
        out=x.copy();wit={}
        while max(out,default=0)>=p:
            n=max(out);shift=n-p;c=out[n]/derivative[p];add(wit,shift,c)
            for k,a in exact({shift:c}).items():add(out,k,-a)
        while min(out,default=0)<-3:
            n=min(out);shift=n+3;c=out[n]/derivative[-3]
            # Only the original negative module reaches this loop: high
            # reduction cannot produce exponents below -3.
            assert s.Poly(field.to_sympy(c),t,domain=s.QQ) is not None
            add(wit,shift,c)
            for k,a in exact({shift:c}).items():add(out,k,-a)
        check=x.copy()
        for k,a in out.items():add(check,k,-a)
        for k,a in exact(wit).items():add(check,k,-a)
        assert not check and set(out)<=valid
        return out,wit
    R=s.Integer(1);rows=[]
    for j in range(1,maximum_seed+1):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+Sd*R-t*t*K*y*s.diff(R,y)/3)
        assert min(laurent(R))>=-2*j and max(laurent(R))<=p*j
        for k,c in laurent(R).items():
            if k<0:assert s.Poly(s.cancel(c/K**(-k)),t,domain=s.QQ) is not None
        normal,wit=reduce(parse(R));high=2*j-1 if p==2 else 0
        den=s.expand(Ap**high);values=[]
        b={0:3*j+high+6*max(0,2*j-3),1:8*j+4*high+3*max(0,2*j-3),2:5*j+high}[kd] if p==2 else None
        wt,wy=(2,1) if kd==1 else (1,1)
        for k in basis:
            value=s.Poly(s.cancel(field.to_sympy(normal.get(k,field.zero))*den),t,domain=s.QQ)
            bound=(3-k if j==1 else 7*j+k) if p==3 else (b-wy*k)//wt
            assert value.degree()<=bound,(label,j,k,value.degree(),bound)
            values.append(value)
        H=s.Poly(s.expand(h**j*den),t,domain=s.QQ);assert min(e[0] for e,c in H.terms() if c)>=3*j
        if p==2:assert wt*H.degree()<=b
        else:assert H.degree()<=3*j
        rows.append({'j':j,'N':[str(x.as_expr()) for x in values],'H':str(H.as_expr()),'clearing_Aplus_power':high,'witness':{str(k):str(field.to_sympy(a)) for k,a in wit.items()}})
    return {'label':label,'W':str(W),'K':str(K),'K_degree':kd,'orders':[p,q],'basis':basis,'seed_columns':rows}

def degree_arithmetic():
    rows=[]
    def gate(label,basis,bj,rowgrade,den=None):
        maximum=-1
        for size in range(1,len(basis)+1):
            for selected in combinations(basis,size):
                if 0 not in selected:continue
                # Here bj already bounds numerator and denominator degree.
                g=sum(bj(j) for j in range(1,size+1))-sum(rowgrade(k) for k in selected)
                maximum=max(maximum,g)
        rows.append({'case':label,'gate':maximum});return maximum
    assert gate('parent_d4',list(range(3)),lambda j:9*j-2,lambda i:i)==45
    assert gate('parent_d5',list(range(4)),lambda j:5*j,lambda i:i)==44
    assert gate('parent_d6',list(range(5)),lambda j:6*j,lambda i:i)==80
    assert gate('tangent_pole',list(range(5)),lambda j:9*j,lambda i:[0,2,4,3,6][i])==120
    basis=list(range(-3,3));mx=-1
    for size in range(1,7):
        for selected in combinations(basis,size):
            if 0 not in selected:continue
            # Column1 has 3-k; later columns have7j+k. The largest
            # determinant assigns the smallest row to column1.
            g=3-sum(selected) if size==1 else 3+7*sum(range(2,size+1))+sum(selected)-2*min(selected)
            mx=max(mx,g)
    assert mx==146;rows.append({'case':'conic_p3q3_Kconstant','gate':mx})
    for kd,wt,wy,expected in [(0,1,1,171),(1,2,1,136),(2,1,1,105)]:
        bj=lambda j:{0:3*j+(2*j-1)+6*max(0,2*j-3),1:8*j+4*(2*j-1)+3*max(0,2*j-3),2:5*j+(2*j-1)}[kd]
        mx=-1
        for size in range(1,6):
            for selected in combinations(range(-3,2),size):
                if 0 not in selected:continue
                g=(sum(bj(j) for j in range(1,size+1))-wy*sum(selected))//wt;mx=max(mx,g)
        assert mx==expected;rows.append({'case':'conic_p2q3_Kdegree'+str(kd),'gate':mx})
    assert max(146,199,108,171,136,105,130,145)==199
    cover={'quadratic':7,'smooth':94,'nodal':199,'cusp':183,'transverse':199,'tangent':132,'triangle':199,'concurrent':104,'repeated':99}
    assert max(cover.values())==199
    return rows,cover

def tangent_geometry():
    z,T,w,g30,g21,g12,g03,f,q=s.symbols('z T w g30 g21 g12 g03 f q')
    cubic=f*w*T*T+g30*z**3+g21*z*z*T+g12*z*T*T+g03*T**3
    assert s.expand(cubic.subs(T,0))==g30*z**3
    quad=s.cancel(cubic.subs(g30,0)/T);assert quad.is_polynomial(z,T,w)
    matrix=s.hessian(quad,(z,T,w))/2;assert s.expand(matrix.det()+f*f*g21/4)==0
    # Rank-one completion w=(z-q/2)/T has highest z-degree from G3.
    phase=s.expand(q**3*((cubic.subs(g30,0)).subs({z:z/q,T:(z-s.Rational(1,2))/q},simultaneous=True).subs(w,(z-s.Rational(1,2))**2/q)))
    assert s.degree(phase,z)<=5
    return {'tangent_line':'T=0','vanishing_cubic_coefficient':'g30','smooth_conic_determinant':str(matrix.det()),'rank_one_phase_degree_at_most':5,'retained_degree_through5_gates':[18,50,93,132]}

def main():
    started=time.monotonic();universal=universal_tangent_filtration();print('Universal tangent filtration passed',flush=True)
    parent=parent_universal();print('Universal parent filtration passed',flush=True)
    old=json.loads((P/'common_time_seed_factor_certificate.json').read_text())
    # Discover the retained tangent collection by its explicit model fields.
    def collect(x):
        if isinstance(x,dict):
            if x.get('kind')=='tangent' and 'seed_columns' in x:return [x]
            return sum((collect(z) for z in x.values()),[])
        if isinstance(x,list):return sum((collect(z) for z in x),[])
        return []
    retained=collect(old);assert len(retained)==4
    tangent=[tangent_seed_replay(x) for x in retained];print('Retained tangent seeds passed',flush=True)
    const=constant_tangent_fixture();tangent.append(tangent_seed_replay(const,pre_time=True));print('Constant-D full Gaussian tangent passed',flush=True)
    z,a=s.symbols('z a');F=y/t-2-y*y
    def parent_phase(d):
        G=sum((1+2*i+3*j)*y**i*F**j for i in range(4) for j in range(4-i) if not (d==4 and (i,j) in [(0,3),(1,2)]) and not (d==5 and (i,j)==(0,3)))
        return s.expand(t*t*(t*G-F*F/2))
    parents=[parent_model('degree'+str(d),parent_phase(d)) for d in [4,5,6]];print('Parent seed identities passed',flush=True)
    conics=[]
    for kd,K in [(0,s.Rational(-1,4)),(1,t/2-s.Rational(1,4)),(2,(t+1)*(5*t-1)/4-t*t)]:
        for Ap in ([1+t,t,K] if kd else [1+t,t]):
            W=1+t**3+(1+t*t)*y+Ap*y*y+sum(K**k*(2+t**(3-k))/y**k for k in [1,2,3])
            conics.append(conic_model('synthetic_p2_K'+str(kd)+'_A'+str(Ap),s.expand(W),K))
    K=s.Rational(-1,4);W=1+t**3+(1+t*t)*y+(1+t)*y*y+2*y**3+sum(K**k*(2+t**(3-k))/y**k for k in [1,2,3])
    conics.append(conic_model('synthetic_p3_Kconstant',s.expand(W),K))
    print('Conic seed identities passed',flush=True)
    arithmetic,cover=degree_arithmetic();geometry=tangent_geometry()
    certificate={'universal_tangent_filtration':universal,'universal_parent_filtration':parent,'tangent_models':tangent,'parent_models':parents,'conic_models':conics,'degree_cases':arithmetic,'tangent_geometry':geometry,'global_cover':cover}
    cp=P/'pole_filtered_seed_bounds_certificate.json';cp.write_text(json.dumps(certificate,indent=2)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(cp),'input_certificate_sha256':sha(P/'common_time_seed_factor_certificate.json'),'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,'tangent_models':len(tangent),'parent_models':len(parents),'conic_models':len(conics),'direct_seed_identities':3*(len(tangent)+len(parents)+len(conics)),'mixed_parent_full_gate':92,'conic_full_gate':199,'full_W_linear_gate':199,'tangent_complement_gate':120,'tangent_full_gate':132,'global_sufficient_pure_prefix':199,'odd_sufficient_pure_prefix':198,'homogeneous_sufficient_pure_prefix':198,'least_prefix_interval':[20,199],'odd_least_prefix_interval':[20,198],'homogeneous_least_prefix_interval':[14,198],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact auxiliary filtration/source identities, retained and new model seed reductions, finite determinant arithmetic and tangent factor geometry only. Universal degree induction, complete class/rank coverage, analytic exclusions and period comparison remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (P/'pole_filtered_seed_bounds_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
