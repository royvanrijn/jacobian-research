#!/usr/bin/env python3
"""Exact weighted connection certificates for104/720/935 class gates."""
from pathlib import Path
from functools import lru_cache
import hashlib,json
import sympy as S

u,v,t,T=S.symbols('u v t T')
CASES={
 'concurrent':{'weights':(1,1),'top_weight':3,'lower_weight':2,'basis':[(0,0),(1,0),(0,1),(1,1)],'max_source':5,'minimum_reducible':2,'gap':1,'steps':4,'q':2,'B':10,'rank':4,'gate':104},
 'cusp':{'weights':(3,4),'top_weight':12,'lower_weight':10,'basis':[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1)],'max_source':22,'minimum_reducible':8,'gap':2,'steps':8,'q':4,'B':34,'rank':6,'gate':720},
 'tangent':{'weights':(2,3),'top_weight':8,'lower_weight':7,'basis':[(0,0),(1,0),(2,0),(0,1),(0,2)],'max_source':14,'minimum_reducible':5,'gap':1,'steps':10,'q':6,'B':62,'rank':5,'gate':935}}


def weight(exponent,weights):return sum(a*b for a,b in zip(exponent,weights))


def monomials(limit,weights):
    return [(i,j) for i in range(limit//weights[0]+1) for j in range(limit//weights[1]+1) if weight((i,j),weights)<=limit]


def principal_data(kind,D):
    a=t*t/2;b=t*D
    if kind=='concurrent':return b*(u**3+v**3),b
    if kind=='cusp':return a*u**4+b*v**3,a*b
    return a*u**4+b*u*v*v,a*b*b


def principal_division(kind,exponent,a,b):
    i,j=exponent;m=u**i*v**j
    if exponent in CASES[kind]['basis']:return m,S.Integer(0),S.Integer(0)
    if kind=='concurrent':
        if i>=2:return S.Integer(0),u**(i-2)*v**j/(3*b),S.Integer(0)
        return S.Integer(0),S.Integer(0),u**i*v**(j-2)/(3*b)
    if kind=='cusp':
        if i>=3:return S.Integer(0),u**(i-3)*v**j/(4*a),S.Integer(0)
        return S.Integer(0),S.Integer(0),u**i*v**(j-2)/(3*b)
    if i and j:return S.Integer(0),S.Integer(0),u**(i-1)*v**(j-1)/(2*b)
    if i==3:return -b*v*v/(4*a),1/(4*a),S.Integer(0)
    if i>=4:return S.Integer(0),u**(i-3)/(4*a),-u**(i-4)*v/(8*a)
    assert i==0 and j>=3
    return S.Integer(0),v**(j-2)/b,-2*a*u*u*v**(j-3)/(b*b)


def check_universal_divisions():
    a,b=S.symbols('a b',nonzero=True);rows=[]
    for kind,data in CASES.items():
        top=b*(u**3+v**3) if kind=='concurrent' else a*u**4+(b*v**3 if kind=='cusp' else b*u*v*v)
        gu,gv=S.diff(top,u),S.diff(top,v)
        count=0
        for exponent in monomials(data['max_source'],data['weights']):
            nf,hu,hv=principal_division(kind,exponent,a,b)
            assert S.expand(u**exponent[0]*v**exponent[1]-nf-hu*gu-hv*gv)==0
            assert all(e in data['basis'] for e,c in S.Poly(nf,u,v).terms() if c)
            count+=1
        rows.append({'kind':kind,'symbolic_divisions':count})
    return rows


def model_connection(kind,D):
    data=CASES[kind];ring=S.QQ.poly_ring(t)
    poly=lambda expr:S.Poly(expr,u,v,domain=ring)
    zero=poly(0);zero_t=S.Poly(0,t,domain=S.QQ)
    top,A=principal_data(kind,D)
    support=monomials(data['lower_weight'],data['weights'])
    R=sum((2+3*i+5*j)*(1+t+t*t)*u**i*v**j for i,j in support)
    U=S.expand(top+R)
    Up=poly(U);Uu=Up.diff(u);Uv=Up.diff(v)
    Ru,Rv=poly(R).diff(u),poly(R).diff(v)
    Ap=poly(S.expand(A));Dp=poly(D)
    At=S.Poly(S.expand(A),t,domain=S.QQ)
    a=t*t/2;b=t*D
    basis=data['basis'];index={e:i for i,e in enumerate(basis)}
    @lru_cache(None)
    def step(exponent):
        nf,hu,hv=principal_division(kind,exponent,a,b)
        # Clear the single common principal divisor before all recursion.
        nfp=poly(S.cancel(A*nf));hup=poly(S.cancel(A*hu));hvp=poly(S.cancel(A*hv))
        assert all(S.Poly(c,t).degree()<=data['q'] for e,c in nfp.terms() if c)
        assert all(S.Poly(c,t).degree()<=data['q']-2 for h in [hup,hvp] for e,c in h.terms() if c)
        remainder=nfp-hup*Ru-hvp*Rv-Dp*(hup.diff(u)+hvp.diff(v))
        assert all(e in index or weight(e,data['weights'])<=weight(exponent,data['weights'])-data['gap'] for e,c in remainder.terms() if c)
        assert all(S.Poly(c,t).degree()<=data['q'] for e,c in remainder.terms() if c)
        return remainder,hup,hvp
    @lru_cache(None)
    def reduce_monomial(exponent,remaining):
        if exponent in index:
            values=[zero_t]*len(basis);values[index[exponent]]=At**remaining
            return tuple(values),zero,zero
        assert remaining>0,(kind,exponent)
        remainder,hu,hv=step(exponent)
        values=[zero_t]*len(basis);Cu=zero;Cv=zero
        for e,c in remainder.terms():
            if not c:continue
            vector,cu,cv=reduce_monomial(e,remaining-1)
            coefficient=S.Poly(c,t,domain=S.QQ)
            values=[x+coefficient*y for x,y in zip(values,vector)]
            scalar=poly(c);Cu+=scalar*cu;Cv+=scalar*cv
        Cu+=Ap**(remaining-1)*hu;Cv+=Ap**(remaining-1)*hv
        return tuple(values),Cu,Cv
    source=poly(S.expand(S.diff(U,t)*D-U*S.diff(D,t)-S.diff(D,t)*D/2))
    H=S.expand(D*D*A**data['steps'])
    columns=[];certificates=[]
    for exponent in basis:
        f=source*poly(u**exponent[0]*v**exponent[1])
        values=[zero_t]*len(basis);Cu=zero;Cv=zero
        for e,c in f.terms():
            if not c:continue
            vector,cu,cv=reduce_monomial(e,data['steps'])
            coefficient=S.Poly(c,t,domain=S.QQ)
            values=[x+coefficient*y for x,y in zip(values,vector)]
            scalar=poly(c);Cu+=scalar*cu;Cv+=scalar*cv
        normal=sum((poly(value.as_expr()*u**e[0]*v**e[1]) for e,value in zip(basis,values)),zero)
        # Independently expand the full twisted integration-by-parts identity.
        defect=Ap**data['steps']*f-normal-Dp*(Cu.diff(u)+Cv.diff(v))-Uu*Cu-Uv*Cv
        assert defect.is_zero
        assert all(value.degree()<=data['B'] for value in values)
        columns.append([str(value.as_expr()) for value in values])
        certificates.append({'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())})
    assert S.Poly(H,t).degree()<=data['B']
    return {'kind':kind,'D':str(D),'D_at_zero':str(D.subs(t,0)),'U':str(U),'basis':[str(u**e[0]*v**e[1]) for e in basis],
            'H':str(H),'N_columns':columns,'integration_by_parts_certificates':certificates,
            'connection_identity':'A^r*source*basis_j-N_j=D*(partial_u Cu+partial_v Cv)+U_u*Cu+U_v*Cv',
            'H_degree':int(S.Poly(H,t).degree()),'N_max_degree':int(max(S.Poly(x,t).degree() for col in columns for x in col)),
            'columns_certified':len(basis),'B':data['B'],'rank':data['rank'],'prefix':data['gate']}


def main():
    p=Path(__file__).resolve().parent
    universal=check_universal_divisions()
    # Common completion with all lower coefficients left symbolic.
    alpha,m=S.symbols('alpha m');D=m-2*alpha*t
    ell=S.Symbol('ell_u')*u+S.Symbol('ell_v')*v+S.Symbol('ell_0')
    J=S.Symbol('j_u')*u+S.Symbol('j_v')*v+S.Symbol('j_0')
    g2=sum(S.Symbol('g'+str(i)+str(j))*u**i*v**j for i,j in monomials(2,(1,1)))
    q2=sum(S.Symbol('q'+str(i)+str(j))*u**i*v**j for i,j in monomials(2,(1,1)))
    for kind,data in CASES.items():
        H2=0 if kind=='concurrent' else u*u
        G3=u**3+v**3 if kind=='concurrent' else (v**3 if kind=='cusp' else u*v*v)
        L=t*(H2+ell)-J
        U=S.expand(D*(t*(G3+g2)-q2)+L*L/2)
        original=-D*T*T/2+L*T+t*(G3+g2)-q2
        assert S.expand(original.subs(T,T+L/D)-(-D*T*T/2+U/D)).simplify()==0
        top,A=principal_data(kind,D)
        remainder=S.Poly(S.expand(U-top),u,v)
        assert all(weight(e,data['weights'])<=data['lower_weight'] for e,c in remainder.terms() if c)
        assert S.Poly(U,t).degree()<=2
        assert data['steps']==(data['max_source']-data['minimum_reducible'])//data['gap']+1
        assert data['B']==2+data['q']*data['steps']
        assert data['gate']==data['B']*data['rank']*(data['rank']+1)//2+data['rank']
    models=[]
    for kind in CASES:
        for denominator in [1-2*t,t]:
            models.append(model_connection(kind,denominator))
            print('Certified',kind,str(denominator),flush=True)
    modelpath=p/'squarefree_singular_connection_models.json'
    modelpath.write_text(json.dumps({'completed':True,'models':models},indent=2)+'\n')
    retained={name:hashlib.sha256((p/name).read_bytes()).hexdigest() for name in ['concurrent_line_profiles_verification.json','cuspidal_affine_profiles_verification.json','affine_line_conic_profiles_verification.json','rank_one_wt2_verification.json','mixed_period_infinity_verification.json','repeated_line_finite_prefix_verification.json']}
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'model_certificate_sha256':hashlib.sha256(modelpath.read_bytes()).hexdigest(),'retained_receipts':retained,
            'universal_principal_divisions':universal,
            'class_bounds':{kind:{key:data[key] for key in ['weights','top_weight','lower_weight','max_source','minimum_reducible','gap','steps','q','B','rank','gate']} for kind,data in CASES.items()},
            'connection_models':[{key:row[key] for key in ['kind','D','D_at_zero','H_degree','N_max_degree','columns_certified','B','rank','prefix']} for row in models],
            'frame_gates':{'concurrent':92,'cusp':425,'tangent':425},
            'global_arbitrary_cubic_pure_gate':None,
            'remaining_unquantified_leading_types':['irreducible nodal','transverse line-smooth-conic','three independent lines'],
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Universal symbolic quadratic completion, principal weighted divisions and lower-weight/parameter-degree bounds; six exact rational connection models with all column integration-by-parts certificates, including D(0)=0;104/720/935 arithmetic. Universal filtered-basis exactness, complete class/field/frame choices, seed bridge and nonzero zeroth coefficients remain written mathematics with retained analytic dependencies. No sharpness, high-order moment expansion, numeric global pure-only gate, independent review or formal verification.'}
    (p/'squarefree_singular_finite_prefixes_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','global_arbitrary_cubic_pure_gate','remaining_unquantified_leading_types','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
