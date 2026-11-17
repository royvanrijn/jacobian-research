"""Auxiliary ordinary Gaussian-square D residues and general125 cofactor cover.
All-order seed/minor induction, complete coefficient/frame cover and analytic
nonzero zeroth arguments remain written. Prior checker mains are not run.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
import verify_ordinary_node_finite_prefixes as ordinary
import verify_euler_seed_bounds as euler
PTH=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()

def universals():
    D,eps,R=s.symbols('D eps R');a,d,Le,Lf,J,L0,N0,K=s.symbols('a d Le Lf J L0 N0 K')
    F=u*v+R;out=[]
    for kind in ['nodal','transverse','triangle']:
        G=(a*u**3 if kind!='triangle' else 0)+(d*v**3 if kind=='nodal' else 0)+Le*u*u+Lf*v*v+J*u*v+L0*u+N0*v+K
        W=D*G+F*F/2;delta=t*t*D;ids=[]
        for i,k in ordinary.monomials(8):
            if i+k==0:continue
            f=u**i*v**k
            if i:
                C=u**(i-1)*v**k;defect=F*f-(delta*s.diff(C,v)+s.diff(W,v)*C)+D*(s.diff(G,v)*C+t*t*s.diff(C,v))
            else:
                C=v**(k-1);defect=F*f-(delta*s.diff(C,u)+s.diff(W,u)*C)+D*(s.diff(G,u)*C+t*t*s.diff(C,u))
            assert s.expand(defect)==0;ids.append([i,k])
        S0=s.expand(-t*eps*(F*F/2-(u*s.diff(F*F/2,u)+v*s.diff(F*F/2,v))/4))
        assert s.expand(S0+t*eps*R*F/2)==0
        bs=ordinary.CASES[kind]['basis'];b=s.Matrix([-t*eps*R*R/2 if ex==(0,0) else -t*eps*R/2 if ex==(1,1) else 0 for ex in bs]);unit=s.zeros(1,len(bs));unit[0,0]=1
        A=b*unit;aug=s.Matrix.hstack(b,A);minors=0
        for rows in combinations(range(len(bs)),2):
            for cols in combinations(range(len(bs)+1),2):assert s.expand(aug.extract(rows,cols).det())==0;minors+=1
        out.append({'kind':kind,'exact_F_monomial_identities':ids,'source_residue':str(S0),'residue_vector':[str(x) for x in b],'residue_outer_product':'b*e_unit^T','augmented_two_by_two_minors':minors,'all_R_zero_specializations_allowed':True,'requires_retained_polynomial_D_normal_module':True})
    return out

def parameters(kind,mode):
    aa=1 if kind!='triangle' else 0;dd=2 if kind=='nodal' else 0
    if mode=='zero_shared':
        den=t;m=0;A=18 if kind=='transverse' else 0 if kind=='triangle' else 11;B=0 if kind!='nodal' else 13
    else:den=1-2*t if mode=='linear' else s.Integer(1);m=1;A=11;B=13
    ju,jv=2,3;b,c=3,5
    F=u*v-ju*jv+t*t-t
    x,y=u+jv,v+ju
    G=aa*x**3+dd*y**3+b*t*x*x+c*t*y*y+7*t*t*x+11*t*t*y+2*t**3-(A*x*x+2*x*y+B*y*y)/2-5*t*x-7*t*y-2*t*t
    W=s.expand(den*G+F*F/2);metric=s.Matrix([[A,1,ju],[1,B,jv],[ju,jv,m]]).det();assert metric!=0
    return W,den,metric

def row_parameters(kind,W,D):
    pv=s.Poly(W,u,v);Le=s.cancel(pv.coeff_monomial(u*u)/D);Lf=s.cancel(pv.coeff_monomial(v*v)/D)
    bs=ordinary.CASES[kind]['basis']
    if kind=='nodal':scales=[s.Integer(1)]*len(bs);grades=[sum(e) for e in bs];clear=lambda j:s.Integer(1)
    elif kind=='transverse':scales=[s.Integer(1)]*len(bs);grades=[sum(e) for e in bs];clear=lambda j:Lf**(2*j-2)
    else:
        scales=[Le**int(s.ceiling(i+s.Rational(k,2)))*Lf**int(s.ceiling(s.Rational(i,2)+k)) for i,k in bs]
        grades=[i+k+int(s.ceiling(i+s.Rational(k,2)))+int(s.ceiling(s.Rational(i,2)+k)) for i,k in bs]
        clear=lambda j:Le**(2*j)*Lf**(2*j)
    return Le,Lf,scales,grades,clear

def det(cols,rank):
    ring=s.QQ.poly_ring(t);mat=DomainMatrix.from_list_sympy(rank,rank,[[x.as_expr() for x in row] for row in zip(*cols)]).convert_to(ring)
    return s.Poly(ring.to_sympy(mat.det()),t,domain=s.QQ)

def model(kind,mode):
    W,D,metric=parameters(kind,mode);eps=s.diff(D,t);bs=ordinary.CASES[kind]['basis'];rank=len(bs);delta=t*t*D;h=t**3*D*D
    source=s.cancel(euler.ordinary_source(W,D)/t);vector=lambda f:s.cancel(delta*s.diff(delta,t)*(u*s.diff(f,u)+v*s.diff(f,v))/(4*t))
    if not eps:assert s.Poly(source,t,u,v).total_degree()<=4
    pv=s.Poly(W,u,v);co=[pv.coeff_monomial(x) for x in [u*u*v*v,u**3,v**3,u*u,v*v,u*v,u,v,1]];prepared=ordinary.operator_data(kind,*co,delta)
    field=s.QQ.frac_field(t);poly=lambda z:s.Poly(z,u,v,domain=field);zero=poly(0);wp,dp=poly(W),poly(delta)
    @lru_cache(None)
    def reduce(ex,depth):
        if ex in bs:
            out=[field.zero]*rank;out[bs.index(ex)]=field.one;return tuple(out),zero,zero
        assert depth>0,(kind,ex)
        rem,Cu,Cv,_,_=ordinary.one_step(kind,ex,co,delta,prepared);rem,Cu,Cv=[s.cancel(x/prepared[-1]) for x in [rem,Cu,Cv]]
        out=[field.zero]*rank;Cu,Cv=poly(Cu),poly(Cv)
        for f,c in poly(rem).terms():
            if c:
                vals,x,z=reduce(f,depth-1);cp=field.convert(c);out=[a+cp*b for a,b in zip(out,vals)];Cu+=poly(c)*x;Cv+=poly(c)*z
        return tuple(out),Cu,Cv
    Le,Lf,scales,grades,clear=row_parameters(kind,W,D);rows=[];columns=[[s.Poly(int(i==0),t,domain=s.QQ) for i in range(rank)]];Hcols=[s.Poly(1,t,domain=s.QQ)];R=s.Integer(1)
    slope={'nodal':5,'transverse':7,'triangle':9}[kind]-(0 if eps else 1);offset=2 if kind=='transverse' else 0
    expected={('nodal',1):123,('nodal',0):101,('transverse',1):122,('transverse',0):106,('triangle',1):114,('triangle',0):103}[(kind,int(bool(eps)))]
    for j in range(1,rank+1):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+vector(R)+source*R);out=[field.zero]*rank;Cu,Cv=zero,zero
        for ex,c in poly(R).terms():
            if c:
                vals,x,z=reduce(ex,3*j-(1 if kind=='triangle' else 2));cp=field.convert(c);out=[a+cp*b for a,b in zip(out,vals)];Cu+=poly(c)*x;Cv+=poly(c)*z
        E=clear(j);values=[s.Poly(s.cancel(field.to_sympy(c)*E/scale),t,domain=s.QQ) for c,scale in zip(out,scales)]
        Cu,Cv=poly(E)*Cu,poly(E)*Cv;normal=sum((poly(z.as_expr()*scale*u**ex[0]*v**ex[1]) for z,scale,ex in zip(values,scales,bs)),zero)
        assert (poly(E*R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        assert all(z.degree()<=slope*j-offset-g for z,g in zip(values,grades)),(kind,mode,j,[z.degree() for z in values],[slope*j-offset-g for g in grades])
        H=s.Poly(s.expand(h**j*E),t,domain=s.QQ);assert min(e[0] for e,c in H.terms() if c)>=3*j
        columns.append(values);Hcols.append(H);rows.append({'j':j,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':min(e[0] for e,c in H.terms() if c),'row_degree_bounds':[slope*j-offset-g for g in grades],'exact_twisted_identity':True})
    cofactors=[(-1)**k*det(columns[:k]+columns[k+1:],rank) for k in range(rank+1)]
    assert all(sum((cofactors[k]*columns[k][i] for k in range(rank+1)),s.Poly(0,t)).is_zero for i in range(rank))
    scalars=[]
    for k,C in enumerate(cofactors):
        if eps:
            C.exquo(s.Poly(D**(rank-1 if k==0 else max(0,rank-2)),t,domain=s.QQ))
            full=(C*Hcols[k]).exquo(s.Poly(D**(rank-1),t,domain=s.QQ))
        else:full=C*Hcols[k]
        if k==0:assert not full.is_zero and full.degree()<=expected,(kind,mode,full.degree(),expected)
        val=min(e[0] for e,c in full.terms() if c) if full else None
        if k and full:assert val>=3*k,(kind,mode,k,val)
        scalars.append({'order':k,'raw_minor_degree':int(C.degree()) if C else -1,'minor_D_factor':rank-1 if k==0 and eps else max(0,rank-2) if eps else 0,'reduced_scalar_degree':int(full.degree()) if full else -1,'reduced_scalar_valuation':val})
    return {'kind':kind,'mode':mode,'D':str(D),'D_degree':int(bool(eps)),'Le':str(Le),'Lf':str(Lf),'metric_determinant':int(metric),'seed_columns':rows,'scalar_cofactors':scalars,'cofactor_gate':expected}

def arithmetic():
    rows=[]
    for kind in ['nodal','transverse','triangle']:
        bs=ordinary.CASES[kind]['basis'];dummyW=s.Poly(0,u,v)
        if kind=='triangle':grades=[i+k+int(s.ceiling(i+s.Rational(k,2)))+int(s.ceiling(s.Rational(i,2)+k)) for i,k in bs]
        else:grades=[sum(e) for e in bs]
        for kd in [1,0]:
            slope={'nodal':5,'transverse':7,'triangle':9}[kind]-(0 if kd else 1);offset=2 if kind=='transverse' else 0;cases=[]
            for size in range(1,len(bs)+1):
                for selected in combinations(range(len(bs)),size):
                    if 0 not in selected:continue
                    before=slope*size*(size+1)//2-offset*size-sum(grades[k] for k in selected);after=before-(size-1 if kd else 0)
                    cases.append({'rank':size,'rows':selected,'old_degree_bound':before,'new_degree_bound':after})
            gate=max(c['new_degree_bound'] for c in cases);expected={'nodal':[101,123],'transverse':[106,122],'triangle':[103,114]}[kind][kd];assert gate==expected
            rows.append({'kind':kind,'D_degree':kd,'gate':gate,'all_rank_row_cases':cases})
    cover={'quadratic':7,'smooth':94,'nodal':123,'cusp':125,'transverse':122,'tangent':120,'triangle':114,'concurrent':104,'repeated':99};assert max(cover.values())==125
    return rows,cover

def main():
    start=time.monotonic();univ=universals();models=[];print('Universal Gaussian-square residue identities passed',flush=True)
    for kind in ['nodal','transverse','triangle']:
        for mode in ['linear','constant','zero_shared']:
            models.append(model(kind,mode));print('Model',kind,mode,'passed',round(time.monotonic()-start,2),flush=True)
    cases,cover=arithmetic();cert={'universal':univ,'models':models,'degree_cases':cases,'global_cover':cover};path=PTH/'ordinary_rank_one_lattice_certificate.json';path.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(path),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'universal_strata':3,'exact_F_monomial_identities':sum(len(c['exact_F_monomial_identities']) for c in univ),'rank_one_augmented_minors':sum(c['augmented_two_by_two_minors'] for c in univ),'models':9,'exact_seed_identities':sum(len(m['seed_columns']) for m in models),'full_scalar_cofactor_relations':9,'linear_D_gates':{'nodal':123,'transverse':122,'triangle':114},'constant_D_gates':{'nodal':101,'transverse':106,'triangle':103},'global_sufficient_pure_prefix':125,'odd_sufficient_pure_prefix_unchanged':98,'homogeneous_sufficient_pure_prefix_unchanged':96,'least_prefix_interval':[20,125],'odd_least_prefix_interval':[20,98],'homogeneous_least_prefix_interval':[14,96],'global_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact Gaussian-square F-monomial identities, source rank-one residues and augmented minors,54 complete original twisted seed identities and nine scalar cofactor relations through full rank, constant-D source/normal degree bounds, common D-factor and scalar valuations including shared axis/D/t roots, all-rank/row arithmetic. Universal all-order polynomial-module/seed/minor induction, complete coefficient/frame cover and retained analytic exclusions remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (PTH/'ordinary_rank_one_lattice_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
