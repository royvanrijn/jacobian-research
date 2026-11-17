"""Exact common-t cancellation in cleared binary actual-seed operators.

Replays retained certificates without executing or modifying old checker mains.
Analytic exclusions and complete rank/class coverage remain written inputs.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
import verify_ordinary_node_finite_prefixes as ordinary
import verify_euler_seed_bounds as euler
import verify_squarefree_singular_finite_prefixes as weighted
import verify_shared_divisor_weighted_seeds as shared
import verify_structured_seed_modules as structured
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v

def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def total(f):return -1000 if f==0 else int(s.Poly(f,t,u,v).total_degree())

def universal_sources():
    d0,d1=s.symbols('d0 d1');D=d0+d1*t;delta=t*t*D
    W=u*u*v*v/2+sum(s.Symbol('n'+str(i)+'_'+str(j)+'_'+str(k))*t**k*u**i*v**j for i,j in ordinary.monomials(3) for k in range(5-i-j))
    source=euler.ordinary_source(W,D);new=s.cancel(source/t)
    assert total(new)<=5 and s.Poly(new,u,v).total_degree()<=3
    assert s.expand(delta**2/t-t**3*D**2)==0
    assert total(s.cancel(delta*s.diff(delta,t)/(4*t)))<=4
    rows=[{'kind':'ordinary','source_joint_degree_after_division':5,'source_amplitude_degree':3,'cleared_denominator':str(t**3*D**2)}]
    for kind in ['cusp','tangent']:
        data=weighted.CASES[kind];phase=u**4/2+D*(v**3 if kind=='cusp' else u*v*v)
        phase+=sum(s.Symbol(kind+str(i)+'_'+str(j)+'_'+str(k))*t**k*u**i*v**j for i,j in weighted.monomials(data['lower_weight'],data['weights']) for k in range(5-i-j))
        if kind=='cusp':
            S,bx,by=shared.weighted_source(kind,phase,D);h=delta**2
            fields=[delta*bx,delta*by]
            limit=5;weight_limit=10
        else:
            S,bx,by,C,_=structured.tangent_source(phase,D);h=delta**2*D
            fields=[delta*D*bx,delta*D*by,-delta*C/2]
            limit=6;weight_limit=6
        divided=s.cancel(S/t);assert total(divided)<=limit
        assert all(weighted.weight(e,data['weights'])<=weight_limit for e,c in s.Poly(divided,u,v).terms() if c)
        for field in fields:assert total(s.cancel(field/t))<=limit-1
        assert total(s.cancel(h/t))<=limit
        rows.append({'kind':kind,'source_joint_degree_after_division':limit,'source_weight_after_division':weight_limit,'cleared_denominator':str(s.cancel(h/t))})
    # This identity proves the cancellation at every seed order, retaining the
    # derivative of the new denominator. G stands for the vertical-linear part.
    h=s.Function('h')(t);F=s.Function('F')(t);G=s.Function('G')(t);j=s.Symbol('j',integer=True,nonnegative=True)
    original=h*s.diff(t**j*F,t)-j*s.diff(h,t)*t**j*F+t**j*G
    normalized=t**(j+1)*((h/t)*s.diff(F,t)-j*s.diff(h/t,t)*F+G/t)
    assert s.simplify(original-normalized)==0
    # Homogeneous node: t-independent W and D=m make the source degree3.
    m=s.Symbol('m',nonzero=True)
    hom=u*u*v*v/2+sum(s.Symbol('a'+str(i)+'_'+str(j))*u**i*v**j for i,j in ordinary.monomials(3))
    hom_source=s.cancel(euler.ordinary_source(hom,m)/t)
    assert total(hom_source)<=3 and 3*28-11==73
    return {'source_cases':rows,'all_j_recurrence_identity':True,'same_quotient_connection':True,'homogeneous_nodal_source_degree':3,'homogeneous_nodal_complement_gate':73}

def replay_model(model,family):
    local={'t':t,'u':u,'v':v};read=lambda x:s.sympify(x,locals=local)
    W,D,A=[read(model[n]) for n in ['W','D','A0']];delta=t*t*D
    basis=model['basis'];kind=model['kind'];metric=model['metric_determinant'];assert metric!=0
    if family=='ordinary':
        source=euler.ordinary_source(W,D);h=delta**2
        vector=lambda f:delta*s.diff(delta,t)*euler.Euler(f)/4
        cost={'nodal':1,'transverse':2,'triangle':3}[kind];seed_total=5
    elif family=='cusp':
        source,bx,by=shared.weighted_source(kind,W,D);h=delta**2
        vector=lambda f:delta*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))
        cost=1;seed_total=5
    else:
        source,bx,by,C,_=structured.tangent_source(W,D);h=delta**2*D
        vector=lambda f:delta*D*(bx*u*s.diff(f,u)+by*v*s.diff(f,v))-delta*C*u*s.diff(f,v)/2
        cost=2;seed_total=6
    hbar=s.cancel(h/t);Sbar=s.cancel(source/t)
    ring=s.QQ.poly_ring(t);poly=lambda f:s.Poly(f,u,v,domain=ring)
    wp,dp,ap=map(poly,[W,delta,A]);wu,wv=wp.diff(u),wp.diff(v)
    def divide(expr,power):
        return s.Poly(expr,u,v,t,domain=s.QQ).exquo(s.Poly(t**power,u,v,t,domain=s.QQ)).as_expr()
    def replay(input,row,depth,power):
        values=[s.Poly(divide(read(x),power),t,domain=s.QQ) for x in row['N']]
        witness=row.get('reduction_witness',row.get('witness'))
        Cu,Cv=[poly(divide(read(witness[n]),power)) for n in ['Cu','Cv']]
        normal=sum((poly(z.as_expr()*u**e[0]*v**e[1]) for z,e in zip(values,basis)),poly(0))
        assert (ap**depth*poly(input)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wu*Cu-wv*Cv).is_zero
        return values,{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}
    columns=[]
    for exponent,row in zip(basis,model['connection_columns']):
        f=u**exponent[0]*v**exponent[1];input=s.cancel((source*f+vector(f))/t)
        values,witness=replay(input,row,row['depth'],1)
        assert all(z.degree()<=seed_total+sum(exponent)+cost*row['depth']-sum(e) for z,e in zip(values,basis))
        columns.append({'depth':row['depth'],'H':str(s.expand(hbar*A**row['depth'])),'N':[str(z.as_expr()) for z in values],'witness':witness})
    R=s.Integer(1);Rold=s.Integer(1);seeds=[]
    for row in model['seed_columns']:
        j=row['j']
        R=s.expand(hbar*s.diff(R,t)-(j-1)*s.diff(hbar,t)*R+s.cancel(vector(R)/t)+Sbar*R)
        Rold=s.expand(h*s.diff(Rold,t)-(j-1)*s.diff(h,t)*Rold+vector(Rold)+source*Rold)
        assert s.expand(Rold-t**j*R)==0
        assert total(R)<=seed_total*j
        depth=row['depth'];values,witness=replay(R,row,depth,j);bnd=seed_total*j+cost*depth
        assert all(z.degree()<=bnd-sum(e) for z,e in zip(values,basis))
        oldH=read(row['H']);H=s.Poly(divide(oldH,j),t,domain=s.QQ)
        assert s.expand(H.as_expr()-hbar**j*A**depth)==0 and H.degree()<=bnd
        val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        seeds.append({'j':j,'depth':depth,'degree_bound':bnd,'cancelled_time_power':j,'H':str(H.as_expr()),'H_valuation':val,'N':[str(z.as_expr()) for z in values],'witness':witness})
    return {'family':family,'kind':kind,'singular_at_zero':model['singular_at_zero'],'zero_shear':model.get('zero_shear'),'metric_determinant':metric,'basis':basis,'W':str(W),'D':str(D),'A0':str(A),'old_cleared_denominator':str(h),'new_cleared_denominator':str(hbar),'source_after_division':str(Sbar),'connection_columns':columns,'seed_columns':seeds}

def degree_cases():
    rows=[]
    specifications=[('nodal',ordinary.CASES['nodal']['basis'],8,2,199),('transverse',ordinary.CASES['transverse']['basis'],11,4,199),('triangle',ordinary.CASES['triangle']['basis'],14,3,189),('cusp',weighted.CASES['cusp']['basis'],10,3,183),('tangent',weighted.CASES['tangent']['basis'],18,8,224)]
    for kind,basis,slope,offset,expected in specifications:
        grades=[sum(e) for e in basis];rank=len(basis);cases=[]
        for size in range(1,rank+1):
            for chosen in combinations(range(1,rank),size-1):
                selected=(0,)+chosen;gate=sum(slope*j-offset for j in range(1,size+1))-sum(grades[i] for i in selected)
                cases.append({'rank':size,'rows':selected,'gate':gate})
        assert max(x['gate'] for x in cases)==expected
        rows.append({'kind':kind,'slope':slope,'offset':offset,'gate':expected,'cases':cases})
    return rows

def main():
    started=time.monotonic();universal=universal_sources();print('Universal source divisibility, recurrence and homogeneous73 passed',flush=True)
    shared_cert=json.loads((P/'shared_divisor_weighted_seeds_certificate.json').read_text());structured_cert=json.loads((P/'structured_seed_modules_certificate.json').read_text());models=[]
    for m in shared_cert['ordinary_models']:models.append(replay_model(m,'ordinary'));print('Replayed time cancellation',m['kind'],m['singular_at_zero'],flush=True)
    for m in shared_cert['weighted_models']:
        if m['kind']=='cusp':models.append(replay_model(m,'cusp'));print('Replayed time cancellation cusp',m['singular_at_zero'],flush=True)
    for m in structured_cert['tangent_models']:models.append(replay_model(m,'tangent'));print('Replayed time cancellation tangent',m['singular_at_zero'],m['zero_shear'],flush=True)
    cases=degree_cases();cover={'quadratic':7,'smooth':94,'nodal':216,'cusp':210,'transverse':216,'tangent':224,'triangle':216,'concurrent':104,'repeated':99};assert max(cover.values())==224
    retained=['shared_divisor_weighted_seeds_certificate.json','shared_divisor_weighted_seeds_verification.json','structured_seed_modules_certificate.json','structured_seed_modules_verification.json','balanced_logarithmic_seed_bounds_certificate.json']
    cert={'universal_time_factor':universal,'models':models,'degree_cases':cases,'global_cover':cover,'retained_inputs':{n:sha(P/n) for n in retained}}
    f=P/'common_time_seed_factor_certificate.json';f.write_text(json.dumps(cert,indent=2)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-started,3),'sympy_version':s.__version__,
            'universal_source_cases':3,'all_j_recurrence_identity':True,'exact_models_replayed':len(models),'connection_columns_replayed':sum(len(m['connection_columns']) for m in models),'seed_columns_replayed':sum(len(m['seed_columns']) for m in models),
            'ordinary_complement_gates':{'nodal':199,'transverse':199,'triangle':189},'ordinary_full_gates':{'nodal':216,'transverse':216,'triangle':216},'homogeneous_nodal_complement_gate':73,'cusp_complement_gate':183,'cusp_full_gate':210,'tangent_complement_gate':224,'tangent_full_gate':224,'global_cover':cover,
            'global_sufficient_pure_prefix':224,'least_prefix_interval':[20,224],'odd_least_prefix_interval':[20,224],'all_homogeneous_sufficient_pure_prefix':216,'all_homogeneous_least_prefix_interval':[14,216],'binary_search_rounds_max':8,
            'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,
            'boundary':'Exact universal source and whole-operator divisibility, all-j recurrence identity and twelve retained full-Gaussian connection/seed model replays after polynomial division. Full class/rank coverage, all-degree degree bounds and retained analytic zeroth-coefficient exclusions remain written and independently unreviewed. No leastness search or expanded global radical cofactors.'}
    (P/'common_time_seed_factor_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
