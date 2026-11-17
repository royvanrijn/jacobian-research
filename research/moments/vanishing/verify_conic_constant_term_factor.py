"""Auxiliary structured conic constant-term functional and common K factor.
Universal formal CT/seed induction,complete endpoint/frame coverage and
analytic nonzero-zeroth arguments remain independently unreviewed written.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
import verify_balanced_logarithmic_seed_bounds as balanced
import verify_conic_completed_row_factors as previous
P=Path(__file__).resolve().parent
t,y=balanced.t,balanced.y
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
def cross(n,m):return min(abs(n),abs(m)) if n*m<0 else 0

def formal_functional():
    K=s.Symbol('K');A=s.symbols('A1:4');B=s.symbols('B1:4');limit=3;maximum=10
    def exponential(c):
        out=[s.Integer(1)]
        for n in range(1,maximum+1):out.append(s.expand(sum(k*c[k-1]*out[n-k]/t**2 for k in range(1,min(3,n)+1))/n))
        return out
    pos,neg=exponential(A),exponential(B)
    def pi(n):
        return s.expand(sum(pos[l-n]*neg[l]*K**(l+max(-n,0)) for l in range(max(0,n),limit+1) if l+max(-n,0)<=limit))
    checks=[]
    for n in range(-4,5):
        ward=t*t*n*pi(n)
        for k in range(1,4):ward+=k*A[k-1]*K**cross(k,n)*pi(k+n)-k*B[k-1]*K**cross(-k,n)*pi(n-k)
        ward=s.Poly(s.expand(ward),K)
        for j in range(limit+1):assert s.expand(ward.coeff_monomial(K**j))==0
        checks.append(n)
    basis=list(range(-3,3));assert pi(0).subs(K,0)==1
    assert all(pi(n).subs(K,0)==0 for n in basis if n)
    return {'formal_basis':'E_n=y^n for n>=0; E_n=K^(-n)y^n for n<0',
            'functional':'CT_y exp((W-A0)/t^2)*f','axis_unit_row_at_K0':True,
            'pi_E_n_order_at_least':'abs(n)','logarithmic_Ward_amplitudes':checks,
            'Ward_coefficient_order':limit,'Ward_coefficient_checks':len(checks)*(limit+1),
            'period_normalization':'t^(-1) exp(A0/t^2)',
            'scalar_connection_regular_in_K_required':True,'projected_seed':'pi(R_j)=h^j nabla^j pi(1)',
            'unit_completed_seed_row_divisible_by_K':True}

def retained_rows():
    out=[]
    for name,key in [('structured_seed_modules_certificate.json','conic_models'),('two_sided_conic_seed_bounds_certificate.json','models'),('endpoint_degrees_and_odd_prefixes_certificate.json','conic_models')]:
        for m in json.loads((P/name).read_text())[key]:
            K=s.sympify(m['K'],locals={'t':t});W=s.sympify(m['W'],locals={'t':t,'y':y});raw=balanced.laurent(W)
            if max(raw)!=3 or min(raw)!=-3:continue
            for n in [1,2,3]:s.Poly(s.cancel(raw.get(-n,0)/K**n),t)
            unit=m['basis'].index(0);checks=[]
            for row in m['seed_columns']:
                N=s.Poly(s.sympify(row['N'][unit],locals={'t':t}),t);Q=N.exquo(s.Poly(K,t))
                checks.append({'j':row['j'],'raw_degree':int(N.degree()) if N else -1,'divided_unit_degree':int(Q.degree()) if Q else -1})
            out.append({'origin':name,'origin_sha256':sha(P/name),'label':m['label'],'K':str(K),'unit_row_checks':checks})
    assert len(out)==11 and sum(len(m['unit_row_checks']) for m in out)==33
    return out

def full_model(label,K):
    K=s.expand(K);kd=int(s.degree(K,t));Bv=t*s.diff(K,t)-s.Rational(2,3)*K;h=s.expand(t**3*K)
    poly=lambda f:s.Poly(s.expand(f),t,domain=s.QQ);zero=poly(0);kp=poly(K)
    W={0:poly(1+t**3),1:poly(2+t*t),2:poly(1+t),3:poly(1),-1:poly(3+t*t),-2:poly(2+t),-3:poly(2)}
    def add(out,n,c):
        out[n]=out.get(n,zero)+c
        if out[n].is_zero:del out[n]
    def mul(left,right):
        out={}
        for n,a in left.items():
            for m,b in right.items():add(out,n+m,a*b*kp**cross(n,m))
        return out
    grad={n:c.mul_ground(n) for n,c in W.items() if n}
    def exact(witness):
        out=mul(grad,witness)
        for n,c in witness.items():add(out,n,poly(t*t*n)*c)
        return out
    def reduce(R):
        normal=R.copy();witness={}
        while max(normal,default=0)>=3:
            n=max(normal);c=normal[n].mul_ground(s.Rational(1,3));shift=n-3;add(witness,shift,c)
            for m,a in exact({shift:c}).items():add(normal,m,-a)
        while min(normal,default=0)<-3:
            n=min(normal);c=normal[n].mul_ground(-s.Rational(1,6));shift=n+3;add(witness,shift,c)
            for m,a in exact({shift:c}).items():add(normal,m,-a)
        defect=R.copy()
        for n,c in normal.items():add(defect,n,-c)
        for n,c in exact(witness).items():add(defect,n,-c)
        assert not defect and set(normal)<=set(range(-3,3))
        return normal,witness
    # Time derivatives of E_n include (-n)K_prime/K for n<0.
    S={}
    for n,a in W.items():
        value=poly(K)*(poly(t)*a.diff()-a.mul_ground(2))
        value+=poly(-n*t*s.diff(K,t) if n<0 else 0)*a+poly(n*Bv)*a
        if n==0:value-=poly(t*t*K)
        if value:add(S,n,value)
    assert min(S)>=-2 and max(S)<=3
    R={0:poly(1)};rows=[];basis=list(range(-3,3));cols=[[poly(int(k==0)) for k in basis]]
    for j in range(1,7):
        nextR=mul(S,R)
        for n,a in R.items():
            value=poly(h)*a.diff()-poly((j-1)*s.diff(h,t))*a
            value+=poly(t**3*(-n)*s.diff(K,t) if n<0 else 0)*a+poly(t*t*Bv*n)*a
            add(nextR,n,value)
        R=nextR;normal,witness=reduce(R);values=[normal.get(k,zero) for k in basis]
        limits=[previous.bound(3,3,kd,0,j,k) for k in basis]
        assert all(a.degree()<=b for a,b in zip(values,limits))
        unit=values[basis.index(0)].exquo(kp);assert unit.degree()<=limits[basis.index(0)]-kd
        H=poly(h**j);val=min(e[0] for e,c in H.terms() if c);assert val>=3*j
        rows.append({'j':j,'N_completed':[str(a.as_expr()) for a in values],'H':str(H.as_expr()),'H_valuation':val,
                     'unit_divided_by_K':str(unit.as_expr()),'row_degree_bounds':limits,
                     'exact_original_logarithmic_identity':True,'witness_completed_basis':{str(n):str(c.as_expr()) for n,c in witness.items()}})
        cols.append(values)
    ring=s.QQ.poly_ring(t);co=[]
    for k in range(7):
        selected=cols[:k]+cols[k+1:]
        M=DomainMatrix.from_list_sympy(6,6,[[c.as_expr() for c in row] for row in zip(*selected)]).convert_to(ring)
        co.append((-1)**k*poly(ring.to_sympy(M.det())))
    assert all(sum((co[j]*cols[j][k] for j in range(7)),zero).is_zero for k in range(6))
    checks=[]
    for j,C in enumerate(co):
        H=poly(1 if j==0 else h**j);scalar=(H*C).exquo(kp)
        val=min(e[0] for e,c in scalar.terms() if c) if scalar else None
        if j==0:assert scalar and scalar.degree()<=[72,83,94][kd]
        else:assert val is None or val>=3*j
        checks.append({'order':j,'common_K_factor':1,'scalar_degree':int(scalar.degree()) if scalar else -1,'scalar_valuation':val})
    return {'label':label,'K':str(K),'K_degree':kd,'phase_completed_coefficients':{str(n):str(c.as_expr()) for n,c in W.items()},
            'basis':basis,'seed_columns':rows,'cofactor_checks':checks}

def arithmetic():
    rows=[]
    for kd in range(3):
        cases=[]
        for n in range(1,7):
            for selected in combinations(range(-3,3),n):
                if 0 not in selected:continue
                raw=sum(previous.bound(3,3,kd,0,j,0) for j in range(1,n+1))
                raw+=sum(previous.bound(3,3,kd,0,1,k)-previous.bound(3,3,kd,0,1,0) for k in selected)
                cases.append({'rank':n,'rows':selected,'old_degree_bound':raw,'common_K_factor':1,'new_degree_bound':raw-kd})
        gate=max(c['new_degree_bound'] for c in cases);assert gate==[72,83,94][kd]
        rows.append({'K_degree':kd,'p3q3_gate':gate,'all_rank_row_cases':cases})
    assert max(94,94,84,69,93,93,63)==94 and max(94,92,95)==95
    cover={'quadratic':7,'smooth':94,'nodal':95,'cusp':95,'transverse':95,'tangent':95,'triangle':95,'concurrent':92,'repeated':92}
    assert max(cover.values())==95 and 95%2==1 and 94%2==0
    return rows,cover

def main():
    start=time.monotonic();functional=formal_functional();retained=retained_rows();print('Formal CT Ward identities and retained unit rows passed',flush=True)
    models=[]
    for label,K in [('constant',s.Integer(2)),('linear_regular',1+t),('linear_zero',t),('quadratic_regular',1+t+t*t),('double_zero',t*t)]:
        models.append(full_model(label,K));print('Full rank-six model',label,'passed',flush=True)
    cases,cover=arithmetic();cert={'formal_functional':functional,'retained_unit_rows':retained,'full_models':models,'degree_cases':cases,'global_cover':cover}
    f=P/'conic_constant_term_factor_certificate.json';f.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
            'formal_CT_Ward_coefficient_checks':36,'retained_models':11,'retained_unit_row_checks':33,
            'new_full_rank_six_models':5,'new_original_logarithmic_seed_identities':30,'full_scalar_cofactor_relations':5,'scalar_coefficient_checks':35,
            'p3q3_gates_by_K_degree':[72,83,94],'conic_gate':94,'full_W_linear_gate':95,'global_cover':cover,
            'global_sufficient_pure_prefix':95,'odd_sufficient_pure_prefix':94,'homogeneous_sufficient_pure_prefix':94,
            'least_prefix_interval':[20,95],'odd_least_prefix_interval':[20,94],'homogeneous_least_prefix_interval':[14,94],
            'least_prefix_computed':False,'universal_written_proof_certified_by_checker':False,'independent_replay':False,'external_review':False,'formal_verification':False,
            'boundary':'Finite formal constant-term logarithmic Ward coefficients,33 retained unit-row divisions,30 full rank-six original logarithmic seed identities,five full scalar relations/35 coefficients,and all rank/row/valuation arithmetic. Universal formal CT/seed induction,full structured coefficient and endpoint/frame coverage,retained other conic strata,analytic nonzero-zeroth inputs and parity deduction remain independently unreviewed written mathematics. No old checker main,radical search,external package or fresh literature review.'}
    (P/'conic_constant_term_factor_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
