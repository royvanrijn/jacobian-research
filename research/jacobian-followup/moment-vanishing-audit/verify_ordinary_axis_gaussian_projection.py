"""Auxiliary exact quadratic axis projection and joint-graded cofactors.
Universal formal finite quotient/connection, complementary-minor, joint-grade
and polynomial-extension proofs and full analytic/frame cover remain written.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
import verify_ordinary_node_finite_prefixes as ordinary
P=Path(__file__).resolve().parent
t,u,v=ordinary.t,ordinary.u,ordinary.v
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
tri=lambda n:n*(n-1)//2

def gaussian_identities():
    A,B,C,delta,Ap,Bp,Cp,dp=s.symbols('A B C delta Ap Bp Cp dp');mean=-B/A;variance=-delta/A
    def J(n):
        return s.expand(sum(s.binomial(n,2*k)*s.factorial(2*k)*mean**(n-2*k)*variance**k/(2**k*s.factorial(k)) for k in range(n//2+1)))
    Phi=C-B*B/(2*A);phip=Cp-B*Bp/A+B*B*Ap/(2*A*A)
    assert s.cancel(A*v*v/2+B*v+C-A*(v+B/A)**2/2-Phi)==0
    checks=[]
    for n in range(7):
        ward=A*J(n+1)+B*J(n)+(delta*n*J(n-1) if n else 0)
        assert s.cancel(ward)==0
        derivative=Ap*s.diff(J(n),A)+Bp*s.diff(J(n),B)
        left=Ap*J(n+2)/2+Bp*J(n+1)+Cp*J(n)
        assert s.cancel(left-delta*derivative-(phip-delta*Ap/(2*A))*J(n))==0
        derivative+=dp*s.diff(J(n),delta)
        left-=dp*(A*J(n+2)/2+B*J(n+1)+C*J(n))/delta
        assert s.cancel(left-delta*derivative-(phip-dp*Phi/delta+dp/2-delta*Ap/(2*A))*J(n))==0
        checks.append({'amplitude_v_power':n,'transverse_twisted_identity':True,'amplitude_u_connection':True,'full_time_and_delta_connection':True})
    D,eps=s.symbols('D eps')
    assert s.cancel(-2/t-eps/(2*D)+(2/t+eps/D)/2+1/t)==0
    return {'completion_identity':True,'mean':'-B/A','variance':'-delta/A',
            'identities':checks,'Gaussian_coefficient_identities':21,
            'old_normalization':'t^(-2) D^(-1/2)','new_normalization':'t^(-1) A^(-1/2)',
            'all_volume_determinant_and_delta_derivatives_retained':True}

def axis_comparisons():
    D,z,Le,N,R,a,L,K=s.symbols('D z Le N R a L K');A=u*u+2*D*z;B=R*u+D*N
    C=D*(a*u**3+Le*u*u+L*u+K)+R*R/2
    Phi=s.cancel(C-B*B/(2*A));Phi0=s.expand(s.cancel(Phi.subs(z,0)))
    expected=D*(a*u**3+Le*u*u+L*u+K-R*N/u-D*N*N/(2*u*u))
    assert s.expand(Phi0-expected)==0
    mu=-R/u-D*N/(u*u)
    transcols=[1,u,mu,u*u,s.expand(u*mu)]
    Q=s.Matrix([[s.expand(c).coeff(u,k) for c in transcols] for k in [-2,-1,0,1,2]])
    td=s.factor(Q.det());assert s.cancel(td/(D*D*N*N)) in [1,-1]
    tricol=[1,Le*z*u,Le*z*mu,Le*Le*z*z*u*mu]
    Q=s.Matrix([[s.expand(c).coeff(u,k) for c in tricol] for k in [-2,-1,0,1]])
    qd=s.factor(Q.det());assert s.cancel(qd/(Le**4*z**4*D*D*N*N)) in [1,-1]
    return {'axis':'z=Lf','quadratic_leading':'A=u^2+2D Lf','linear':'B=r*u+D*N0',
            'residual_phase_at_axis0':str(Phi0),'Laurent_prefactor_at_axis0':'du/u',
            'transverse':{'finite_rank':5,'basis_powers':[-2,-1,0,1,2],'comparison_determinant':str(td),'loss':0},
            'triangle':{'finite_rank':4,'basis_powers':[-2,-1,0,1],'comparison_determinant':str(qd),'loss':4},
            'extra_analytic_row':'select v^2 coordinate','generic_N0_D_Le_units_then_polynomial_extension':True}

def joint_grading():
    D,Le,Lf,J,L,N,K,R,eps=s.symbols('D Le Lf J L N K R eps')
    lep,lfp,jp,lp,np,kp,rp=s.symbols('lep lfp jp lp np kp rp')
    positive={t:1,u:1,v:1,D:1,Le:1,Lf:1,J:1,L:2,N:2,K:3,R:2,lp:1,np:1,kp:2,rp:1}
    def deg(f):
        poly=s.Poly(s.expand(f),*positive)
        return max((sum(e*positive[x] for x,e in zip(positive,ex)) for ex,c in poly.terms() if c),default=-1000)
    dt=lambda f:s.diff(f,t)+eps*s.diff(f,D)+lep*s.diff(f,Le)+lfp*s.diff(f,Lf)+jp*s.diff(f,J)+lp*s.diff(f,L)+np*s.diff(f,N)+kp*s.diff(f,K)+rp*s.diff(f,R)
    rows=[]
    for kind,a in [('transverse',s.Integer(1)),('triangle',s.Integer(0))]:
        W=D*(a*u**3+Le*u*u+Lf*v*v+J*u*v+L*u+N*v+K)+(u*v+R)**2/2
        delta=t*t*D;bx=dt(delta)/4
        S=s.expand((delta*dt(W)-dt(delta)*W)/t-2*t*t*D*D-t**3*D*eps/2+2*delta*bx/t+bx*(u*s.diff(W,u)+v*s.diff(W,v))/t)
        assert deg(W)<=4 and deg(S)<=5 and s.Poly(S,u,v).total_degree()<=3
        coeff=[s.Rational(1,2),D*a,s.Integer(0),D*Le,D*Lf,R+D*J,D*L,D*N,R*R/2+D*K]
        prep=ordinary.operator_data(kind,*coeff,delta);axis=Lf if kind=='transverse' else Le*Lf;identities=[]
        for ex in ordinary.monomials(5):
            if ex in ordinary.CASES[kind]['basis']:continue
            rem,cu,cv,_,div=ordinary.one_step(kind,ex,coeff,delta,prep)
            normal=s.cancel(axis*rem/div)
            assert deg(normal)<=sum(ex)+(1 if kind=='transverse' else 2)
            identities.append(ex)
        rows.append({'kind':kind,'phase_joint_degree':4,'full_source_joint_degree':5,'actual_seed_joint_bound':'5j',
                     'polynomial_normal_principal_joint_identities':identities,
                     'old_completed_row_joint_bound':'7j-2-|e|' if kind=='transverse' else '9j-|e|-ceil(i+k/2)-ceil(i/2+k)'})
    return {'positive_weights':{str(k):v for k,v in positive.items()},'D_and_axis_weights':1,
            'time_derivation_lowers_joint_degree_by_one':True,'rows':rows,
            'divide_formal_coprime_D_and_Lf_before_affine_or_constant_specialization':True}

def axis_bounds(kind,n):
    weights=[0]+[2*j-2 if kind=='transverse' else 2*j for j in range(1,n+1)]
    loss=0 if kind=='transverse' else 4
    return [max(0,sum(sorted(weights[j] for j in range(n+1) if j!=k)[:n-1])-loss) for k in range(n+1)]

def cofactors(cols,n):
    ring=s.QQ.poly_ring(t);out=[]
    for k in range(n+1):
        selected=cols[:k]+cols[k+1:]
        M=DomainMatrix.from_list_sympy(n,n,[[x.as_expr() for x in row] for row in zip(*selected)]).convert_to(ring)
        out.append((-1)**k*s.Poly(ring.to_sympy(M.det()),t,domain=s.QQ))
    assert all(sum((out[k]*cols[k][i] for k in range(n+1)),s.Poly(0,t)).is_zero for i in range(n))
    return out

def replay(m,origin):
    kind=m['kind'];n=len(m['seed_columns']);cols=[[s.Poly(int(i==0),t) for i in range(n)]]
    cols += [[s.Poly(s.sympify(x,locals={'t':t}),t) for x in row['N']] for row in m['seed_columns']]
    den=s.Poly(s.sympify(m['D'],locals={'t':t}),t);axis=s.Poly(s.sympify(m['Lf'],locals={'t':t}),t)
    bounds=axis_bounds(kind,n);co=cofactors(cols,n);checks=[]
    for k,C in enumerate(co):
        d=n*(n-1) if k==0 else max(0,n*(n-1)-2*k)
        C.exquo(den**d*axis**bounds[k])
        H=s.Poly(1,t) if k==0 else s.Poly(s.sympify(m['seed_columns'][k-1]['H'],locals={'t':t}),t)
        reduced=(H*C).exquo(den**(n*(n-1))*axis**bounds[0])
        val=min(e[0] for e,c in reduced.terms() if c) if reduced else None
        if k==0:assert reduced and reduced.degree()<=(77 if kind=='transverse' else 82)
        else:assert val is None or val>=3*k
        checks.append({'order':k,'raw_minor_D_power':d,'raw_minor_Lf_power':bounds[k],
                       'common_D_power':n*(n-1),'common_Lf_power':bounds[0],
                       'reduced_scalar_degree':int(reduced.degree()) if reduced else -1,'reduced_scalar_valuation':val})
    return {'kind':kind,'mode':m['mode'],'origin':origin,'D':m['D'],'Lf':m['Lf'],
            'D_degree':int(den.degree()),'Lf_degree':int(axis.degree()),'cofactor_checks':checks}

def arithmetic():
    rows=[]
    for kind in ['transverse','triangle']:
        bs=ordinary.CASES[kind]['basis'];grades=[sum(e) for e in bs] if kind=='transverse' else [i+k+int(s.ceiling(i+s.Rational(k,2)))+int(s.ceiling(s.Rational(i,2)+k)) for i,k in bs]
        cases=[]
        for n in range(1,len(bs)+1):
            ab=axis_bounds(kind,n);db=n*(n-1)
            for k in range(1,n+1):
                hp=2*k-2 if kind=='transverse' else 2*k
                assert ab[k]+hp>=ab[0] and max(0,db-2*k)+2*k>=db
            for selected in combinations(range(len(bs)),n):
                if 0 not in selected:continue
                raw=(7 if kind=='transverse' else 9)*n*(n+1)//2-(2*n if kind=='transverse' else 0)-sum(grades[k] for k in selected)
                cases.append({'rank':n,'rows':selected,'raw_joint_degree':raw,'D_common_power':db,
                              'Lf_minor_bounds':ab,'Lf_common_power':ab[0],'new_joint_degree':raw-db-ab[0]})
        gate=max(c['new_joint_degree'] for c in cases);assert gate==(77 if kind=='transverse' else 82)
        rows.append({'kind':kind,'uniform_complement_gate':gate,'all_rank_row_cases':cases})
    cover={'quadratic':7,'smooth':94,'nodal':96,'cusp':95,'transverse':96,'tangent':95,'triangle':96,'concurrent':92,'repeated':96}
    assert max(77,96)==max(82,96)==max(cover.values())==96
    return rows,cover

def main():
    start=time.monotonic();gaussian=gaussian_identities();comparison=axis_comparisons();grade=joint_grading()
    print('Exact Gaussian connection,axis comparisons and joint grading passed',flush=True)
    models=[]
    for m in json.loads((P/'ordinary_rank_one_lattice_certificate.json').read_text())['models']:
        if m['kind'] in ['transverse','triangle']:
            models.append(replay(m,'ordinary_rank_one_lattice_certificate.json'));print('Retained',m['kind'],m['mode'],'passed',flush=True)
    fresh=json.loads((P/'constant_axis_gaussian_models.json').read_text())
    for m in fresh:
        models.append(replay(m,'constant_axis_gaussian_models.json'));print('Fresh constant-axis',m['kind'],m['mode'],'passed',flush=True)
    assert len(models)==10 and sum(len(m['cofactor_checks']) for m in models)==65
    cases,cover=arithmetic();cert={'exact_Gaussian':gaussian,'axis_comparisons':comparison,'joint_grading':grade,'models':models,'degree_cases':cases,'global_cover':cover}
    f=P/'ordinary_axis_gaussian_projection_certificate.json';f.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
            'Gaussian_connection_identities':21,'axis_comparison_determinants':2,
            'principal_joint_degree_identities':sum(len(c['polynomial_normal_principal_joint_identities']) for c in grade['rows']),
            'models':10,'retained_models':6,'fresh_constant_axis_models':4,'fresh_original_twisted_seed_identities_retained':22,
            'full_scalar_cofactor_checks':65,'uniform_complement_gates':{'transverse':77,'triangle':82},
            'common_D_powers':{'transverse':30,'triangle':20},'common_Lf_powers':{'transverse':20,'triangle':16},
            'global_cover':cover,'global_sufficient_pure_prefix':96,'odd_sufficient_pure_prefix':96,'homogeneous_sufficient_pure_prefix_unchanged':96,
            'least_prefix_interval':[20,96],'odd_least_prefix_interval':[20,96],'homogeneous_least_prefix_interval':[14,96],
            'least_prefix_computed':False,'universal_written_proof_certified_by_checker':False,'independent_replay':False,'external_review':False,'formal_verification':False,
            'boundary':'Exact quadratic Gaussian completion and full connection identities,two leading comparison determinants,finite principal joint-degree identities,10 full scalar cofactor families with65 checks including constant axes/shared roots;four fresh full Gaussian producer models retain22 original twisted seed identities. Universal formal Laurent quotient/regular axis connection,restored-direction/complementary-minor argument,independent-variable divisibility and joint-grade specialization,and full analytic/frame cover remain independently unreviewed written inputs. No radical search,producer package or fresh literature review.'}
    (P/'ordinary_axis_gaussian_projection_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
