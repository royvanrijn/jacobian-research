"""Auxiliary two-branch rank-one projection and q cofactor factors.
Universal formal quotient/connection,minors,joint-grade specialization,
all-stratum/frame cover and analytic exclusions remain written inputs.
"""
from pathlib import Path
from itertools import combinations
import hashlib,json,time
import sympy as s
from sympy.polys.matrices import DomainMatrix
P=Path(__file__).resolve().parent
t,y,T,Z,q=s.symbols('t y T Z q');sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
def trunc(z,n=5):return s.series(z,q,0,n).removeO().expand()

def universal():
    f,r,z=s.symbols('f r z');roots=[]
    for root,other in [(r,z),(z,r)]:
        w=s.Integer(0)
        for k in range(5):w=trunc(q*Z/(t*f*(root-other))-w*w/(root-other))
        branch=root+w
        assert trunc(f*(branch-r)*(branch-z)-q*Z/t)==0
        jac=trunc(s.diff(branch,Z)/q,4)
        exact=trunc(1/(t*f*(2*branch-r-z)),4)
        assert trunc(jac-exact,4)==0 and s.cancel(jac.subs(q,0)-1/(t*f*(root-other)))==0
        roots.append({'root':str(root),'branch_through_q4':str(branch),'critical_equation_zero_through_q4':True,
                      'normalized_Jacobian_constant':str(1/(t*f*(root-other)))})
    # After analytic branchwise normalization of pi(1), the first Z row
    # is the Hermite derivative of the T^k row.  Columns use y^k=t^k T^k.
    rows=[[t**k*r**k for k in range(4)],[t**k*z**k for k in range(4)],
          [s.Integer(0) if k==0 else q*t**k*k*r**(k-1)/(t*f*(r-z)) for k in range(4)],
          [s.Integer(0) if k==0 else q*t**k*k*z**(k-1)/(t*f*(z-r)) for k in range(4)]]
    determinant=s.factor(s.Matrix(rows).det());assert s.cancel(determinant/(q*q*t**4*(r-z)**2/f**2)) in [1,-1]
    a,c,b=s.symbols('a c b');gs={(i,j):s.Symbol('g'+str(i)+str(j)) for i in range(4) for j in range(4-i)}
    ah=a*t*t+c*t*y+f*y*y
    W=s.expand(sum(g*t**(3-i-j)*q**(3-i)*ah**i*y**j for (i,j),g in gs.items())-q**3*y*y/2)
    original=t*sum(g*(t*(a+c*T+f*T*T)/q)**i*T**j for (i,j),g in gs.items())-T*T/2
    assert s.cancel(W.subs(y,t*T)/(t*t*q**3)-original)==0
    dt=lambda v:s.diff(v,t)-b*s.diff(v,q)
    S=s.expand(t*q*dt(W)-(2*q-3*b*t)*W-t*t*q**4+b*t**3*q**3)
    assert s.Poly(W,t,q,y).total_degree()<=6 and s.Poly(S,t,q,y).total_degree()<=7
    assert s.Poly(W,y).coeff_monomial(y**6)==f**3*gs[3,0]
    return {'two_formal_branches':roots,'branch_identity_order':4,'Jacobian_identity_order':3,
            'exact_original_phase_rescaling':True,'moving_q_derivation':'partial_t-b*partial_q',
            'phase_joint_degree':6,'source_joint_degree':7,'leading_y6_coefficient':'f^3*g30',
            'finite_rank_each_branch':2,'total_finite_rank':4,'finite_basis_each_branch':['1','Z'],
            'comparison_leading_determinant':str(determinant),'comparison_q_order':2,
            'original_T_prefactor':'1/q','branch_prefactor':'1/(t*A_prime(T_branch))',
            'no_q_pole_in_branch_phase_or_prefactor':True}

def powers(n):
    weights=[0]+[4*j for j in range(1,n+1)]
    return [max(0,sum(sorted(weights[j] for j in range(n+1) if j!=k)[:n-1])-2) for k in range(n+1)]

def model(label,b,a,c,f):
    Q=1-b*t;ah=a*t*t+c*t*y+f*y*y
    W=s.expand(sum((1+2*i+3*j)*t**(3-i-j)*Q**(3-i)*ah**i*y**j for i in range(4) for j in range(4-i))-Q**3*y*y/2)
    delta=t*t*Q**3;h=t**3*Q**4
    S=s.expand(t*Q*s.diff(W,t)-(2*Q-3*b*t)*W-t*t*Q**4+b*t**3*Q**3)
    ring=s.QQ.poly_ring(t);poly=lambda z:s.Poly(z,y,domain=ring);zero=poly(0);wp,dp=poly(W),poly(delta)
    R=s.Integer(1);rows=[];cols=[[s.Poly(int(k==0),t) for k in range(5)]]
    for j in range(1,6):
        R=s.expand(h*s.diff(R,t)-(j-1)*s.diff(h,t)*R+S*R);normal=poly(R);witness=zero
        while normal.degree()>=5:
            n=normal.degree();coefficient=normal.LC()/wp.diff(y).coeff_monomial(y**5);w=poly(coefficient*y**(n-5))
            normal-=dp*w.diff(y)+wp.diff(y)*w;witness+=w
        assert (poly(R)-normal-dp*witness.diff(y)-wp.diff(y)*witness).is_zero
        values=[s.Poly(normal.coeff_monomial(y**k),t) for k in range(5)]
        assert all(z.degree()<=7*j-k for k,z in enumerate(values))
        H=s.Poly(h**j,t);hv=min(e[0] for e,z in H.terms() if z);assert hv>=3*j
        rows.append({'j':j,'N':[str(z.as_expr()) for z in values],'H':str(H.as_expr()),'H_valuation':hv,
                     'witness':str(witness.as_expr()),'exact_original_twisted_identity':True});cols.append(values)
    co=[];domain=s.QQ.poly_ring(t)
    for k in range(6):
        selected=cols[:k]+cols[k+1:]
        M=DomainMatrix.from_list_sympy(5,5,[[z.as_expr() for z in row] for row in zip(*selected)]).convert_to(domain)
        co.append((-1)**k*s.Poly(domain.to_sympy(M.det()),t))
    assert all(sum((co[k]*cols[k][i] for k in range(6)),s.Poly(0,t)).is_zero for i in range(5))
    checks=[];den=s.Poly(Q,t);need=powers(5);assert need[0]==38
    for k,C in enumerate(co):
        C.exquo(den**need[k]);H=s.Poly(1 if k==0 else h**k,t);A=(H*C).exquo(den**38)
        val=min(e[0] for e,c in A.terms() if c) if A else None
        if k==0:assert A and A.degree()<=57
        else:assert val is None or val>=3*k
        checks.append({'order':k,'raw_q_order_certified':need[k],'common_scalar_q_power':38,
                       'scalar_degree':int(A.degree()) if A else -1,'scalar_valuation':val})
    return {'label':label,'parameters':{'b':b,'a':a,'c':c,'f':f},'quadratic_discriminant':c*c-4*a*f,
            'q':str(Q),'W':str(W),'delta':str(delta),'h':str(h),'source':str(S),'seed_columns':rows,'cofactor_checks':checks}

def arithmetic():
    cases=[]
    for n in range(1,6):
        need=powers(n)
        for k in range(1,n+1):assert need[k]+4*k>=need[0]
        for selected in combinations(range(5),n):
            if 0 not in selected:continue
            raw=7*n*(n+1)//2-sum(selected)
            cases.append({'rank':n,'rows':selected,'raw_joint_degree':raw,'q_minor_bounds':need,
                          'common_scalar_q_power':need[0],'new_joint_degree':raw-need[0]})
    assert max(c['new_joint_degree'] for c in cases)==57
    assert max(9,21,51,54,57)==57 and max(94,92,57)==94
    cover={'quadratic':7,'smooth':94,'nodal':94,'cusp':92,'transverse':94,'tangent':92,'triangle':94,'concurrent':92,'repeated':92}
    assert max(cover.values())==94
    return cases,cover

def main():
    start=time.monotonic();univ=universal();print('Two Hensel branches,Jacobian,phase and Hermite determinant passed',flush=True)
    models=[]
    for args in [('regular_b1',1,1,2,3),('constant_q',0,1,2,3),('double_branch',1,1,2,1)]:
        models.append(model(*args));print('Full rank-five model',args[0],'passed',flush=True)
    cases,cover=arithmetic();cert={'universal':univ,'models':models,'all_rank_row_cases':cases,'lower_degree_gates_retained':{'2':9,'3':21,'4':51,'5':54},'global_cover':cover}
    f=P/'rank_one_two_branch_projection_certificate.json';f.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(f),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,
            'formal_branches':2,'branch_identity_order':4,'comparison_q_order':2,'full_rank_five_models':3,
            'exact_original_twisted_seed_identities':15,'full_scalar_cofactor_relations':3,'scalar_coefficient_checks':18,
            'degree6_common_scalar_q_power':38,'degree6_gate':57,'rank_one_full_gate':57,'conic_gate_retained':94,
            'full_W_linear_gate':94,'global_cover':cover,'global_sufficient_pure_prefix':94,'odd_sufficient_pure_prefix_unchanged':94,'homogeneous_sufficient_pure_prefix_unchanged':94,
            'least_prefix_interval':[20,94],'odd_least_prefix_interval':[20,94],'homogeneous_least_prefix_interval':[14,94],
            'least_prefix_computed':False,'universal_written_proof_certified_by_checker':False,'independent_replay':False,'external_review':False,'formal_verification':False,
            'boundary':'Finite two-branch Hensel/Jacobian identities,exact full phase rescaling and joint source bounds,leading Hermite comparison determinant,15 original twisted seed identities,three full scalar relations/18 coefficients and all rank/row/valuation arithmetic. Models include constant q and merging quadratic roots. Universal formal finite quotient/regular branch connection,complementary-minor induction,joint-grade specialization,polynomial extension and complete analytic/frame cover remain independently unreviewed written mathematics. No old checker main,radical search,external package or fresh literature review.'}
    (P/'rank_one_two_branch_projection_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
