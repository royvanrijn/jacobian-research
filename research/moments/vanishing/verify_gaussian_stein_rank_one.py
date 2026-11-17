"""Focused exact controls for Stein action and faithful rank-one presentations."""
from pathlib import Path
import hashlib,json
import sympy as s
import gaussian_stein_rank_one as G
import gaussian_generic_cutoff_producer as W

p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()


def main():
    x,y,z,w,T,e=s.symbols('x y z w T e');t=G.t
    pair=s.Matrix([[0,1],[1,0]])
    cases=[(x**3+y,(x,y),pair,[1,x,y,x*y,x*x*y*y]),
           (T*T/2,(T,),s.eye(1),[1,T,T*T,T**3,T**4]),
           (w*(1+z)-(1+s.Rational(3,2)*z+z*z/2)*T*T,(z,w,T),s.diag(pair,s.eye(1)),
            [1,z,w,T*T,z*w])]
    actions=[];identities=0
    for P,X,C,markers in cases:
        wick=W.wick_engine(C)
        expectation=lambda A:sum(c*wick(alpha) for alpha,c in s.Poly(s.expand(A),*X).terms())
        for A in markers:
            A=s.sympify(A);gamma=G.stein(P,A,X,C);inv=G.inverse_number(A,X,C)
            vector=C*s.Matrix([s.diff(inv,v) for v in X])
            divergence=(s.Matrix(X).T*C.inv()*vector)[0]-sum(s.diff(vector[i],X[i]) for i in range(len(X)))
            assert s.expand(divergence-A+expectation(A))==0
            for m in range(1,6):
                assert s.expand(expectation(A*P**m)-expectation(A)*expectation(P**m)-
                                m*expectation(gamma*P**(m-1)))==0
                identities+=1
            actions.append(dict(phase=str(P),variables=[str(v) for v in X],covariance=[[str(c) for c in row] for row in C.tolist()],
                                marker=str(A),Stein_image=str(gamma),divergence_identity=True))
    # An ordinary nilpotent control, with no probability or field assumption.
    P=e*T*T;A=T*T;gamma=G.stein(P,A,(T,),s.eye(1));wick=W.wick_engine(s.eye(1))
    E=lambda B:sum(c*wick(alpha) for alpha,c in s.Poly(s.expand(B),T).terms())
    assert s.expand(E(A*P)-E(A)*E(P)-E(gamma))==0
    nil=dict(ring='Q[e]/(e^2)',phase=str(P),Stein_image=str(gamma),unit_period='1+e*t')
    models=[]
    data=[('pure safe',0,[t*t],True,True),('pure pole',0,[1/(1-t)],False,True),
          ('exponential nilpotent Stein',2*t,[t*t],False,False),
          ('Gaussian quadratic',1/(2*(1-t)),[1/(1-t)],False,False),
          ('abstract nonpure safe seed',2/(t-1),[1],True,False),
          ('Gaussian polynomial seed with hidden pole',2/(t-1),[1/(1-t)**3],False,False),
          ('mixed finite and infinite poles',2/(t-1)+3/(t-2),[1/(t-1)**3],False,False),
          ('abstract regular transcendental',1/(1+t*t),[1],False,False),
          ('irregular finite pole',1/(1-t)**2,[1],False,False)]
    for name,r,ratios,safe,pure in data:
        model=G.classify_rank_one(r,ratios)
        assert model['safe']==safe and model['pure']==pure,name
        assert not model['Gaussian_realization_inferred']
        models.append(dict(name=name,model=model))
    assert s.expand(s.sympify(models[4]['model']['finite_pole_polynomial'])-(1-t)**2)==0
    assert models[7]['model']['globally_regular'] and not models[8]['model']['globally_regular']
    matrices=[]
    for f in [1-t,(1-t)**2,1+t+t*t]:
        matrix=G.seed_matrix(f);matrices.append(dict(seed=str(f),**matrix))
        m=matrix['dimension'];V=s.Matrix([[s.sympify(a) for a in row] for row in matrix['matrix']])
        if s.expand(f-(1-t)**2)==0:assert (V-s.eye(m))**2==s.zeros(m) and V!=s.eye(m)
    poles=[]
    for j in [1,2,4]:
        value=(1-t)**(-j)
        for k in range(j-1):value=s.cancel(G.shift(value)-value)
        assert value!=0 and s.cancel(G.shift(value)-value)==0
        poles.append(dict(pole_order=j,Jordan_chain_length=j,eigenvalue=1))
    # Published-v2 Gaussian right-ideal telescoping index correction.
    source_controls=[]
    for k,m in [(0,1),(1,2),(2,3)]:
        operator={}
        for r in range(m):
            h=s.hermite_prob(k+r,x);j=m-1-r
            operator[j+1]=operator.get(j+1,0)+h
            operator[j]=operator.get(j,0)+s.diff(h,x)-x*h
        operator={j:s.expand(a) for j,a in operator.items() if s.expand(a)!=0}
        assert operator=={m:s.hermite_prob(k,x),0:-s.hermite_prob(k+m,x)}
        source_controls.append(dict(k=k,m=m,corrected_telescoping_identity=True))
    assert s.diff(x,x)-x*x != -x # Printed m1,k0 RHS uses H1 instead of H0.
    rejects=[]
    for fn in [lambda:G.stein(T,T,(T,),s.zeros(1)),lambda:G.stein(1/T,T,(T,),s.eye(1)),
               lambda:G.seed_matrix(2+t),lambda:G.classify_rank_one(1/t,[1]),
               lambda:G.shift(1/t)]:
        try:fn()
        except (ValueError,s.PolynomialError):rejects.append(True)
        else:raise AssertionError('Invalid input accepted')
    packet=dict(original_Stein_controls=actions,original_IBP_identities=identities,nilpotent_control=nil,
                rank_one_models=models,finite_seed_matrices=matrices,Jordan_controls=poles,
                source_index_controls=source_controls,rejections=rejects,
                universal_claims_use_written_proofs=True,actual_rank_inferred_from_finite_jets=False,
                nonpure_safe_Gaussian_counterexample_claimed=False,regular_transcendental_Gaussian_realization_claimed=False)
    cert=p/'gaussian_stein_rank_one_certificate.json';cert.write_text(json.dumps(packet,indent=2)+'\n')
    imports=['RATIONAL_CORES_AND_BOUNDARY_SUBSPACES.md','TIME_ZERO_SATURATION_AND_PRINCIPAL_PART_RECONSTRUCTION.md',
             'POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md','FILTERED_SAFETY_WITHOUT_A_PURITY_PREMISE.md',
             'BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md','gaussian_generic_cutoff_producer.py']
    runtime=dict(python=__import__('sys').version,sympy=s.__version__,sealed_imports_sha256={n:sha(p/n) for n in imports})
    run=p/'gaussian_stein_rank_one_runtime.json';run.write_text(json.dumps(runtime,indent=2)+'\n')
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'STEIN_ACTION_AND_RANK_ONE_LATTICES.md'),
        producer_sha256=sha(p/'gaussian_stein_rank_one.py'),checker_sha256=sha(Path(__file__)),
        certificate_sha256=sha(cert),runtime_sha256=sha(run),
        exact_controls=dict(original_Stein_markers=len(actions),original_IBP_identities=identities,
                            rank_one_models=len(models),finite_seed_matrices=len(matrices),
                            Jordan_controls=len(poles),source_index_controls=len(source_controls),rejections=len(rejects)),
        independent_replay=False,formal_verification=False,external_review=False,
        universal_claims_use_written_proofs=True,nonpure_safe_Gaussian_counterexample_claimed=False)
    (p/'gaussian_stein_rank_one_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
