"""Small exact controls; universal lifting and imported theorems are written proofs."""
from pathlib import Path
import hashlib,itertools,json,platform
import sympy as s


def wick(f,z,w):
    return s.expand(sum(c*s.factorial(i) for (i,j),c in s.Poly(s.expand(f),z,w).terms() if i==j))


def run():
    a,b,c,z,w,e,x=s.symbols('a b c Z W e X')
    cs=(a,b,c);basis=(z,w,z*w);p=sum(u*v for u,v in zip(cs,basis))
    identities=0
    for m in range(1,5):
        for marker in (s.Integer(1),z,z*w*w):
            mu=wick(marker*p**m,z,w)
            for k in range(3):
                for inds in itertools.combinations_with_replacement(range(3),k):
                    derivative=mu;insert=marker
                    for i in inds:derivative=s.diff(derivative,cs[i]);insert*=basis[i]
                    expected=0 if k>m else s.prod(m-j for j in range(k))*wick(insert*p**(m-k),z,w)
                    assert s.expand(derivative-expected)==0
                    identities+=1

    # Differential order two does not imply ordinary square membership.
    generators=(a*b,a*c,b*c)
    gj=s.groebner(generators,a,b,c);gj2=s.groebner([u*v for u in generators for v in generators],a,b,c)
    f=a*b*c
    assert all(gj.reduce(h)[1]==0 for h in [f]+[s.diff(f,v) for v in cs])
    assert gj2.reduce(f)[1]==f

    # Noetherian nil-ideal bound for three independently nilpotent generators.
    exponents=(2,3,2);nilbound=1+sum(k-1 for k in exponents)
    gn=s.groebner([v**k for v,k in zip(cs,exponents)],*cs)
    nilchecks=0
    for aa in itertools.product(range(nilbound+1),repeat=3):
        if sum(aa)!=nilbound:continue
        assert gn.reduce(s.prod(v**k for v,k in zip(cs,aa)))[1]==0;nilchecks+=1
    assert gn.reduce(a*b*b*c)[1]!=0

    # Pure-null nilpotent lift aX over Q[a]/a^2: its marked m=1 moment survives.
    def normal_gaussian(f):
        out=0
        for (i,),coef in s.Poly(s.expand(f),x).terms():
            if i%2==0:out+=coef*(s.factorial2(i-1) if i else 1)
        return s.expand(out)
    liftchecks=0
    for m in range(1,7):
        assert s.rem(normal_gaussian((a*x)**m),a*a,a)==0
        for delta in range(1,5):
            val=s.rem(normal_gaussian(x**delta*(a*x)**m),a*a,a)
            if m>delta:assert val==0
            liftchecks+=1
    assert normal_gaussian(x*a*x)==a

    # The constant pure equation cannot be omitted from the completion quotient.
    assert s.groebner([e**m for m in range(1,6)],e).reduce(e)[1]==0
    assert s.rem(e,e*e,e)==e
    bounds=[]
    for C,d,q,L in [(1,1,2,1),(1,1,2,2),(3,4,35,3),(2,3,20,4)]:
        B=C+(C*d+1)*(q*L-1)
        for delta in range(1,6):
            m=B*delta+1
            assert m>C*delta+(C*d+1)*(q*L-1)
        bounds.append(dict(C=C,d=d,q=q,L=L,sufficient_integer_slope=B))
    return dict(completed=True,python_version=platform.python_version(),sympy_version=s.__version__,
        coefficient_derivative_identities=identities,differential_ordinary_power_mutation_detected=True,
        nil_ideal_bound=nilbound,nil_ideal_monomial_controls=nilchecks,
        pure_null_nilpotent_lift_controls=liftchecks,marked_first_moment='a',
        pure_quotient_omission_detected=True,affine_to_linear_bound_controls=bounds,
        universal_proof_inferred_from_samples=False,imported_symbolic_theorems_reproved=False,
        general_safe_ideal_or_completion_solver_executed=False,numeric_universal_field_cutoff_computed=False,
        individual_cutoff_ideals_identified_with_safe_ideal_powers=False,
        ordinary_cutoff_stabilization_claimed=False,all_nonreduced_uniformity_without_thickness_claimed=False,
        independent_review=False,formal_verification=False,novelty_claim=False)


if __name__=='__main__':
    path=Path(__file__).resolve();result=run()
    result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('nilpotent_safe_completion_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
