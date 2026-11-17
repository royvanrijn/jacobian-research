"""Exact finite-seed controls; no universal symbolic-power proof from samples."""
from pathlib import Path
import hashlib,itertools,json,platform
import sympy as s


def apply_symbol(symbol,polynomial,xi,x):
    return s.expand(sum(c*s.diff(polynomial,x,j) for (j,),c in s.Poly(symbol,xi).terms()))


def run():
    a0,a1,a2,b0,b1,b2,x,xi=s.symbols('a0 a1 a2 b0 b1 b2 X xi')
    variables=(a0,a1,a2,b0,b1,b2)
    p=b0+b1*x+b2*x*x;lam=a0+a1*xi+a2*xi*xi
    first=apply_symbol(lam,p,xi,x);gens=s.Poly(first,x).all_coeffs()
    gb=s.groebner(gens,*variables)
    second=apply_symbol(lam**2,p*p,xi,x)
    remainders=[s.expand(gb.reduce(c)[1]) for c in s.Poly(second,x).all_coeffs()]
    remainder=a1*a1*b1*b1-4*a2*a2*b2*b2
    assert remainders==[0,0,0,0,remainder] and remainder!=0
    assert gb.reduce(s.expand(remainder**2))[1]==0

    # The entire affine seed radical and its exact nilpotence index.
    af=(a0,a1,b0,b1);seed=(a0*b1,a0*b0+a1*b1)
    ga=s.groebner(seed,*af);rad=(a0*b0,a0*b1,a1*b1)
    assert any(ga.reduce(f)[1]!=0 for f in rad)
    products=[s.expand(f*g) for f in rad for g in rad]
    assert all(ga.reduce(f)[1]==0 for f in products)
    components=[(b0,b1),(a0,a1),(a0,b1)]
    # Monomial intersection has precisely the displayed squarefree generators.
    intersection=[s.Integer(1)]
    for comp in components:
        intersection=[s.lcm(f,g) for f in intersection for g in comp]
        intersection=list(set(intersection))
    assert list(s.groebner(intersection,*af))==list(s.groebner(rad,*af))

    # General univariate radical formula: intersections of coordinate primes.
    monomial_controls=0
    for degree in range(5):
        intersections={0}
        for split in range(-1,degree+1):
            component=[1<<i for i in range(split+1)]+[1<<(degree+1+j) for j in range(split+1,degree+1)]
            proposed={f|g for f in intersections for g in component}
            minimal=[]
            for f in sorted(proposed,key=lambda v:(v.bit_count(),v)):
                if not any((f&g)==g for g in minimal):minimal.append(f)
            intersections=set(minimal)
        expected={(1<<i)|(1<<(degree+1+j)) for i in range(degree+1) for j in range(i,degree+1)}
        assert intersections==expected;monomial_controls+=1

    # Polynomial inputs only inspect a finite symbol jet, even for high-order tails.
    truncations=[]
    for degree,delta,m in [(1,2,2),(2,1,2),(3,2,1)]:
        cutoff=degree*m+delta
        phase=(x+1)**degree;marker=x**delta+1
        low=xi+2*xi*xi
        high=low+3*xi**(cutoff+1)+5*xi**(cutoff+4)
        left=apply_symbol(high**m,marker*phase**m,xi,x)
        right=apply_symbol(low**m,marker*phase**m,xi,x)
        assert s.expand(left-right)==0
        truncations.append(dict(D=degree,delta=delta,m=m,symbol_jet_bound=cutoff))

    # Univariate sharp field slope one in several degrees; not transferred to rings.
    slopes=[]
    for degree in range(1,5):
        for delta in range(1,4):
            symbol=xi**(degree+1);phase=x**degree
            boundary=apply_symbol(symbol**delta,x**delta*phase**delta,xi,x)
            assert boundary==s.factorial((degree+1)*delta)
            after=apply_symbol(symbol**(delta+1),x**delta*phase**(delta+1),xi,x)
            assert after==0
            slopes.append(dict(D=degree,delta=delta,boundary=str(boundary)))
    bounds=[]
    for n,nu in [(1,1),(1,2),(2,1),(2,3)]:
        C=3*(n+1);h=2*s.binomial(n+2,2);E=(2*C+1)*(h*nu-1)
        K=max(n+1,int(E));N=C+int(E)
        assert K>=n+1 and K>=E and N==C+E
        bounds.append(dict(n=n,nu=nu,h=int(h),E=int(E),K=K,N=N))
    return dict(completed=True,python_version=platform.python_version(),sympy_version=s.__version__,
       ordinary_field_gate_counterexample=dict(generators=[str(f) for f in gens],groebner_basis=[str(f) for f in gb],
           second_output_normal_forms=[str(f) for f in remainders],nonzero_constant_normal_form=str(remainder),its_square_reduces_to_zero=True),
       affine_seed_radical=[str(f) for f in rad],affine_seed_nilradical_index=2,
       affine_seed_square_products_checked=len(products),affine_seed_component_intersection_verified=True,
       univariate_seed_radical_monomial_intersection_controls=monomial_controls,
       finite_symbol_truncation_controls=truncations,univariate_sharp_field_controls=slopes,
       quadratic_bound_formula_controls=bounds,
       universal_proof_inferred_from_samples=False,infinite_coefficient_ring_assumed_Noetherian=False,
       general_seed_radical_or_finite_jet_producer_executed=False,
       universal_quadratic_or_all_degree_univariate_ordinary_ideals_computed=False,
       numerical_general_nilradical_indices_or_ring_cutoffs_computed=False,
       least_formal_mixed_cutoff_algorithm_claimed=False,sharp_field_slopes_transferred_to_rings=False,
       nonpolynomial_phases_or_multipliers_covered=False,independent_review=False,
       formal_verification=False,novelty_claim=False)


if __name__=='__main__':
    path=Path(__file__).resolve();result=run();result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('formal_operator_finite_seed_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
