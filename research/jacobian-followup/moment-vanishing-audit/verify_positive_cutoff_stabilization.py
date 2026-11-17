"""Small exact controls for written stabilization and coefficient-power proofs."""
from pathlib import Path
import hashlib,itertools,json,platform
import sympy as s


def apply_symbol(symbol,f,xi,x):
    return s.expand(sum(c*s.diff(f,x,i) for (i,),c in s.Poly(symbol,xi).terms()))


def run():
    a,b,c,d,e,f,x,xi,t=s.symbols('a b c d e f X xi t')
    av=(a,b,c);bv=(d,e,f);ops=(1,xi,xi**2);pbs=(1,x,x**2)
    lam=a+b*xi+c*xi**2;p=d+e*x+f*x*x
    identities=0;chains=0
    for m in range(1,4):
        Q=x+1;out=apply_symbol(lam**m,Q*p**m,xi,x)
        for opinds,phaseinds in [((),()),((0,),()),((1,2),()),((),(2,)),((),(0,2)),((1,),(2,)),((2,),(0,1)),((0,1),(2,))]:
            u=len(opinds);v=len(phaseinds);k=max(u,v)
            observed=out
            for i in opinds:observed=s.diff(observed,av[i])
            for i in phaseinds:observed=s.diff(observed,bv[i])
            if u>m or v>m:expected=0
            else:
                marker=Q*s.prod(pbs[i] for i in phaseinds)*p**(k-v)
                inner=apply_symbol(lam**(m-k),marker*p**(m-k),xi,x)
                outer=s.prod(ops[i] for i in opinds)*lam**(k-u)
                expected=s.prod(m-j for j in range(u))*s.prod(m-j for j in range(v))*apply_symbol(outer,inner,xi,x)
            assert s.expand(observed-expected)==0;identities+=1
        # A nonlinear parameterization uses coefficient derivatives through chain rule.
        change={a:t*t,b:t+1,c:t**3,d:t,e:t*t+2,f:t+3}
        direct=s.diff(out.subs(change),t)
        reconstructed=sum(s.diff(out,vv).subs(change)*s.diff(change[vv],t) for vv in av+bv)
        assert s.expand(direct-reconstructed)==0;chains+=1
    # Omitting P^(k-v) when an operator coefficient is differentiated is false.
    true=s.diff(apply_symbol((a*xi)**2,x*(b*x*x)**2,xi,x),a)
    wrong=2*apply_symbol(xi*(a*xi),x*(b*x*x),xi,x)
    assert s.expand(true-wrong)!=0

    # A complete universal one-dimensional linear family: pure ideal (a^2).
    # m>=2 outputs are divisible by a^2; at N>=1 no positive marker permits m=1.
    def E(expr):
        return s.expand(sum(co*(s.factorial2(i-1) if i else 1) for (i,),co in s.Poly(expr,x).terms() if i%2==0))
    pure=[E((a*x)**m) for m in range(1,7)]
    gp=s.groebner(pure,a);assert list(gp)==[a*a]
    bounded=[]
    for N in range(4):
        vals=[E(x**delta*(a*x)**m) for delta in range(5) for m in range(1,8) if m>N*delta]
        g=s.groebner(vals,a);expected=[a] if N==0 else [a*a]
        assert list(g)==expected;bounded.append([str(v) for v in g])
    assert s.rem(a,a*a,a)!=0

    # Failure-line ideals remain strict, but an idempotent restriction disconnects them.
    idem=s.groebner([e*e-e],e);scales=[]
    for N in range(6):
        assert s.rem(e**(N+1),e**(N+2),e)==e**(N+1)
        assert idem.reduce(e**(N+1))[1]==e
        scales.append(N+1)
    # A nonreduced safe component survives in the clopen stabilized representative.
    # Synthetic topology control, not asserted to be a new moment family.
    nu=3;quot=s.groebner([e*e-e,(1-e)*t**nu],e,t)
    def equal_mod(gens,other):
        left=s.groebner([e*e-e,(1-e)*t**nu]+list(gens),e,t)
        right=s.groebner([e*e-e,(1-e)*t**nu]+list(other),e,t)
        return list(left)==list(right)
    assert equal_mod([e,t**nu],[e])
    assert not equal_mod([e,t**(nu-1)],[e])
    assert s.groebner([e,t**nu],e,t).reduce(t)[1]==t
    bounds=[]
    for length in range(1,5):
        ternary=3+(3*3+1)*(20*length-1)
        assert ternary==200*length-7
        bounds.append(dict(nu=length,ternary_sufficient_slope=ternary))
    return dict(completed=True,python_version=platform.python_version(),sympy_version=s.__version__,
      direct_GVC_power_alignment_identities=identities,nonlinear_parameter_chain_controls=chains,
      missing_power_insertion_mutation_detected=True,linear_Gaussian_pure_ideal=['a**2'],
      bounded_linear_Gaussian_cutoff_ideals=bounded,linear_nilpotent_first_marked_moment_survives=True,
      failure_scale_strict_power_controls=scales,idempotent_failure_restriction_stabilizes=True,
      nonreduced_clopen_component_control=True,positive_ternary_bound_controls=bounds,
      universal_proof_inferred_from_samples=False,general_ideal_or_stabilization_solver_executed=False,
      universal_ternary_or_quadratic_GVC_ordinary_ideal_computed=False,
      numerical_universal_nilradical_index_or_ring_cutoff_computed=False,
      sharp_field_slope_transferred_to_nonreduced_targets=False,
      infinite_formal_symbol_universal_ring_covered=False,protected_frontiers_improved=False,
      independent_review=False,formal_verification=False,novelty_claim=False)


if __name__=='__main__':
    path=Path(__file__).resolve();result=run();result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('positive_cutoff_stabilization_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
