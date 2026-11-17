"""Explicit sufficient-bound calculator and small exact source/scope controls."""
from pathlib import Path
import hashlib,itertools,json,math,platform
import sympy as s


def repair_degree(degree):
    if not isinstance(degree,int) or isinstance(degree,bool) or degree<1:
        raise ValueError('Nonzero homogeneous generators of positive degree required')
    return 3 if degree==1 else 4 if degree==2 else degree


def kollar_exponent(degrees,coefficient_variables):
    """Only the published homogeneous case p=q-1>=2, all degrees>=3."""
    ds=sorted(degrees,reverse=True);p=coefficient_variables-1
    if p<2 or not ds or any(d<3 for d in ds):
        raise ValueError('Requires q>=3 and all homogeneous generator degrees>=3')
    return math.prod(ds) if len(ds)<=p else math.prod(ds[:p-1])*ds[-1]


def _positive_ints(*values):
    if any(not isinstance(v,int) or isinstance(v,bool) or v<1 for v in values):
        raise ValueError('Positive integers required')


def _descriptor(*,n,D,k0,C,q,h,base,exponent,model,max_bits=4096):
    _positive_ints(n,D,k0,C,q,h,base,exponent,max_bits)
    nu=pow(base,exponent) if base.bit_length()*exponent<=max_bits else None
    E=(C*D+1)*(h*nu-1) if nu is not None else None
    K=max(k0,E) if E is not None else None
    N=C+E if E is not None else None
    if E is not None and max(E.bit_length(),K.bit_length(),N.bit_length())>max_bits:
        nu=E=K=N=None
    power=f'{base}^{exponent}'
    Eformula=f'({C}*{D}+1)*({h}*({power})-1)'
    return dict(model=model,n=n,phase_degree=D,field_pure_gate=k0,field_mixed_slope=C,
       seed_coefficient_variables=q,component_codimension_bound=h,
       nilradical_exponent_upper_bound=dict(base=base,exponent=exponent,value=nu),
       sufficient_E=E,sufficient_ordinary_pure_prefix=K,sufficient_mixed_slope=N,
       sufficient_symbol_jet=D*K if K is not None else None,
       exact_formulas=dict(E=Eformula,K=f'max({k0},{Eformula})',N=f'{C}+{Eformula}',symbol_jet=f'{D}*max({k0},{Eformula})'),
       integer_evaluation_bit_budget=max_bits,all_Q_algebras=True,formal_constant_coefficient_symbols=model!='ternary_Gaussian_cubics',
       least_values_computed=False,ordinary_ideal_computed=False)


def conditional_gvc_seed_bounds(n,D,k0,C,h=None,max_bits=4096):
    """Caller supplies the proved order-independent field gate/cutoff hypotheses."""
    _positive_ints(n,D,k0,C)
    q=math.comb(n+D,D)+math.comb(n+D*k0,D*k0)
    count=sum(math.comb(n+D*j,n) for j in range(1,k0+1))
    result=_descriptor(n=n,D=D,k0=k0,C=C,q=q,h=q if h is None else h,
       base=max(4,2*k0),exponent=min(count,q-1),model='conditional_GVC_field_seed',max_bits=max_bits)
    result.update(seed_generator_count_upper_bound=count,
       requires_proved_field_hypotheses=True,requires_proved_codimension_bound=h is not None)
    return result


def quadratic_formal_bounds(n,max_bits=4096):
    _positive_ints(n)
    q=math.comb(n+2,2)+math.comb(3*n+2,n)
    return _descriptor(n=n,D=2,k0=n+1,C=3*(n+1),q=q,h=2*math.comb(n+2,2),
       base=2*(n+1),exponent=q-1,model='quadratic_GVC_formal_symbols',max_bits=max_bits)


def univariate_formal_bounds(D,max_bits=4096):
    _positive_ints(D)
    return _descriptor(n=1,D=D,k0=1,C=1,q=2*(D+1),h=D+1,
       base=4,exponent=D+1,model='univariate_GVC_formal_symbols',max_bits=max_bits)


def ternary_gaussian_cubic_bounds(max_bits=4096):
    result=_descriptor(n=3,D=3,k0=94,C=3,q=20,h=20,
       base=94,exponent=19,model='ternary_Gaussian_cubics',max_bits=max_bits)
    result['sufficient_symbol_jet']=None;result['exact_formulas']['symbol_jet']=None
    return result


def run():
    # Exact homogeneous degree checks for small pure seed polynomials.
    x,xi=s.symbols('X xi');av=s.symbols('a0:7');bv=s.symbols('b0:3')
    symbol=sum(v*xi**i for i,v in enumerate(av));phase=sum(v*x**i for i,v in enumerate(bv))
    degree_controls=[]
    for m in range(1,4):
        out=0
        for (j,),co in s.Poly(s.expand(symbol**m),xi).terms():
            if j<=2*m:out+=co*s.diff(phase**m,x,j)
        coefficients=s.Poly(s.expand(out),x).all_coeffs()
        assert all(set(sum(monomial) for monomial,co in s.Poly(c,*av,*bv).terms())=={2*m} for c in coefficients)
        degree_controls.append(dict(m=m,coefficient_count=len(coefficients),raw_degree=2*m,repaired_degree=repair_degree(2*m)))
    assert [repair_degree(d) for d in [1,2,3,4,94]]==[3,4,3,4,94]
    # The raw degree-two case and the p=1 source exception are rejected.
    for degrees,q in [([2,3],4),([3,3],2)]:
        try:kollar_exponent(degrees,q)
        except ValueError:pass
        else:raise AssertionError('Source restriction was silently lost')

    # Known full ordinary exponents for coordinate-power ideals.
    monomial=[]
    for degrees in [(3,3,3),(4,3,3),(4,4,4,3)]:
        q=len(degrees);least=1+sum(d-1 for d in degrees)
        top=tuple(d-1 for d in degrees);assert sum(top)==least-1
        for powers in itertools.product(range(least+1),repeat=q):
            if sum(powers)==least:assert any(v>=d for v,d in zip(powers,degrees))
        bound=kollar_exponent(degrees,q);assert bound>=least
        monomial.append(dict(degrees=list(degrees),exact_ordinary_exponent=least,published_bound=bound))
    # Empty projective zero set does not give ordinary containment at the vertex.
    cone_degrees=(3,3);cone_last=tuple(d-1 for d in cone_degrees)
    assert sum(cone_last)==4 and not any(v>=d for v,d in zip(cone_last,cone_degrees))
    assert all(any(v>=d for v,d in zip((i,5-i),cone_degrees)) for i in range(6))
    assert not any(v>=d for v,d in zip((2,1),cone_degrees))

    quadratic=[quadratic_formal_bounds(n) for n in [1,2,3]]
    assert [b['seed_coefficient_variables'] for b in quadratic]==[8,34,175]
    assert quadratic[0]['nilradical_exponent_upper_bound']['value']==16384
    assert quadratic[0]['sufficient_ordinary_pure_prefix']==1277939
    assert quadratic[0]['sufficient_mixed_slope']==1277945
    univariate=[univariate_formal_bounds(D) for D in range(1,5)]
    assert [b['sufficient_mixed_slope'] for b in univariate]==[63,574,4093,25596]
    generic=conditional_gvc_seed_bounds(1,2,2,6,h=6)
    assert generic['sufficient_mixed_slope']==quadratic[0]['sufficient_mixed_slope']
    ternary=ternary_gaussian_cubic_bounds();power=94**19
    assert ternary['sufficient_ordinary_pure_prefix']==200*power-10
    assert ternary['sufficient_mixed_slope']==200*power-7
    assert kollar_exponent([3,4]+list(range(3,95)),20)<=power
    large=quadratic_formal_bounds(10,max_bits=256)
    assert large['nilradical_exponent_upper_bound']['value'] is None
    assert large['sufficient_mixed_slope'] is None and 'max(' in large['exact_formulas']['K']
    return dict(completed=True,python_version=platform.python_version(),sympy_version=s.__version__,
       homogeneous_seed_degree_controls=degree_controls,degree_two_source_mutation_detected=True,
       projective_dimension_exception_detected=True,coordinate_power_ideal_controls=monomial,
       projective_sheaf_vs_cone_ordinary_exponent=5,quadratic_examples=quadratic,
       univariate_examples=univariate,ternary_cubic_example=ternary,
       conditional_seed_specialization_agrees=True,large_integer_guard_verified=True,
       bound_formula_producer_implemented=True,explicit_sufficient_general_bounds_produced=True,
       universal_proof_inferred_from_samples=False,exact_general_nilradical_indices_computed=False,
       least_general_pure_prefix_or_mixed_cutoff_computed=False,huge_ordinary_pure_ideals_expanded=False,
       full_seed_radical_or_ordinary_ideal_solver_implemented=False,
       recent_preprint_theorem_imported=False,sharp_field_slopes_transferred_to_rings=False,
       protected_field_prefix_or_operator_frontiers_improved=False,
       independent_review=False,formal_verification=False,novelty_claim=False)


if __name__=='__main__':
    path=Path(__file__).resolve();result=run();result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('effective_ring_bounds_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(completed=result['completed'],checker_sha256=result['checker_sha256'],
       quadratic_n1_slope=result['quadratic_examples'][0]['sufficient_mixed_slope'],
       univariate_slopes=[r['sufficient_mixed_slope'] for r in result['univariate_examples']],
       ternary_slope=result['ternary_cubic_example']['sufficient_mixed_slope'],
       large_integer_guard_verified=result['large_integer_guard_verified']),indent=2))
