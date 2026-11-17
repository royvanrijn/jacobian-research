"""Exact tail-cutoff arithmetic and finite controls for the written proof."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,platform
import sympy as s


def tail_cutoff(D,C,tail_order):
    if any(not isinstance(v,int) or isinstance(v,bool) or v<1 for v in [D,C,tail_order]):
        raise ValueError('Positive integer degree/cutoff/order data required')
    a=C*D+1;denominator=tail_order-D*a
    if denominator<=0:raise ValueError('Tail order must be strictly greater than D(CD+1)')
    ratio=Fraction(tail_order*C+a,denominator)
    return dict(phase_degree=D,base_cutoff=C,tail_order=tail_order,
        strict_order_threshold=D*a,ratio_numerator=ratio.numerator,ratio_denominator=ratio.denominator,
        sufficient_integer_cutoff=(ratio.numerator+ratio.denominator-1)//ratio.denominator,
        baseline_cutoff_is_a_required_proof_input=True,pure_outputs_unchanged_under_hypotheses=True)


def finite_core_bound(D,C):
    if any(not isinstance(v,int) or isinstance(v,bool) or v<1 for v in [D,C]):raise ValueError('Positive integers required')
    core_cutoff=2*C+1;r=D*(core_cutoff*D+1)
    assert r+1>=2*D*(C*D+1)
    assert r+1>D*(core_cutoff*D+1)
    return dict(phase_degree=D,formal_pair_cutoff=C,jet_order=r,truncation_sufficient_cutoff=core_cutoff,
        tail_transfer=tail_cutoff(D,core_cutoff,r+1))


def apply_symbol(symbol,f,derivatives,variables):
    answer=0
    for powers,coefficient in s.Poly(s.expand(symbol),*derivatives).terms():
        term=f
        for x,k in zip(variables,powers):term=s.diff(term,x,k)
        answer+=coefficient*term
    return s.expand(answer)


def run():
    x,y,dx,dy,e=s.symbols('X Y dx dy e');vs=(x,y);ds=(dx,dy)
    p=x*y**3;base=dx**2;bad_tail=dy**5
    assert apply_symbol(base+bad_tail,p,ds,vs)==0
    second=apply_symbol((base+bad_tail)**2,p*p,ds,vs);assert second==2880*y
    nil=apply_symbol((base+e*bad_tail)**2,p*p,ds,vs)
    assert s.rem(nil,e*e,e)==2880*e*y and s.rem(nil,e*e,e)!=0
    try:tail_cutoff(4,1,5)
    except ValueError:pass
    else:raise AssertionError('Higher than phase degree was substituted for the actual threshold')
    # Every strict boundary must be rejected.
    rejected=0;ineq=0;cores=[];rows=[]
    for D in range(1,4):
        for C in range(1,4):
            a=C*D+1;threshold=D*a
            try:tail_cutoff(D,C,threshold)
            except ValueError:rejected+=1
            else:raise AssertionError('Strict order boundary lost')
            cores.append(finite_core_bound(D,C))
            for order in [threshold+1,2*threshold,2*threshold+1]:
                row=tail_cutoff(D,C,order);N=row['sufficient_integer_cutoff']
                if order>=2*threshold:assert N<=2*C+1
                for delta in range(4):
                    m=N*delta+1
                    low=max(0,(m-C*delta+a-1)//a)
                    high=min(m,(delta+D*m)//order)
                    assert low>high;ineq+=1
                rows.append(row)
    # Binomial identities with nilpotent/transcendental tail coefficients retained.
    identities=0
    for m in range(1,4):
        for marker in [s.Integer(1),x+y,x*y]:
            observed=apply_symbol((base+e*bad_tail)**m,marker*p**m,ds,vs)
            expected=0
            for l in range(m+1):
                inner=apply_symbol(base**(m-l),marker*p**l*p**(m-l),ds,vs)
                expected+=s.binomial(m,l)*apply_symbol((e*bad_tail)**l,inner,ds,vs)
            assert s.expand(observed-expected)==0;identities+=1
    # Admitted tail: pure outputs agree in small controls, marked outputs can change.
    admitted=e*dy**21
    for m in range(1,4):assert apply_symbol((base+admitted)**m,p**m,ds,vs)==0
    marked=apply_symbol(base+admitted,y**21*p,ds,vs)
    assert marked==e*s.factorial(24)/s.factorial(3)*x*y**3 and marked!=0
    assert apply_symbol(base,y**21*p,ds,vs)==0
    return dict(completed=True,python_version=platform.python_version(),sympy_version=s.__version__,
       insufficient_tail_order_second_output=str(second),nilpotent_tail_second_output=str(s.rem(nil,e*e,e)),
       merely_above_phase_degree_mutation_detected=True,strict_boundary_rejections=rejected,
       cutoff_rows=rows,necessary_inequality_controls=ineq,finite_core_controls=cores,
       binomial_identity_controls=identities,admitted_tail_pure_controls=3,
       admitted_tail_changes_marked_output=True,admitted_tail_marked_output=str(marked),
       tail_bound_calculator_implemented=True,finite_core_bound_calculator_implemented=True,
       universal_proof_inferred_from_samples=False,baseline_linear_cutoff_verified_by_calculator=False,
       full_formal_eventuality_implies_linearity_claimed=False,general_formal_safety_decision_claimed=False,
       finite_operator_oracle_or_semidecision_implemented=False,ordinary_cutoff_ideal_equality_claimed=False,
       phase_perturbations_or_variable_coefficients_covered=False,independent_review=False,
       formal_verification=False,novelty_claim=False)


if __name__=='__main__':
    path=Path(__file__).resolve();result=run();result['checker_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name('operator_tail_safety_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(completed=result['completed'],checker_sha256=result['checker_sha256'],
        necessary_inequality_controls=result['necessary_inequality_controls'],
        binomial_identity_controls=result['binomial_identity_controls'],
        insufficient_tail_order_second_output=result['insufficient_tail_order_second_output']),indent=2))
