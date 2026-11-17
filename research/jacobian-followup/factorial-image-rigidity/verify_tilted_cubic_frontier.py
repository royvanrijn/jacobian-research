"""Exact replay of prime-power tilted-cubic arithmetic and fold data.

The all-parameter proof is TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex.
Finite replays check bookkeeping and exact failed sufficient hypotheses.
"""
from pathlib import Path
from math import comb,factorial,gcd
import hashlib,json
import sympy as s
import tilted_cubic_frontier as f
from prime_power_windows import newton_residuals,factorial_valuation as vf
from newton_slope_coverage import slope_square_divisor_degrees
HERE=Path(__file__).resolve().parent
counts={}
def check(key,condition):
    assert condition,key
    counts[key]=counts.get(key,0)+1

# The uniform single-level inequality is proved by the cases in the note.
# These exhaustive small base replays include both special endpoints.
for b in [5,7,11,25,49]:
    for i in range(min(3,(b-1)//3)+1):
        for c in range(b):
            v=(0 if c==0 else 2+(2*c-1)//b) if i==0 else (i-1+2*c)//b+int(c>i)
            check('single_level_inequality',s.Rational(v)>=s.Rational(2*(c-i),b))

cases=[(1,5,1,0),(1,5,2,0),(1,5,2,1),(1,11,1,3),(2,7,1,1),
       (2,7,2,2),(3,11,1,0),(3,11,2,1),(4,13,1,0),(2,11,1,3)]
for a,p,r,i in cases:
    n=a*p**r+i;actual=newton_residuals(f.coefficients(n),p)
    expected=f.expected_shifted_edges(a,p,r,i)
    for row in actual:row.pop('residual_squarefree')
    check('shifted_polygon_and_full_residuals',actual==expected)
    check('independent_valuation_polygon',f.valuation_edges(n,p)==[{k:v for k,v in row.items() if k!='residual_coefficients'} for row in expected])
    P=p**r;S=(P-1)//(p-1)
    check('reduced_slope_denominator',gcd(2*S,P)==1)

# Exact normalized factorial units at positive and negative digit offsets.
for p,r in [(5,1),(5,2),(7,1),(7,2)]:
    P=p**r;S=(P-1)//(p-1)
    for c in range(1,p):
        for h in [0,1,p-2]:
            n=c*P+h;unit=(factorial(n)//p**vf(n,p))%p
            check('factorial_digit_positive_offset',unit==((-1)**(c*S)*factorial(c)*factorial(h))%p)
        n=c*P-1;unit=(factorial(n)//p**vf(n,p))%p
        check('factorial_digit_minus_one',unit==((-1)**(c*S)*factorial(c-1))%p)

for p,rmax in [(2,7),(3,4)]:
    for r in range(1,rmax+1):
        P=p**r;top=2*P-1 if p==2 else P-1
        check('small_prime_primitive_polygon',f.valuation_edges(P,p)==[{'left':[0,0],'right':[P,top],'e':P,'h':top}])
        for j in range(1,P+1):
            t=0;value=j
            while value%p==0:t+=1;value//=p
            actual=vf(P,p)-vf(j,p)-vf(P-j,p)+vf(P+2*j-1,p)-vf(P-1,p)
            predicted=2*r-t+vf(2*j-1,p)+(int(j>P//2) if p==2 else 0)
            check('small_prime_exact_coefficient_valuation',actual==predicted)
            if j<P:check('small_prime_interior_above_chord',actual>=(2*j if p==2 else j))
        cert=f.prime_power_certificate(p,r)
        check('small_prime_structural_certificate',cert['irreducible_over_Qp'] and not cert['all_starting_indices_claimed'])
check('unbounded_power_certificate_without_expansion',f.prime_power_certificate(5,1000000)['irreducible_over_Qp'])

for a,i in sorted(f.FAMILIES):
    patterns=[(u,v) for u in range(i//2+1) for v in range(a//2+1) if (u or v) and (u+v)%2==0]
    remaining=[(u,v) for u,v in patterns if not (2*u==i and 2*v==a)]
    check('all_r_parity_and_full_square_patterns',not remaining)
    p=int(s.nextprime(3*max(a,i)))
    check('all_r_family_scope',f.uniform_shift_family(a,p,i)['all_terminals_in_family_certified'])
    # The retained leading coefficient argument is universal by Bertrand;
    # these literals check that the correct factorial ratio is being used.
    n=a*p+i
    leading=factorial(3*n-1)//factorial(n-1)
    ell=int(s.nextprime((3*n-1)//2))
    check('leading_coefficient_odd_prime_valuation',ell<=3*n-1 and ell>n-1 and vf(3*n-1,ell)-vf(n-1,ell)==1)
check('unproved_shift_family_not_claimed',not f.uniform_shift_family(3,11,2)['all_terminals_in_family_certified'])

# Universal exact fold geometry and nondegeneracy, not a saddle asymptotic.
x,w=s.symbols('x w');F=w*x**3-3*w*x**2+x-1
Phi=s.log(x)+s.log(1+w*x*x)-x
check('saddle_polynomial_identity',s.cancel(s.diff(Phi,x)+F/(x*(1+w*x*x)))==0)
check('saddle_discriminant',s.factor(s.discriminant(F,x)+4*w*(1-9*w+27*w*w))==0)
for sign in [1,-1]:
    xc=(3+sign*s.I*s.sqrt(3))/2;wc=s.Rational(1,6)+sign*s.I*s.sqrt(3)/18
    at=lambda expr:s.simplify(s.expand_complex(expr.subs({x:xc,w:wc})))
    check('fold_zero',at(F)==0 and at(s.diff(F,x))==0)
    check('fold_phase_first_second_zero',at(s.diff(Phi,x))==0 and at(s.diff(Phi,x,2))==0)
    third=at(s.diff(Phi,x,3));cross=at(s.diff(Phi,x,w))
    check('fold_third_derivative',s.simplify(third+3*(xc-1)**2/xc**3)==0)
    check('fold_third_modulus',s.simplify(third*s.conjugate(third))==s.Rational(1,3))
    check('fold_cross_modulus',s.simplify(cross*s.conjugate(cross))==s.Rational(27,4))
    check('fold_third_saddle',s.simplify(F.subs({x:3-2*xc,w:wc}))==0)
    check('fold_value_modulus',s.simplify(wc*s.conjugate(wc))==s.Rational(1,27))

small=f.slope_certificate(732,[2,3,5,7,11]);extended=f.slope_certificate(732,[2,3,5,7,11,13])
check('fixed_small_prime_set_failure',small['surviving_even_degrees']==[6] and not small['squarefree_at_this_index_certified'])
check('additional_prime_excludes_survivor',extended['surviving_even_degrees']==[] and extended['excluding_prime_by_degree']['6']==13)
failures=[]
for n,k,p,ratio in [(51,4,11,2),(118,4,13,3)]:
    edges=f.valuation_edges(n,p);allowed=slope_square_divisor_degrees(edges)
    check('candidate_prime_shortcut_failure',p>ratio*k+(1 if ratio==2 else 0) and n%p<2*k and k in allowed)
    failures.append({'n':n,'candidate_square_degree':k,'prime':p,'n_mod_prime':n%p,
        'edges':edges,'allowed_degrees':allowed,'actual_counterexample_to_SFC_claimed':False})
kernel11=s.Poly(sum((-1)**j*11*factorial(10+j)//(factorial(j)*factorial(11-2*j))*f.Z**j for j in range(6)),f.Z)
check('kernel_real_rooted_shortcut_false',kernel11.degree()==5 and kernel11.count_roots(-s.oo,s.oo)==3 and kernel11.count_roots(0,s.oo)==3)
check('nonreal_kernel_not_a_multiple_root_counterexample',s.gcd(kernel11,kernel11.diff()).degree()==0)
diagnostics={'fixed_set_failure_at_index732':small,'prime13_certificate_at_same_index':extended,
    'failed_large_prime_sufficient_guesses':failures,
    'all_primes_limitation_claimed':False,'all_indices_coverage_claimed':False,
    'kernel11_real_rooted_shortcut_failure':{'coefficients_in_increasing_degree':[int(kernel11.nth(j)) for j in range(6)],'degree':5,'exact_real_root_count':3,'exact_positive_real_root_count':3,'squarefree':True},
    'numerical_root_evidence_used':False,'counterexample_to_SFC_claimed':False}
(HERE/'tilted_cubic_frontier_diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
audit={'proof':'TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex',
    'shifted_scope':'ALL r>=1,a>=1,i>=0,prime p>3max(a,i),n=a*p^r+i: horizontal width i plus slope2*(p^r-1)/(p^r*(p-1)), reduced denominator p^r',
    'two_residual_types':'i0 uses T_a with Gamma shape a; i>0 uses the raw a-th cubic moment normalized by a!, Gamma shape a+1. Horizontal residual T_i.',
    'prime_power_scope':'Every prime p,r>=1: T_(p^r,2) irreducible over Qp and Q. Small primes2,3 have separate primitive endpoint proofs.',
    'completed_index_families':[[a,i] for a,i in sorted(f.FAMILIES)],
    'analytic_scope':'Exact scaled integral, saddle cubic/discriminant, two nondegenerate folds and their derivatives; no contour coverage or asymptotic theorem',
    'failed_guesses':'K11 has only3 real roots at degree5, exact Sturm count, so universal kernel real-rootedness fails. Specific fixed set2..11 fails at732; prime13 works. p>2k+1 or p>3k and nmodp<2k are insufficient slope exclusions at retained examples.',
    'finite_extrapolation_used':False,'uniform_prime_coverage_proved':False,'effective_simplicity_cutoff_proved':False,
    'all_index_cubic_SFC_closed':False,'general_FC_closed':False,'unrestricted_SFC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'tilted_cubic_frontier_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
files=['TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex','tilted_cubic_frontier.py','verify_tilted_cubic_frontier.py','tilted_cubic_frontier_diagnostics.json','tilted_cubic_frontier_audit.json','prime_power_windows.py','newton_slope_coverage.py','DYADIC_BINOMIAL_SFC.tex','dyadic_binomial_sfc_verification.json']
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
    'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
    'uniform_shifted_polygon_and_residuals_proved':True,'all_prime_power_indices_irreducible':True,
    'ten_uniform_shifted_index_families_proved':True,'exact_analytic_fold_data_proved':True,
    'finite_extrapolation_used':False,'uniform_prime_coverage_proved':False,'effective_simplicity_cutoff_proved':False,
    'all_index_cubic_SFC_closed':False,'general_FC_closed':False,'unrestricted_SFC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'tilted_cubic_frontier_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'all_prime_power_indices_irreducible':True,'uniform_shifted_families':10,'all_index_cubic_SFC_closed':False}))
