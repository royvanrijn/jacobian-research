"""Exact algebra and boundary checks for the dyadic binomial transfer.

The universal proof is the companion note. Finite instances below are
bookkeeping checks, not all-index or all-dilation extrapolation.
"""
from pathlib import Path
from math import comb,factorial
import hashlib,json
import sympy as s
import dyadic_binomial_sfc as p
from prime_power_windows import newton_residuals
HERE=Path(__file__).resolve().parent
z=p.Z
counts={}
def check(category,condition):
    assert condition,category
    counts[category]=counts.get(category,0)+1

n,k,q,L,delta,h=s.symbols('n k q L delta h',integer=True,positive=True)
check('universal_Bessel_coefficient_identity',s.expand((n+k)*(n+k-1)-2*k*(2*n-1)-(n-k)*(n-k-1))==0)
W=lambda qq,ll,dd:s.factorial(qq*ll+dd)/(2**((qq-1)*ll)*s.factorial(ll)*s.factorial(dd))
odd_factor=s.factorial(2*q*L)/(2**(q*L)*s.factorial(q*L))
check('universal_doubling_factorial_identity',s.simplify(W(2*q,L,0)/(W(q,L,0)*odd_factor))==1)
check('universal_offset_binomial_identity',s.combsimp(W(q,L,delta)/(W(q,L,0)*s.binomial(q*L+delta,delta)))==1)
check('universal_same_variable_scaling',s.expand((q-1)*(n+k+h)-(q-1)*h-(q-1)*n-q*k+k)==0)
# Formal differentiation of the recurrence proves the shifted relation
# for arbitrary n and arbitrary predecessor polynomials.
Y0=s.Function('Y0')(z);Y1=s.Function('Y1')(z);Yn=(2*n-1)*z*Y1+Y0
D=lambda expr,nn:(nn+1)*expr+z*s.diff(expr,z)
check('universal_shifted_Bessel_recurrence',s.expand(D(Yn,n)-(2*n-1)*z*D(Y1,n-1)-D(Y0,n-2)-2*Yn)==0)

for nn in range(0,14):
    for hh in [0,1]:
        P=p.bessel_polynomial(nn,hh)
        check('Bessel_coefficient_integrality',all(c.q==1 for c in P.all_coeffs()))
        check('Bessel_leading_unit',int(P.LC())%2==1 and P.degree()==nn)
        if nn>=2:
            previous=p.bessel_polynomial(nn-1,hh);before=p.bessel_polynomial(nn-2,hh)
            correction=2*p.bessel_polynomial(nn,0).as_expr() if hh else 0
            check('Bessel_literal_recurrences',s.expand(P.as_expr()-(2*nn-1)*z*previous.as_expr()-before.as_expr()-correction)==0)
        if nn>=1:
            a=s.Poly(P,z,modulus=2);bb=s.Poly(p.bessel_polynomial(nn-1,hh),z,modulus=2)
            check('Bessel_reduction_coprimality',s.gcd(a,bb).degree()==0)

for qq in [1,2,4,8]:
    for dd in range(qq):
        for ll in range(0,7):
            unit=p.factorial_unit(qq,ll,dd)
            literal=s.Rational(factorial(qq*ll+dd),2**((qq-1)*ll)*factorial(ll)*factorial(dd))
            check('literal_factorial_unit',literal==unit and unit%2==1)
            if dd:
                for j in range(1,dd+1):
                    v=lambda x:(x&-x).bit_length()-1
                    check('offset_valuation_identity',v(qq*ll+j)==v(j))
    for beta in range(1,2*qq+1):
        hh,dd=divmod(beta-1,qq)
        for nn in range(0,7):
            G=p.normalized_moment(qq,beta,nn)
            for kk in range(nn+1):
                actual=s.Rational(factorial(hh)*p.factorial_unit(qq,hh,dd)*comb(nn,kk)*factorial(beta-1+qq*(nn+kk)),factorial(beta-1)*factorial(nn)*2**((qq-1)*nn+qq*kk))
                check('Gamma_normalization_coefficient_identity',G.nth(kk)==actual)
            check('Gamma_reduction_identity',s.Poly(G,z,modulus=2)==s.Poly(p.bessel_polynomial(nn,hh),z,modulus=2))
            check('Gamma_leading_degree_preserved',G.degree()==nn and int(G.LC())%2)
            if nn>=1:
                other=p.normalized_moment(qq,beta,nn-1)
                check('Gamma_reduction_consecutive_coprimality',s.gcd(s.Poly(G,z,modulus=2),s.Poly(other,z,modulus=2)).degree()==0)

# Exact raw moment resultants independently verify the scaling valuation.
for qq in [1,2,4]:
    for beta in range(1,2*qq+1):
        for nn in range(1,6):
            raw=lambda d:s.Poly(sum(s.binomial(d,j)*s.rf(beta,qq*(d+j))*z**j for j in range(d+1)),z)
            value=int(s.resultant(raw(nn),raw(nn-1)))
            check('raw_resultant_nonzero',value!=0)
            value=abs(value);valuation=(value&-value).bit_length()-1
            check('raw_resultant_exact_valuation',valuation==p.exact_resultant_valuation(qq,nn))

# The general unit theorem does not require the exact Bessel recurrence
# in characteristic zero. Nonconstant rational unit weights check this.
for hh in [0,1]:
    for nn in range(1,8):
        G=s.Poly(sum(s.Rational(2*(nn+kk)+1,4*(nn+kk)+1)*c*z**kk for kk,c in enumerate(p.bessel_polynomial(nn,hh).all_coeffs()[::-1])),z)
        check('rational_unit_denominators',all(int(c.q)%2 for c in G.all_coeffs()))
        reduced=s.Poly(sum(int(c.p)%2*z**kk for kk,c in enumerate(G.all_coeffs()[::-1])),z,modulus=2)
        check('rational_unit_same_reduction',reduced==s.Poly(p.bessel_polynomial(nn,hh),z,modulus=2))

for rr,beta in [(0,1),(0,2),(1,1),(1,4),(3,16),(1000000,1)]:
    cert=p.structural_certificate(rr,beta)
    check('unbounded_structural_input',cert['accepted'] and cert['all_starting_indices'] and not cert['all_positive_Gamma_shapes'])
try:p.structural_certificate(1,5)
except ValueError:check('outside_shape_range_rejected',True)
else:raise AssertionError('shape range silently extended')

for nn in range(1,7):
    for gap in range(1,5):
        T=p.lower_one_tilted(nn,gap).as_expr();theta=z*s.diff(T,z)
        raw=lambda j:sum(s.binomial(j,k)*s.Rational(factorial(j+gap*k),factorial(j))*z**k for k in range(j+1))
        check('lower_one_current_moment_identity',s.expand(raw(nn)-T-s.Rational(gap,nn)*theta)==0)
        check('lower_one_previous_moment_identity',s.expand(raw(nn-1)-T+s.Rational(1,nn)*theta)==0)
        check('lower_one_nonzero_constant',p.lower_one_tilted(nn,gap).nth(0)==1)
check('lower_one_moment_matrix_determinant',s.cancel(s.Matrix([[1,q/n],[1,-1/n]]).det()+(q+1)/n)==0)
for offset in range(1,13):
    A=p.cubic_residual_kernel(offset).as_expr();B=p.cubic_residual_kernel(offset+1).as_expr();K=p.cubic_square_kernel(offset).as_expr()
    check('cubic_kernel_linear_identity',s.expand(3*K-A-2*B)==0)
    check('cubic_kernel_derivative_identity',s.expand(3*z*s.diff(K,z)-offset*(B-A))==0)

# Universal arithmetic for the conditional cubic prime window. Positivity
# and slope separation on the full domain are proved in the note.
hh,dd,qq=s.symbols('hh dd qq',integer=True,positive=True)
check('cubic_even_cross_width_gap',s.expand((qq-1)/2-(3*hh-qq+1)-3*(dd-1)/2).subs(qq,2*hh+dd).expand()==0)
check('cubic_odd_cross_width_gap',s.expand((qq-1)/2-(3*hh-qq+2)-(3*dd-2)/2).subs(qq,2*hh+1+dd).expand()==0)
check('cubic_even_width_total',s.expand((qq-2*hh-1)/2+(qq-1)/2+(3*hh-qq+1)-2*hh)==0)
check('cubic_odd_width_total',s.expand((qq-2*hh-3)/2+(qq+1)/2+(3*hh-qq+2)-(2*hh+1))==0)
for j in range(0,7):
    check('cubic_symbolic_even_falling_identity',s.expand(s.rf(1-dd,2*j)-s.ff(dd-1,2*j))==0)
    check('cubic_symbolic_kernel_coefficient',s.combsimp(s.rf(dd,j)*s.factorial(dd-1)/(s.factorial(j)*s.factorial(dd-1-2*j))-s.factorial(dd-1+j)/(s.factorial(j)*s.factorial(dd-1-2*j)))==0)
windows=[]
for nn in range(6,41):
    C=lambda index:[comb(index,k)*factorial(index+2*k) for k in range(index+1)]
    for qq in s.primerange(nn+1,(3*(nn-1))//2+1):
        qq=int(qq);row=p.cubic_prime_window(nn,qq)
        E=[newton_residuals(C(index),qq) for index in [nn,nn-1]]
        slopes=lambda edges:{s.Rational(e['h'],e['e']) for e in edges if e['h']}
        check('cubic_no_common_nonhorizontal_slopes',not slopes(E[0])&slopes(E[1]))
        for index in [nn,nn-1]:
            raw=s.Poly(sum(s.Rational(c,factorial(index))*z**j for j,c in enumerate(C(index))),z)
            kernel=p.cubic_residual_kernel(qq-index)
            check('cubic_residual_kernel_specialization',s.Poly(raw,z,modulus=qq)==s.Poly(kernel,z,modulus=qq))
        T=p.lower_one_tilted(nn,2)
        edges=newton_residuals(list(reversed([int(c) for c in T.all_coeffs()])),qq)
        check('cubic_tilted_nonhorizontal_primitive',all(e['right'][1]-e['left'][1]==1 for e in edges if e['h']))
        horizontal=next((e['right'][0]-e['left'][0] for e in edges if not e['h']),0)
        check('cubic_tilted_horizontal_width',horizontal==(qq-nn)//2)
        check('cubic_tilted_kernel_specialization',s.Poly(T,z,modulus=qq)==s.Poly(p.cubic_square_kernel(qq-nn),z,modulus=qq))
        if row['pair_certified_coprime']:
            check('cubic_window_exact_coprimality',s.gcd(s.Poly(sum(c*z**j for j,c in enumerate(C(nn))),z),s.Poly(sum(c*z**j for j,c in enumerate(C(nn-1))),z)).degree()==0)
        windows.append(row)

# Retained exact limitations of the first residue-field certificate.
D2=s.Poly(sum(comb(2+kk+2,2)*c*z**kk for kk,c in enumerate(p.bessel_coefficients(2))),z)
D3=s.Poly(sum(comb(3+kk+2,2)*c*z**kk for kk,c in enumerate(p.bessel_coefficients(3))),z)
check('shape_extension_shared_reduction',s.Poly(D2,z,modulus=2)==s.Poly(z**2,z,modulus=2)==s.Poly(D3,z,modulus=2))
check('shape_extension_not_counterexample',s.gcd(D2,D3).degree()==0)
C2=[comb(2,kk)*factorial(2+2*kk) for kk in range(3)]
C3=[comb(3,kk)*factorial(3+2*kk) for kk in range(4)]
E2=newton_residuals(C2,2);E3=newton_residuals(C3,2)
check('gap_two_common_slope',E2[0]['e']==E3[0]['e']==2 and E2[0]['h']==E3[0]['h']==3)
check('gap_two_common_residual',E2[0]['residual_coefficients']==E3[0]['residual_coefficients']==[1,1])
check('gap_two_not_counterexample',s.gcd(s.Poly(sum(c*z**j for j,c in enumerate(C2)),z),s.Poly(sum(c*z**j for j,c in enumerate(C3)),z)).degree()==0)
frontier={'first_uncovered_shape_base_shift_h':2,'base_reductions_at_n2_n3':'both z^2; exact rational polynomials are coprime',
    'gap2_lower_exponent1':{'terminal_index':3,'moment_coefficients':[C2,C3],'Newton_residuals':[E2,E3],
        'common_slope':'3/2','common_residual':'1+Y','characteristic_zero_coprime':True},
    'prime_window_diagnostics':windows,'all_prime_windows_covered_claimed':False,
    'counterexample_claimed':False,'uniform_all_gaps_inferred':False}
(HERE/'dyadic_binomial_frontier_diagnostics.json').write_text(json.dumps(frontier,indent=2)+'\n')
audit={'proof':'DYADIC_BINOMIAL_SFC.tex','lower_one_equivalence':'For all positive n,e, SFC adjacent moments for X*(A+B*X^e) are equivalent to simplicity of T_(n,e)=sum binom(n,k)*(n)_(e*k)*t^k; identity only, universal simplicity unproved.',
    'conditional_cubic_scope':'For prime n<q<=3(n-1)/2,n>=6, every common factor of moments n-1,n of X*(A+B*X^2) has even degree<=floor((q-n-1)/2); degree-preserving reduction divides gcd(H_(q-n),H_(q-n+1)). No uniform all-index coverage proved.',
    'all_index_scope':'q=2^r,r>=0; X^q(A+B*X^q); integer Gamma shape1<=beta<=2q; all complex A,B giving a nonzero polynomial',
    'uniform_input':'arbitrary rational 2-adic unit sequences multiplying (L+h)!,h0 or1',
    'resultant_valuation':'3q*n*(n-1)-(n-1)*s2(n)-n*s2(n-1), independent of beta in the proved integer range',
    'mechanism':'one common variable scaling; integral normalized coefficient polynomials; odd leading coefficients; Bessel/Fibonacci consecutive reductions; resultant odd',
    'source_scope':'Rocks2015 author draft Lemma3.8 for unshifted factorial units compared; every input rederived; no full final-publication or priority audit',
    'independent_of':['prior contour/simplicity results','parabola separation','MR/density','prime gaps','finite-prefix searches'],
    'finite_parameter_extrapolation':False,'other_lower_exponents_or_gaps_proved':False,'all_positive_Gamma_shapes':False,
    'general_FC_closed':False,'unrestricted_SFC_closed':False,'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'dyadic_binomial_sfc_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
files=['DYADIC_BINOMIAL_SFC.tex','dyadic_binomial_sfc.py','verify_dyadic_binomial_sfc.py','dyadic_binomial_frontier_diagnostics.json','dyadic_binomial_sfc_audit.json','prime_power_windows.py','sources/rocks-2015-author-preprint.pdf','sources/rocks-2015-author-preprint.txt']
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
    'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
    'dyadic_equal_lower_exponent_and_gap_all_indices':True,'integer_Gamma_shapes_through_2q':True,
    'arbitrary_rational_unit_sequence_transfer':True,'same_variable_scaling_at_consecutive_orders':True,'exact_resultant_valuation_proved':True,'conditional_cubic_prime_confinement_proved':True,'cubic_residual_transfer_proved':True,'cubic_all_index_SFC_closed':False,'lower_one_multiple_root_equivalence_proved':True,'cubic_quarter_gap_square_degree_bound_proved':True,
    'finite_parameter_extrapolation':False,'other_lower_exponents_or_gaps_proved':False,'all_positive_Gamma_shapes':False,
    'general_FC_closed':False,'unrestricted_SFC_closed':False,'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'dyadic_binomial_sfc_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'all_indices_dyadic_gap_family':True,'Gamma_shape_range':'integer1..2q','unrestricted_SFC_closed':False}))
