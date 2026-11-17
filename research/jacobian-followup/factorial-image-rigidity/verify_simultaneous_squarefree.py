"""Exact parameters, finite common-offset controls and the scoped MR II source audit."""
from pathlib import Path
from fractions import Fraction as Q
from pypdf import PdfReader
from html.parser import HTMLParser
import hashlib,json,re
import sympy as s
import simultaneous_squarefree as a
import quantitative_squarefree as q
import prime_power_windows as p

HERE=Path(__file__).resolve().parent
counts={};examples=[]
def check(name,value):
    if not value:raise AssertionError(name)
    counts[name]=counts.get(name,0)+1

def unpack(row):return Q(row['numerator'],row['denominator'])
cert=a.parameter_certificate();get=lambda key:unpack(cert[key])
nu1,nu2,eta,rho,delta=(get(key) for key in ('nu1','nu2','eta','rho','delta'))
mp,hp=get('M_power'),get('H_power')
check('dense_source_nu_conditions',0<nu1<nu2<Q(1,6))
check('dense_source_eta_condition',0<eta<Q(1,6)-nu2/3)
check('dense_source_rho_condition_from_pi_greater_than_three',0<rho<Q(1,9)==get('rho1_strict_lower'))
check('dense_source_delta_condition',0<delta<Q(1,1000))
check('source_prime_interval_powers_are_strictly_ordered',
      0<get('P1_power')<nu1 and nu1**2<get('middle_interval_power_squared')<nu2**2)
check('source_interval_spacing_is_vacuous_at_J_one',cert['J_MR']==1)
check('power_P1_eventually_exceeds_required_log_power',get('P1_power')>0 and get('source_P1_log_power')==300)
check('source_H_range_and_Q1_inclusion_eventually_hold',nu1<hp<Q(1,2))
check('dense_error_log_exponent_is_one_over_120',get('dense_error_log_power')==Q(1,120))
check('dense_first_saving_exactly_matches_displayed_error',get('dense_first_saving')==Q(33,10**10))
check('dense_second_saving_exactly_matches_displayed_error',get('dense_second_saving')==Q(5,10**18))
check('final_power_saving_has_strict_slack',get('final_epsilon')==Q(1,10**18)<get('dense_second_saving')<get('dense_first_saving'))
check('sieve_square_root_comparison_is_exact',get('sqrt_ratio_strict_upper')**2>get('sqrt_nu1_over_nu2_squared')==Q(1,15000))
check('ordinary_sieve_complement_upper_is_exact',get('sieve_complement_upper')==Q(13,625))
check('positivity_counts_the_complement_once',get('total_loss_upper')==Q(13,625)+2*delta==Q(109,5000))
check('exponential_rational_certificate_proves_Delta_margin',
      get('e_strict_upper')<get('e_comparison') and get('e_comparison')*Q(20,19)<get('margin_log_argument')==Q(1024,357))
check('usable_offset_margin_retains_strict_slack',get('witness_margin_strict_lower')==Q(141,5000)>get('usable_offset_fraction')==Q(1,50))
zp=(mp+hp)/2;yp=(1+mp)/2
check('common_offset_count_logarithmic_ratio',hp/zp==get('J_density_log_argument')==Q(32,21))
check('common_smooth_mean_logarithmic_ratio',1/yp==get('smooth_complement_log_argument')==Q(32,17))
check('small_offset_and_gap_deletions_are_lower_order',2*mp<hp and 2*hp/3<hp)
check('uniform_local_Y_exceeds_MH_eventually',yp>mp+hp and 2*mp<1)
check('incidence_division_by_H_is_essential_and_sufficient',
      Q(1,2)+Q(3,2)*mp+2*hp==get('incidence_error_power')==Q(159,160)<1
      and Q(1,2)+Q(3,2)*mp+3*hp==Q(191,160)>1)
check('published_fallback_polynomial_saving',get('published_2016_first_saving')==Q(1,12500))
check('analytic_exception_is_common_to_all_exponents',not cert['analytic_exception_multiplied_by_M'])
check('sparse_branch_and_advertised_exponent_are_unused',
      not cert['Section10_sparse_sieve_used'] and not cert['preprint_Theorem1_7_advertised_exponent_used'])
def has_float(value):
    if isinstance(value,float):return True
    if isinstance(value,dict):return any(has_float(v) for v in value.values())
    if isinstance(value,list):return any(has_float(v) for v in value)
    return False
check('tiny_power_savings_are_never_rounded_to_floats',not has_float(cert))
for X in (2,10**8,2**80):
    data=a.dyadic_parameters(X);M,H=data['M'],data['H']
    check('dyadic_parameters_use_exact_integer_roots',M**16<=2*X<(M+1)**16 and H**5<=X<(H+1)**5)
    check('finite_dyadic_parameters_do_not_claim_asymptotic_conditions_at_small_X',not data['analytic_asymptotic_conditions_claimed_at_this_X'])

for H,M,i in ((200,3,94),(200,3,194),(500,10,202)):
    item=a.common_offset_certificate(H,M,i);prime=item['largest_prime_factor'];prev=item['preceding_prime'];gap=i-prev
    check('common_offset_criteria_are_literal_integer_inequalities',
          i>(M+1)**2 and gap*gap<i and prime*prime>M*H and i % prime==0 and s.isprime(prime))
    check('common_large_factor_cannot_have_positive_complement',
          prev>M*max(1,gap) and gap<prime and [degree for degree in range(prev,i+1) if degree % prime==0]==[i])
    for m in (1,M):
        coeffs=q.coefficients(m,i)
        at_prime=p.newton_residuals(coeffs,prime)
        at_prev=p.newton_residuals(coeffs,prev)
        check('literal_common_offset_hull_has_single_divisibility_edge',len(at_prime)==1 and at_prime[0]['e']==prime)
        check('literal_preceding_hull_supplies_one_primitive_large_edge',
              at_prev[-1]['e']==prev and at_prev[-1]['right'][0]-at_prev[-1]['left'][0]==prev)
    examples.append(item)
for m,i in ((1,1),(3,2),(6,3)):
    item=q.resultant_size_certificate(m,i);bound=int(item['hadamard_absolute_bound'])
    check('absolute_resultant_log_bound_has_no_hidden_m_constant',bound<=(2*m*i)**(5*m*i*i))
for operation in (lambda:a.common_offset_certificate(200,3,13),lambda:a.common_offset_certificate(200,3,49),
                  lambda:a.common_offset_certificate(200,3,126),lambda:a.common_offset_certificate(200,3,201),
                  lambda:a.common_offset_certificate(200,True,94),lambda:a.dyadic_parameters(True)):
    try:operation()
    except ValueError:check('outside_finite_common_offset_domain_fails_closed',True)
    else:raise AssertionError('invalid input accepted')

mr_pdf=HERE/'sources/mr-ii-2007.04290v1.pdf';reader=PdfReader(mr_pdf)
plain='\n'.join(page.extract_text() or '' for page in reader.pages)
start=re.search(r'Proof of Theorem\s+9\.2\.',plain).start()
end=start+re.search(r'10\.\s*Sieve\s+estimates',plain[start:]).start()
dense_proof=plain[start:end]
check('dense_proof_uses_earlier_comparison_proposition','Proposition 8.3' in dense_proof)
check('dense_part_two_replacements_are_explicit_in_source',all(token in dense_proof for token in ('5.4','5.5','(29)')))
check('affected_sparse_lemma_and_bad_book_reference_absent_before_section_ten',
      '10.3' not in plain[:end] and '6.18' not in plain[:end])
linear=reader.pages[20].extract_text() or ''
check('separate_ordinary_linear_sieve_is_present','Lemma 4.3.' in linear and 'F (s)' in linear)
warning=(HERE/'sources/matomaki-2022-almost-primes-accepted-v3.txt').read_text()
check('pinned_accepted_warning_retains_the_published_caveat',
      'incorrect' in warning and 'version of Lemma 6.1' in warning and 'not yet on arXiv' in warning)
class HTMLText(HTMLParser):
    def __init__(self):super().__init__();self.parts=[]
    def handle_data(self,data):self.parts.append(data)
parser=HTMLText();parser.feed((HERE/'sources/matomaki-publications-20261117.html').read_text())
author_text=' '.join(' '.join(parser.parts).split())
idx=author_text.index('Multiplicative functions in short intervals II')
check('author_catalogue_still_labels_MR_II_as_preprint','Pre-print' in author_text[idx:idx+120])
references=sorted(set(re.findall(r'(?:Lemma|Proposition|Theorem|Corollary)\s+\d+\.\d+',dense_proof)))
check('selected_dense_proof_has_no_section_ten_reference',not any(re.search(r'\s10\.',name) for name in references))

source_files=('sources/mr-ii-2007.04290v1.pdf','sources/mr-ii-2007.04290v1.txt',
    'sources/matomaki-2022-almost-primes-accepted-v3.pdf','sources/matomaki-2022-almost-primes-accepted-v3.txt',
    'sources/matomaki-publications-20261117.html','sources/filaseta-rocks-2016-published.pdf',
    'sources/filaseta-rocks-2016-published.txt','sources/matomaki-radziwill-2016-annals.pdf',
    'sources/matomaki-radziwill-2016-annals.txt')
pins={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in source_files}
audit={'date':'2026-11-17','canonical_proof':'SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex',
    'power_claim':'FIRSIMULTPOWER20261117','published_fallback_claim':'FIRSIMULTLOG20261117',
    'analytic_input_kind':'MR II preprint dense Theorem9.2(ii), arXiv2007.04290v1,8July2020',
    'theorem_pinpoint':'p44; Section9 conditions pp42-43; dense proof pp45-51; ordinary sieve Lemma4.3 p21 and equation78 p62',
    'author_status':'publication catalogue checked17November2026 lists MR II as Pre-print; arXiv has v1 only',
    'warning_source':{'accepted_version':'arXiv2012.11565v3, accepted author copy, not publisher typesetting',
        'published_url':'https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.12592',
        'publisher_warning_checked_online':True,'pinpoint':'Section6 footnote before Lemma6.1 and Remark6.2',
        'scope':'first MR II arXiv sparse-sieve calculation uses an incorrect version; missing coprimality creates off-diagonal terms'},
    'selected_dependency_route':['dense Theorem9.2(ii)','Proposition8.3 and earlier inputs',
        'dense replacements Lemmas5.4(ii),5.5(ii),equation29','ordinary linear sieve Lemma4.3 as in equation78'],
    'dense_proof_direct_named_references':references,
    'sparse_Section10_including_Lemma10_3_used':False,'Theorem1_7_advertised_exponent_used':False,
    'complete_preprint_reconstruction_or_validation_claimed':False,
    'corrected_MR_II_version_obtained':False,'published_MR_II_claimed':False,
    'long_mean_and_common_offset_density':'PNT at fixed power ranges; no fixed-m density theorem applied with growing m',
    'common_J_density':'log(32/21)','smooth_mean':'1-log(32/17)','margin':'log(1024/357)-1>1/20',
    'complement_loss':'once by positivity; upper13/625+2/2000=109/5000+o(1)',
    'usable_offsets':'at least H/50 outside one common analytic exception; no factor M',
    'absolute_resultant_bound':'log abs R_(i,m)<=5*m*i^2*log(2mi); sum<=5*M^2*H^3*log(2MH)',
    'incidence_error':'divide witness count by H/50, yielding X^(159/160)',
    'dense_errors':['X^(1-33/10^10)*(log X)^2','X^(1-5*10^-18)'],
    'power_rate':'D(X)<< X^(1-10^-18), absolute implied constant',
    'published_2016_independent_fallback':'D(X)<< X/(log X)^(1/50), absolute constant, delta1/100,H=X^(1/5)',
    'reciprocal_sum_converges_using_power_route':True,
    'fixed_m_reciprocal_sum':'finite initial part n<m^16; later indices lie in the common exceptional set',
    'numeric_implied_constant_or_starting_threshold_extracted':False,
    'all_index_SFC_closed':False,'finite_controls_are_not_numerical_density_evidence':True,
    'source_pins':pins,'provenance':'User-supplied simultaneous argument and explicit parameters checked; no priority claim',
    'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'simultaneous_squarefree_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
(HERE/'simultaneous_squarefree_parameters.json').write_text(json.dumps(cert,indent=2)+'\n')
(HERE/'simultaneous_squarefree_examples.json').write_text(json.dumps({'common_offsets':examples,
    'finite_controls_are_not_numerical_density_evidence':True},indent=2)+'\n')
files=('SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex','simultaneous_squarefree.py','verify_simultaneous_squarefree.py',
    'simultaneous_squarefree_audit.json','simultaneous_squarefree_parameters.json','simultaneous_squarefree_examples.json',
    'quantitative_squarefree.py','prime_power_windows.py',*source_files)
receipt={'status':'pass','checks':counts,'total':sum(counts.values()),
    'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
    'universal_evidence':'SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex with explicitly qualified analytic imports; not finite tests',
    'simultaneous_range':'all integer m<=n^(1/16), one common exception',
    'published_2016_fallback_proved':True,'preprint_dense_power_application_checked':True,
    'final_epsilon':{'numerator':1,'denominator':10**18},'implied_constant_absolute':True,
    'MR_II_preprint_source_qualification_retained':True,'complete_preprint_validation_claimed':False,
    'sparse_Section10_used':False,'Theorem1_7_advertised_exponent_used':False,
    'reciprocal_sum_converges_using_power_route':True,'analytic_exception_multiplied_by_M':False,
    'numeric_implied_constant_or_starting_threshold_extracted':False,
    'finite_controls_are_not_numerical_density_evidence':True,'all_index_SFC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'simultaneous_squarefree_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'range':receipt['simultaneous_range'],
    'epsilon':'1/10^18','power_input':'MR II dense preprint theorem; source-qualified',
    'published_log_fallback':True,'all_index_SFC_closed':False}))
