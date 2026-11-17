"""Exact constants and complete finite prefix for all ordinary A+B Y^m.

Universal evidence is ALL_UNIVARIATE_TWO_TERM_SFC.tex. Arithmetic controls
are not formal verification or an independent external proof review.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib,json
import sympy as S
import uniform_two_term_sfc as u
HERE=Path(__file__).resolve().parent
B=F(11,5);R=F(21,10);D=F(1,25)
qt=B*F(35,81)**6;qc=B*F(5,9)**6;qf=B*F(5,18)**6
rl=F(101,150)*F(141,100);rv=F(38,45)
ineq={
 'e_upper_49_over18':(F(65,24)+F(1,100),F(49,18)),
 'e_lower27_over10':(F(27,10),sum(F(1,factorial(j)) for j in range(5))),
 'e_lower8_over3':(F(8,3),sum(F(1,factorial(j)) for j in range(5))),
 'high_k_coefficient':(F(300,299)*F(125,123),F(21,20)),
 'low_k_coefficient':(F(54,10**5),F(21,20)),
 'log10000_upper':(F(10000),F(27,10)**10),
 'principal_saddle_Rouche':(B,D*(1-D)*F(27,10)**6/6),
 'curvature_real':(F(9,10),1-D),
 'curvature_modulus':((1+D)/(1-D/6),F(11,10)),
 'phase_height_loss':(D/6+D**2/(12*(1-D/6)),F(1,100)),
 'local_radius':(F(6,7),F(7,8)-D/6),
 'third_derivative':(F(7,6)**3*(2+36*qt*(1+qt)/(1-qt)**3+18*qt/(1-qt)**2+2*qt/(1-qt)),F(8)),
 'third_derivative_monotonicity':(F(35,81)*F(7,6)**2,F(1)),
 'horizontal_curvature':(F(1,3),F(9999,10000)/F(10001,10000)**2-6*qc/(1-qc)**2-qc/(1-qc)),
 'curvature_monotonicity':(F(5,9)*F(7,6),F(1)),
 'far_derivative':(F(3,4)/(1-qf)-1,-F(1,5)),
 'far_height_loss':(F(1,100),F(3,32)*F(49,150)**2),
 'compact_height_loss':(F(1,1000),F(3,2048)),
 'low_envelope_increases':(B,F(27,20)**6),
 'low_base_decreases':(rl*F(50,49),F(1)),
 'low_envelope_m6':(B*F(10,27)**3+rl**6/F(22,25),F(1)),
 'vertical_base_decreases':(rv*F(50,49),F(1)),
 'vertical_envelope_m6':(B*F(10,27)**2+rv**6/F(22,25),F(1)),
 'exp_one_third':(F(49,18),F(141,100)**3),
 'sqrt_e':(F(49,18),F(25,9)),
 'inverse_sqrt_e':(F(64,25),F(27,10)),
 'horizontal_log_chart':(B*F(20,27)**6,F(1)),
 'local_Gaussian_error':(F(4,3)/F(17,60)**2,F(17)),
 'Gaussian_tail_constant':(F(160,9),F(18)),
 'total_error_constant':(40+F(6*1000**3,36*10000**2)+F(2*100**2,10000)*F(6+D,36),F(44)),
 'relative_error_bound_squared':(F(22**2*100),F(6*10000)),
 'relative_derivative_bound_squared':(F(220**2),F(6*10000)),
 'saddle_derivative':(6*F(10,27)**6/(1-D)**2,F(1)),
 'prefactor_log_derivative':(F(5,9)/(1-D/6)**2,F(1)),
 'endpoint_decay_m6':(1+F(5,2)*F(49,36)**6,sum(F(57,20)**j/factorial(j) for j in range(13))),
 'endpoint_phase_base':(F(200,117),F(2)),
 'endpoint_h_positive':(F(19,40),F(5,6)),
 'endpoint_relative_error':(150*6*10000*F(1,500)**6+F(8,6**2*10000**2),F(1,1000)),
 'endpoint_relative_derivative_error':(375*6*10000*F(1,500)**6,F(1,1000)),
 'endpoint_log_derivative_error':(F(1,1000)/(1-F(1,1000)),F(1)),
 'inner_envelope':(F(5,8)**6+2*F(5,6)**6,F(1)),
 'inner_comparison':(10000*F(50,99)**10000,F(3)),
 'inner_comparison_monotonicity':(F(3,2)*F(50,99),F(1)),
 'log_derivative_gap':(F(1,3),(1-D/6)/R),
 'large_index_gap':(F(4),F(10000,3)),
 'small_exponent_overlap':(F(2)+4*7,F(100)),
 'finite_exponent_cap':(F(4*9999**2),F(1<<30)),
}
for name,(lo,hi) in ineq.items():assert lo<hi,name
m=S.symbols('m',integer=True,positive=True)
s=S.symbols('s',positive=True);w,v,theta=S.symbols('w v theta')
phase=S.log(w+s**m)-s;sw=s**(m-1)*(m-s)
identities={
 'entropy_derivative':S.diff((m-1)*S.log(theta)-(1-theta)/theta*S.log(1-theta),theta)-(m*theta+S.log(1-theta))/theta**2,
 'phase_stationarity':S.diff(phase,s).subs(w,sw),
 'phase_curvature':S.diff(phase,s,2).subs(w,sw)-(m-1-s)/s,
 'saddle_map_derivative':S.diff(sw,s)-m*s**(m-2)*(m-1-s),
 'phase_derivative_gap':(m-s)/(m*v)-1/v+s/(m*v),
 'local_Gaussian_exponent':S.Rational(9,20)-S.Rational(1,6)-S.Rational(17,60),
 'local_Gaussian_integral':S.Rational(4,3)/S.Rational(17,60)**2-S.Rational(4800,289),
 'high_k_factor_constant':10*(1+S.Rational(3,2)*10)/10000-S.Rational(2,125),
}
for name,value in identities.items():assert S.simplify(value)==0,name
prefix=json.loads((HERE/'all_exponent_prefix_certificates.json').read_text())
assert prefix['range_n']==[2,9999]
expected={(n,m) for m in range(1,30) if m!=6 and not u.is_prime_power_or_one(m)
          for n in range(2,10000) if (1<<m)<4*n*n}
assert u.finite_exponents()==[10,12,14,15,18,20,21,22,24,26,28]
seen=set();modes={};gcds={};by_m={};max_kernel=0;max_full=0
for row in prefix['records']:
    n,m,p,d,mode,g,a,b=row
    assert (n,m) not in seen and (n,m) in expected
    u.verify_record(row)
    seen.add((n,m));modes[mode]=modes.get(mode,0)+1;gcds[len(g)-1]=gcds.get(len(g)-1,0)+1;by_m[m]=by_m.get(m,0)+1
    if mode=='full_degree_modular':max_full=max(max_full,n)
    else:max_kernel=max(max_kernel,d)
assert seen==expected and len(seen)==92935
assert modes=={'preceding_prime_kernel':92900,'full_degree_modular':35}
assert gcds=={0:92914,1:21} and max_kernel==35 and max_full==126
# Verify source receipts before using the retained m6 and prime-power/large-m
# proofs; do not reconstruct unrelated calculations.
for receipt_name in ['exponent_six_simplicity_verification.json','all_exponent_squarefree_verification.json']:
    source=json.loads((HERE/receipt_name).read_text());assert source['status']=='pass'
    for name,digest in source['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,(receipt_name,name)
audit={'uniform_analytic_cutoff':10000,'analytic_exponent_range':'6<=m<=sqrt(n)',
       'normalization':'z=n^m*(m/e)^m*v','root_annulus_in_analytic_range':['1/2','21/10'],
       'tail_disk':'|v|<=11/5','relative_error':'22/sqrt(m*n)','log_derivative_gap':'n/3-4',
       'large_exponent_bridge':'2^m>=4n^2 excludes the complementary range; if it fails, m<sqrt(n) for n>=10000',
       'finite_prefix_pairs':len(seen),'finite_exponents':u.finite_exponents(),'by_exponent':by_m,
       'certificate_modes':modes,'gcd_degree_counts':gcds,'max_kernel_degree':max_kernel,'max_full_degree_index':max_full,
       'retained_cases':'all prime powers and m1; all-index m6 proof reused',
       'all_finite_prefix_Bezout_identities_checked':True,'complete_pair_domain_checked':True,
       'analytic_inequalities':list(ineq),'symbolic_identities':list(identities),
       'universal_evidence':'ALL_UNIVARIATE_TWO_TERM_SFC.tex: exact interval split and full uniform analytic bounds',
       'all_positive_exponents_all_indices_squarefree':True,'ordinary_A_plus_B_Ym_SFC_all_starts':True,
       'ordinary_exceptional_set':'EMPTY; E_m(X)=D_all(X)=0',
       'unrestricted_SFC_closed':False,'general_FC_closed':False,'other_Gamma_shapes_claimed':False,
       'numerical_root_scan':False,'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'uniform_two_term_sfc_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
files=('ALL_UNIVARIATE_TWO_TERM_SFC.tex','uniform_two_term_sfc.py','verify_uniform_two_term_sfc.py',
       'all_exponent_prefix_certificates.json','uniform_two_term_sfc_audit.json',
       'EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex','exponent_six_simplicity_verification.json',
       'ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex','all_exponent_squarefree_verification.json',
       'PRIME_POWER_TWO_TERM_SFC.tex','QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex',
       'MIXED_PRIME_RESIDUAL_WINDOWS.tex')
receipt={'status':'pass','total':len(ineq)+len(identities)+len(seen),
         'checks':{'uniform_exact_inequalities':len(ineq),'symbolic_identities':len(identities),'finite_prefix_Bezout_certificates':len(seen)},
         'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
         'effective_uniform_simplicity_cutoff':10000,'all_positive_exponents_all_indices_squarefree':True,
         'ordinary_A_plus_B_Ym_SFC_all_starts':True,'complete_finite_prefix':True,
         'ordinary_exceptional_set_empty':True,'unrestricted_SFC_closed':False,'general_FC_closed':False,
         'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'uniform_two_term_sfc_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'uniform_simplicity_cutoff':10000,
                  'finite_prefix_pairs':len(seen),'all_positive_exponents_all_indices_squarefree':True,
                  'unrestricted_SFC_closed':False,'general_FC_closed':False}))
