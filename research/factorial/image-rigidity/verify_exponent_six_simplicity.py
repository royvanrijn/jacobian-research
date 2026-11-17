"""Exact constants and complete finite-prefix certificates for exponent six.

The large-index universal evidence is EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex.
This script verifies constants and every modular Bezout identity in the prefix;
it does not formally verify the analytic proof or solve unrestricted SFC.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb,isqrt
from bisect import bisect_right
import hashlib,json
import sympy as S
HERE=Path(__file__).resolve().parent
D=F(1,25);q=F(250,4**6);r=F(250,5**6)
ineq={
 'e_lower_19_over7':(F(19,7),sum(F(1,factorial(j)) for j in range(6))),
 'e_lower_27_over10':(F(27,10),sum(F(1,factorial(j)) for j in range(5))),
 'large_k_factor':(F(18000),F(101,100)**1000),
 'entropy_critical_location':(F(300),F(27,10)**5*(1+F(49,50)+F(49,50)**2/2)),
 'outer_large_coefficients':(F(42,19)**6*F(300,299)*F(101,100),F(120)),
 'outer_small_k':(F(177147,128)*F(1,10)**5,F(120)),
 'outer_k1':(F(720,3**5),F(120)),
 'outer_k2_squared':(F(factorial(12),2),F(120**2*3**10)),
 'saddle_Rouche_disk':(F(250),(6-D)**5*D),
 'curvature_real':(F(1,7),1-F(5)/(6-D)),
 'curvature_modulus':((1+D)/(6-D),F(1,5)),
 'third_derivative':((12+216*r*(1+r)/(1-r)**3+108*r/(1-r)**2+12*r/(1-r))/125,F(1,4)),
 'phase_loss':(D/6+5*D**2/(72*(1-D/6)),F(1,100)),
 'log6_lower':(F(179,100),2*sum(F(5,7)**(2*j+1)/F(2*j+1) for j in range(12))),
 'phase_height':(F(47,10),6*F(179,100)-6-F(1,100)),
 'horizontal_curvature':(F(3),6*(1-F(1,10000))/(1+F(1,10000))**2-(36*q/(1-q)**2+6*q/(1-q))),
 'vertical_at_three':(F(250)+(3+D)**6,F(27,10)**7),
 'horizontal_three_to_four':(F(250)+(4+D)**6,F(27,10)**8*F(41,25)),
 'sqrt_e_lower':(F(41,25)**2,F(27,10)),
 'far_derivative':(6/(8*(1-F(250,8**6)))-1,-F(1,5)),
 'far_height_gap':(F(1,20),F(3,128)*F(49,25)**2),
 'medium_height_gap':(F(1,100),F(3,128)*F(3,4)**2),
 'local_Gaussian_error':(F(1,24)/F(9,224)**2,F(26)),
 'Gaussian_tail_constant':(F(56,3),F(19)),
 'complete_error_constant':(F(50)+F(80000+50,10000)+F(8,25*10000),F(60)),
 'endpoint_decay':(F(779,50),F(27,10)**2*(1+F(17,20)+F(17,20)**2/2)),
 'principal_exponential_above100':(F(100),F(27,10)**4*(1+F(7,10)+F(7,10)**2/2)),
 'prefactor_log_derivative':(F(35,2)/((6-D)**2*6*(6-D)**4*(1-D)),F(1,1000)),
 'principal_log_derivative_error':(F(1,1000)+F(12,1000)/(1-F(12,100)),F(1)),
 'endpoint_relative_error':(F(720*20**7,50*10000**5)+F(2,9*10000**2),F(1,1000)),
 'endpoint_relative_derivative_error':(F(720*20**7,2500*10000**5),F(1,1000)),
 'endpoint_log_derivative_error':(F(1,1000)/(1-F(1,1000)),F(1)),
 'inner_root_exclusion':(3*F(3,5)**10000,F(4,100)),
 'inner_comparison_monotone_squared':(2*F(3,5)**2,F(1)),
 'log_derivative_gap':(F(2),F(10000,288))}
for name,(lo,hi) in ineq.items():assert lo<hi,name
s,w,t,theta=S.symbols('s w t theta');phase=S.log(w+s**6)-s
entropy=5*S.log(theta)-(1-theta)/theta*S.log(1-theta)
identities={
 'entropy_derivative':S.diff(entropy,theta)-(6*theta+S.log(1-theta))/theta**2,
 'phase_stationary':S.diff(phase,s).subs(w,s**5*(6-s)),
 'phase_curvature':S.diff(phase,s,2).subs(w,s**5*(6-s))-(5-s)/s,
 'phase_endpoint_derivative_gap':(6-s)/(6*w)-1/w+s/(6*w),
 'saddle_inverse_derivative':S.diff(s**5*(6-s),s)-6*s**4*(5-s),
 'local_Gaussian_exponent':S.Rational(1,14)-S.Rational(3,4)/24-S.Rational(9,224)}
for name,value in identities.items():assert S.simplify(value)==0,name
# Polynomials below use increasing degree. Every Bezout identity is checked
# separately by multiplication and addition, rather than trusting gcd output.
def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a

def add(a,b,p):
 out=[0]*max(len(a),len(b))
 for j in range(len(a)):out[j]=a[j]
 for j in range(len(b)):out[j]=(out[j]+b[j])%p
 return trim(out)

def neg(a,p):return [(-c)%p for c in a]

def mul(a,b,p):
 out=[0]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):out[i+j]=(out[i+j]+u*v)%p
 return trim(out)

def divrem(a,b,p):
 a=a[:];q=[0]*max(1,len(a)-len(b)+1);inv=pow(b[-1],-1,p)
 while len(a)>=len(b) and a!=[0]:
  c=a[-1]*inv%p;shift=len(a)-len(b);q[shift]=c
  for j,v in enumerate(b):a[j+shift]=(a[j+shift]-c*v)%p
  trim(a)
 return trim(q),a

def bezout(f,g,p):
 r0,r1=f,g;u0,u1=[1],[0];v0,v1=[0],[1]
 while r1!=[0]:
  quotient,remainder=divrem(r0,r1,p)
  r0,r1=r1,remainder
  u0,u1=u1,add(u0,neg(mul(quotient,u1,p),p),p)
  v0,v1=v1,add(v0,neg(mul(quotient,v1,p),p),p)
 scale=pow(r0[-1],-1,p)
 return [c*scale%p for c in r0],[c*scale%p for c in u0],[c*scale%p for c in v0]

N=10000;sieve=bytearray([1])*N;sieve[:2]=b'\0\0'
for p in range(2,isqrt(N-1)+1):
 if sieve[p]:sieve[p*p:N:p]=b'\0'*len(range(p*p,N,p))
primes=[j for j in range(N) if sieve[j]]
# Independently verify each prime actually selected, by trial division.
trial_prime=lambda p:p>=2 and all(p%v for v in range(2,isqrt(p)+1))
facts=[factorial(6*j) for j in range(36)]
records=[];modes={'full_degree_modular':0,'preceding_prime_kernel':0};gcds={0:0,1:0};max_offset=0
for n in range(2,N):
 if n<30:p=211;d=n;mode='full_degree_modular';assert p>6*n
 else:
  p=primes[bisect_right(primes,n)-1];d=n-p;mode='preceding_prime_kernel'
  assert p>6*max(1,d)
  max_offset=max(max_offset,d)
 assert trial_prime(p)
 f=[comb(d,j)*(facts[j]%p)%p for j in range(d+1)]
 assert f[-1]!=0
 derivative=trim([j*f[j]%p for j in range(1,len(f))] or [0])
 g,u,v=bezout(f,derivative,p)
 assert len(g)<=2 and g[-1]==1,(n,p,d,g)
 assert add(mul(u,f,p),mul(v,derivative,p),p)==g,(n,'Bezout identity')
 gcds[len(g)-1]+=1;modes[mode]+=1
 records.append([n,p,d,mode,g,u,v])
assert [r[0] for r in records]==list(range(2,N)) and max_offset==35
prefix={'range_inclusive':[2,N-1],'n1':'linear',
 'columns':['n','prime','kernel_degree_or_full_degree','mode','G','U','V'],
 'identity':'U*F_d+V*F_d_prime=G modulo prime; degree G<=1',
 'polynomial_coefficient_order':'increasing degree','records':records}
(HERE/'exponent_six_prefix_certificates.json').write_text(json.dumps(prefix,separators=(',',':'))+'\n')
# Keep audit short: the enormous rational powers above are checked exactly,
# but do not need to be printed as decimal strings in the source ledger.
audit={'analytic_cutoff':N,'root_annulus_for_n_ge_cutoff':[60,240],
 'analytic_constant_checks':list(ineq),'symbolic_checks':list(identities),
 'finite_prefix_count':len(records),'certificate_modes':modes,'gcd_degree_counts':gcds,
 'largest_preceding_prime_offset':max_offset,'max_polynomial_degree_in_prefix_certificates':35,
 'complete_finite_prefix':True,'all_selected_primes_trial_division_checked':True,
 'all_Bezout_identities_reconstructed':True,'numerical_root_scan':False,
 'universal_evidence':'EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex: exact endpoint/tail split, full saddle contour, differentiated analytic errors and log-derivative gap',
 'eventual_simplicity_proved':True,'all_index_exponent_six_two_term_SFC':True,
 'unrestricted_SFC_closed':False,'general_FC_closed':False,'other_Gamma_shapes_claimed':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,'novelty_claim':False}
(HERE/'exponent_six_simplicity_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
files=('EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex','verify_exponent_six_simplicity.py',
       'exponent_six_prefix_certificates.json','exponent_six_simplicity_audit.json',
       'EFFECTIVE_SADDLE_ANNULUS.tex','PRIME_GAPS_AND_SMALL_KERNELS.tex',
       'MIXED_PRIME_RESIDUAL_WINDOWS.tex','QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex')
receipt={'status':'pass','total':len(ineq)+len(identities)+len(records),
 'checks':{'exact_analytic_inequalities':len(ineq),'symbolic_identities':len(identities),'finite_prefix_Bezout_certificates':len(records)},
 'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
 'effective_simplicity_cutoff':N,'finite_prefix_complete':True,
 'all_index_exponent_six_two_term_SFC':True,'unrestricted_SFC_closed':False,'general_FC_closed':False,
 'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'exponent_six_simplicity_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'effective_simplicity_cutoff':N,
                 'finite_prefix':len(records),'modes':modes,'gcd_degrees':gcds,
                 'all_index_exponent_six_two_term_SFC':True,'unrestricted_SFC_closed':False}))
