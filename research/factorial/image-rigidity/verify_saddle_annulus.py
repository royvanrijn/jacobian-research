"""Check exact constants in the analytic saddle-annulus proof, without a root scan.

The universal evidence is the accompanying contour proof. These arithmetic
and algebraic controls are not formal verification or a finite-prefix audit.
"""
from pathlib import Path
from fractions import Fraction as F
from math import factorial
import hashlib,json
import sympy as S
HERE=Path(__file__).resolve().parent
D=F(1,1000);x=F(1,5**6)
values={
 'outer_k_ge3':(4*6**6*F(3,8)**5,F(1400)),
 'outer_k2_square':(F(factorial(12),2),F((500*32)**2)),
 'saddle_Rouche':(F(1),(6-D)**5*D),
 'curvature_real':(F(1,7),1-F(5)/(6-D)),
 'curvature_modulus':((1+D)/(6-D),F(1,5)),
 'third_derivative':((12+216*x*(1+x)/(1-x)**3+108*x/(1-x)**2+12*x/(1-x))/125,F(1)),
 'saddle_phase_loss':(5*D/(6-D)+D,F(1,500)),
 'log6_lower':(F(179,100),2*sum(F(5,7)**(2*j+1)/F(2*j+1) for j in range(12))),
 'saddle_phase_positive':(F(47,10),6*F(179,100)-6-F(1,500)),
 'low_tail_absolute':(F(1000,999)+(4+D)**6*F(3,8)**4,F(83)),
 'low_tail_logarithm':(F(83),F(8,3)**4*sum(F(3,5)**k/F(factorial(k)) for k in range(5))),
 'large_x_perturbation':(3*D**2/16+F(1,4**6),F(1,4000)),
 'medium_tail_gap':(F(1,2500),F(3,64)*F(6,25)**2-F(1,4000)-F(1,500)),
 'far_tail_gap':(F(1,10),F(3,16)-F(1,4000)-F(1,500)),
 'local_remainder':(F(1,6)/F(5,168)**2,F(189)),
 'final_error_scaled':(F(245,1000)+F(8*2500**2,10**9)+F(4,1000)+F(800,10**9)+F(1,4*10**6),F(1)),
 'principal_saddle_circle':(F(2800),F(21,4)**5*F(3,4)),
 'small_saddle_circle':(F(2800),F(23,5)**5*F(7,5)),
 'saddle_lower_radius':(F(1,2)**5*F(13,2),F(1)),
 'saddle_pair_separation':(F(14**28),F(6**6*325*10**26)),
 'saddle_nondegeneracy':(F(1,4000),F(325,7*6*7**4*12)),
 'saddle_phase_derivative_gap':(F(1,10**18),F(1,10**13*6*2800)),
 'endpoint_phase_derivative_gap':(F(1,40000),F(1,2*6*2800)),
 'e_lower_from_partial_series':(F(8,3),sum(F(1,factorial(j)) for j in range(5)))}
for name,(lower,upper) in values.items():assert lower<upper,name
n,k=S.symbols('n k',integer=True,positive=True)
s,w,t=S.symbols('s w t');phi=S.log(w+s**6)-s;ws=s**5*(6-s)
identities={
 'Appell_derivative_coefficients':S.simplify((n+1-k)*S.binomial(n+1,k)-(n+1)*S.binomial(n,k)),
 'saddle_map_derivative':S.diff(ws,s)-6*s**4*(5-s),
 'phase_stationarity':S.diff(phi,s).subs(w,ws),
 'phase_curvature':S.diff(phi,s,2).subs(w,ws)-(5-s)/s,
 'saddle_phase_w_derivative':S.diff(phi,w).subs(w,ws)-(6-s)/(6*ws),
 'saddle_phase_derivative_difference':((6-s)-(6-t))/(6*w)-(t-s)/(6*w),
 'endpoint_phase_derivative_difference':1/w-(6-s)/(6*w)-s/(6*w),
 'saddle_discriminant':S.discriminant(s**6-6*s**5+w,s)-6**6*w**4*(3125-w),
 'local_integral_constant':S.Rational(1,6)/S.Rational(5,168)**2-S.Rational(4704,25),
 'factorial_ratio_constant':4*6**6*S.Rational(3,8)**5-S.Rational(177147,128)}
for name,value in identities.items():assert S.simplify(value)==0,name
# Certify logarithmic and exponential bounds by exact finite expressions;
# no floating-point evaluation or polynomial-zero sampling enters these checks.
audit={'inequalities':{name:{'strict_lower':str(lo),'strict_upper':str(hi)} for name,(lo,hi) in values.items()},
       'symbolic_identities':list(identities),'n_min_inner_disk':10**6,'exterior_radius':2800,
       'inner_radius':1,'uniform_saddle_distance':str(F(1,10**13)),
       'universal_evidence':'EFFECTIVE_SADDLE_ANNULUS.tex: complete shifted contour and all tails',
       'root_scan_performed':False,'finite_prefix_verified':False,'eventual_simplicity_proved':False,
       'all_index_SFC_closed':False,'formal_verification':False,'independent_replay':False,'external_review':False}
(HERE/'saddle_annulus_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
files=('EFFECTIVE_SADDLE_ANNULUS.tex','verify_saddle_annulus.py','saddle_annulus_audit.json')
receipt={'status':'pass','total':len(values)+len(identities),
 'checks':{'exact_rational_inequalities':len(values),'symbolic_identities':len(identities)},
 'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
 'inner_disk_cutoff':10**6,'cutoff_is_for_zero_exclusion_only':True,
 'all_roots_in_annulus_above_cutoff':True,'saddle_data_uniformly_separated':True,
 'eventual_simplicity_proved':False,'simplicity_cutoff_extracted':False,'finite_prefix_verified':False,
 'all_index_SFC_closed':False,'formal_verification':False,'independent_replay':False,'external_review':False}
(HERE/'saddle_annulus_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'pass','checks':receipt['total'],'zero_exclusion_cutoff':10**6,
                 'annulus':'1<|z/n^6|<2800','eventual_simplicity_proved':False}))
