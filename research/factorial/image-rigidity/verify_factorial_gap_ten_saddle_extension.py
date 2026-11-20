#!/usr/bin/env python3
if not __debug__: raise SystemExit('Refusing -O/-OO: exact assertions required')
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse,json
HERE=Path(__file__).resolve().parent
PROOF='FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md'
CONSTANTS='factorial_gap_ten_saddle_constants.json'
LEDGER='factorial_gap_ten_saddle_source_ledger.json'
RECEIPT='factorial_gap_ten_saddle_verification.json'
def digest(p):return sha256(p.read_bytes()).hexdigest()
def exp_bounds(x,M=80):
 x=Q(x);t=Q(1);s=t
 for j in range(1,M+1):t*=x/j;s+=t
 ratio=x/(M+2);assert 0<=ratio<1
 return s,s+t*x/(M+1)/(1-ratio)
def find_root(arg):
 if arg:return arg.resolve()
 for p in [HERE,*HERE.parents]:
  if(p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Use --repo-root for scratch replay')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--write-receipt',action='store_true');args=ap.parse_args();root=find_root(args.repo_root)
 inp=json.loads((HERE/CONSTANTS).read_text());assert inp=={'minimum_gap':10,'split_numerator':3,'split_denominator':8,'threshold_factor':10**12,'threshold_power':10,'finite_box_maximum_gap':511}
 d=10;D=11;q=Q(3,8);e_lo=Q(163,60);e_hi=Q(87,32);checks=[]
 def lt(name,a,b):
  a,b=Q(a),Q(b);assert a<b,(name,a,b);checks.append({'name':name,'left':str(a),'relation':'<','right':str(b)})
 def le(name,a,b):
  a,b=Q(a),Q(b);assert a<=b,(name,a,b);checks.append({'name':name,'left':str(a),'relation':'<=','right':str(b)})
 for row in [
 ('Euler_lower_series',e_lo,exp_bounds(1,6)[0]),('Euler_upper_series',exp_bounds(1,6)[1],e_hi),('W_base',e_hi,Q(11,4)),
 ('W_threshold',400*D,4**d),('W_growth_ratio',Q(D+1,4*D),1),
 ('q_polynomial_ratio',Q(12,11)**3*q,1),('half_polynomial_ratio',Q(12,11)**3/2,1),('quarter_polynomial_ratio',Q(12,11)**3/4,1),('lower_power_ratio',Q(12,11)**2*Q(2,5),1),
 ('upper_Rouche',40*D**2*q**d,1),('upper_local_log_third',14*D**3*q**d,2),('upper_separation',300*D*q**d,1),
 ('lower_curvature',800*D*Q(2,5)**d,1),('lower_third',8*D**2*Q(2,5)**d,1),
 ('upper_local_y_envelope',42*D*q**d,1),('upper_local_y_small',7*D*q**d,Q(1,100)),('upper_local_modulus',Q(47,50),1-Q(3,5*D)),
 ('upper_third_correction_factor',(Q(12,11)+Q(3,100)+Q(2,10000))*Q(50,47)**3,2),('upper_third_base',2*Q(50,47)**3,3),('upper_third_remainder',5,18),
 ('upper_local_height',Q(3,500*D)+5*D*q**d,Q(1,100)),('lower_local_height',Q(1,10000*D**2)+Q(5,400*D**2),Q(1,100)),
 ('upper_saddle_curvature_error',Q(1,10)*(1+Q(1,10*D))/(1-Q(1,10*D))**2,Q(1,5)),('lower_saddle_curvature_error',Q(1,100)/(1-Q(1,100*D))**2,Q(1,10)),
 ('upper_local_exp',exp_bounds(Q(3,5))[1],2),('upper_Rouche_exp',exp_bounds(Q(1,10))[1],Q(10,9)),('lower_Rouche_exp',exp_bounds(Q(1,100))[1],Q(4,3)),
 ('upper_leading_lower_squared',4,Q(110,111)**2*5),('upper_leading_upper_squared',Q(110,109)**2*Q(55,7),9),
 ('lower_leading_lower_squared',4,Q(1100,1101)**2*Q(60,11)),('lower_leading_upper_squared',Q(1100,1099)**2*Q(440,63),9),
 ('upper_central_imaginary_ratio',Q(2,165),Q(1,80)),('upper_inverse_square',Q(99,100),(1-Q(1,80**2))/(1+Q(1,80**2))**2),
 ('lower_imaginary_ratio',Q(1,100*D)/Q(1,2),Q(1,100)),('lower_inverse_square',Q(99,100),(1-Q(1,10000))/(1+Q(1,10000))**2),
 ('upper_central_curvature',Q(1,3),Q(59,100)*Q(16,25)),('upper_far_y',3*D*Q(3,10)**d/(1-3*D*Q(3,10)**d),Q(1,100)),
 ('upper_left_imaginary',Q(8,225),Q(1,25)),
 ('upper_left_first_high_squared',(e_hi/4)**11,Q(1,64)),('upper_left_second_high_fourth',Q(3,4)**44*e_hi**11,Q(2,3)**4),
 ('upper_left_first_low_exp',22,e_lo**3*Q(9,8)),('upper_left_second_low_exp_squared',88**2,e_lo**9),
 ('first_left_low_monotonicity',Q(12,11)/(1+Q(3,8)),1),('second_left_low_monotonicity',Q(12,11)/(1+Q(1,2)),1),
 ('first_left_sum',Q(1,8)+Q(9,16),Q(7,8)),('second_left_sum',Q(2,3)+Q(3,16),Q(7,8)),
 ('left_height_loss_exp',Q(44,39),Q(8,7)),('left_integrated_amplitude',Q(3,4)+Q(4,15*D),1),('left_log2',2**4,Q(8,3)**3),
 ('lower_local_ratio',Q(1101,4000),Q(2,5)),('lower_central_ratio',Q(1501,4000),Q(2,5)),
 ('lower_log_third_factor',(1+Q(3,100)+Q(2,10000))/Q(89,100)**3,2),('lower_base_third',2/Q(89,100)**3,3),
 ('lower_far_first_exp',Q(40,9),exp_bounds(Q(3,2))[0]),
 ('lower_far_second_base_eighth',Q(3,8)**8*e_hi**5,1),
 ('lower_far_second_eighth',Q(9,8)**8*Q(3,8)**88*e_hi**55,Q(9,400)**8),
 ('lower_far_margin',e_hi*Q(1000,989),Q(25,9)),('lower_far_imaginary',Q(1,45000*D),Q(1,10000)),
 ('lower_origin_base',Q(1,2)*Q(101,100)/Q(8,5),Q(8,25)),('lower_vertical_base',Q(501,1000)*Q(101,100)/Q(8,5),Q(33,100)),
 ('lower_endpoint_loss',Q(11,4)/Q(49,50),Q(100,33)),('lower_vertical_integrated',Q(2,100*D),1),
 ('upper_a_derivative',1/(1-Q(1,10*D))**2+Q(2,D)/(1-Q(1,10*D))**3,2),
 ('lower_a_derivative',Q(1100,1099)**2+2*Q(1100,1099)**3,4),('upper_prefactor_s_derivative',Q(10,109)+Q(5,4),Q(3,2)),
 ('lower_prefactor_s_derivative',2+Q(20*D,9),3*D),('main_derivative_gap',Q(1,4),Q(1,3)-Q(1,100))]:lt(*row)
 r=Q(3*D,2**d);y=r/(1-r);assert D*(y+y*y)==Q(371712,982081);lt('upper_central_factor',D*(y+y*y),Q(2,5))
 le('upper_local_decay',Q(1,5),Q(2,5)-Q(3,2*D));le('lower_local_decay',Q(1,3),Q(9,20)-Q(1,10));le('left_height_exponent',Q(11,100)+Q(1,25*D),Q(5,44))
 v=Q(D)-Q(1,10);rho=Q(1,2);cc=Q(1,5*D);B=Q(3,D**2)
 lt('upper_local_error',B/((v-rho)*cc**2)+1/(v*(v-rho)*cc)+1/(v*cc*rho),20)
 v=Q(999,1000);rho=Q(1,10);cc=Q(1,3)
 lt('lower_local_error',1/((v-rho)*cc**2)+1/(v*(v-rho)*cc)+1/(v*cc*rho),46)
 N=10**12*D**10
 le('upper_tail_threshold',96*24**4*D**4,N**3);le('lower_tail_threshold',96*2000**4*D,N**3)
 lt('relative_error',Q(1,10000*D**3),Q(1,100));lt('error_Cauchy',Q(13,10000*D**3),1);le('zero_free_threshold',12*D,N);lt('zero_separation_threshold',16,N)
 # Deliberate changes fail; no floating root/contour test is used.
 bad_r=Q(4*D,2**d);bad_y=bad_r/(1-bad_r);assert not(D*(bad_y+bad_y*bad_y)<Q(2,5))
 assert not(800*10*Q(2,5)**9<1)
 assert not(Q(6,5)**8*Q(3,8)**88*e_hi**55<Q(9,400)**8)
 data=json.loads((HERE/LEDGER).read_text());mandatory={('packet',PROOF),('packet',Path(__file__).name),('packet',CONSTANTS),('repository','research/factorial/image-rigidity/FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md'),('repository','research/factorial/image-rigidity/ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex'),('repository','research/factorial/image-rigidity/GAP_SIX_RESULTANT_STRICT_POSITIVITY.md'),('repository','research/factorial/image-rigidity/BINOMIAL_SHAPE_RESULTANT.tex')}
 keys={(v['location'],v['path'])for v in data['files']};assert mandatory<=keys and len(keys)==len(data['files'])
 pins=[]
 for item in data['files']:
  p=(root if item['location']=='repository'else HERE)/item['path'];h=digest(p);assert h==item['sha256'],str(p);pins.append({'path':item['path'],'sha256':h})
 result={'status':'PASS','scope':'Exact changed-bound/monotonicity ledger for written all-contour d>=10 proof; finite all-terminal box remains unevaluated.','minimum_gap':10,'terminal_threshold':'10^12*(d+1)^10','finite_all_terminal_box':'10<=d<=511,2<=N<10^12*(d+1)^10','rational_inequality_count':len(checks),'inequalities':checks,'negative_controls':['loosened_central_y','unsupported_d9_K_bound','oversized_far_K_constant'],'proof_sha256':digest(HERE/PROOF),'checker_sha256':digest(Path(__file__).resolve()),'constants_sha256':digest(HERE/CONSTANTS),'ledger_sha256':digest(HERE/LEDGER),'verified_source_pins':pins}
 out=args.output or HERE/RECEIPT;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print('PASS',len(checks),'exact rational comparisons;3 negative controls');print(out)
if __name__=='__main__':main()
