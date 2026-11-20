#!/usr/bin/env python3
if not __debug__:raise SystemExit('Refusing -O/-OO: exact assertions required')
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
from collections import defaultdict
import argparse,json
HERE=Path(__file__).resolve().parent
PROOF='FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md'
INPUTS='factorial_uniform_large_gap_saddle_constants.json'
LEDGER='factorial_uniform_large_gap_saddle_source_ledger.json'
RECEIPT='factorial_uniform_large_gap_saddle_verification.json'
def exp_bounds(x,M=80):
 x=Q(x);term=Q(1);total=term
 for j in range(1,M+1):term*=x/j;total+=term
 ratio=x/(M+2);assert 0<=ratio<1
 return total,total+term*x/(M+1)/(1-ratio)
def clean(p):return {e:Q(c)for e,c in p.items()if c}
def add(*ps):
 d=defaultdict(Q)
 for p in ps:
  for e,c in p.items():d[e]+=c
 return clean(d)
def scale(p,c):return clean({e:v*c for e,v in p.items()})
def mul(p,q):
 d=defaultdict(Q)
 for e,c in p.items():
  for f,v in q.items():d[tuple(a+b for a,b in zip(e,f))]+=c*v
 return clean(d)
def var(i):
 e=[0,0];e[i]=1;return {tuple(e):Q(1)}
ONE={(0,0):Q(1)}
def digest(p):return sha256(p.read_bytes()).hexdigest()
def root_path(arg):
 if arg:return arg.resolve()
 for p in [HERE,*HERE.parents]:
  if(p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Scratch replay requires --repo-root')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--write-receipt',action='store_true');args=ap.parse_args();root=root_path(args.repo_root)
 inp=json.loads((HERE/INPUTS).read_text());assert inp=={'all_terminal_gap_cutoff':512,'minimum_gap':20,'threshold_factor':10**12,'threshold_power':10}
 d=20;D=d+1;x=Q(1,D);checks=[]
 def lt(name,a,b):
  a,b=Q(a),Q(b);assert a<b,(name,a,b);checks.append({'name':name,'left':str(a),'relation':'<','right':str(b)})
 def le(name,a,b):
  a,b=Q(a),Q(b);assert a<=b,(name,a,b);checks.append({'name':name,'left':str(a),'relation':'<=','right':str(b)})
 for label,a in [('upper_Rouche',40*D**2),('upper_local_third',14*D**3),('upper_real_curvature',160*D**2),('small_log_ratios',700*D),('main_derivative_separation',300*D)]:lt(label,a,2**d)
 for label,a in [('lower_real_curvature',800*D),('lower_local_third',8*D**2)]:lt(label,a,5**d)
 for row in [
 ('two_power_ratio_j<=3',Q(22,21)**3/2,1),('five_power_ratio_j<=2',Q(22,21)**2/5,1),
 ('W_lower_bound',400*D,7**d),('seven_power_growth',Q(D+1,7*D),1),
 ('Euler_lower',Q(8,3),exp_bounds(1,6)[0]),('Euler_upper',exp_bounds(1,6)[1],Q(87,32)),('Euler_weakening',Q(87,32),Q(11,4)),
 ('local_exp_3_5',exp_bounds(Q(3,5))[1],2),('left_upper_envelope',exp_bounds(Q(9,32))[1],Q(4,3)),
 ('small_log_radius',exp_bounds(Q(1,10))[1],Q(10,9)),('low_saddle_radius',exp_bounds(Q(1,100))[1],Q(4,3)),
 ('lower_far_first',Q(40,9),exp_bounds(Q(3,2))[0]),('lower_far_second_squared',2500*11**21,16**21),
 ('lower_far_margin',Q(87,32)*Q(1000,989),Q(25,9)),
 ('real_inverse_square',Q(99,100),(1-Q(1,10000))/(1+Q(1,10000))**2),
 ('upper_imaginary_ratio',Q(1,10)/Q(3*D,4),Q(1,100)),('lower_imaginary_ratio',Q(1,100*D)/Q(1,2),Q(1,100)),
 ('upper_saddle_curvature_error',Q(1,10)*(1+x/10)/(1-x/10)**2,Q(1,5)),
 ('lower_saddle_curvature_error',Q(1,100)/(1-Q(1,2100))**2,Q(1,10)),
 ('upper_a_derivative',1/(1-x/10)**2+2*x/(1-x/10)**3,2),
 ('lower_a_derivative',Q(2100,2099)**2+2*Q(2100,2099)**3,4),
 ('upper_prefactor_s_derivative',Q(10,209)+Q(5,4),Q(3,2)),
 ('upper_global_curvature_margin',Q(1,3),Q(396,625)-Q(1,10)),('lower_global_curvature_margin',Q(1,3),Q(98,225)),
 ('upper_local_cube_denominator',Q(9,10),Q(97,100)**3),
 ('upper_log_third_factor',Q(10,9)*(1+Q(3,100)+Q(2,10000)),2),
 ('lower_log_third_factor',(1+Q(3,100)+Q(2,10000))/Q(89,100)**3,2),('lower_base_third',2/Q(89,100)**3,3),
 ('upper_leading_lower_squared',4,Q(210,211)**2*5),('upper_leading_upper_squared',Q(210,209)**2*Q(55,7),9),
 ('lower_leading_lower_squared',4,Q(2100,2101)**2*Q(60,11)),('lower_leading_upper_squared',Q(2100,2099)**2*Q(440,63),9),
 ('upper_height',Q(3,500*D)+Q(5*D,2**d),Q(1,100)),('lower_height',Q(1,10000*D**2)+Q(5,400*D**2),Q(1,100)),
 ('upper_left_low_term_base',Q(3,2)*Q(87,32)*D,exp_bounds(Q(D,4))[0]),('upper_left_monotonicity',Q(1,D),Q(1,4)),
 ('lower_central_ratio',Q(1501,1000)*Q(11,84),Q(1,5)),('lower_local_ratio',Q(1101,1000)*Q(11,84),Q(1,5)),
 ('lower_origin_base',Q(1,2)*Q(101,100)/Q(8,5),Q(8,25)),('lower_vertical_base',Q(501,1000)*Q(101,100)/Q(8,5),Q(33,100)),
 ('lower_origin_endpoint_loss',Q(11,4)/Q(49,50),Q(100,33)),('lower_far_imaginary',Q(1,45000*D),Q(1,10000)),
 ('lower_left_integrated_amplitude',Q(2,100*D),1),
 ('overlap_decimal',10**12,2**40),('overlap_power_two',513,2**10),('overlap_log_two',2**4,Q(8,3)**3),
 ('overlap_log_margin',Q(435,4),Q(511,4)),('overlap_derivative',Q(10,513),Q(1,4))]:lt(*row)
 for row in [
 ('lower_prefactor_s_derivative',2+Q(20*D,9),3*D),('upper_far_derivative',-Q(1,5)+Q(1,100),-Q(1,10)),
 ('upper_local_decay',Q(1,3),Q(2,5)-Q(1,D)),('lower_local_decay',Q(1,3),Q(9,20)-Q(1,10)),
 ('upper_left_log_sum',1/(1+Q(7*D,32)+Q(7*D,32)**2/2),Q(D,160)),
 ('upper_left_endpoint_loss',Q(1,100)+Q(1,50*D),Q(D,200))]:le(*row)
 v=Q(D)-Q(1,10);rho=Q(1,2);cc=Q(1,3*D);B=Q(2,D**2)
 lt('upper_local_error_bracket',B/((v-rho)*cc**2)+1/(v*(v-rho)*cc)+1/(v*cc*rho),9)
 v=Q(999,1000);rho=Q(1,10);cc=Q(1,3)
 lt('lower_local_error_bracket',1/((v-rho)*cc**2)+1/(v*(v-rho)*cc)+1/(v*cc*rho),46)
 N0=10**12*D**10
 le('upper_tail_threshold',96*24**4*D**4,N0**3);le('lower_tail_threshold',96*2000**4*D,N0**3)
 lt('common_relative_error',Q(1,10000*D**3),Q(1,100));lt('log_error_Cauchy',Q(13,10000*D**3),1)
 le('zero_free_threshold',12*D,N0);lt('main_derivative_gap',Q(1,4),Q(1,3)-Q(1,100));lt('main_derivative_order',16,N0)
 DD,ss=var(0),var(1);dm=add(DD,scale(ONE,-1))
 stationary=add(mul(add(dm,scale(ONE,-1)),add(DD,scale(ss,-1))),scale(add(ss,scale(ONE,-1)),-1),scale(mul(dm,add(dm,scale(ss,-1))),-1));assert stationary=={}
 lhs=add(mul(DD,add(mul(ss,ss),scale(mul(DD,ss),-1),DD)),scale(mul(ss,ss),-1))
 rhs=mul(add(ss,scale(DD,-1)),add(mul(add(DD,scale(ONE,-1)),ss),scale(DD,-1)));assert lhs==rhs
 assert add(mul(ss,ss),scale(mul(DD,ss),-1),DD,scale(mul(ss,ss),-1))==scale(mul(DD,add(ss,scale(ONE,-1))),-1)
 assert not(Q(396,625)-Q(1,2)>Q(1,3));assert not(2500*11**20<16**20);assert not(Q(9,20)-Q(1,5)>=Q(1,3))
 source=json.loads((HERE/LEDGER).read_text());pins=[]
 mandatory={('packet',PROOF),('packet',Path(__file__).name),('packet',INPUTS),
 ('repository','research/factorial/image-rigidity/FACTORIAL_HIGH_GAP_ROUCHE_REGION.md'),
 ('repository','research/factorial/image-rigidity/GAP_SIX_RESULTANT_STRICT_POSITIVITY.md'),
 ('repository','research/factorial/image-rigidity/BINOMIAL_SHAPE_RESULTANT.tex')}
 keys={(v['location'],v['path'])for v in source['files']}
 assert mandatory<=keys and len(keys)==len(source['files']),'missing/duplicate essential proof or source pin'
 for item in source['files']:
  p=(root if item['location']=='repository'else HERE)/item['path'];h=digest(p);assert h==item['sha256'],str(p);pins.append({'path':item['path'],'sha256':h})
 result={'status':'PASS','scope':'Exact base and monotonicity ledger supports the written uniform contour proof; no numerical root/contour scan or missing small-gap completion inferred.','minimum_gap':20,'terminal_threshold':'10^12*(d+1)^10','all_terminal_gap_cutoff':512,'rational_inequality_count':len(checks),'phase_identity_count':3,'inequalities':checks,'negative_controls':['weakened_upper_curvature','D20_far_envelope_not_certified','oversized_lower_local_radius'],'verified_source_pins':pins,'proof_sha256':digest(HERE/PROOF),'checker_sha256':digest(Path(__file__).resolve()),'inputs_sha256':digest(HERE/INPUTS),'ledger_sha256':digest(HERE/LEDGER)}
 out=args.output or HERE/RECEIPT;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print('PASS',len(checks),'rational inequalities;3 phase identities;3 corruptions rejected');print(out)
if __name__=='__main__':main()
