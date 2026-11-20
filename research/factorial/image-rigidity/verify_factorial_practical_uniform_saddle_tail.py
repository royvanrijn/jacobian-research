#!/usr/bin/env python3
if not __debug__:raise SystemExit('Refusing -O/-OO: exact assertions required')
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse,json
HERE=Path(__file__).resolve().parent
PROOF='ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md'
INPUTS='factorial_practical_uniform_saddle_constants.json'
LEDGER='factorial_practical_uniform_saddle_source_ledger.json'
RECEIPT='factorial_practical_uniform_saddle_verification.json'
def digest(p):return sha256(p.read_bytes()).hexdigest()
def exp_bounds(x,M=100):
 x=Q(x);t=Q(1);s=t
 for j in range(1,M+1):t*=x/j;s+=t
 r=x/(M+2);assert 0<=r<1
 return s,s+t*x/(M+1)/(1-r)
def add(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,Q(0))+c
 return {e:c for e,c in out.items()if c}
def mul(p,q):
 out={}
 for e,c in p.items():
  for f,v in q.items():out[e+f]=out.get(e+f,Q(0))+c*v
 return {e:c for e,c in out.items()if c}
def scale(p,c):return {e:v*c for e,v in p.items()if v*c}
def find_root(arg):
 if arg:return arg.resolve()
 for p in [HERE,*HERE.parents]:
  if(p/'research/MATH_STATUS.json').is_file():return p
 raise SystemExit('Use --repo-root for scratch replay')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path);ap.add_argument('--write-receipt',action='store_true');args=ap.parse_args();root=find_root(args.repo_root)
 inp=json.loads((HERE/INPUTS).read_text());assert inp=={'minimum_gap':10,'absolute_terminal_threshold':5000,'linear_terminal_factor':128,'all_terminal_gap_cutoff':44,'finite_box_maximum_gap':43,'finite_box_pair_count':171812,'finite_box_terminal_ceiling_exclusive':5632}
 checks=[]
 def lt(name,a,b):
  a,b=Q(a),Q(b);assert a<b,(name,a,b);checks.append({'name':name,'left':str(a),'relation':'<','right':str(b)})
 def le(name,a,b):
  a,b=Q(a),Q(b);assert a<=b,(name,a,b);checks.append({'name':name,'left':str(a),'relation':'<=','right':str(b)})
 D=11;N=5000;delta=Q(1,3);q=Q(3,8)
 for row in [
 ('third_derivative_raw_upper',14*D**3*q**10,2),('local_remainder_budget',5,6),('local_decay',Q(1,3),Q(2,5)-Q(1,2*D)),
 ('bracket_monotone',0,5*D-3),('upper_bracket',Q(375,52),Q(29,4)),('lower_bracket',Q(46),Q(47)),
 ('central_upper_gaussian_amplitude',Q(4,3*D),Q(2,D)),
 ('lower_far_loss',Q(87,32)*Q(1000,979),Q(25,9)),('lower_imaginary_allowance',Q(1,45000*D),Q(1,1000)),
 ('lower_far_log_amplitude',Q(D,4),D),('lower_total_remaining_amplitude',2+D+1,4*D),
 ('sqrt128_lower',100,128),('sqrt5000_upper',N,71**2),('sqrt5000_above60',60**2,N),
 ('upper_main_relative_squared',(Q(29,8)+8*q**5)**2,128*Q(33,100)**2),
 ('upper_left_relative',Q(12000,N**3),Q(1,1000)),('upper_far_relative',Q(400000,(N*D)**2),Q(1,1000)),
 ('upper_total_relative',Q(33,100)+Q(2,1000),delta),
 ('lower_local_relative_squared',Q(23**2,N),Q(33,100)**2),('lower_central_relative',q**8,Q(1,2000)),
 ('lower_far_decreasing',Q(3,2*N),Q(1,100)),('lower_far_relative',Q(N*71,64)*q**50,Q(1,1000)),
 ('lower_total_relative',Q(33,100)+Q(1,2000)+Q(1,1000),delta),
 ('leading_magnitude_lower',1,2*(1-delta)),
 ('zero_free_exponential',Q(25,49),Q(2,3)),('zero_free_base',2*Q(2,3)**2,1),
 ('separation_threshold',56,N),('main_reciprocal_margin',Q(1,4),Q(1,3)-Q(1,100)),
 ('small_exp_series_bound',exp_bounds(Q(1,50))[1],Q(50,49)),
 ('sparse_overlap',46080,exp_bounds(Q(1936,180))[0]),
 ('overlap_derivative',0,Q(1,4)-Q(1,4*45**2)-Q(1,45)),
 ('overlap_derivative_increasing',4*45+1,45**2),('overlap_linear_regime',5000,128*45),
 ('canonical_large_N_not_needed',N,10**12*D**10)]:lt(*row)
 for row in [
 ('lower_far_exp_geometric',Q(21,1000),Q(1)-Q(979,1000)),
 ('upper_central_exponent',5,Q(128,24)),('lower_central_exponent',8,Q(N,600)),
 ('coarse_upper_magnitude',3*(1+delta),4),('Cauchy_logerror',12*delta/(1-delta),6),
 ('logderivative_allowance',6+1,7),('zero_free_monotone_ratio',Q(3,2)*Q(2,3),1)]:le(*row)
 # Exact local bracket identity after clearing denominators, with D an indeterminate.
 z={1:Q(1)};one={0:Q(1)};v=add(z,{0:-Q(1,10)});vm=add(z,{0:-Q(3,5)})
 lhs=add(scale(v,9),scale(z,3),scale(mul(z,vm),6));rhs=mul(v,add(scale(z,6),{0:Q(9)}));assert lhs==rhs
 # Same numerator is (6+63/(5D-3)) after cancelling v.
 assert add(scale(vm,6),{0:Q(63,5)})==add(scale(z,6),{0:Q(9)})
 # The Gaussian-tail coefficients retain the amplitude1/s.
 assert Q(4,3*D)*Q(12*D)==16
 assert 2*Q(6)/Q(1,10)==120
 # Exact overlap derivative numerator D²-4D-1.
 assert add(mul(add(z,{0:-1}),add(z,{0:1})),scale(z,-4))=={2:1,1:-4,0:-1}
 bracket=Q(1,D**2)/((Q(D)-Q(3,5))*Q(1,3*D)**2)+1/((Q(D)-Q(1,10))*(Q(D)-Q(3,5))*Q(1,3*D))+1/((Q(D)-Q(1,10))*Q(1,3*D)*Q(1,2));assert bracket==Q(375,52)
 n0=lambda d:max(5000,128*(d+1))
 assert n0(44)==5760
 assert sum(n0(d)-2 for d in range(10,44))==171812
 assert max(n0(d)for d in range(10,44))==5632
 # Rejected proof budgets (not numerical-root counterexamples).
 # Old central length-only estimate is too weak at N=128D,D40.
 assert exp_bounds(Q(16,3))[1]<300
 assert Q(20*40,300)>delta
 # A coarse e<11/4 no longer certifies the new1/100 lower far loss.
 assert not(Q(11,4)*Q(1000,979)<Q(25,9))
 # Delta1/2 does not retain the prefactor upper bound4 used for zero-free regions.
 assert not(3*(1+Q(1,2))<=4)
 data=json.loads((HERE/LEDGER).read_text());mandatory={('packet',PROOF),('packet',Path(__file__).name),('packet',INPUTS),('repository','research/factorial/image-rigidity/FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md'),('repository','research/factorial/image-rigidity/FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md'),('repository','research/factorial/image-rigidity/lower_exponent_initial_and_boundary.tex'),('repository','research/factorial/image-rigidity/GAP_SIX_RESULTANT_STRICT_POSITIVITY.md'),('repository','research/factorial/image-rigidity/ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex')}
 keys={(v['location'],v['path'])for v in data['files']};assert mandatory<=keys and len(keys)==len(data['files'])
 pins=[]
 for item in data['files']:
  p=(root if item['location']=='repository'else HERE)/item['path'];h=digest(p);assert h==item['sha256'],str(p);pins.append({'path':item['path'],'sha256':h})
 result={'status':'PASS','scope':'Exact proof-ledger support for the full-contour practical tail; the actual finite global discriminant/resultant box is not checked here.','minimum_gap':10,'terminal_onset':'max(5000,128*(d+1))','all_terminal_gap_cutoff':44,'finite_box_pair_count':171812,'finite_box_terminal_ceiling_exclusive':5632,'rational_comparison_count':len(checks),'inequalities':checks,'symbolic_identity_count':5,'negative_controls':['old_length_only_upper_central_bound','coarse_Euler_far_loss_budget','oversized_relative_error'],'proof_sha256':digest(HERE/PROOF),'checker_sha256':digest(Path(__file__).resolve()),'inputs_sha256':digest(HERE/INPUTS),'ledger_sha256':digest(HERE/LEDGER),'verified_source_pins':pins}
 out=args.output or HERE/RECEIPT;out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print('PASS',len(checks),'exact comparisons;5 symbolic identities;3 rejected budgets');print(out)
if __name__=='__main__':main()
