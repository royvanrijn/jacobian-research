#!/usr/bin/env python3
"""Exact supporting corner/Young/monotonicity guards for the relative-gap tail."""
from __future__ import annotations
import argparse,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
if not __debug__:raise SystemExit('Refusing -O/-OO: exact assertions are required.')
HERE=Path(__file__).resolve().parent
ROWS=[]
def text(x):
 x=F(x);return str(x.numerator)if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def check(label,l,r,kind='<'):
 l,r=F(l),F(r);ok=l<r if kind=='<'else l==r
 assert ok,(label,text(l),kind,text(r));ROWS.append({'label':label,'left':text(l),'relation':kind,'right':text(r)})
def exp_bounds(x):
 x=F(x)
 if x<0:
  if -x>100:return F(0),1/exp_bounds(F(100))[0]
  lo,hi=exp_bounds(-x);return 1/hi,1/lo
 term=total=F(1)
 for j in range(1,257):term*=x/j;total+=term
 q=x/258;assert q<1
 return total,total+term*x/257/(1-q)
def eu(x):return exp_bounds(x)[1]
def el(x):return exp_bounds(x)[0]
def root_bounds(x):
 x=F(x);m=10**12;k=math.isqrt(x.numerator*m*m//x.denominator)
 assert F(k*k,m*m)<=x<F((k+1)**2,m*m)
 return F(k,m),F(k+1,m)
def sl(x):return root_bounds(x)[0]
def su(x):return root_bounds(x)[1]
def bernstein(p,l,h):
 n=len(p)-1;a=[F(0)]*(n+1)
 for k,c in enumerate(p):
  for j in range(k+1):a[j]+=c*math.comb(k,j)*l**(k-j)*(h-l)**j
 return [sum((a[i]*F(math.comb(k,i),math.comb(n,i))for i in range(k+1)),F(0))for k in range(n+1)]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'relative_gap_three_weighted_tail_verification.json');args=ap.parse_args();repo=args.repo_root.resolve()
 # Exact logarithm bounds used by partial-derivative arguments, avoiding floats.
 check('log4<3/2',4,el(F(3,2)))
 check('log4>4/3',eu(F(4,3)),4)
 check('1+log5/8>1/2',F(8,5),el(F(1,2)))
 check('upper left F negative',4,F(8,3)**3*F(5,8)**3)
 check('upper left q<1 at crude corner',F(6,5)*16*F(3,5)**6,1)
 check('upper-left h rho derivative',F(1,12)-F(1,4),0)
 check('upper central a derivative margin',F(1,2)+F(3,2)-F(18,7),0)
 check('upper central rho derivative margin',F(2,3)-F(17,14),0)
 check('lower strip F negative',F(7,4)**3*F(87,32)**3,F(4)**4)
 check('lower strip log7/16<-1/2',eu(F(1,2)),F(16,7))
 check('lower strip h rho derivative',F(1,3)-F(4,5),0)
 check('lower local a derivative margin',1-F(47,24),0)
 check('lower local rho derivative margin',1-F(29,12),0)
 check('lower central a derivative margin',F(3,5)-F(941,600),0)
 check('lower central rho derivative margin',F(11,15)-F(647,300),0)
 # Analytic branches; no finite parameter grid is substituted for monotonicity.
 check('upper branch corner',F(85761906966528,93776847890045),1)
 check('upper branch a derivative',F(1,2)+F(3,2)-3,0)
 check('upper branch rho derivative',F(1,3)-F(3,2)+F(2,15),0)
 lowerR=F(5,2)*4*F(49,8)*F(5,3)*F(87,32)**6/F(4)**8
 check('lower branch uniform cap',lowerR,1)
 check('lower branch a derivative',F(1,2)+3-F(16,3),0)
 check('lower branch rho derivative',F(1,4)+F(1,3)-F(8,3),0)
 q=32*F(60,163)**6
 check('upper uniform stationary height',q/(1-q),F(7,80))
 check('lower uniform stationary height',F(1,240)+F(1,46),F(1,32))
 check('lower curvature norm difference',F(18,25),1)
 # New local majorant and its polynomial-growth derivatives.
 q=F(96,5)*F(60,163)**5;y=q/(1-q)
 check('uniform upper local third derivative',2*F(8,7)**3+64*(y+3*y*y+2*y**3),18)
 check('D2q upper a derivative',1+F(3,2)-3,0)
 check('D2q upper rho derivative',F(1,2)+F(1,2)-2,0)
 # Gaussian-amplitude Young constants, valid for arbitrarily large a.
 Cp=F(24,23)**2*F(511,475);Ca=F(24,23)*F(511,475)
 check('Young lower exponential',eu(F(36,511)),F(511,475))
 for a in [F(2),F(125,16)]:
  e=48*Cp/sl(3*3999*a)+6*Ca*su(a/F(3*3999))
  check('lower weighted local endpoint a='+text(a),e,F(9,20))
 check('upper Young phase constant',F(24,23)**2*F(127968,85968),F(5,3))
 check('upper Young amplitude constant',F(24,23)*F(127968,85968),F(8,5))
 check('upper Young complete local budget',F(5,3)*F(1000,999)*F(33,40)+F(8,5)*F(5,8),F(5,2))
 # Growing left contour lengths are controlled by the same onset.
 check('left length prefactor',su(F(3,4*512)),F(1,16))
 check('unweighted lower left',F(4000,16)*eu(-40),F(1,1000))
 check('weighted lower left exponent',-F(4000,100)+F(4000,512)+F(1,8)+F(1,100),-F(4000,128))
 check('weighted lower left',F(4000,16)*eu(-F(4000,128)),F(1,1000))
 check('new lower weighted complete sum',F(9,20)+F(5,1000)+F(165,1000)+F(3,1000),F(18,25))
 check('old upper weighted complete sum',F(5,2)+F(1,50)+F(2,5)+F(1,20)+F(2,1000),3)
 check('sparse overlap d64',F(4096*64,3),F(8,3)**12)
 check('sparse overlap derivative',0,F(3,16)-F(1,64))
 # The three corner certificates are re-expanded independently coefficientwise.
 oldconst=repo/'research/factorial/image-rigidity/lower_two_weighted_saddle_constants.json'
 old=json.loads(oldconst.read_text());count=0
 for c in old['bernstein_certificates']:
  p=[F(v)for v in c['power_coefficients']];l,h=map(F,c['interval']);bs=bernstein(p,l,h)
  assert bs==[F(v)for v in c['bernstein_coefficients']] and all(v>0 for v in bs)
  count+=len(bs)
 assert count==47
 # Finite-domain statements here are counts only, not coprimality claims.
 supports=[(a,d)for d in range(9,64)for a in range(3,d//3+1)]
 assert len(supports)==532
 assert sum(max(4000,128*(a+d))-2 for a,d in supports)==3694040
 assert max(max(4000,128*(a+d))-1 for a,d in supports)==10751
 assert max((a+d)*(max(4000,128*(a+d))-1)for a,d in supports)==903084
 # Negative controls: the naive fixed-radius amplitude cap is not uniform;
 # omitting the -1 in d(d-1)'s derivative changes a strict exact inequality.
 assert el(100/F(16))>100
 actual_log_derivative=F(9,9*2-3)
 assert actual_log_derivative==F(3,5)>F(1,2)
 proof=HERE/'RELATIVE_GAP_THREE_LOWER_WEIGHTED_SADDLE_TAIL.md';assert proof.is_file()
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 files=[proof,Path(__file__).resolve(),oldconst,repo/'research/factorial/image-rigidity/ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md']
 out={'status':'PASS','scope':'Exact support for uniform analytical tail and finite-domain counts only; no execution of the3694040 rational coprimality obligation.','counts':{'rational_guards':len(ROWS),'reexpanded_Bernstein_coefficients':47,'support_count':532,'terminal_count':3694040,'negative_controls':2},'sources':[{'path':str(p.relative_to(repo))if p.is_relative_to(repo)else p.name,'sha256':sha(p)}for p in files],'comparisons':ROWS}
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','counts':out['counts']}))
if __name__=='__main__':main()
