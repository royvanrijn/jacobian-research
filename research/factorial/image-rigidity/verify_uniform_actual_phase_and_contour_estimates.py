#!/usr/bin/env python3
"""Exact support for the written uniform actual-phase bridge, not sampled proof."""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO: exact checks must execute.')
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse,json,math
HERE=Path(__file__).resolve().parent
ROWS=[]
def ck(name,left,right,op='<'):
 left,right=Q(left),Q(right)
 assert(left<right if op=='<'else left==right),(name,left,op,right)
 ROWS.append(dict(name=name,left=str(left),relation=op,right=str(right)))
def eu(x):
 x=Q(x);assert 0<=x<65
 term=total=Q(1)
 for k in range(1,65):term*=x/k;total+=term
 return total+term*x/65/(1-x/66)
def add(a,b):
 out=dict(a)
 for k,v in b.items():out[k]=out.get(k,Q(0))+v
 return{k:v for k,v in out.items()if v}
def scale(a,c):return{k:v*c for k,v in a.items()if v*c}
def mul(a,b):
 out={}
 for i,u in a.items():
  for j,v in b.items():
   k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,Q(0))+u*v
 return{k:v for k,v in out.items()if v}
def var(n,i):return{tuple(int(i==j)for j in range(n)):Q(1)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'uniform_actual_phase_and_contour_verification.json');args=ap.parse_args()
 outer=Q(1,10**5);actual=Q(1,10**10);R=Q(2)
 ck('outer stationary perturbation coefficient',(R+outer)*R*R*eu(outer*R)/2+R,Q(601,100))
 R=Q(3,2);ck('selected circle perturbation coefficient',(R+outer)*R*R*eu(outer*R)/2+R,Q(16,5))
 ck('outer selected two-root count',Q(16,5)*outer,Q(17,200))
 nu=Q(601,100)*outer/Q(39,50);ck('outer log ratio',nu,Q(78,10**6))
 logL=2*nu/(1-nu);derL=4*nu/(1-nu)
 ck('projected L modulus',logL/(1-logL),Q(161,10**6));ck('projected log derivative',derL,Q(32,10**5))
 R0=Q(61,100)+Q(9,8)*outer*eu(Q(3,2)*outer);ck('selected logarithm disk',R0,Q(610012,10**6))
 merr=(Q(161,10**6)+Q(3,4)*outer)/(1-Q(3,4)*outer)
 ck('uniform K real half-plane',Q(3,5),1/(1+R0)-merr/(1-R0))
 ck('uniform K modulus',(1+merr)/(1-R0),Q(13,5))
 ck('uniform K logarithmic derivative',Q(32,10**5)+(outer/2)/(1-Q(3,4)*outer)+eu(Q(3,2)*outer)/(1-R0),3)
 ck('cube-root J lower',Q(4,5)**3,Q(3,5));ck('cube-root J upper',Q(13,5),Q(7,5)**3)
 ck('stationary H modulus squared',Q(7,3),4);ck('stationary B quotient bound',3*Q(5,4),4)
 ck('outer zeta bound',Q(1,4)*Q(7,5)**2,Q(49,100),'=')
 ratio=actual/(outer-actual);ck('outer-to-actual Cauchy ratio',ratio,Q(1,99999),'=')
 ck('uniform actual A lower',Q(98,100),Q(989,1000)-2*ratio)
 ck('uniform actual B upper',Q(500,989)+4*ratio,Q(52,100))
 ck('uniform actual zeta',Q(1,9)*Q(118,117)**2+Q(49,100)*ratio,Q(1,8))
 ck('uniform actual zeta derivative',Q(3,100)+Q(490,9)*ratio,Q(1,4))
 ck('selected roots physical inclusion',Q(1,2)+Q(501,10000),Q(47,50));ck('physical disk P inclusion',Q(47,50)+Q(501,10000),1)
 delta=actual;R=Q(2);den=Q(199,20000)
 point=R**3/(6*(1-delta*R/2)**2)+(R*R/2)*eu(delta*R)/(den-delta*R*R*eu(delta*R)/2)
 ck('uniform physical phase delta derivative',point,203)
 R=Q(3,2);critical=R**3/(6*(1-delta*R/2)**2)+(R**3/2)*eu(delta*R)/(1-delta*R)
 ck('uniform critical phase derivative',critical,Q(23,10));ck('actual mean and physical allowance',203+Q(23,10),206)
 ck('uniform physical boundary',Q(447,5000),Q(9,100)-Q(52,100000)-206*actual)
 ck('uniform cubic radius',Q(91,200)**3/3+Q(1,8)*Q(91,200),Q(447,5000))
 ck('first-log physical analytic',actual,1);ck('second-log physical analytic',2*actual*eu(2*actual),Q(199,20000))
 ep=R**3/(6*(1-delta*R/2)**2)+(R*R/2)*eu(delta*R)/(Q(18999,20000)-Q(2,5)-delta*R*R*eu(delta*R)/2)
 ck('uniform endpoint point derivative',ep,Q(27,10))
 ck('uniform endpoint W inclusion',2*Q(1,10**5)+5*actual,2*Q(1,10**5)+4*Q(1,10**8))
 V=Q(751,1000);pert=Q(1,10000)*V+actual*((V+actual)*V*V*eu(actual*V)/2+V)
 ck('all-height stationary root disk',pert,Q(1,1000)*(Q(7,12)-Q(1,1000)))
 Vm=Q(16,5);q=actual*Vm/2;denH=Q(21,10)-V-actual*Vm*Vm*eu(actual*Vm)/2
 ck('all-height horizontal denominator',Q(4,3),denH)
 curvature=-1+(2*q+q*q)/(1-q)**2+actual*eu(actual*Vm)/Q(4,3)+eu(2*actual*Vm)/Q(4,3)**2
 ck('all-height whole curvature',curvature,-Q(2,5))
 absA=1/(1-q)**2+actual*eu(actual*Vm)/Q(4,3)+eu(2*actual*Vm)/Q(4,3)**2
 ck('all-height curvature modulus',absA,Q(8,5))
 loc=Q(801,1000);denloc=Q(21,10)-loc-actual*loc*loc*eu(actual*loc)/2
 ck('all-height local denominator',Q(5,4),denloc)
 third=actual/(1-actual*loc/2)**3+actual**2*eu(actual*loc)/Q(5,4)+3*actual*eu(2*actual*loc)/Q(5,4)**2+2*eu(3*actual*loc)/Q(5,4)**3
 ck('all-height third derivative',third,Q(6,5))
 ck('stable U lower allowance',V*V/2+actual*V**3/(6*(1-actual*V/2)),Q(283,1000))
 ck('stable log ratio',Q(752,1000),Q(2,5)*Q(21,10));x=Q(717,1000);ck('stable peak log comparison',Q(5,3),1+x+x*x/2)
 ck('real endpoint U',-Q(9,2)+actual*27/(6*(1-3*actual/2)),-Q(449,100))
 h3=3+actual*9*eu(3*actual)/2;ck('real endpoint logy ratio',1+((h3+Q(1,10000))/Q(21,10))**2,Q(7,4)**2)
 ck('far vertical U',-(9-V*V)/2+actual*Vm**3/(6*(1-actual*Vm/2)),-Q(21,5))
 ck('far vertical logy ratio',1+Vm/Q(21,10),Q(8,3));ck('exp1 exceeds 8/3',Q(8,3),sum(Q(1,math.factorial(k))for k in range(6)))
 ck('uniform right tail derivative',(-3+Q(11,10)*(Q(1,3)+actual))/(1+3*actual/2),-2)
 ck('uniform left tail derivative',2,3-1/(3/(1+3*actual)-Q(1,10000)))
 r=Q(1,10**6);R=1+r;dx=Q(2,10**7)
 ck('bulk stationary perturbation coefficient',(R+actual)*R*R*eu(actual*R)/2+R,Q(8,5))
 ck('bulk root isolation',dx*R+Q(8,5)*actual,r*(Q(74,125)-r))
 t=dx+r+actual*R*R*eu(actual*R)/2;q=actual*R/2;e=actual*R*eu(actual*R)
 Adiff=(2*q+q*q)/(1-q)**2+actual*eu(actual*R)/(1-t)+(2*e+e*e)/(1-t)**2+(2*t+t*t)/(1-t)**2
 ck('bulk curvature difference',Adiff,Q(3,10**6));ck('bulk real curvature',Q(17,100),2-Q(191,100)**2/2-Q(3,10**6))
 ck('bulk local denominator',Q(47,50),1-t-Q(1,20)*eu(actual*Q(1051,1000)))
 loc=Q(1051,1000);den=Q(47,50)
 third=actual/(1-actual*loc/2)**3+actual**2*eu(actual*loc)/den+3*actual*eu(2*actual*loc)/den**2+2*eu(3*actual*loc)/den**3
 ck('bulk third derivative',third,Q(5,2));ck('bulk horizontal local decay',Q(3,50),Q(17,200)-Q(5,2)*Q(1,20)/6)
 ck('bulk compact U perturbation coefficient',Vm**3/(6*(1-actual*Vm/2)),6);ck('bulk compact h perturbation coefficient',Vm**2*eu(actual*Vm)/2,6)
 derr=r+dx+6*actual;phase=6*actual+Vm*r+r*r/2+derr/(Q(1,10)-derr)
 ck('bulk whole phase comparison',phase,Q(1,50000))
 near0=(Q(191,200)+r)**2/2+6*actual+Q(3,2);ck('near log-zero profile',(Q(1,10)+derr)*eu(near0),1)
 ck('bulk nonlocal horizontal margin',Q(1,10000),Q(3,50)*Q(49,1000)**2-Q(2,50000))
 ck('bulk centre margin',Q(1,1000),Q(1,600)-Q(2,50000));ck('bulk far vertical phase',-(9-(Q(191,200)+r)**2)/2+6*actual+2,-2)
 R=Vm;D0=Q(18999,20000);eh=actual*R*R*eu(actual*R)/2
 hor=R**3/(6*(1-actual*R/2)**2)+(R*R/2)*eu(actual*R)/(D0-eh)
 ck('uniform nonlocal horizontal delta derivative',hor,11);ck('uniform normalized nonlocal perturbation',Q(11,5)*Q(1,10**5)+14*actual,Q(3,100000))
 R=Q(3,2);jv=(1+actual*R)*eu(Q(3,2)*actual*R/(1-actual*R/2))
 ck('uniform critical action z derivative',jv,Q(1001,1000));ck('uniform individual critical action derivative',Q(2,5)*jv,Q(401,1000))
 ck('uniform individual critical action perturbation',Q(401,1000)*Q(1,10**5)+Q(46,10)*actual,Q(1,10**5))
 R=Vm;vert=R**3/(6*(1-actual*R/2)**2)+(R*R/2)*eu(actual*R)/(Q(8998,10000)-actual*R*R*eu(actual*R)/2)
 ck('uniform nonlocal vertical delta derivative',vert,12);ck('uniform nonlocal vertical phase',-Q(5,2)+12*actual,-Q(19,10));ck('uniform nonlocal real endpoint',-3+12*actual,-2)
 ck('uniform matched mean lower',Q(2,5),D0**2-Q(1,10**5)**2/4-Q(1,2)-Q(23,10)*actual)
 V,Z,H,D=[var(4,j)for j in range(4)];one={(0,0,0,0):Q(1)}
 P=add(add(mul(Z,V),mul(add(V,scale(D,-1)),H)),scale(one,-1))
 fvnum=add(scale(mul(V,add(Z,H)),-1),add(one,mul(D,H)));assert fvnum==scale(P,-1)
 radial=add(add(mul(D,Z),scale(one,-1)),scale(mul(add(D,scale(V,-1)),add(Z,H)),-1));assert radial==P
 G,E,V=[var(3,j)for j in range(3)];one3={(0,0,0):Q(1)};p=add(one3,mul(E,V));jn=add(p,scale(mul(mul(G,E),V),-1))
 assert jn==add(one3,scale(mul(mul(add(G,scale(one3,-1)),E),V),-1))
 count=0
 for d in range(2,65):
  c=Q(1)
  for k in range(1,d+1):
   c*=Q(d-k+1,d);assert 0<c<=1 and c/Q(math.factorial(k))==Q(math.comb(d,k),d**k);count+=1
 assert Q(3*2*1,3**3)==Q(2,9) # Generic h_third for d3 is nonzero.
 assert Q(2*3,2**2)>1 # Rising coefficients violate the majorant.
 proof=HERE/'UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md'
 out=dict(schema=1,status='PASS_UNIFORM_PHASE_AND_CONTOUR_NUMERIC_SUPPORT',scope='Exact support for the written uniform local, matched-data, inverse and complete-contour estimates; analytic arguments supply the all-parameter scope, not the finite controls.',proof_sha256=sha256(proof.read_bytes()).hexdigest(),checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),rational_guards=ROWS,universal_identity_count=3,falling_coefficient_controls=count,negative_controls=2)
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='PASS',rational_guards=len(ROWS),universal_identities=3,falling_controls=count,negative_controls=2)))
if __name__=='__main__':main()
