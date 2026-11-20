#!/usr/bin/env python3
"""Outward-ball cubic endpoint certificates; no whole-period/SFC assertion."""
if not __debug__:raise SystemExit('Unoptimized Python required')
from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
import argparse,json,cmath
import flint
from flint import arb,acb,ctx
ctx.prec=128
I=acb(0,1)
def av(q):return arb(str(q))
def ball(lo,hi):return arb(str((lo+hi)/2),str((hi-lo)/2))
def ac(re,im):return acb(av(re),av(im))
R=Q(1,1000);NODE=Q(1,10**9)
EW=Q(2,10**5)+Q(4,10**8)
EZ=Q(5,4)*Q(1,10**5)+Q(49,100)/99999

def data(side,lo,hi):
 p=ball(lo,hi)
 if side=='left':D=p;u=acb(av(Q(2,5)))
 else:
  D=(1+p*p).sqrt();X=(av(Q(4,25))-p*p).sqrt();u=acb(X,p)
 e=D-1;A=arb(1);c=Q(1)
 for k in range(1,33):
  c*=Q(3-2*k,2*k);A+=av(c*Q(3,2*k+3))*(e/2)**k
 # Absolute binomial coefficient and3/(2k+3) are at most1.
 q=Q(51,2000);tail=q**33/(1-q)
 A+=arb(0,str(tail))
 assert A>av(Q(99,100)) and A<av(Q(101,100))
 zeta=2*e*A.root(3)**2
 v=u-I*D
 W=-v*v/2+(u+I*D).log()-D*D+av(Q(1,2))-I*arb.pi()/2
 return zeta,W,u,D

def cmid(a):return complex(float(a.real.mid()),float(a.imag.mid()))
def centre(side,lo,hi):
 z,W,_,_=data(side,(lo+hi)/2,(lo+hi)/2)
 z=cmid(acb(z));W=cmid(W);t=.04-.4j
 for j in range(40):
  change=(t*t*t/3-z*t-W)/(t*t-z);t-=change
  if abs(change)<1e-16:break
 # These machine values propose only a rational center. Every conclusion
 # below follows from outward-ball inequalities using the exact decimals.
 return Q(str(t.real)),Q(str(t.imag))

def certify(side,lo,hi,re,im,radius=R,actual=True):
 z,W,u,D=data(side,lo,hi);t=ac(re,im);r=av(radius)
 residual=abs(t**3/3-z*t-W)
 if actual:residual+=av(EW)+abs(t)*av(EZ)
 linear=abs(t*t-z)
 if actual:linear-=av(EZ)
 lhs=residual+abs(t)*r*r+r**3/3;rhs=linear*r
 return lhs<rhs,(lhs,rhs,z,W,t)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path(__file__).with_name('gap_two_horizontal_endpoint_cubic_balls_verification.json'));args=ap.parse_args()
 leaves=[];attempts=0
 domains={'left':(Q(18999,20000),Q(1)),'right':(Q(0),Q(801,2500))}
 for side,(start,end) in domains.items():
  stack=[(start,end,0)];sideleaves=[]
  while stack:
   lo,hi,depth=stack.pop();attempts+=1;re,im=centre(side,lo,hi)
   ok,info=certify(side,lo,hi,re,im)
   if ok:
    t=info[-1]
    assert abs(t)-av(R)>av(Q(379,1000))
    assert abs(t)+av(R)<av(Q(403,1000))
    assert t.real-av(R)>av(Q(1,100))
    assert t.imag+av(R)<-av(Q(1,5))
    angle=(-t.imag-av(R))**2-(t.real+av(R))**2/3
    assert angle>av(Q(1,1000))
    assert 3*angle+info[2]-av(EZ)>av(Q(1,10))
    sideleaves.append(dict(side=side,lo=str(lo),hi=str(hi),centre=[str(re),str(im)],radius=str(R),left=str(info[0]),right=str(info[1]),depth=depth))
   else:
    assert depth<24,('unresolved cell',side,lo,hi,info)
    mid=(lo+hi)/2;stack.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
  sideleaves.sort(key=lambda q:Q(q['lo']));assert Q(sideleaves[0]['lo'])==start and Q(sideleaves[-1]['hi'])==end
  for a,b in zip(sideleaves,sideleaves[1:]):assert Q(a['hi'])==Q(b['lo'])
  # Point-root subdiscs certify the SAME continued cubic root across every
  # adjacent parameter cell, rather than inferring it from overlapping discs.
  for j in range(len(sideleaves)+1):
   point=start if j==0 else Q(sideleaves[j-1]['hi']);re,im=centre(side,point,point)
   ok,_=certify(side,point,point,re,im,NODE,False);assert ok
   for k in [j-1,j]:
    if 0<=k<len(sideleaves):
     leaf=sideleaves[k];other=ac(*map(Q,leaf['centre']))
     assert abs(ac(re,im)-other)+av(NODE)<av(R)
  leaves+=sideleaves
 # Known critical base normalization t=-i*u*(3I(u))^(1/3).
 z,W,u,D=data('left',Q(1),Q(1));base=-I*u*(3*W/(I*u**3)).root(3)
 for side,point in [('left',Q(1)),('right',Q(0))]:
  re,im=centre(side,point,point);assert abs(base-ac(re,im))<av(NODE)
 # Exact parameter and perturbation guards.
 assert Q(801,2500)**2>Q(21001,20000)**2-1
 assert Q(801,2500)**2<Q(4,25)
 assert (av(1)+av(Q(801,2500))**2).sqrt()<av(Q(2101,2000))
 assert Q(51,2000)>Q(101,2000)/2
 assert (Q(2,5)+Q(1051,1000))/2+1/(2*(Q(18999,20000)-Q(2,5)))<2
 assert Q(3,2)**3/(6*(1-Q(1,10**8)*Q(3,2)/2)**2)+Q(3,2)**2/(4*(Q(18999,20000)-Q(2,5)-Q(1,10**8)))<2
 assert EW==Q(501,25000000) and EZ<Q(9,500000)
 assert Q(1,1000)-EZ*Q(403,1000)>Q(99,100000)
 assert Q(403,1000)<Q(41,100)<Q(91,200)
 # Strong controls: coarse whole-domain ball is NOT a certificate, wrong
 # global radius .38 is rejected, and the wrong conjugate branch fails the actual cubic.
 re,im=centre('right',Q(0),Q(801,2500));assert not certify('right',Q(0),Q(801,2500),re,im)[0]
 assert abs(base)>av(Q(38,100))
 bre,bim=centre('left',Q(1),Q(1));assert not certify('left',Q(1),Q(1),bre,-bim)[0]
 assert EW>0 and EZ>0
 stream=sha256(json.dumps(leaves,separators=(',',':'),sort_keys=True).encode()).hexdigest()
 result=dict(status='PASS_CUBIC_ENDPOINT_ROOTS_ONLY',scope='All Hermite horizontal endpoints and their actual |Re z|<=1e-5,|delta|<=1e-8 perturbations: selected cubic roots satisfy .379<|T|<.403,Re T>.01 and opposite Im signs of modulus>.2. This does not certify complete-period connections or an SFC threshold.',source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),arithmetic='outward-rounded Arb/Acb balls; exact rational subdivision and centers',precision_bits=ctx.prec,versions=dict(python_flint=flint.__version__,FLINT=flint.__FLINT_VERSION__),root_disc_radius=str(R),angle_guard='1/1000',radial_bracket_guard='1/10',actual_W_perturbation=str(EW),actual_zeta_perturbation=str(EZ),attempted_cells=attempts,certified_cells=len(leaves),shared_endpoint_subdiscs=len(leaves)+2,maximum_depth=max(q['depth'] for q in leaves),stream_sha256=stream,leaves=leaves,negative_controls=3)
 args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in result.items()if k!='leaves'}))
if __name__=='__main__':main()
