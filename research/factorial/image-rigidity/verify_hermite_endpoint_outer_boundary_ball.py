#!/usr/bin/env python3
"""Rigorous ball enclosure ONLY of a Hermite boundary; no CFU or SFC claim."""
if not __debug__:raise SystemExit('Unoptimized Python required')
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as Q
import argparse,json,time
import flint
from flint import arb,acb,ctx
ctx.prec=100
I=acb(0,1);rho=arb('47/50');target=arb('9/100')
def ball(lo,hi):
 mid=(lo+hi)/2;rad=(hi-lo)/2
 return arb(str(mid),str(rad))
def phase_box(alo,ahi,ylo,yhi):
 angle=arb.pi()*ball(alo,ahi);y=ball(ylo,yhi);D=acb(0,y/2)
 u=acb(rho*angle.cos(),rho*angle.sin());v=u-D
 assert (D+u).imag>arb('99/10000')
 return -v*v/2+(D+u).log()+D*D+arb('1/2')-I*arb.pi()/2

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path('/tmp/hermite_endpoint_outer_boundary_ball_verification.json'));args=ap.parse_args()
 start=time.monotonic();leaves=[];stack=[];attempts=0;maxdepth=0
 # theta/pi covers[0,2], y covers[1.8999,2.1001]; all cells closed.
 ymin=Q(18999,10000);ymax=Q(21001,10000)
 for j in range(8):
  yl=ymin+(ymax-ymin)*j/8;yh=ymin+(ymax-ymin)*(j+1)/8
  for k in range(64):stack.append((Q(k,32),Q(k+1,32),yl,yh,0))
 while stack:
  al,ah,yl,yh,depth=stack.pop();attempts+=1;s=phase_box(al,ah,yl,yh)
  if abs(s)>target:
   leaves.append([str(al),str(ah),str(yl),str(yh),str(s.abs_lower())]);maxdepth=max(maxdepth,depth);continue
  assert depth<28,('UNRESOLVED_BOUNDARY_CELL',al,ah,yl,yh,s)
  # Balance angular and y widths in the original dimensional coordinates.
  if (ah-al)*3>yh-yl:
   md=(al+ah)/2;stack.extend([(al,md,yl,yh,depth+1),(md,ah,yl,yh,depth+1)])
  else:
   md=(yl+yh)/2;stack.extend([(al,ah,yl,md,depth+1),(al,ah,md,yh,depth+1)])
 assert leaves
 # Exact pole-domain and target consistency controls.
 assert Q(18999,20000)-Q(47,50)>Q(99,10000)
 control=phase_box(Q(1,2),Q(1,2),Q(21,10),Q(21,10));assert control.real>target and control.imag.contains(0)
 assert not abs(control)>arb('1/10') # A false .1 global lower bound must be rejected.
 assert not (arb('9/100')-arb('52/100000')-arb('110/100000000'))>arb('9/100')
 result={'status':'PASS_RIGOROUS_BOUNDARY_ONLY','scope':'For every y in[1.8999,2.1001] and theta in[0,2pi], with z=iy,v=-z/2+.94exp(i theta), |F0(v,z)-Phi0(z)|>.09. No uniform inverse, period remainder or SFC cutoff follows from this verifier alone.',
 'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'arithmetic':'python-flint Arb/Acb outward-rounded complex balls','precision_bits':ctx.prec,'versions':{'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__},'root_boxes':512,'attempted_boxes':attempts,'certified_leaf_boxes':len(leaves),'maximum_bisection_depth':maxdepth,'leaf_stream_sha256':sha256(json.dumps(leaves,separators=(',',':')).encode()).hexdigest(),'leaves':leaves,'negative_controls':2,'seconds':time.monotonic()-start}
 args.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k!='leaves'}))
if __name__=='__main__':main()
