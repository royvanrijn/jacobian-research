#!/usr/bin/env python3
"""Source-guarded outward-ball endpoint roots and local connector margins."""
from __future__ import annotations
if not __debug__:raise SystemExit('Refusing -O/-OO: exact certification must execute.')
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import argparse,importlib.util,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
ROWS=[]
def sha(p):return sha256(p.read_bytes()).hexdigest()
def ck(name,left,right,op='<'):
 left,right=Q(left),Q(right)
 assert (left<right if op=='<' else left==right),(name,left,op,right)
 ROWS.append(dict(name=name,left=str(left),relation=op,right=str(right)))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'gap_two_horizontal_endpoint_location_verification.json');args=ap.parse_args();repo=args.repo_root.resolve()
 ledgerpath=HERE/'gap_two_horizontal_endpoint_source_ledger.json';ledger=json.loads(ledgerpath.read_text());assert ledger.get('schema')==1
 local={'GAP_TWO_HORIZONTAL_ENDPOINT_CUBIC_LOCATION.md',Path(__file__).name,'gap_two_horizontal_endpoint_inputs.json','gap_two_horizontal_endpoint_primary_audit.md','verify_gap_two_horizontal_endpoint_cubic_balls.py','gap_two_horizontal_endpoint_cubic_balls_verification.json'}
 canonical={'QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md','HERMITE_CRITICAL_CUBIC_INVERSE_DISK.md','GAP_TWO_FULL_ENDPOINT_UNIFORM_CUBIC_INVERSE.md','verify_gap_two_full_endpoint_uniform_inverse.py','gap_two_full_endpoint_inverse_source_ledger.json','sources/dlmf-1-10-rouche.html'}
 mandatory={('local',s)for s in local}|{('repo','research/factorial/image-rigidity/'+s)for s in canonical}
 assert len(ledger['files'])==len(mandatory) and {(v['root'],v['path'])for v in ledger['files']}==mandatory
 for row in ledger['files']:
  assert row['root']in{'local','repo'};p=(HERE if row['root']=='local'else repo)/row['path'];assert p.is_file()and sha(p)==row['sha256'],row['path']
 inputs=json.loads((HERE/'gap_two_horizontal_endpoint_inputs.json').read_text());assert inputs['schema']==1
 assert [r['declaration']for r in inputs['coverage']]==['Theorem HorizontalEndpointCubicLocation','Corollary EndpointRadialCubicDecay']
 assert inputs['root_radius']=='1/1000'and inputs['root_modulus_bounds']==['379/1000','403/1000']and inputs['root_real_lower']=='1/100'and inputs['root_imaginary_abs_lower']=='1/5'
 assert inputs['angle_guard']=='1/1000'and inputs['radial_bracket_guard']=='1/10'and inputs['common_connector_radius']=='41/100'and inputs['radial_derivative_bound']=='-99/100000'
 source=HERE/'verify_gap_two_horizontal_endpoint_cubic_balls.py';baseline=json.loads((HERE/'gap_two_horizontal_endpoint_cubic_balls_verification.json').read_text());expected=inputs['expected_execution']
 for name,value in expected.items():assert baseline[name]==value,(name,baseline[name],value)
 assert sha(source)==baseline['source_sha256'];leaves=baseline['leaves'];assert len(leaves)==7262
 assert sha256(json.dumps(leaves,separators=(',',':'),sort_keys=True).encode()).hexdigest()==baseline['stream_sha256']
 spec=importlib.util.spec_from_file_location('executed_cubic_endpoint',source);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 domains={'left':(Q(18999,20000),Q(1)),'right':(Q(0),Q(801,2500))}
 independent_shared=0
 for side,(start,end)in domains.items():
  cells=sorted((r for r in leaves if r['side']==side),key=lambda r:Q(r['lo']))
  assert cells and Q(cells[0]['lo'])==start and Q(cells[-1]['hi'])==end
  for a,b in zip(cells,cells[1:]):assert Q(a['hi'])==Q(b['lo'])
  for row in cells:
   lo,hi=Q(row['lo']),Q(row['hi']);assert lo<hi and Q(row['radius'])==Q(1,1000)
   re,im=map(Q,row['centre']);z,W,u,D=mod.data(side,lo,hi);t=mod.ac(re,im);r=mod.av(Q(1,1000))
   # Independent expansion of the full actual-family Rouché inequality.
   lhs=abs(t**3/3-z*t-W)+mod.av(mod.EW)+abs(t)*mod.av(mod.EZ)+abs(t)*r*r+r**3/3
   rhs=(abs(t*t-z)-mod.av(mod.EZ))*r
   assert lhs<rhs
   assert abs(t)-r>mod.av(Q(379,1000))and abs(t)+r<mod.av(Q(403,1000))
   assert t.real-r>mod.av(Q(1,100))and t.imag+r<-mod.av(Q(1,5))
   angle=(-t.imag-r)**2-(t.real+r)**2/3
   assert angle>mod.av(Q(1,1000))and 3*angle+z-mod.av(mod.EZ)>mod.av(Q(1,10))
  for j in range(len(cells)+1):
   point=start if j==0 else Q(cells[j-1]['hi']);re,im=mod.centre(side,point,point)
   assert mod.certify(side,point,point,re,im,mod.NODE,False)[0]
   for k in[j-1,j]:
    if 0<=k<len(cells):assert abs(mod.ac(re,im)-mod.ac(*map(Q,cells[k]['centre'])))+mod.av(mod.NODE)<mod.av(mod.R)
   independent_shared+=1
 assert independent_shared==7264
 with tempfile.TemporaryDirectory(prefix='horizontal_endpoint_replay_')as td:
  out=Path(td)/'executed.json';run=subprocess.run([sys.executable,str(source),'--output',str(out)],capture_output=True,text=True)
  assert run.returncode==0,(run.stdout,run.stderr);fresh=json.loads(out.read_text())
 assert fresh==baseline
 # All scalar margins are exact rationals, not floating estimates.
 ew=Q(501,25000000);ez=Q(139199,7999920000);r=Q(41,100);outer=Q(17,40)
 ck('Rouché perturbation W',mod.EW,ew,'=');ck('Rouché perturbation zeta',mod.EZ,ez,'=')
 ck('radial derivative strict loss',Q(99,100000),Q(1,1000)-ez*Q(403,1000))
 ck('endpoint radius below connector',Q(403,1000),r);ck('connector below inverse',r,Q(91,200))
 ck('actual zeta real lower',Q(1001,10000)*Q(101,100)+ez,Q(102,1000))
 ck('fixed Airy exit strict phase',-r**3/3+Q(102,1000)*r/2+ez*r,-Q(1,500))
 ck('larger inverse Jacobian',Q(47,50)/(Q(9,20)-outer),38)
 ck('annular G',outer*38/((outer-r)*(outer*outer-Q(1,8))),20000)
 ck('annular G derivative',outer*38/((outer-r)**2*(outer*outer-Q(1,8))),1300000)
 ck('circle convexity on endpoint-to-ray arc',Q(79,1000),8*Q(1,100)-ez*r*r/Q(1,5)**3)
 # The larger circle-arc range used by other compatible canonical cuts.
 ck('circle large-p branch',ez*r*r,8*Q(1,5)*Q(1,50)**3)
 ck('circle small-p ordinate',(Q(7,20))**2,r*r-Q(1,5)**2)
 ck('circle small-p branch',ez*r*r,8*Q(1,100)*Q(7,20)**3)
 # Universal cubic radial/circle derivative identities via exact coefficient
 # comparison. For t=x+iy, Re(t^3)/3=x^3/3-x*y^2.
 # At radius r, y²=r²-x², hence Re phi=4x³/3-(r²+alpha)x+beta*y.
 # Radial derivative is -x[3lambda²(y²-x²/3)+alpha]+beta*y.
 # Independent exact sparse polynomial multiplication checks the identity.
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
 def var(i):return{tuple(int(i==j)for j in range(5)):Q(1)}
 x,y,lam,alpha,beta=[var(i)for i in range(5)]
 lhs=add(add(mul(mul(lam,lam),add(mul(x,mul(x,x)),scale(mul(x,mul(y,y)),-3))),scale(mul(alpha,x),-1)),mul(beta,y))
 rhs=add(scale(mul(x,add(scale(mul(mul(lam,lam),add(mul(y,y),scale(mul(x,x),-Q(1,3)))),3),alpha)),-1),mul(beta,y))
 assert lhs==rhs
 # Explicit negative controls beyond the original computation's three:
 # a removed cell cannot cover the full interval, and the unsafe conjugate
 # branch fails the actual plus-endpoint cubic at the known critical base.
 left=[q for q in leaves if q['side']=='left'];assert Q(left[1]['lo'])>domains['left'][0]
 re,im=mod.centre('left',Q(1),Q(1));assert not mod.certify('left',Q(1),Q(1),re,-im)[0]
 result=dict(schema=1,status='PASS_HORIZONTAL_ENDPOINT_LOCATION_SUPPORT',scope=inputs['scope'],coverage=inputs['coverage'],independent_ball_leaf_checks=len(leaves),independent_shared_endpoint_discs=independent_shared,complete_source_replay={k:v for k,v in fresh.items()if k!='leaves'},rational_guards=ROWS,universal_radial_identity_count=1,negative_controls=5,pins=ledger['files'],source_ledger_sha256=sha(ledgerpath))
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(status='PASS',ball_cells=len(leaves),shared_endpoint_discs=independent_shared,rational_guards=len(ROWS),universal_identities=1,negative_controls=5)))
if __name__=='__main__':main()
