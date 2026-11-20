#!/usr/bin/env python3
"""Exact margins for the written actual original-period nonlocal estimate."""
if not __debug__:raise SystemExit('Unoptimized Python required')
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,hashlib
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--output',type=Path,default=HERE/'actual_shifted_hermite_nonlocal_verification.json')
args=ap.parse_args();repo=args.repo_root.resolve()
lp=HERE/'actual_shifted_hermite_nonlocal_source_ledger.json'
ledger=json.loads(lp.read_text())
mandatory={('local','ACTUAL_SHIFTED_HERMITE_NONLOCAL_BOUND.md'),
 ('local',Path(__file__).name),('local','actual_shifted_hermite_nonlocal_inputs.json'),
 ('repo','research/factorial/image-rigidity/QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md'),
 ('repo','research/factorial/image-rigidity/EXPLICIT_AIRY_CONTOUR_ENVELOPE.md'),
 ('repo','research/factorial/image-rigidity/sources/dlmf-9-2-endpoint.html'),
 ('repo','research/factorial/image-rigidity/sources/dlmf-9-7-endpoint.html'),
 ('repo','research/factorial/image-rigidity/sources/dlmf-1-10-rouche.html')}
assert ledger.get('schema')==1 and len(ledger['files'])==len(mandatory)
assert {(v['root'],v['path'])for v in ledger['files']}==mandatory
for v in ledger['files']:
 assert v['root']in{'local','repo'}
 path=(HERE if v['root']=='local'else repo)/v['path']
 assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==v['sha256'],v['path']
inp=json.loads((HERE/'actual_shifted_hermite_nonlocal_inputs.json').read_text())
assert inp['schema']==1
assert inp['domain']=={'delta_upper':'1/100000000','real_z_upper':'1/100000','D_lower':'18999/20000','D_upper':'21001/20000','m_lower':1,'actual_shape':'a=4/delta^2 is a positive integer'}
assert inp['conclusion']=={'coefficient':60000,'m_power':'1/2','exponential_gap':'33/50000','normalization':'m^(1/3)*exp(-m*Phi)/(2*pi*i)','weak_coefficient':10**15,'weak_exponential_gap':'1/100000'}
assert inp['expected_rational_guards']==32 and inp['negative_controls']==2
ROWS=[]
def ck(name,l,r):
 l,r=Q(l),Q(r);assert l<r,(name,l,r);ROWS.append({'name':name,'left':str(l),'right':str(r)})
dx=Q(1,10**5);dd=Q(1,10**8);D0=Q(18999,20000);D1=Q(21001,20000);V=Q(16,5)
ck('cut above .057',Q(57,1000),Q(4,25)-(D1*D1-1))
ck('bulk q above .0624',Q(39,625),Q(4,25)+D0*D0-1)
q=Q(39,625);ck('bulk horizontal loss',Q(9,10000),q*q/(4*(1+q)))
q=Q(57,1000);ck('stable horizontal loss',Q(7,10000),q*q/(4*(1+q)))
ck('stable shift s below .321',D1*D1-1,Q(321,1000)**2)
ck('D+s below1.372',D1+Q(321,1000),Q(343,250))
ck('horizontal denominator upper',(Q(300001,100000))**2+Q(343,250)**2,Q(33,10)**2)
ck('horizontal z derivative',Q(33,20)+1/(2*D0),Q(11,5))
der=V**3/(6*(1-dd*V/2)**2)+V*V/(4*(D0-dd*V*V/4))
ck('horizontal delta phase derivative',der,9)
ck('horizontal actual perturbation',Q(11,5)*dx+11*dd,Q(3,100000))
D2=max(1-D0*D0,D1*D1-1)+dx*D1+dx*dx/4
ck('critical quadratic centre',D2,Q(103,1000))
vp=Q(3,2);criterr=dd*(vp+vp**3/4)+dd*dd*vp*vp/4
ck('critical numerator perturbation',criterr,Q(3,10**8))
ck('critical small circle',D2+criterr,Q(2,5)**2)
zder=Q(2,5)*(1+dd*vp/2)*vp*vp/(1-dd*vp/2)**4
ck('critical z action derivative',zder,Q(901,1000))
pder=vp**3/(6*(1-dd*vp/2)**2)+vp**3/(4*(1-dd*vp/2)**2)
ck('critical delta phase derivative',pder,Q(141,100))
ck('critical delta action derivative',Q(141,100)+2,Q(7,2))
ck('critical action perturbation',Q(901,1000)*dx+Q(7,2)*dd,Q(1,10**5))
ck('kappa squared',Q(4,9)*Q(1,8)**3,Q(3,100)**2)
ck('vertical norm',9+1,V*V)
ck('vertical logarithm norm',(3+dx)**2+(2*D1)**2,14)
E3=sum(Q(3)**j/Q(__import__('math').factorial(j))for j in range(5));ck('exp3 lower polynomial',14,E3)
derv=V**3/(6*(1-dd*V/2)**2)+V*V/(4*(Q(8998,10000)-dd*V*V/4))
ck('vertical delta derivative',derv,9)
ck('vertical actual phase',-Q(5,2)+9*dd,-Q(19,10))
ck('real endpoint actual phase',-3+9*dd,-2)
ck('right log derivative numerator difference',0,1-dx-Q(3,4)*dd)
ck('right tail phase derivative',-3/(1+Q(3,2)*dd)+Q(1,2),-2)
ck('left tail phase derivative',2,3-1/(2*D0))
ck('mean phase lower',Q(2,5),D0*D0-dx*dx/4-Q(1,2)-2*dd)
ck('joint lower constant',Q(1,36000),Q(3,16))
ck('far gap dominates horizontal',Q(33,50000),Q(23,10)-Q(3,100)-Q(1,10**5))
ck('weak target coefficient',60000,10**15)
ck('weak target exponential',Q(1,100000),Q(33,50000))
# The preceding endpoint convention is essential: t endpoints differ from +/-3.
assert (Q(-3)+dx/2,Q(3)+dx/2)!=(Q(-3),Q(3))
# A modulus-only lower bound exp(-.03m) cannot cover near-cut horizontal.
assert Q(3,100)>Q(7,10000)
o={'status':'PASS','scope':'Exact rational support for the written shifted-Hermite ORIGINAL-period nonlocal bound only; no numerical SFC or complete-period conclusion.','rational_guards':ROWS,'negative_controls':2,
'proof_sha256':hashlib.sha256((HERE/'ACTUAL_SHIFTED_HERMITE_NONLOCAL_BOUND.md').read_bytes()).hexdigest()}
o.update(source_ledger_sha256=hashlib.sha256(lp.read_bytes()).hexdigest(),pins=ledger['files'],input_sha256=hashlib.sha256((HERE/'actual_shifted_hermite_nonlocal_inputs.json').read_bytes()).hexdigest())
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','guards':len(ROWS),'negative_controls':2}))

