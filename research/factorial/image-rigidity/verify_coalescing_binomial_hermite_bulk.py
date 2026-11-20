#!/usr/bin/env python3
"""Exact support for LOCAL Hermite bulk lemmas; no global SFC certificate."""
from fractions import Fraction as F
from pathlib import Path
from math import comb,isqrt
import argparse,hashlib,json
if not __debug__:raise SystemExit('Refusing -O/-OO: exact checks required.')
HERE=Path(__file__).resolve().parent
RECORDS=[]
def st(x):
 x=F(x);return str(x.numerator)if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def ck(name,l,r,relation='<'):
 l,r=F(l),F(r);assert (l<r if relation=='<' else l==r),(name,st(l),relation,st(r));RECORDS.append({'name':name,'left':st(l),'relation':relation,'right':st(r)})
def eu(x):
 x=F(x)
 if x<0:
  term=total=F(1)
  for k in range(1,257):term*=(-x)/k;total+=term
  return 1/total
 term=total=F(1)
 for k in range(1,257):term*=x/k;total+=term
 assert x/258<1
 return total+term*x/257/(1-x/258)
def sqrtbound(x,upper):
 x=F(x);m=10**12;k=isqrt(x.numerator*m*m//x.denominator)
 assert F(k*k,m*m)<=x<F((k+1)**2,m*m)
 return F(k+int(upper),m)
def padd(a,b):return[(a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0)for i in range(max(len(a),len(b)))]
def pscale(a,c):return[c*v for v in a]
def pmul(a,b):
 out=[F(0)]*(len(a)+len(b)-1)
 for i,u in enumerate(a):
  for j,v in enumerate(b):out[i+j]+=u*v
 return out
def deriv(a):return[(i+1)*v for i,v in enumerate(a[1:])]
def equal(a,b):return not any(padd(a,pscale(b,-1)))

MANDATORY_SOURCES = {'airy_proof': ('packet', 'COALESCING_BINOMIAL_SADDLE_BRIDGE.md'), 'bulk_proof': ('packet', 'COALESCING_BINOMIAL_HERMITE_BULK.md'), 'airy_checker': ('packet', 'verify_coalescing_binomial_saddle_bridge.py'), 'bulk_checker': ('packet', 'verify_coalescing_binomial_hermite_bulk.py'), 'large_lower': ('repository', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'), 'lower_four': ('repository', 'research/factorial/image-rigidity/LOWER_FOUR_QUADRATIC_SECTOR_SFC.tex'), 'lower_five': ('repository', 'research/factorial/image-rigidity/LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex'), 'gamma_obstruction': ('repository', 'research/factorial/image-rigidity/GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex'), 'etna_context': ('repository', 'research/factorial/image-rigidity/sources/dunster_gil_ruiz_segura_etna65_2026_140_155.pdf'), 'wolfs_context': ('repository', 'research/factorial/image-rigidity/sources/wolfs-2024-multiple-orthogonal-2401.08312v1.pdf')}

def verify_sources(repo_arg):
    repo=repo_arg.resolve() if repo_arg else next((p for p in (HERE,*HERE.parents,Path.cwd(),*Path.cwd().parents) if (p/'research/MATH_STATUS.json').is_file()),None)
    if repo is None:raise SystemExit('Use --repo-root outside the checkout.')
    ledger_path=HERE/'coalescing_binomial_local_source_ledger.json'
    ledger=json.loads(ledger_path.read_text())
    assert ledger['schema_version']==1 and set(ledger['files'])==set(MANDATORY_SOURCES),'Missing or extra mandatory source record'
    for key,(scope,name) in MANDATORY_SOURCES.items():
        rec=ledger['files'][key]
        assert rec['scope']==scope and rec['path']==name
        p=(HERE if scope=='packet' else repo)/name
        assert p.is_file() and hashlib.sha256(p.read_bytes()).hexdigest()==rec['sha256'],'Stale source: '+key
    return ledger_path

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path);ap.add_argument('--output',type=Path,default=HERE/'coalescing_binomial_hermite_bulk_verification.json');args=ap.parse_args();ledger_path=verify_sources(args.repo_root)
 ck('quadratic discriminant real-part gap',F(13,10)**2,F(7,4))
 ck('quadratic-root ellipse',F(3,11)**2+F(45,61)**2,1)
 R=F(123,100);delta=F(1,100)
 ck('small-circle stationary remainder',R+(R+delta)*R*R/2*F(10000,9877),F(11,5))
 ck('small-circle Rouche margin',F(11,5),F(381,100))
 ck('big-circle stationary remainder',delta*(3+F(301,100)*F(9,2)*F(100,97)),F(1,5))
 ck('big-circle Rouche margin',F(1,5),F(16,5))
 ck('real-part branch separation',F(6,5),F(13,10)-F(6,100))
 re=1+F(57,100)**2-F(1089,1300)**2-F(123,10000)
 ck('curvature numerator real part',F(61,100),re)
 t=F(123,10000);M=F(6313,2500)
 ck('curvature real part',F(1,2),F(61,100)-M*(2*t+t*t)/(1-t)**2)
 ck('curvature norm',M/(1-t)**2,3)
 H=F(625,617)
 ck('local logarithm denominator',F(3,4),F(9877,12300)-F(125,2468))
 third=F(2,100)/(1-F(32,2500))**3+H/F(10000)/F(3,4)+3*(H/100)*H/F(3,4)**2+2*H**3/F(3,4)**3
 ck('local third derivative',third,6)
 ck('Taylor decay algebra',F(1,4)-F(1,20),F(1,5),'=')
 ck('local Gaussian relative error',(25+100*eu(-F(9999,2000)))/sqrtbound(2*9999,False),F(1,5))
 ck('inner z norm',F(1,20)**2+F(29,20)**2,F(3,2)**2)
 ck('radial parameter shift',F(203,9998),F(1,40))
 ck('inner-to-outer real margin',F(1,20)+F(1,40),F(1,10))
 ck('inner-to-outer imaginary margin',F(29,20)+F(1,40),F(3,2))
 ck('hard-domain order lower bound',9999,F(10000*2*2,2*2)+F(1,2))
 ck('action drop order margin',F(5,6),F(9999,10000))
 ck('conditional packet prefactor ratio',F(3,2)*sqrtbound(6,True),4)
 ck('conditional action separation',16*eu(-200),1)
 # Coefficient-level identities, with rational epsilon on exact square-a fixtures.
 identity_count=0
 for d in [2,3,7]:
  eps=F(1,100*d);de=d*eps
  A=[F(comb(d,k))*eps**k for k in range(d+1)]
  h=[F(0)]+[A[k]/de for k in range(1,d+1)]
  assert equal(A,padd([F(1)],pscale(h,de)));identity_count+=1
  vh=pmul([F(0),F(1)],h)
  left=padd(vh,pscale(A,-1));right=padd(pmul([-de,F(1)],h),[F(-1)])
  assert equal(left,right);identity_count+=1
  assert equal(deriv(h),[F(comb(d-1,k))*eps**k for k in range(d)]);identity_count+=1
  B=[F(1),eps];B2=pmul(B,B)
  assert equal(pmul(deriv(deriv(h)),B2),pscale(A,(d-1)*eps));identity_count+=1
  assert equal(pmul(pmul(deriv(h),deriv(h)),B2),pmul(A,A));identity_count+=1
  # -Fvv at a stationary point after clearing(1+eps*v)^2:
  # 1 -(d-1)eps*v +v², exactly1+v²+(eps-delta)v.
  assert equal([F(1),-(d-1)*eps,F(1)],[F(1),eps-de,F(1)]);identity_count+=1
  if d==3:
   P=padd(pmul([-de,F(1)],h),[F(-1)])
   assert len(P)-1==d+1>2 # two LOCAL roots cannot mean only two GLOBAL roots.
 # Genuine corruptions of the exact normalization are distinguishable.
 assert not equal(padd(vh,pscale(A,-1)),padd(vh,[F(-1)]))
 assert not equal([F(1),-de,F(1)],[F(1),eps-de,F(1)])
 proof=HERE/'COALESCING_BINOMIAL_HERMITE_BULK.md';assert proof.is_file()
 def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
 out={'status':'PASS','scope':'Local uniform stationary branches/curvature/Gaussian and exact conditional action-drop support only. No global root localization, full contour, or SFC theorem.','counts':{'rational_guards':len(RECORDS),'coefficient_identities':identity_count,'negative_controls':3},'files':[{'path':p.name,'sha256':sha(p)}for p in[proof,Path(__file__).resolve(),ledger_path]],'comparisons':RECORDS}
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','counts':out['counts']}))
if __name__=='__main__':main()
