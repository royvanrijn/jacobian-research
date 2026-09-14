#!/usr/bin/env sage-python
"""Retain three disjoint Hensel balls at each discovered local split fibre."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import PolynomialRing, Qp, QQ
from sage.version import version

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-common-quartic-splitting-local-gate-v1'
SOURCE='artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
CHECKER='elkies-k3/scripts/verify_common_quartic_local_splitting.py'
PRIMES=[2,3,5,7,11,13,17,19,23,29,31,37,41,43]

def digest(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def write_new(path,value):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('x') as f:json.dump(value,f,indent=2,sort_keys=True);f.write('\n')

def freeze(out):
 source=read(ROOT/SOURCE)
 initial=read(out/'hensel-preview.json')['rows']
 completed={(r['parent'],r['p']):r for r in read(out/'panel-completion-preview.json')['rows']}
 choices=[]
 for row in initial:
  selected=row if row['status']=='ROOT' else completed[(row['parent'],row['p'])]
  assert selected['status'] in ['ROOT','SPLITS']
  choices.append({k:selected[k] for k in ['parent','p','chart','coordinate']})
 assert len(choices)==56
 write_new(out/'input.json',{'schema':'common-quartic-local-splitting-input-v1',
  'source':SOURCE,'source_sha256':digest(ROOT/SOURCE),'parents':source['parents'],
  'primes':PRIMES,'choices':choices,'sage_version':version,'padic_precision':40,
  'producer_sha256':digest(Path(__file__)),'checker_sha256':digest(ROOT/CHECKER),
  'preflight_bindings':{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file()},
  'cpu_seconds':40,'address_space_bytes':4*1024**3,
  'selection':'The fixed four-parent fourteen-prime local panel. Retained discovery precedes this frozen proof; no global branch point, cover or section search.',
  'boundary':'Local split fibres over each Qp separately, not a global split fibre, low-degree rational divisor, coefficient-system solution or positive rank19 construction.'})
 print('Frozen56 local fibres',flush=True)

def run(out):
 packet=read(out/'input.json')
 assert packet['source_sha256']==digest(ROOT/SOURCE)
 assert packet['producer_sha256']==digest(Path(__file__))
 assert packet['checker_sha256']==digest(ROOT/CHECKER)
 assert packet['parents']==read(ROOT/SOURCE)['parents']
 for rel,h in packet['preflight_bindings'].items():assert digest(out/rel)==h
 start=time.process_time();R=PolynomialRing(QQ,'t');records=[]
 for parent in packet['parents']:
  rows=[]
  for choice in [c for c in packet['choices'] if c['parent']==parent['name']]:
   p=choice['p'];c=choice['coordinate'];reverse=choice['chart']=='v'
   A=R(parent['A'][::-1] if reverse else parent['A'])
   B=R(parent['B'][::-1] if reverse else parent['B'])
   a,b=A(c),B(c)
   assert -4*a**3-27*b**2
   scale=min(int(a.valuation(p))//2,int(b.valuation(p))//3)
   a,b=a/QQ(p)**(2*scale),b/QQ(p)**(3*scale)
   assert a.valuation(p)>=0 and b.valuation(p)>=0
   K=Qp(p,prec=packet['padic_precision']);S=PolynomialRing(K,'z')
   roots=S([b,a,0,1]).roots()
   assert len(roots)==3 and all(n==1 for _,n in roots)
   centres=sorted(QQ(z.lift()) for z,_ in roots)
   balls=[]
   for z in centres:
    f=z**3+a*z+b;d=3*z*z+a
    assert d
    fv=None if not f else int(f.valuation(p));dv=int(d.valuation(p))
    assert fv is None or fv>2*dv
    balls.append({'centre':str(z),'f_valuation':fv,'derivative_valuation':dv})
   for i in range(3):
    for j in range(i):
     distance=int((centres[i]-centres[j]).valuation(p))
     for ball in [balls[i],balls[j]]:
      assert ball['f_valuation'] is None or distance<ball['f_valuation']-ball['derivative_valuation']
   rows.append({**choice,'root_scale':scale,'balls':balls})
  name=parent['name']+'-checkpoint.json'
  write_new(out/name,{'parent':parent['name'],'rows':rows})
  records.append({'path':name,'sha256':digest(out/name)})
  print(parent['name'],len(rows),'local split witnesses',flush=True)
 write_new(out/'result.json',{'schema':'common-quartic-local-splitting-result-v1','status':'PASS',
  'input_sha256':digest(out/'input.json'),'records':records,'witness_count':56,
  'local_splitting_at_all_tested_primes':True,'global_splitting_divisor':'UNKNOWN',
  'positive_mw17_target_complete':False,'cpu_seconds':round(time.process_time()-start,6)})

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('mode',choices=['freeze','run']);parser.add_argument('--output',type=Path,default=DEFAULT)
 args=parser.parse_args()
 resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
 (freeze if args.mode=='freeze' else run)(args.output)
