#!/usr/bin/env python3
"""Independent exact-rational replay of local complete-splitting witnesses."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-common-quartic-splitting-local-gate-v1'
SOURCE='artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
PRIMES=[2,3,5,7,11,13,17,19,23,29,31,37,41,43]

def require(test,message):
 if not test:raise ValueError(message)
def digest(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def valuation(q,p):
 q=F(q)
 if not q:return None
 a,b=q.numerator,q.denominator;v=0
 while a%p==0:a//=p;v+=1
 while b%p==0:b//=p;v-=1
 return v
def integral(q,p):
 v=valuation(q,p);return v is None or v>=0
def ev(coefficients,t):
 answer=F(0)
 for c in reversed(coefficients):answer=answer*t+F(c)
 return answer

def check_balls(a,b,p,balls):
 require(integral(a,p) and integral(b,p),'monic integral cubic required')
 require(-4*a**3-27*b*b!=0,'singular fibre')
 require(len(balls)==3,'three root balls required')
 checked=[]
 for ball in balls:
  z=F(ball['centre']);require(integral(z,p),'integral root centre required')
  f=z*z*z+a*z+b;d=3*z*z+a
  fv,dv=valuation(f,p),valuation(d,p)
  require(dv is not None,'zero derivative')
  require((fv,dv)==(ball['f_valuation'],ball['derivative_valuation']),'valuation attachment')
  require(fv is None or fv>2*dv,'Hensel inequality')
  checked.append((z,None if fv is None else fv-dv))
 for i in range(3):
  for j in range(i):
   distance=valuation(checked[i][0]-checked[j][0],p)
   require(distance is not None,'duplicate root centre')
   require(all(radius is None or distance<radius for _,radius in [checked[i],checked[j]]),'overlapping root balls')

def check_witness(parent,row):
 require(row['parent']==parent['name'],'parent attachment')
 p=row['p'];c=F(row['coordinate'])
 require(row['chart'] in ['t','v'],'unknown chart')
 if row['chart']=='t':a,b=ev(parent['A'],c),ev(parent['B'],c)
 elif c:
  # Check the reciprocal chart through the original equation and weights.
  a,b=c**8*ev(parent['A'],1/c),c**12*ev(parent['B'],1/c)
 else:a,b=F(parent['A'][8]),F(parent['B'][12])
 scale=row['root_scale'];a,b=a/F(p)**(2*scale),b/F(p)**(3*scale)
 check_balls(a,b,p,row['balls'])

def verify(out):
 start=time.monotonic();packet,result=read(out/'input.json'),read(out/'result.json')
 require(packet['source']==SOURCE and packet['source_sha256']==digest(ROOT/SOURCE),'source binding')
 require(packet['parents']==read(ROOT/SOURCE)['parents'] and len(packet['parents'])==4,'parent scope')
 require(packet['primes']==PRIMES,'prime scope')
 require(packet['producer_sha256']==digest(ROOT/'elkies-k3/scripts/certify_common_quartic_local_splitting.sage'),'producer binding')
 require(packet['checker_sha256']==digest(Path(__file__)),'checker binding')
 for path,h in packet['preflight_bindings'].items():require(digest(out/path)==h,'preflight binding')
 expected=[(c['name'],p) for c in packet['parents'] for p in PRIMES]
 require([(c['parent'],c['p']) for c in packet['choices']]==expected,'choice coverage')
 require(result['input_sha256']==digest(out/'input.json') and result['status']=='PASS','result binding')
 require(len(result['records'])==4,'checkpoint coverage')
 count=0;summary=[]
 for parent,record in zip(packet['parents'],result['records']):
  require(record['path']==parent['name']+'-checkpoint.json' and digest(out/record['path'])==record['sha256'],'checkpoint binding')
  checkpoint=read(out/record['path']);rows=checkpoint['rows']
  require(checkpoint['parent']==parent['name'] and [r['p'] for r in rows]==PRIMES,'prime coverage')
  choices=[c for c in packet['choices'] if c['parent']==parent['name']]
  for row,choice in zip(rows,choices):
   require(all(row[k]==v for k,v in choice.items()),'chosen local fibre')
   check_witness(parent,row);count+=1
  summary.append({'parent':parent['name'],'certified_local_split_fibres':len(rows)})
 require(result['witness_count']==count==56 and result['local_splitting_at_all_tested_primes'] is True,'complete result')
 require(result['global_splitting_divisor']=='UNKNOWN' and result['positive_mw17_target_complete'] is False,'local evidence promoted globally')
 return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
         'witness_count':count,'root_balls_checked':3*count,'parents':summary,
         'independent_exact_arithmetic':True,'global_splitting_divisor':'UNKNOWN',
         'positive_mw17_target_complete':False,'elapsed_seconds':round(time.monotonic()-start,6)}

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT)
 parser.add_argument('--record',type=Path);args=parser.parse_args();result=verify(args.input)
 if args.record:
  with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps(result,sort_keys=True),flush=True)
