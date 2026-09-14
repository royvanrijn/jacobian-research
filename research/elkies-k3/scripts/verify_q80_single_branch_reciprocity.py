#!/usr/bin/env python3
"""Independent rational identities and finite root-character certificate."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-single-branch-reciprocity-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
PARENT='artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'

def require(test,message):
 if not test:raise ValueError(message)
def digest(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,v in enumerate(a):c[i]+=v
 for i,v in enumerate(b):c[i]+=v
 return c
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def trim(a):
 while a and not a[-1]:a.pop()
 return a
def rem(a,b,p):
 a=trim(a[:])
 while len(a)>=len(b):
  v=a[-1]*pow(b[-1],-1,p)%p;shift=len(a)-len(b)
  for i,c in enumerate(b):a[shift+i]=(a[shift+i]-v*c)%p
  trim(a)
 return a
def squarefree_degree24(d,p):
 require(len(d)==25 and d[-1]%p,'good infinity')
 a=[c%p for c in d];b=trim([i*d[i]%p for i in range(1,25)])
 while b:a,b=b,rem(a,b,p)
 require(len(a)==1,'bad curve reduction')
def mod(q,p):
 q=F(q);require(q.denominator%p,'nonintegral coefficient')
 return q.numerator*pow(q.denominator,-1,p)%p
def ev(a,t,p):
 v=0
 for c in reversed(a):v=(v*t+c)%p
 return v
def rational_at(n,d,t,weight,p):
 if t is not None:
  a,b=ev(n,t,p),ev(d,t,p)
 else:
  n,d=trim(n[:]),trim(d[:]);require(bool(d),'zero section denominator')
  shift=len(d)-1+weight
  require(len(n)-1<=shift,'section pole at infinity')
  a,b=(n[shift] if len(n)>shift else 0),d[-1]
 require(b!=0,'section pole')
 return a*pow(b,-1,p)%p

def check_character(a,b,troot,x,y,p):
 require((4*a**3+27*b*b)%p!=0,'singular parent fibre')
 require((troot**3+a*troot+b)%p==0,'not a cubic root')
 require((y*y-x**3-a*x-b)%p==0,'section point identity')
 require(y%p!=0,'section meets two-torsion')
 value=(x-troot)%p
 require(value and pow(value,(p-1)//2,p)==p-1,'character is not a nonzero nonsquare')
 return value

def section_data(section,A,B,p):
 values=[]
 for key in ['X','Y']:
  f=section[key]
  values.extend([list(map(F,f['numerator_coefficients_low_to_high'])),list(map(F,f['denominator_coefficients_low_to_high']))])
 xn,xd,yn,yd=values
 require(any(xd) and any(yd),'zero rational denominator')
 xd2=mul(xd,xd);xd3=mul(xd2,xd)
 rhs=add(add(mul(xn,mul(xn,xn)),mul(A,mul(xn,xd2))),mul(B,xd3))
 lhs=mul(mul(yn,yn),xd3);rhs=mul(mul(yd,yd),rhs)
 require(not any(add(lhs,[-x for x in rhs])),'characteristic-zero section identity')
 return [[mod(c,p) for c in a] for a in values]

def verify(out):
 start=time.monotonic();packet,result=read(out/'input.json'),read(out/'result.json')
 require(packet['source']==SOURCE and packet['source_sha256']==digest(ROOT/SOURCE),'source binding')
 require(packet['parent']==PARENT and packet['parent_sha256']==digest(ROOT/PARENT),'parent binding')
 require(packet['producer_sha256']==digest(ROOT/'elkies-k3/scripts/certify_q80_single_branch_reciprocity.sage'),'producer binding')
 require(packet['checker_sha256']==digest(Path(__file__)),'checker binding')
 require(packet['prime']==131 and packet['section_indices']==list(range(17)),'frozen scope')
 require(result['input_sha256']==digest(out/'input.json') and result['status']=='PASS','result binding')
 require(set(result['records'])=={'finite-model.json','witnesses.json'},'record scope')
 for path,h in result['records'].items():require(digest(out/path)==h,'record binding')
 p=131;source=read(ROOT/SOURCE);model=source['weierstrass_model']
 A0=list(map(F,model['A_coefficients_low_to_high']));B0=list(map(F,model['B_coefficients_low_to_high']))
 parent=next(c for c in read(ROOT/PARENT)['parents'] if c['name']=='alternate-q80')
 require(A0==list(map(F,parent['A'])) and B0==list(map(F,parent['B'])),'literal fixed parent')
 A,B=[[mod(c,p) for c in a] for a in [A0,B0]]
 delta=[c%p for c in add([4*c for c in mul(A,mul(A,A))],[27*c for c in mul(B,B)])]
 squarefree_degree24(delta,p)
 finite=read(out/'finite-model.json');require((A,B,delta)==(finite['A'],finite['B'],finite['discriminant']),'finite model identity')
 fibres=[];root_points=[]
 for t in list(range(p))+[None]:
  a,b=(A[8],B[12]) if t is None else (ev(A,t,p),ev(B,t,p))
  require((4*a**3+27*b*b)%p,'singular projective residue fibre')
  roots=[x for x in range(p) if (x**3+a*x+b)%p==0]
  fibres.append({'t':t,'A':a,'B':b,'roots':roots})
  root_points.extend((t,x) for x in roots)
 require(finite['fibres']==fibres and finite['all_parent_fibres_smooth_over_Fp'] is True,'complete projective roster')
 sections=source['sections']['records'];require([s['basis_index'] for s in sections]==list(range(17)),'section source indexing')
 reduced=[section_data(s,A0,B0,p) for s in sections]
 witness=read(out/'witnesses.json');rows=witness['rows']
 require(witness['prime']==p and [(w['t'],w['root']) for w in rows]==root_points,'root-point coverage')
 for w in rows:
  t=w['t'];i=w['section_index'];require(i in range(17),'unknown inherited section')
  xn,xd,yn,yd=reduced[i]
  x,y=rational_at(xn,xd,t,4,p),rational_at(yn,yd,t,6,p)
  a,b=(A[8],B[12]) if t is None else (ev(A,t,p),ev(B,t,p))
  value=check_character(a,b,w['root'],x,y,p)
  require((x,y,value)==(w['section_x'],w['section_y'],w['difference']),'selected character attachment')
 require(result['root_points']==len(rows)==110 and result['projective_base_values']==132,'result coverage')
 require(result['fully_split_fibres']==sum(len(f['roots'])==3 for f in fibres)==16,'full-splitting roster')
 require(result['used_sections']==sorted({w['section_index'] for w in rows}),'used section list')
 require(result['rational_single_complement_ramification']==result['q80_genus0_k1']==result['q80_genus1_k3']=='EXCLUDED','result theorem boundary')
 require(result['positive_mw17_target_complete'] is False,'false positive endpoint')
 return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
  'exact_rational_section_identities':17,'projective_fibres':132,'root_points':110,'character_witnesses':len(rows),
  'independent_arithmetic_replay':True,'written_reciprocity_formally_verified':False,
  'positive_mw17_target_complete':False,'elapsed_seconds':round(time.monotonic()-start,6)}

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT)
 parser.add_argument('--record',type=Path);args=parser.parse_args();result=verify(args.input)
 if args.record:
  with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps(result,sort_keys=True),flush=True)
