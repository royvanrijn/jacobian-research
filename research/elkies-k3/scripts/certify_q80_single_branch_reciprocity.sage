#!/usr/bin/env sage-python
"""Finite character witnesses for a single-branch Weil-reciprocity gate."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import GF, PolynomialRing, QQ
from sage.version import version

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-single-branch-reciprocity-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
PARENT='artifacts/generated-results/elkies-k3-mestre-parent-branch-cancellation-v1/input.json'
CHECKER='elkies-k3/scripts/verify_q80_single_branch_reciprocity.py'

def digest(path):return sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text())
def write_new(path,data):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')

def freeze(out):
 write_new(out/'input.json',{'schema':'q80-single-branch-reciprocity-input-v1',
  'source':SOURCE,'source_sha256':digest(ROOT/SOURCE),'parent':PARENT,'parent_sha256':digest(ROOT/PARENT),
  'producer_sha256':digest(Path(__file__)),'checker_sha256':digest(ROOT/CHECKER),
  'prime':131,'section_indices':list(range(17)),'sage_version':version,
  'cpu_seconds':40,'address_space_bytes':4*1024**3,
  'selection':'One retained good prime and the original17 equation-side generic sections. The full110-root preview discovered this proof before freezing; no coefficient, point or branch-field search.',
  'boundary':'Finite inherited-section characters, plus a written Weil-reciprocity deduction. No genus0/1 cover, global root-curve divisor or rank19 source is constructed.'})
 print('Frozen Q80, prime131 and17 inherited sections',flush=True)

def run(out):
 packet=read(out/'input.json');start=time.process_time()
 assert packet['source_sha256']==digest(ROOT/SOURCE) and packet['parent_sha256']==digest(ROOT/PARENT)
 assert packet['producer_sha256']==digest(Path(__file__)) and packet['checker_sha256']==digest(ROOT/CHECKER)
 data=read(ROOT/SOURCE);model=data['weierstrass_model'];parent=next(p for p in read(ROOT/PARENT)['parents'] if p['name']=='alternate-q80')
 A0=list(map(QQ,model['A_coefficients_low_to_high']));B0=list(map(QQ,model['B_coefficients_low_to_high']))
 assert A0==list(map(QQ,parent['A'])) and B0==list(map(QQ,parent['B']))
 p=packet['prime'];R0=PolynomialRing(QQ,'t');K0=R0.fraction_field()
 R=PolynomialRing(GF(p),'t');t=R.gen();K=R.fraction_field();A,B=R(A0),R(B0)
 delta=4*A**3+27*B**2
 assert delta.degree()==24 and delta.gcd(delta.derivative())==1
 sections=[]
 for s in data['sections']['records']:
  assert s['basis_index'] in packet['section_indices'];xy=[]
  for key in ['X','Y']:
   q=s[key];xy.append(K0(R0(list(map(QQ,q['numerator_coefficients_low_to_high']))))/R0(list(map(QQ,q['denominator_coefficients_low_to_high']))))
  x,y=xy;assert y*y==x**3+R0(A0)*x+R0(B0)
  xx,yy=[K(R(z.numerator()))/R(z.denominator()) for z in xy]
  assert yy*yy==xx**3+A*xx+B
  sections.append((s['basis_index'],xx,yy))
 def at(f,u,n):
  if u is not None:return f(u)
  return K(t**n*f(1/t))(0)
 fibres=[];witnesses=[]
 for u in list(GF(p))+[None]:
  a,b=(A[8],B[12]) if u is None else (A(u),B(u))
  assert 4*a**3+27*b*b
  roots=[r for r in GF(p) if r**3+a*r+b==0]
  fibres.append({'t':None if u is None else int(u),'A':int(a),'B':int(b),'roots':list(map(int,roots))})
  for root in roots:
   found=None
   for index,x,y in sections:
    try:xx,yy=at(x,u,4),at(y,u,6)
    except ZeroDivisionError:continue
    if yy and not (xx-root).is_square():
     found={'t':None if u is None else int(u),'root':int(root),'section_index':index,
            'section_x':int(xx),'section_y':int(yy),'difference':int(xx-root)};break
   assert found is not None,(u,root)
   witnesses.append(found)
 assert len(witnesses)==110
 finite={'A':list(map(int,A)),'B':list(map(int,B)),'discriminant':list(map(int,delta)),
         'fibres':fibres,'all_parent_fibres_smooth_over_Fp':True}
 write_new(out/'finite-model.json',finite);write_new(out/'witnesses.json',{'prime':p,'rows':witnesses})
 write_new(out/'result.json',{'schema':'q80-single-branch-reciprocity-result-v1','status':'PASS',
  'input_sha256':digest(out/'input.json'),
  'records':{name:digest(out/name) for name in ['finite-model.json','witnesses.json']},
  'projective_base_values':p+1,'root_points':len(witnesses),'fully_split_fibres':sum(len(f['roots'])==3 for f in fibres),
  'used_sections':sorted({w['section_index'] for w in witnesses}),
  'rational_single_complement_ramification':'EXCLUDED',
  'q80_genus0_k1':'EXCLUDED','q80_genus1_k3':'EXCLUDED',
  'positive_mw17_target_complete':False,'cpu_seconds':round(time.process_time()-start,6)})
 print('PASS:110 root points over132 smooth fibres, with exact17-section identities',flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['freeze','run'])
 parser.add_argument('--output',type=Path,default=DEFAULT);args=parser.parse_args()
 resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
 (freeze if args.mode=='freeze' else run)(args.output)
