#!/usr/bin/env sage-python
"""Exact structural marking gates; general geometric implications are written inputs."""
import argparse, hashlib, itertools, json
from pathlib import Path
from collections import Counter
from sage.all import *
ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/'artifacts/generated-results/elkies-k3-rank7-t-arithmetic-v1.json'
PERIOD=ROOT/'artifacts/generated-results/elkies-k3-real-rank19-period-v1.json'
DET800=ROOT/'artifacts/generated-results/elkies-k3-det800-full-marking-v1.json'
OUTPUT=ROOT/'artifacts/generated-results/elkies-k3-marking-structural-gates-v1.json'
cat=json.loads(INPUT.read_text())['surfaces']
period=json.loads(PERIOD.read_text())
det800=json.loads(DET800.read_text())
assert det800['rational_non_CM_marking_exists'] is False
frontier=set(period['surviving_queue_ids'])-{det800['surface_id']}
assert len(frontier)==721
rows=[]
for row in cat:
 if row['similarity_normalization']['literal_content']!=2 or not row['rational_isotropy']['isotropic']:continue
 G=matrix(ZZ,row['literal_transcendental_gram']);g=G/2
 if any(g[i,i]%2 for i in range(3)):continue
 e=vector(ZZ,row['rational_isotropy']['primitive_isotropic_vector']);z=g*e
 if gcd(z)!=1:continue
 h,a,b=xgcd(z[0],z[1]);one,c,d=xgcd(h,z[2]);assert one==1
 w=vector(ZZ,[c*a,c*b,d]);assert e*g*w==1 and e*g*e==0
 f=w-(w*g*w/2)*e
 k=matrix(ZZ,[list(e*g),list(f*g)]).right_kernel().basis()[0]
 B=matrix(ZZ,[e,k,f]).transpose();assert abs(B.det())==1
 standard=B.transpose()*G*B
 N=int(standard[1,1]/4)
 assert standard==matrix(ZZ,[[0,0,2],[0,4*N,0],[2,0,0]])
 assert N>=5
 rows.append({'surface_id':row['surface_id'],'determinant':row['determinant'],'N':N,'basis':[[int(x) for x in r] for r in B.rows()],'gram':row['literal_transcendental_gram']})

assert len(rows)==84
split_ids={r['surface_id'] for r in rows}
assert len(split_ids & frontier)==83
assert split_ids-frontier=={det800['surface_id']}
local_rows=[]
eligible=0
moduli=[8,16,3,5,7,11,13]
for r in cat:
    if r['surface_id'] not in frontier or r['rational_isotropy']['isotropic']:
        continue
    eligible+=1
    g=r['literal_transcendental_gram']
    a,b,c=g[0][0],g[1][1],g[2][2]
    d,e,f=2*g[0][1],2*g[0][2],2*g[1][2]
    for m in moduli:
        witness=next(((x,y,z) for x,y,z in itertools.product(range(m),repeat=3)
                      if (a*x*x+b*y*y+c*z*z+d*x*y+e*x*z+f*y*z-2)%m==0),None)
        if witness is not None:
            continue
        q=QuadraticForm(QQ,matrix(QQ,g));dd=q.rational_diagonal_form()
        aa,bb,cc=[dd[i,i] for i in range(3)]
        prime=int(r['rational_isotropy']['pari_qfsolve_obstruction_prime'])
        assert hilbert_symbol(-aa*bb,-aa*cc,prime)==-1
        local_rows.append({'surface_id':r['surface_id'],'determinant':r['determinant'],
                           'gram':g,'obstructing_modulus':m,'residues_exhausted':m**3,
                           'anisotropy_prime':prime,'hilbert_parameters':[str(-aa*bb),str(-aa*cc)],
                           'hilbert_symbol':-1})
        break
local_ids={r['surface_id'] for r in local_rows}
assert eligible==218 and len(local_ids)==36
assert not local_ids & split_ids
control_rows=[]
for N in [1,2,3,4,5,6,7,8,10,12,35,41,50]:
 G=matrix(QQ,[[0,0,2],[0,4*N,0],[2,0,0]]);Gi=G.inverse();I=identity_matrix(QQ,3)
 def action(g):
  a,b,c,d=g.list();delta=g.det()
  A=matrix(QQ,[[a*a,-2*a*b,-b*b/N],[-a*c,a*d+b*c,b*d/N],[-N*c*c,2*N*c*d,d*d]])/delta
  assert A.transpose()*G*A==G
  return A
 W=I
 if N%2==0:
  n=2**valuation(N,2);_,v,u=xgcd(n,N//n)
  W0=matrix(QQ,[[n,-u],[N,n*v]])
  J=matrix(QQ,[[0,-1],[1,0]]);W=J.inverse()*W0*J
 counts=Counter()
 for a,b,c,d in itertools.product(range(4),repeat=4):
  g=matrix(QQ,[[a,N*b],[c,d]])
  if g.det()%2!=1:continue
  for coset in range(2 if N%2==0 else 1):
   A=action(g if not coset else W*g)
   for sign in [1,-1]:
    stable=all(x==0 or x.valuation(2)>=0 for x in ((sign*A-I)*Gi).list())
    expected=coset==0 and sign==1 and b%2==c%2==0
    assert stable==expected,(N,a,b,c,d,coset,sign)
    counts[(coset,sign)]+=stable
 control_rows.append({'N':N,'counts':{str(k):int(v) for k,v in counts.items()}})

remaining=[x for x in period['surviving_queue_ids'] if x!=det800['surface_id'] and x not in split_ids|local_ids]
assert len(remaining)==602
result={'schema':'elkies-k3.marking-structural-gates.v1',
        'status':'PASS_EXACT_AUDIT_WITH_WRITTEN_GENERAL_THEOREMS',
        'inputs':[{'path':str(p.relative_to(ROOT)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
                  for p in [INPUT,PERIOD,DET800]],
        'scaled_split_rows':rows,
        'local_square2_rows':local_rows,
        'local2_controls':control_rows,
        'moduli':moduli,
        'counts':{'starting_frontier':len(frontier),'scaled_split_total':len(rows),
                  'scaled_split_new_exclusions':len(split_ids&frontier),
                  'anisotropic_rows_tested':eligible,'new_local_square2_exclusions':len(local_ids),
                  'unresolved_remaining':len(remaining)},
        'surviving_queue_ids':remaining,'new_positive_handoff':[],
        'boundary':'General parity/normalizer/period arguments are written theorem inputs; surviving residues or rows are not marking certificates.'}
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--check',action='store_true')
args=parser.parse_args()
if args.check:
    assert json.loads(OUTPUT.read_text())==result
else:
    OUTPUT.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result['counts'],sort_keys=True))
