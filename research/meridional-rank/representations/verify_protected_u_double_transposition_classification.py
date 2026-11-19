#!/usr/bin/env python3
"""Complete protected-U DT classification and component-slot controls.
The written note proves the unrestricted tail consequence; this program
replays the finite classification exactly and checks all local certificates."""
if not __debug__:
    raise RuntimeError("This exact classification requires Python assertions enabled.")

import argparse,hashlib,importlib.util,itertools,json,platform
from collections import Counter
from pathlib import Path
import numpy as np
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3]
def compose(a,b):return tuple(a[b[i]]for i in range(len(a)))
def cs(a,b):return compose(compose(a,b),a)
def make(n,edges):
 p=list(range(n))
 for a,b in edges:p[a]=b;p[b]=a
 return tuple(p)
def colors(n,size):
 if size==0:return [tuple(range(n))]
 if size==1:return [make(n,[(a,b)])for a,b in itertools.combinations(range(n),2)]
 out=[]
 for a,b,c,d in itertools.combinations(range(n),4):
  out.extend([make(n,[(a,b),(c,d)]),make(n,[(a,c),(b,d)]),make(n,[(a,d),(b,c)])])
 return out
def connected(t):
 n=len(t[0]);seen={0};q=[0]
 while q:
  v=q.pop()
  for p in t:
   if p[v]not in seen:seen.add(p[v]);q.append(p[v])
 return len(seen)==n
def action(q,table,word=U):
 q=[a.copy()for a in q]
 for letter in word:
  j=abs(letter)-1;a,b=q[j:j+2]
  q[j:j+2]=[b,table[b,a]]if letter>0 else[table[a,b],a]
 return q
def count(p):return sum(p[i]!=i for i in range(len(p)))//2
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2] if HERE.name=='representations' else None
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--write-receipt',type=Path)
parser.add_argument('--forest-checker',type=Path,default=HERE/'verify_protected_block_forest_colorings.py')
args=parser.parse_args()
spec=importlib.util.spec_from_file_location('protected_u_forest_input',args.forest_checker)
forest=importlib.util.module_from_spec(spec);spec.loader.exec_module(forest)
forest_census,tree_witnesses=forest.cycle_normalized_census()
assert set(tree_witnesses)=={(0,0,1,1)}
expected_census=[{'vertices': 2, 'first_color_edges': 0, 'tested': 4, 'outer_returns': 1, 'connected': 1, 'orbits': 1}, {'vertices': 2, 'first_color_edges': 1, 'tested': 7, 'outer_returns': 2, 'connected': 2, 'orbits': 2}, {'vertices': 3, 'first_color_edges': 0, 'tested': 27, 'outer_returns': 0, 'connected': 0, 'orbits': 0}, {'vertices': 3, 'first_color_edges': 1, 'tested': 54, 'outer_returns': 3, 'connected': 2, 'orbits': 1}, {'vertices': 4, 'first_color_edges': 0, 'tested': 540, 'outer_returns': 3, 'connected': 0, 'orbits': 0}, {'vertices': 4, 'first_color_edges': 1, 'tested': 864, 'outer_returns': 9, 'connected': 2, 'orbits': 1}, {'vertices': 4, 'first_color_edges': 2, 'tested': 981, 'outer_returns': 18, 'connected': 14, 'orbits': 3}, {'vertices': 5, 'first_color_edges': 0, 'tested': 10125, 'outer_returns': 0, 'connected': 0, 'orbits': 0}, {'vertices': 5, 'first_color_edges': 1, 'tested': 15300, 'outer_returns': 45, 'connected': 24, 'orbits': 2}, {'vertices': 5, 'first_color_edges': 2, 'tested': 17200, 'outer_returns': 49, 'connected': 24, 'orbits': 3}, {'vertices': 6, 'first_color_edges': 0, 'tested': 91125, 'outer_returns': 0, 'connected': 0, 'orbits': 0}, {'vertices': 6, 'first_color_edges': 1, 'tested': 182250, 'outer_returns': 237, 'connected': 24, 'orbits': 1}, {'vertices': 6, 'first_color_edges': 2, 'tested': 218700, 'outer_returns': 196, 'connected': 88, 'orbits': 6}, {'vertices': 7, 'first_color_edges': 1, 'tested': 1157625, 'outer_returns': 0, 'connected': 0, 'orbits': 0}, {'vertices': 7, 'first_color_edges': 2, 'tested': 1852200, 'outer_returns': 573, 'connected': 48, 'orbits': 1}, {'vertices': 8, 'first_color_edges': 2, 'tested': 9261000, 'outer_returns': 1906, 'connected': 192, 'orbits': 1}]
components=[];census=[]
for n in range(2,9):
 groups=[colors(n,i)for i in range(3)];ps=[p for g in groups for p in g];lookup={p:i for i,p in enumerate(ps)}
 weights=np.array([count(p)for p in ps],dtype=np.int8)
 table=np.array([[lookup[cs(p,q)]for q in ps]for p in ps],dtype=np.int16)
 # First-color class capacity rules are exact bounds, not sampling.
 for first in range(3):
  if not groups[first]or first+6<n:continue
  aa=lookup[make(n,[]if first==0 else[(0,1)]if first==1 else[(0,1),(2,3)])]
  # Direct all-triples vectorization, with only the total edge budget filter.
  grid=np.indices((len(ps),)*3,dtype=np.int16).reshape(3,-1)
  e=first+weights[grid[0]]+weights[grid[1]]+weights[grid[2]]
  keep=(e>=n)&(e<=8);grid=grid[:,keep];e=e[keep]
  q=[np.full(grid.shape[1],aa,dtype=np.int16),*grid]
  after=action(q,table)
  selected=np.flatnonzero((after[0]==aa)&(after[3]==grid[2]))
  good=[tuple(ps[int(q[k][i])]for k in range(4))for i in selected if connected(tuple(ps[int(q[k][i])]for k in range(4)))]
  # Rebuild the whole centralizer by brute S_n; independent of C++ indexing.
  central=[]
  if good:
   for g in itertools.permutations(range(n)):
    gi=[0]*n
    for i,j in enumerate(g):gi[j]=i
    if tuple(g[ps[aa][gi[i]]]for i in range(n))==ps[aa]:central.append((g,gi))
  def key(t):return min(tuple(tuple(g[p[gi[i]]]for i in range(n))for p in t)for g,gi in central)
  keys={key(t)for t in good}
  for t in sorted(keys):
   counts=[count(p)for p in t]
   components.append({'vertices':n,'counts':counts,'cycle_rank':sum(counts)-n+1,'colors':t})
  row={'vertices':n,'first_color_edges':first,'tested':q[0].size,'outer_returns':len(selected),'connected':len(good),'orbits':len(keys)}
  census.append(row)
assert census==expected_census
assert len(components)==22
# A connected forest component is the sole imported finite forest input;
# it has just been independently replayed by its complete existing checker.
tree=tree_witnesses[(0,0,1,1)]
assert len(tree[0])==3
components.append({'cycle_rank':0,'vertices':3,'counts':[0,0,1,1],'colors':tree})
def braid_tuple(t,word=U):
 t=list(t)
 for v in word:
  i=abs(v)-1;a,b=t[i:i+2]
  t[i:i+2]=[b,cs(b,a)]if v>0 else[cs(a,b),a]
 return t
def edges(p):return [(i,j)for i,j in enumerate(p)if i<j]
def mobile_components(edges4,n):
 adj=[set()for _ in range(n)]
 for a,b in edges4:adj[a].add(b);adj[b].add(a)
 loc={};parts=[]
 for a in range(n):
  if not adj[a]or a in loc:continue
  seen={a};q=[a]
  while q:
   v=q.pop()
   for w in adj[v]-seen:seen.add(w);q.append(w)
  for w in seen:loc[w]=len(parts)
  parts.append(sorted(seen))
 return loc,parts
def cycle_count(p):
 seen=set();c=0
 for i in range(len(p)):
  if i in seen:continue
  c+=1;j=i
  while j not in seen:seen.add(j);j=p[j]
 return c
full=[]
for num in range(1,5):
 for ids in itertools.combinations_with_replacement(range(len(components)),num):
  rs=[components[i]for i in ids]
  if [sum(r['counts'][j]for r in rs)for j in range(4)]!=[2]*4:continue
  n=sum(r['vertices']for r in rs);colors4=[list(range(n))for _ in range(4)];offset=0
  for r in rs:
   for j in range(4):
    for i,a in enumerate(r['colors'][j]):colors4[j][offset+i]=offset+a
   offset+=r['vertices']
  colors4=tuple(map(tuple,colors4));after=braid_tuple(colors4)
  assert after[0]==colors4[0]and after[3]==colors4[3]
  assert all(count(p)==2 for p in colors4)
  inp=[edges(p)for p in colors4];oup=[edges(p)for p in after]
  ie=inp[1]+inp[2];oe=oup[1]+oup[2]
  loc,parts=mobile_components(ie,n)
  input_labels=[loc[a]for a,b in ie];output_labels=[]
  for a,b in oe:
   assert a in loc and b in loc and loc[a]==loc[b]
   output_labels.append(loc[a])
  cost=5;perm=None;flips=None
  for ff in itertools.product(range(2),repeat=2):
   target=[output_labels[2*j+(e^ff[j])]for j in range(2)for e in range(2)]
   for pp in itertools.permutations(range(4)):
    if [input_labels[a]for a in pp]!=target:continue
    cc=4-cycle_count(pp)
    if cc<cost:cost=cc;perm=pp;flips=ff
  missing=n-len(loc)
  assert cost+missing<=2
  full.append({'component_ids':ids,'vertices':n,'cycle_rank':sum(r['cycle_rank']for r in rs),
               'input_colors':colors4,'prefix_output_colors':after,
               'local_mobile_components':parts,'input_component_slots':input_labels,
               'output_component_slots':output_labels,'output_edge_order_flips':flips,
               'component_slot_permutation':perm,'component_slot_cost':cost,
               'omitted_mobile_symbols':missing,'cost_plus_omission':cost+missing,
               'same_geometric_mobile_edge_multiset':Counter(ie)==Counter(oe)})
assert len(full)==31
assert Counter(r['cycle_rank']for r in full)=={1:2,2:4,3:11,4:10,5:3,6:1}
assert sum(not r['same_geometric_mobile_edge_multiset']for r in full)==7
assert Counter(r['component_slot_cost']for r in full)=={0:28,1:3}
assert Counter((r['component_slot_cost'],r['omitted_mobile_symbols'])for r in full)=={(0,0):17,(0,1):9,(0,2):2,(1,0):1,(1,1):2}
# Positive and crossing-sign mutation controls use an explicit geometric tuple.
old=tuple(make(8,e)for e in [[(0,1),(2,3)],[(1,4),(2,3)],[(0,5),(1,6)],[(1,2),(6,7)]])
assert braid_tuple(old)[0]==old[0]and braid_tuple(old)[3]==old[3]
mutations=0
for i in range(len(U)):
 bad=U.copy();bad[i]*=-1;after=braid_tuple(old,bad)
 mutations+=after[0]!=old[0]or after[3]!=old[3]
assert mutations>0
# Exhaustive small controls of the unrestricted slot-cycle inequality.
slot_controls=0
for k in range(1,4):
 for extra in range(4-k):
  m=2*k+extra
  for active in itertools.product(range(2),repeat=k):
   tau=list(range(2*m));cost=sum(active)
   for j,on in enumerate(active):
    if on:tau[4*j],tau[4*j+2]=tau[4*j+2],tau[4*j]
   for rest in itertools.permutations(range(1,m)):
    order=(0,)+rest;rho=[0]*m
    for a,b in zip(order,order[1:]+order[:1]):rho[a]=b
    for flips in itertools.product(range(2),repeat=m):
     pi=[2*rho[j]+(e^flips[j])for j in range(m)for e in range(2)]
     assert cycle_count(pi)<=2
     assert cycle_count([pi[tau[i]]for i in range(2*m)])<=cost+2
     slot_controls+=1
# Without the one-cycle knot condition the component inequality can fail.
tau=list(range(12))
for j in range(3):tau[4*j],tau[4*j+2]=tau[4*j+2],tau[4*j]
assert cycle_count(tau)==9>3+2
sources=[Path(__file__),args.forest_checker]
for name in ['PROTECTED_U_DOUBLE_TRANSPOSITION_CLASSIFICATION.md','MULTIPLE_PROTECTED_OBSTRUCTIONS.md',
             'FREE_MERIDIAN_STABILIZERS.md','TWO_SIDED_PROTECTED_OBSTRUCTIONS.md']:
 p=HERE/name
 if p.exists():sources.append(p)
def label(p):
 p=p.resolve()
 return p.relative_to(ROOT).as_posix()if ROOT is not None else p.name
receipt={'status':'passed','python':platform.python_version(),'numpy':np.__version__,
         'scope':'complete local DT classification; exact local component-slot certificates; finite controls supporting the unrestricted written all-k tail barrier',
         'source_sha256':{label(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in sources},
         'complete_forest_census':forest_census,'complete_cyclic_component_census':census,
         'connected_component_representatives':components,'full_local_type_count':len(full),
         'full_cycle_rank_distribution':dict(sorted(Counter(r['cycle_rank']for r in full).items())),
         'local_component_slot_certificates':full,'slot_control_count':slot_controls,
         'crossing_sign_mutations_rejected':mutations,
         'assurance':{'internal_replay':True,'formal_verification':False,'external_review':False},
         'does_not_establish':['actual knot gap upper bound','other finite meridian classes','other meridional invariants','overlapping obstruction blocks','prime or hyperbolic gap-two knot']}
if args.write_receipt:args.write_receipt.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'passed','cyclic_component_types':22,'local_DT_types':31,'component_cost_histogram':dict(Counter(r['component_slot_cost']for r in full)),'slot_controls':slot_controls}))

