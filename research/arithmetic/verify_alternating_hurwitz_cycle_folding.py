#!/usr/bin/env python3
"""Construct and verify braid words reducing DT tuples to doubled-edge trees.

Uniform correctness is proved in ALTERNATING_HURWITZ_CYCLE_FOLDING.md.
Default controls exhaust every transitive degree-eight fixed-product tuple,
then test deterministic braid translates in degrees10,12,26.
Only Python's standard library is used.
"""
from collections import Counter,defaultdict,deque
from itertools import combinations
from pathlib import Path
import argparse,hashlib,importlib.util,json,platform,random

def compose(p,q):return tuple(p[q[i]]for i in range(len(p)))
def inverse(p):
 out=[0]*len(p)
 for i,j in enumerate(p):out[j]=i
 return tuple(out)
def permutation(n,edges):
 p=list(range(n))
 for a,b in edges:p[a],p[b]=p[b],p[a]
 return tuple(p)
def edges(p):return [(i,j)for i,j in enumerate(p)if i<j and p[j]==i]
def ordered_product(tup):
 p=tuple(range(len(tup[0])))
 for q in tup:p=compose(p,q)
 return p
def hurwitz(tup,index,sign):
 a,b=tup[index:index+2]
 pair=(b,compose(compose(b,a),b))if sign==1 else(compose(compose(a,b),a),a)
 return tup[:index]+pair+tup[index+2:]
def graph(tup):
 n=len(tup[0]);occ=[(i,e)for i,p in enumerate(tup)for e in edges(p)]
 assert len(tup)==n//2 and len(occ)==n
 for p in tup:
  assert compose(p,p)==tuple(range(n)) and len(edges(p))==2
 degrees=[0]*n;incident=[[]for _ in range(n)]
 for k,(_,e)in enumerate(occ):
  for v in e:degrees[v]+=1;incident[v].append(k)
 seen={0};queue=[0]
 for v in queue:
  for k in incident[v]:
   for w in occ[k][1]:
    if w not in seen:seen.add(w);queue.append(w)
 assert len(seen)==n
 removed=set();leaves=deque(i for i,d in enumerate(degrees)if d==1)
 while leaves:
  v=leaves.popleft()
  for k in incident[v]:
   if k in removed:continue
   removed.add(k)
   for w in occ[k][1]:
    degrees[w]-=1
    if degrees[w]==1:leaves.append(w)
 core=[(k,i,e)for k,(i,e)in enumerate(occ)if k not in removed]
 assert len(core)>=2 and 0 in {v for _,_,e in core for v in e}
 assert len({i for _,i,_ in core})==len(core)
 # The trajectory of the fixed sheet uses exactly the unique cycle.
 trajectory=[];v=0
 for k in reversed(range(len(occ))):
  _,e=occ[k]
  if v in e:v=e[1]if v==e[0]else e[0];trajectory.append(k)
 assert v==0 and set(trajectory)=={k for k,_,_ in core}
 return occ,core

def normalize(tup):
 initial=tup;product=ordered_product(tup);n=len(product);b=n//2
 assert product[0]==0
 assert len({product[i]for i in range(1,n)})==n-1
 orbit={1};v=product[1]
 while v not in orbit:orbit.add(v);v=product[v]
 assert len(orbit)==n-1
 word=[];folds=0
 def move(index,sign):
  nonlocal tup
  tup=hurwitz(tup,index,sign);word.append(sign*(index+1))
  assert ordered_product(tup)==product
 while True:
  _,core=graph(tup);k=len(core)
  if k==2:break
  colors=sorted(i for _,i,_ in core)
  left,right=colors[:2]
  old_cycle=Counter(e for _,_,e in core)
  for index in range(left,right-1):
   move(index,-1)
   assert Counter(e for _,_,e in graph(tup)[1])==old_cycle
  left=right-1
  _,core=graph(tup)
  core_by_color={i:e for _,i,e in core}
  e=core_by_color[left];f=core_by_color[right]
  common=set(e)&set(f);assert len(common)==1
  v=next(iter(common));assert v!=0
  a=next(x for x in e if x!=v);c=next(x for x in f if x!=v)
  other_a=next(edge for edge in edges(tup[left])if edge!=e)
  other_b=next(edge for edge in edges(tup[right])if edge!=f)
  if a not in other_b:move(left,1)
  elif c not in other_a:move(left,-1)
  else:move(left,1);move(left,1)
  assert len(graph(tup)[1])==k-1
  folds+=1
 # Push the two cycle blocks adjacent, retaining the doubled edge.
 _,core=graph(tup);left,right=sorted(i for _,i,_ in core)
 aux=core[0][2];assert core[1][2]==aux and 0 in aux
 for index in range(left,right-1):move(index,-1)
 occ,core=graph(tup)
 assert len(core)==2 and core[0][2]==core[1][2]==aux
 assert sum(0 in e for _,e in occ)==2
 colors=sorted(i for _,i,_ in core);assert colors[1]==colors[0]+1
 # Independent word replay and complete final graph assertions.
 replay=initial
 for letter in word:replay=hurwitz(replay,abs(letter)-1,1 if letter>0 else-1)
 assert replay==tup
 assert len(word)<=(b+1)*(b-2)
 return tup,word,folds,len(graph(initial)[1])

def coalesce(tup):
 n=len(tup[0]);_,core=graph(tup)
 colors=sorted(i for _,i,_ in core);left,right=colors
 assert right==left+1
 combined=compose(tup[left],tup[right])
 assert combined[0]==0
 old=tup[:left]+(combined,)+tup[right+1:]
 assert all(p[0]==0 for p in old)
 smaller=tuple(tuple(p[i]-1 for i in range(1,n))for p in old)
 assert ordered_product(smaller)==tuple(ordered_product(tup)[i]-1 for i in range(1,n))
 types=[]
 for p in smaller:
  unseen=set(range(n-1));lengths=[]
  while unseen:
   v=next(iter(unseen));cycle=[];w=v
   while w not in cycle:cycle.append(w);unseen.remove(w);w=p[w]
   if len(cycle)>1:lengths.append(len(cycle))
  kind=tuple(sorted(lengths))
  assert kind in((3,),(2,2))
  types.append(kind)
 assert types.count((3,))<=1
 return smaller,types[left]

def seed(n):
 root=Path(__file__).resolve().parent
 spec=importlib.util.spec_from_file_location("matching",root/"verify_real_paired_laurent_matching.py")
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 return tuple(permutation(n,step)for step in module.symmetric_peak_six(n))

def all_transitive_eight():
 n=8;q=3
 A=[(2*i+1,2*i+2)for i in range(q)]
 B=[(2*i+2,2*i+3)for i in range(q)]
 P=compose(permutation(n,A),permutation(n,B))
 C=[permutation(n,es)for a,b,c,d in combinations(range(n),4)
    for es in [((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))]]
 pairs=defaultdict(list)
 for a in C:
  for b in C:pairs[compose(a,b)].append((a,b))
 for ab,left in pairs.items():
  for a,b in left:
   for c,d in pairs.get(compose(inverse(ab),P),()):
    tup=(a,b,c,d)
    if any(p[0]!=0 for p in tup):yield tup

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--write-receipt",type=Path)
 args=parser.parse_args()
 counts=Counter();lengths=Counter();coalesced_types=Counter();fold_total=0;max_word=0;digest=hashlib.sha256()
 for tup in all_transitive_eight():
  out,word,folds,cycle=normalize(tup)
  counts[cycle]+=1;lengths[len(word)]+=1;fold_total+=folds;max_word=max(max_word,len(word))
  _,kind=coalesce(out);coalesced_types[str(kind)]+=1
  digest.update(bytes(x for p in tup for x in p))
  digest.update(bytes((letter+10)for letter in word))
 assert sum(counts.values())==49392
 rng=random.Random(20261118);higher={}
 for n in(10,12,26):
  tup=seed(n);distribution=Counter();moves=0
  for _ in range(1000):
   for __ in range(10):
    tup=hurwitz(tup,rng.randrange(n//2-1),rng.choice((-1,1)))
   out,word,folds,cycle=normalize(tup)
   distribution[cycle]+=1;moves=max(moves,len(word))
   coalesce(out)
  higher[str(n)]={"exact_braid_translate_controls":1000,
                 "initial_cycle_lengths":dict(distribution),"maximum_word_length":moves}
 # Corrupted tuple control: identity is not a double transposition.
 bad=(tuple(range(8)),)+seed(8)[1:]
 try:normalize(bad)
 except AssertionError:pass
 else:raise AssertionError("Invalid DT domain accepted")
 matching=Path(__file__).with_name("verify_real_paired_laurent_matching.py")
 report={
  "schema_version":1,"kind":"exact_braid_cycle_normalization_controls",
  "python":platform.python_version(),
  "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  "matching_seed_checker_sha256":hashlib.sha256(matching.read_bytes()).hexdigest(),
  "degree8_complete_transitive_tuple_count":sum(counts.values()),
  "degree8_initial_cycle_lengths":dict(counts),
  "degree8_braid_word_length_distribution":dict(lengths),
  "degree8_maximum_braid_word_length":max_word,
  "degree8_total_cycle_folds":fold_total,
  "degree8_coalesced_branch_types":dict(coalesced_types),
  "degree8_ordered_certificate_digest":digest.hexdigest(),
  "higher_finite_controls":higher,
  "verified_target":"Unique doubled edge incident fixed sheet0, in adjacent DT blocks; no other edge incident0; deleting0 gives a tree",
  "uniform_word_length_bound":"(b+1)(b-2)",
  "scope":"Written proof gives all-degree normal form. Exhaustive degree8 and higher deterministic exact controls do not prove uniform complex braid connectivity or two-leaf reduction.",
  "independent_replay":False,"external_review":False,"formal_verification":False,
  "result":"PASS"
 }
 output=json.dumps(report,indent=2,sort_keys=True)+"\n"
 if args.write_receipt:args.write_receipt.write_text(output)
 print(output,end="")
if __name__=="__main__":main()

