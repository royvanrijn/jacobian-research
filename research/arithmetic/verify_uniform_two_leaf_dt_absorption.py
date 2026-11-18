#!/usr/bin/env python3
"""Exact six-core absorption and sharp degree-seven extension boundary.

Finite N9/N11/N13 controller endpoints are substituted; no unrecovered
controller word replay is claimed. The written proof supplies all-r scope.
"""
import argparse,hashlib,itertools,json,platform
from collections import deque
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
import verify_uniform_hurwitz_attachment as h
import verify_d5_lobe_core_absorption as u
import verify_imprimitive_leaf_normal_form as normal

BASE=[[(0,1),(2,10)],[(0,4),(1,5)],[(0,6),(1,7)],
      [(0,8),(1,9)],[(0,2),(1,3)]]
CASES={
 "G":{"eta":[(8,11),(6,12)],"anchors":[(11,9),(2,3)],"leaf":8},
 "E1":{"eta":[(2,11),(8,12)],"anchors":[(11,3)],"leaf":0},
 "E2":{"eta":[(10,11),(8,12)],"anchors":[(2,3)],"leaf":0},
 "E3":{"eta":[(0,11),(8,12)],"anchors":[(2,3)],"leaf":0},
 "E4":{"eta":[(0,11),(6,12)],"anchors":[(2,3)],"leaf":0},
 "E5":{"eta":[(0,11),(4,12)],"anchors":[(2,3)],"leaf":0},
}
def power(p,k):
 if k<0:return power(h.inverse(p),-k)
 out=tuple(range(len(p)))
 while k:
  if k&1:out=h.compose(p,out)
  p=h.compose(p,p);k//=2
 return out
def invword(w):return [-a for a in reversed(w)]
def support(g):return {i for i,a in enumerate(g)if a!=i}
def initial(r,a,b):
 n=2*r+3
 return tuple(g+(2*r+1,2*r+2)for g in u.wr(r))+(h.permutation(n,[(a,2*r+1),(b,2*r+2)]),)
def cycle_target(core,leaf):
 n=len(core[0]);m=(n-1)//2;std=u.wr(m);P=h.product(std);C=h.product(core)
 assert len(h.cycles(C))==1
 q=[None]*n;a,b=2*m,leaf
 for _ in range(n):q[a]=b;a=P[a];b=C[b]
 assert sorted(q)==list(range(n))
 target=tuple(h.conjugate(tuple(q),g)for g in std)
 assert h.product(target)==C
 return target
def case_data(r,d):
 if d<=r-5:return "G",r-d-2,(2*(r-d-2),2*(r-2))
 if d==r-4:return "E1",r-3,(2,2*(r-3))
 if d==r-3:return "E2",r-3,(2*r,2*(r-3))
 return "E"+str(d-(r-2)+3),r-3,(0,2*(d-1))
def emap(r,j):
 return {0:0,1:1,2:2,3:3,4:2*r-2,5:2*r-1,
         6:2*r-4,7:2*r-3,8:2*j,9:2*j+1,
         10:2*r,11:2*r+1,12:2*r+2}
def extract(t,selected):
 colors=list(range(len(t)));word=[];done=0
 for color in sorted(selected):
  pos=colors.index(color)
  for i in range(pos-1,done-1,-1):
   t=h.hurwitz(t,i,1);colors[i],colors[i+1]=colors[i+1],colors[i]
   word.append(i+1)
  done+=1
 return t,word
def six_core_controls():
 ans={}
 for name,row in CASES.items():
  core=tuple(h.permutation(13,cs)for cs in BASE+[row["eta"]])
  P=h.product(core);leaf=row["leaf"]
  pairs=u.antipodal_pairs(P,leaf)
  assert all(frozenset(a)in pairs for a in row["anchors"])
  target=cycle_target(core,leaf)
  assert normal.lower_class(target,leaf)[0]=="W"
  ans[name]={"core":core,"target":target,"product":h.cycles(P)[0]}
 return ans
def normalize_prefix(r,d,cores,do_normal=True,wrong_leaf=None):
 name,j,(a,b)=case_data(r,d)
 t=initial(r,a,b);n=2*r+3;em=emap(r,j)
 selected={0,1,2,r-j,r-1,r}
 assert len(selected)==6
 t,word=extract(t,selected)
 assert len(word)<=5*(r-5)
 core=tuple(u.embed(g,em,n)for g in cores[name]["core"])
 assert t[:6]==core
 outside=t[6:];S=set(em.values())
 expected_pairs={frozenset(em[a]for a in pair)for pair in CASES[name]["anchors"]}
 for g in outside:
  pair=support(g)&S
  assert frozenset(pair)in expected_pairs and len(pair)==2
  assert all(len(set(e)&S)==1 for e in h.edges(g))
 leaf=em[CASES[name]["leaf"]]
 target=cores[name]["target"]
 if wrong_leaf is not None:
  target=cycle_target(cores[name]["core"],wrong_leaf);leaf=em[wrong_leaf]
 substituted=tuple(u.embed(g,em,n)for g in target)+outside
 assert h.product(substituted)==h.product(t)
 assert sum(g[leaf]!=leaf for g in substituted)==1
 if wrong_leaf is not None:return normal.lower_class(substituted,leaf)[0]
 assert normal.lower_class(substituted,leaf)[0]=="W"
 length=0
 if do_normal:
  ww,sigma,push=normal.normal(substituted,leaf)
  assert len(ww)<=2*r*r+1
  endpoint=tuple(h.conjugate(sigma,g)for g in u.replay(substituted,ww))
  assert endpoint==u.wr(r+1)
  length=len(ww)
 return {"r":r,"distance":d,"case":name,"selected_index_j":j,
         "prefix_length":len(word),"Wr_normalization_length":length,
         "outside_anchor_pairs":[sorted(a)for a in expected_pairs],
         "eligible_leaf":leaf}
def rotate_pair(r,a,b):
 t=initial(r,a,b);old=t[:r];P=h.product(old);C=h.cycles(P)[0];pos={a:i for i,a in enumerate(C)}
 dist=(pos[b]-pos[a])%(2*r+1)
 if dist>r:
  swap=h.permutation(2*r+3,[(2*r+1,2*r+2)])
  t=tuple(h.conjugate(swap,g)for g in t);a,b=b,a;dist=2*r+1-dist
 name,j,(aa,bb)=case_data(r,dist)
 k=(pos[a]-pos[aa])%(2*r+1)
 if k>r:k-=2*r+1
 twist=list(range(1,r))*r
 word=(invword(twist)if k>0 else twist)*abs(k)
 t=u.replay(t,word)
 label=power(P,-k)
 t=tuple(h.conjugate(label,g)for g in t)
 assert t==initial(r,aa,bb)
 assert len(word)<=r*r*(r-1)
 return dist,len(word)
def small_base(r,a,b):
 t=initial(r,a,b)
 target=cycle_target(t,0)
 assert h.product(target)==h.product(t)
 assert normal.lower_class(target,0)[0]=="W"
 word,sigma,push=normal.normal(target,0)
 assert tuple(h.conjugate(sigma,g)for g in u.replay(target,word))==u.wr(r+1)
 return True
def group_order(t):
 identity=tuple(range(len(t[0])));seen={identity};queue=deque([identity])
 while queue:
  g=queue.popleft()
  for a in t:
   x=h.compose(a,g)
   if x not in seen:seen.add(x);queue.append(x)
 return len(seen)

def boundary():
 rows=[];normalized={}
 for a,b in itertools.combinations(range(5),2):
  t=initial(2,a,b);target=cycle_target(t,0);todo=deque([t]);paths={t:[]}
  edges=0
  while todo:
   z=todo.popleft()
   for j in range(2):
    for sign in(1,-1):
     new=h.hurwitz(z,j,sign)
     assert h.product(new)==h.product(t)
     edges+=1
     if new not in paths:paths[new]=paths[z]+[sign*(j+1)];todo.append(new)
  assert len(paths)==147 and edges==588
  assert group_order(t)==2520
  P=h.product(t[:2]);C=next(c for c in h.cycles(P)if len(c)>1)
  dist=(C.index(b)-C.index(a))%5;dist=min(dist,5-dist)
  assert (target in paths)==(dist==2)
  negative_leaf_choices=0
  if dist==1:
   for z in paths:
    for x in range(7):
     if normal.f.owner(z,x)is not None:
      assert normal.lower_class(z,x)[0]=="S"
      negative_leaf_choices+=1
  else:
   assert u.replay(t,paths[target])==target
  C=h.product(t);q=[None]*7;x=0
  for i in range(7):q[x]=i;x=C[x]
  canonical={tuple(h.conjugate(tuple(q),g)for g in z)for z in paths}
  if dist in normalized:assert normalized[dist]==canonical
  else:normalized[dist]=canonical
  rows.append({"anchors":[a,b],"old_cycle_distance":dist,"literal_orbit_size":147,
               "signed_edges_checked":edges,"upper_group_order":2520,
               "Wr_reachable":target in paths,"literal_word_if_reachable":paths.get(target),
               "all_negative_orbit_leaf_choices_checked":negative_leaf_choices})
 assert len(normalized[1])==len(normalized[2])==147 and normalized[1].isdisjoint(normalized[2])
 return {"all_ten_extensions":rows,"two_normalized_orbits_disjoint":True,
         "normalized_orbit_hashes":{str(k):hashlib.sha256(json.dumps(sorted(v),separators=(",",":")).encode()).hexdigest()
                                    for k,v in normalized.items()}}

def imported_counts():
 receipt=json.loads((HERE/"degree_fourteen_hurwitz_verification.json").read_text())
 result=[]
 for n,V in ((9,6561),(11,1362944),(13,446265625)):
  case=next(a for a in receipt["cases"]if a["kind"]=="DT" and a["polynomial_degree"]==n)
  chars=next(a for a in receipt["character_bounds"]if a["polynomial_degree"]==n)
  assert case["complete"] and case["absolute_orbit"]==case["absolute_total_bound"]==V
  assert chars["literal_DT_count"]==n*V and chars["absolute_DT_count"]==V
  r=(n-1)//2
  bound=V-1+(n-1)*r*(r-1)
  assert bound<=446265984
  result.append({"degree":n,"absolute_total_and_connected_orbit":V,
                 "independent_literal_total":n*V,"literal_controller_bound":bound})
 assert [a["literal_controller_bound"]for a in result]==[6656,1363143,446265984]
 return result

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument("--write-receipt",type=Path,default=HERE/"uniform_two_leaf_dt_absorption_verification.json")
 args=parser.parse_args()
 cores=six_core_controls()
 normalized=[normalize_prefix(r,d,cores)for r in range(6,13)for d in range(1,r+1)]
 normalized+=[normalize_prefix(25,d,cores)for d in (1,2,10,20,21,22,23,24,25)]
 rotation=[]
 for r in range(6,10):
  for a,b in itertools.combinations(range(2*r+1),2):
   d,L=rotate_pair(r,a,b);rotation.append({"r":r,"old_pair":[a,b],"distance":d,"rotation_length":L})
 small=sum(small_base(r,a,b)for r in (3,4,5)for a,b in itertools.combinations(range(2*r+1),2))
 assert normalize_prefix(8,2,cores,wrong_leaf=0)=="S"
 bd=boundary();imports=imported_counts()
 paths=[Path(__file__).resolve(),HERE/"UNIFORM_TWO_LEAF_DT_ABSORPTION.md",
        HERE/"verify_uniform_hurwitz_attachment.py",
        HERE/"verify_d5_lobe_core_absorption.py",
        HERE/"IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md",
        HERE/"verify_imprimitive_leaf_normal_form.py",
        HERE/"imprimitive_leaf_normal_form_verification.json",
        HERE/"DEGREE_FOURTEEN_HURWITZ_ORBIT.md",
        HERE/"verify_degree_fourteen_hurwitz_orbit.py",
        HERE/"degree_fourteen_hurwitz_verification.json",
        HERE/"FIXED_LEAF_HURWITZ_OBSTRUCTION.md",
        HERE/"sources/jones-primitive-cycle-1209.5169v1.pdf"]
 for module in (normal.f,normal.m,normal.leaf):
  p=Path(module.__file__).resolve()
  if p not in paths:paths.append(p)
 result={"schema_version":1,"result":"PASS","python":platform.python_version(),"six_cores":{name:{"product":row["product"],"leaf":CASES[name]["leaf"],
          "outside_anchor_pairs":CASES[name]["anchors"]}for name,row in cores.items()},
         "normalized_case_controls":normalized,"all_pair_rotation_controls":rotation,
         "small_base_extension_cases":small,"degree7_sharp_boundary":bd,
         "imported_finite_controller_bounds":imports,
         "ordinary_move_bound":"446265984+r^3+r^2+5r+1, r>=3",
         "global_Wr_normalizer_bound":"2r^2+1, for r+1 factors",
         "wrong_generic_core_leaf0_rejected_as_symmetric":True,
         "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths},
         "scope":"Arbitrary distinct old anchors and two fresh leaf sheets over canonicalWr_r, everyr>=3. Sharp r2cyclic-distance split: all10extensions haveA7, only distance2 reachesWr. Fixed finiteN9/N11/N13controller endpoints are substituted, not unrecovered words replayed. No arbitrary-DT connectivity or next-degree full classification.",
         "formal_verification":False,"independent_external_replay":False,"external_human_review":False}
 args.write_receipt.write_text(json.dumps(result,indent=2)+"\n")
 print("PASS","normalized",len(normalized),"allpairrotations",len(rotation),"smallbase",small,"sharpboundary",len(bd["all_ten_extensions"]),flush=True)

if __name__=="__main__":main()
