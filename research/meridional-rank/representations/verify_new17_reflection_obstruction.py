#!/usr/bin/env python3
"""Independent integer sparse-polynomial replay of a reflection unit certificate."""
from pathlib import Path
from math import gcd,lcm
from copy import deepcopy
import json,hashlib,time,sys
D=Path(__file__).resolve().parent;Z=(0,0,0,0,0,0)
def const(a):return {Z:a} if a else {}
def add(a,b):
 r=a.copy()
 for e,c in b.items():
  v=r.get(e,0)+c
  if v:r[e]=v
  else:r.pop(e,None)
 return r
def scale(a,c):return {e:v*c for e,v in a.items()} if c else {}
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
 r={}
 for e,x in a.items():
  for f,y in b.items():
   g=(e[0]+f[0],e[1]+f[1],e[2]+f[2],e[3]+f[3],e[4]+f[4],e[5]+f[5])
   r[g]=r.get(g,0)+x*y
 return {e:c for e,c in r.items() if c}
def shift(a,j):
 return {tuple(x+int(i==j) for i,x in enumerate(e)):c for e,c in a.items()}
off={(0,1):0,(0,2):1,(0,3):2,(1,2):3,(1,3):4,(2,3):5}
def dot(u,v):
 w=[]
 for i in range(4):
  t=scale(v[i],2)
  for j in range(4):
   if i!=j:t=add(t,shift(v[j],off[tuple(sorted((i,j)))]))
  w.append(t)
 out={}
 for i in range(4):out=add(out,mul(u[i],w[i]))
 return out
def propagate(word,q):
 q=deepcopy(q)
 for k in word:
  i=abs(k)-1;u,v=q[i],q[i+1];t=dot(u,v)
  if k>0:q[i],q[i+1]=[sub(vj,mul(t,uj)) for uj,vj in zip(u,v)],u
  else:q[i],q[i+1]=v,[sub(uj,mul(t,vj)) for uj,vj in zip(u,v)]
 return q
def primitive(p):
 c=0
 for a in p.values():c=gcd(c,a)
 return {e:a//c for e,a in p.items()} if c else {}
def decode(a,factor=1):
 r={}
 for t in a:
  e=tuple(t["exponents"]);assert len(e)==6 and all(isinstance(x,int) and x>=0 for x in e)
  assert t["denominator"]>0 and factor%t["denominator"]==0 and e not in r
  c=t["numerator"]*(factor//t["denominator"]);assert c;r[e]=c
 return r
def evaluate(p,values):
 return sum(c*prod_int(v**k for v,k in zip(values,e)) for e,c in p.items())
def prod_int(values):
 out=1
 for x in values:out*=x
 return out
def numeric_propagate(word,q,G):
 q=deepcopy(q)
 for k in word:
  i=abs(k)-1;u,v=q[i],q[i+1];t=sum(u[a]*G[a][b]*v[b] for a in range(4) for b in range(4))
  if k>0:q[i],q[i+1]=[y-t*x for x,y in zip(u,v)],u
  else:q[i],q[i+1]=v,[x-t*y for x,y in zip(u,v)]
 return q
def identity(polys,qs):
 total={}
 for p,q in zip(polys,qs):total=add(total,mul(p,q))
 return total

HERE=Path(__file__).resolve().parent
INPUT=HERE/"new17_reflection_inputs.json"
RECEIPT=HERE/"new17_reflection_verification.json"
BETA=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3,2,2,2]
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def diagram(pd):
 ends={}
 for i,row in enumerate(pd):
  assert len(row)==len(set(row))==4
  for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
 assert len(ends)==2*len(pd) and all(len(x)==2 for x in ends.values())
 alpha={a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
 unseen=set(alpha);faces=0
 while unseen:
  start=x=min(unseen)
  while True:
   unseen.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 assert faces==len(pd)+2
 x=0;seen=set()
 for _ in range(2*len(pd)):
  assert x not in seen;seen.add(x);x=alpha[4*(x//4)+(x+2)%4]
 assert x==0 and len(seen)==2*len(pd)
 parent={e:e for e in ends}
 def find(e):
  while parent[e]!=e:parent[e]=parent[parent[e]];e=parent[e]
  return e
 for a,b,c,d in pd:parent[find(b)]=find(d)
 cc={}
 for e in ends:cc.setdefault(find(e),[]).append(e)
 groups=sorted(sorted(g) for g in cc.values());arc={e:i for i,g in enumerate(groups) for e in g}
 assert len(groups)==len(pd)
 return [(arc[b],arc[a],arc[c]) for a,b,c,d in pd],alpha,arc
def braid_projection(word):
 top=list(range(4));cur=top[:];rows=[];nxt=4
 for g in word:
  assert 1<=abs(g)<=3
  i=abs(g)-1;a,b=cur[i:i+2];left,right=nxt,nxt+1;nxt+=2
  rows.append([a,left,right,b] if g>0 else [b,a,left,right])
  cur[i:i+2]=[left,right]
 parent=list(range(nxt))
 def find(e):
  while parent[e]!=e:parent[e]=parent[parent[e]];e=parent[e]
  return e
 for a,b in zip(top,cur):parent[find(a)]=find(b)
 return [[find(e) for e in row] for row in rows],[find(e) for e in top]
def check_braid(case):
 assert case["literal_braid_word"]==BETA and case["hurwitz_word"]==[-g for g in BETA]
 pd=case["pd"];assert len(pd)==17
 crossings,alpha,arc=diagram(pd)
 raw,top=braid_projection(BETA);_,generated_alpha,_=diagram(raw)
 assert alpha==generated_alpha
 edge_map={}
 for row,saved in zip(raw,pd):
  for x,y in zip(row,saved):
   if x in edge_map:assert edge_map[x]==y
   edge_map[x]=y
 assert len(edge_map)==len(set(edge_map.values()))==34
 seeds=[arc[edge_map[e]] for e in top];assert seeds==case["seed_arcs"]==[1,0,8,13]
 known=set(seeds)
 for o,a,b in case["seed_trace"]:
  assert o in known and a in known and b not in known and ((o,a,b) in crossings or (o,b,a) in crossings);known.add(b)
 assert len(case["seed_trace"])==13 and known==set(range(17))
 # A nonconstant involutory colouring tests the literal local crossing convention.
 ps=[tuple(x) for x in case["arc_permutations"]];I=tuple(range(8))
 for p in ps:assert sorted(p)==list(range(8)) and compose(p,p)==I and sum(p[i]!=i for i in range(8))==4
 for o,a,b in crossings:assert compose(compose(ps[o],ps[a]),ps[o])==ps[b]
 initial=[ps[s] for s in seeds];current=initial[:];edge_colours={}
 for g,row in zip(BETA,pd):
  i=abs(g)-1;a,b=current[i:i+2]
  if g>0:
   c=compose(compose(b,a),b);ports=[a,b,c,b];out=[b,c]
  else:
   c=compose(compose(a,b),a);ports=[b,a,c,a];out=[c,a]
  for e,v in zip(row,ports):
   if e in edge_colours:assert edge_colours[e]==v
   edge_colours[e]=v;assert ps[arc[e]]==v
  current[i:i+2]=out
 assert current==initial and len(edge_colours)==34
 return {"literal_braid_to_PD_dart_gluings":68,"PD_sphere_faces":19,"one_component":True,
         "top_meridian_arcs":seeds,"literal_involutory_relations":17,"nonconstant_colour_convention_control":True}
def transports(word,eps):
 from itertools import product
 origin=list(range(4))
 for k in word:i=abs(k)-1;origin[i],origin[i+1]=origin[i+1],origin[i]
 orbit={0};j=origin[0]
 while j!=0:orbit.add(j);j=origin[j]
 assert len(orbit)==4
 out=[]
 for target in product((1,-1),repeat=4):
  if prod_int(target)!=prod_int(eps):continue
  choices=[d for d in product((1,-1),repeat=4) if tuple(d[origin[i]]*d[i]*eps[i] for i in range(4))==target]
  assert len(choices)==2
  out.append({"closing_signs":target,"initial_root_signs":choices[0]})
 assert len(out)==8
 return origin,out
def positive_controls(data):
 ctrl=data["positive_T45_control"];word=ctrl["word"];G=ctrl["Gram"]
 assert ctrl["literal_braid_word"]==[-k for k in word]
 assert word==[1,2,3]*5 and G==[[2 if i==j else 1 for j in range(4)] for i in range(4)]
 E=[[int(i==j) for i in range(4)] for j in range(4)];eps=[-1,-1,-1,1]
 def check(word,G,perms):
  assert numeric_propagate(word,E,G)==[[eps[j]*x for x in row] for j,row in enumerate(E)]
  assert numeric_propagate([-k for k in reversed(word)],[[eps[j]*x for x in row] for j,row in enumerate(E)],G)==E
  cur=[tuple(x) for x in perms];initial=cur[:]
  for k in word:
   i=abs(k)-1;a,b=cur[i:i+2]
   if k>0:cur[i:i+2]=[compose(compose(a,b),a),a]
   else:cur[i:i+2]=[b,compose(compose(b,a),b)]
  assert cur==initial
  _,all_signs=transports(word,eps)
  for item in all_signs:
   ds=item["initial_root_signs"];H=[[ds[i]*ds[j]*G[i][j] for j in range(4)] for i in range(4)]
   expected=[[item["closing_signs"][j]*x for x in row] for j,row in enumerate(E)]
   assert numeric_propagate(word,E,H)==expected
 star=[tuple(j if i==0 else 0 if i==j else i for i in range(5)) for j in range(1,5)]
 assert ctrl["root_final_signs"]==eps and ctrl["S5_seed_permutations"]==[list(p) for p in star] and ctrl["Gram_determinant"]==5
 check(word,G,star)
 chainG=[[2 if i==j else -1 if abs(i-j)==1 else 0 for j in range(4)] for i in range(4)]
 chain=[]
 for j in range(4):
  p=list(range(5));p[j],p[j+1]=p[j+1],p[j];chain.append(tuple(p))
 check([-k for k in word],chainG,chain)
 return {"actual_T45_star_frame":True,"mirror_T45_chain_frame":True,
         "inverse_root_moves":True,"positive_gauge_specializations":16}
def preserve_receipt(path,out):
 # Reviews retain their own source/PDF identities; replay does not reapprove a new edition.
 for attempt in range(5):
  prior=path.read_bytes() if path.exists() else None
  previous=json.loads(prior) if prior is not None else {}
  merged=dict(out)
  if "internal_pdf_reviews" in previous:merged["internal_pdf_reviews"]=previous["internal_pdf_reviews"]
  temp=path.with_name(path.name+".tmp");temp.write_text(json.dumps(merged,indent=2)+"\n")
  if (path.read_bytes() if path.exists() else None)==prior:temp.replace(path);return
 raise RuntimeError("Concurrent receipt changes; review metadata was not overwritten")
def main():
 start=time.monotonic();data=json.loads(INPUT.read_text());case=data["case"]
 assert case["symbol"]=="K_*" and data["field_scope"]=="Characteristic-zero fields only in the written theorem."
 binding=check_braid(case);positive=positive_controls(data)
 cp=HERE/data["certificate"]["path"];raw=cp.read_bytes();cert=json.loads(raw)
 assert hashlib.sha256(raw).hexdigest()==data["certificate"]["sha256"] and len(raw)==data["certificate"]["bytes"]
 word=cert["word"];cut=cert["cut"];eps=cert["signs"]
 assert word==case["hurwitz_word"] and cut==data["certificate"]["midpoint_cut"]==7
 assert eps==data["certificate"]["base_signs"]==[-1,1,1,1]
 assert cert["variables"]==["a","b","c","d","e","f"] and cert["determinant_localized"] is False
 assert prod_int(eps)==(-1)**len(word)
 seeds=[[const(int(i==j)) for i in range(4)] for j in range(4)]
 left=propagate(word[:cut],seeds)
 right=propagate([-k for k in reversed(word[cut:])],[[scale(x,eps[j]) for x in u] for j,u in enumerate(seeds)])
 polys=[primitive(sub(x,y)) for u,v in zip(left,right) for x,y in zip(u,v) if x!=y]
 assert len(polys)==16 and polys==[decode(a) for a in cert["polynomials"]]
 denom=lcm(*(t["denominator"] for a in cert["multipliers"] for t in a))
 assert denom==data["certificate"]["denominator_clearing_factor"]==16
 assert data["certificate"]["identity_constant"]==16 and data["certificate"]["rational_unit_identity_constant"]==1
 qs=[decode(a,denom) for a in cert["multipliers"]];assert len(qs)==16 and sum(map(len,qs))==75136
 total=identity(polys,qs);assert total==const(16)
 # Change a genuine nonzero coefficient, reusing only unaffected verified summands.
 i=min((j for j,q in enumerate(qs) if q),key=lambda j:len(qs[j])*len(polys[j]))
 bad=qs[i].copy();e=next(iter(bad));bad[e]+=1 if bad[e]!=-1 else 2
 altered=add(sub(total,mul(polys[i],qs[i])),mul(polys[i],bad));assert altered!=const(16)
 rejected=["altered_multiplier_coefficient"]
 badcase=deepcopy(case);badcase["literal_braid_word"][0]*=-1
 try:check_braid(badcase)
 except AssertionError:rejected.append("braid_sign_mismatch")
 else:raise AssertionError("Corrupt braid sign accepted")
 badcase=deepcopy(case);badcase["pd"][0][0]=badcase["pd"][0][1]
 try:check_braid(badcase)
 except AssertionError:rejected.append("PD_edge_gluing_corruption")
 else:raise AssertionError("Corrupt PD accepted")
 badterm=deepcopy(cert["multipliers"][i]);badterm[0]["exponents"][0]=-1
 try:decode(badterm,denom)
 except AssertionError:rejected.append("negative_exponent")
 else:raise AssertionError("Nonpolynomial multiplier accepted")
 origin,signs=transports(word,eps)
 assert sorted(tuple(x["closing_signs"]) for x in signs)==sorted(map(tuple,data["closing_sign_patterns"]))
 G=data["positive_T45_control"]["Gram"]
 values=[G[i][j] for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]]
 E=[[int(i==j) for i in range(4)] for j in range(4)]
 numeric_left=numeric_propagate(word[:cut],E,G)
 numeric_right=numeric_propagate([-k for k in reversed(word[cut:])],[[eps[j]*x for x in row] for j,row in enumerate(E)],G)
 assert [[evaluate(x,values) for x in row] for row in left]==numeric_left
 assert [[evaluate(x,values) for x in row] for row in right]==numeric_right
 files=[INPUT,cp,Path(__file__),HERE/"UNIVERSAL_17_REFLECTION_OBSTRUCTION.tex"]
 out={"status":"PASS exact characteristic-zero symmetric-reflection obstruction for K_*",
      "scope":data["description"],"braid_binding":binding,"Hurwitz_word":word,"midpoint_cut":cut,
      "closing_sign_product":-1,"strand_origins":origin,"closing_sign_patterns_covered":8,"sign_transports":signs,
      "constraint_count":len(polys),"constraint_terms":sum(map(len,polys)),
      "nonzero_multiplier_count":sum(bool(q) for q in qs),"integer_multiplier_terms":sum(map(len,qs)),
      "literal_polynomial_identity_constant":16,"rational_unit_identity_constant":1,"common_denominator":16,"Gram_determinant_localized":False,
      "positive_controls":positive,"numeric_polynomial_specialization":True,"corruptions_rejected":rejected,
      "elapsed_seconds":round(time.monotonic()-start,3),"python_version":sys.version.split()[0],
      "input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
      "assurance":"Exact internal standard-library polynomial and literal-braid replay. The sign-gauge, reflection and Coxeter implications are written proofs; no all-characteristics scope, augmentation gap, formal verification or external review is asserted."}
 preserve_receipt(RECEIPT,out);print(json.dumps(out,indent=2),flush=True)
if __name__=="__main__":main()
