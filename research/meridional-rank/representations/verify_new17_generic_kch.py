#!/usr/bin/env python3
"""Exact signed17 folds,392-only core transfer and ordinary3D KCH witness.
Default: Python standard library. Optional --topology: fresh explicit-braid
Spherogram PD/DT and its actual _pieces/knot_group sign convention.
The old392 core contradiction is an imported written premise, not re-proved.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product,permutations,combinations
from copy import deepcopy
import argparse,json,hashlib,sys,time
HERE=Path(__file__).resolve().parent
INPUT=HERE/"new17_generic_kch_inputs.json"
RECEIPT=HERE/"new17_generic_kch_verification.json"
BETA=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3,2,2,2]
class Ring:
 def __init__(self,names):self.names=tuple(names);self.zero=(0,)*len(names)
 def c(self,n):return {self.zero:n} if n else {}
 def v(self,name,power=1):
  e=[0]*len(self.names);e[self.names.index(name)]=power;return {tuple(e):1}
 def add(self,*ps):
  out={}
  for p in ps:
   for e,c in p.items():out[e]=out.get(e,0)+c
  return {e:c for e,c in out.items() if c}
 def neg(self,p):return {e:-c for e,c in p.items()}
 def sub(self,p,q):return self.add(p,self.neg(q))
 def mul(self,*ps):
  out=self.c(1)
  for p in ps:
   q={}
   for e,a in out.items():
    for f,b in p.items():
     g=tuple(x+y for x,y in zip(e,f));q[g]=q.get(g,0)+a*b
   out={e:c for e,c in q.items() if c}
  return out
def pair_engine(R,B):
 add,sub,mul=R.add,R.sub,R.mul;mi=R.v("mu",-1)
 def dot(f,u):return add(*(mul(f[i],B[i][j],u[j]) for i in range(len(B)) for j in range(len(B))))
 def conjugate(o,a,e):
  u,f=o;v,g=a;x=dot(f,v);y=dot(g,u)
  if e==1:return ([sub(v[i],mul(x,u[i])) for i in range(len(u))],[add(g[i],mul(mi,y,f[i])) for i in range(len(f))])
  return ([add(v[i],mul(mi,x,u[i])) for i in range(len(u))],[sub(g[i],mul(y,f[i])) for i in range(len(f))])
 return dot,conjugate
def outer(R,p):return [[R.mul(x,y) for y in p[1]] for x in p[0]]
def reduce_xy(R,p):
 im=R.names.index("mu");ix=R.names.index("x");iy=R.names.index("y");out={}
 for ex,c in p.items():
  e=list(ex);assert e[ix]>=0 and e[iy]>=0
  n=min(e[ix],e[iy]);e[ix]-=n;e[iy]-=n;e[im]+=n
  e=tuple(e);out[e]=out.get(e,0)+((-1)**n)*c
 return {e:c for e,c in out.items() if c}
def free_reduce(word):
 out=[]
 for x in word:
  if out and out[-1]==-x:out.pop()
  else:out.append(x)
 return out
def inverse_word(w):return [-x for x in reversed(w)]
def fold_algebra():
 # First fold: a nonzero root/covector comparison gives
 # t=-lambda,s=mu/lambda,x=k*a,xi=f/k. Necessity/nonzero
 # deductions are written; the exact normalized formulas are checked here.
 R=Ring(("mu","lambda","k"));c,mul,sub=R.c,R.mul,R.sub
 mu,mi=R.v("mu"),R.v("mu",-1);lam,li=R.v("lambda"),R.v("lambda",-1);k,ki=R.v("k"),R.v("k",-1)
 d=sub(c(1),mu);t=R.neg(lam);s=mul(mu,li)
 B=[[d,s],[t,d]];A=([c(1),c(0)],[c(1),c(0)]);D=([c(0),c(1)],[c(0),c(1)])
 X=([k,c(0)],[ki,c(0)])
 dot,conj=pair_engine(R,B)
 left=conj(D,A,1);right=conj(X,D,-1)
 assert left[0]==[mul(lam,z) for z in right[0]]
 assert left[1]==[mul(li,z) for z in right[1]]
 assert outer(R,X)==outer(R,A) and mul(t,s)==R.neg(mu)
 assert dot(left[1],left[0])==d and dot(right[1],right[0])==d
 # Arc1=A^-1 D A, and arc11=arc1^-1 A arc1.
 arc1=[-1,2,1];arc11=free_reduce(inverse_word(arc1)+[1]+arc1)
 assert arc11==[-1,-2,1,2,1]
 assert free_reduce(arc11[:2]+[2,1,2])==[2] # replace ADA by DAD
 # Three negative Hurwitz moves on (D,E), normalized by d=1-mu.
 R=Ring(("mu","x","y"));c,mul,add,sub=R.c,R.mul,R.add,R.sub
 mu,mi,x,y=R.v("mu"),R.v("mu",-1),R.v("x"),R.v("y");d=sub(c(1),mu);h=add(c(1),mul(x,y,mi))
 B=[[d,y],[x,d]];D=([c(1),c(0)],[c(1),c(0)]);E=([c(0),c(1)],[c(0),c(1)])
 dot,conj=pair_engine(R,B)
 N=conj(E,D,-1);Q=conj(N,E,-1);C=conj(Q,N,-1)
 assert Q==([y,h],[R.neg(mul(x,mi)),h])
 assert C==([h,mul(x,mi,add(c(1),h))],[h,R.neg(mul(y,add(c(1),h)))])
 for pair in [N,Q,C]:assert dot(pair[1],pair[0])==d
 assert not reduce_xy(R,h)
 for p,q in [(Q,D),(C,E)]:
  for row,expected in zip(outer(R,p),outer(R,q)):
   for a,b in zip(row,expected):assert not reduce_xy(R,sub(a,b))
 # Third fold is group-theoretic product preservation:
 # N=Y^-1 D Y; D=N^-1 Y N; B=D^-1 N D.
 assert free_reduce([2]+[-2,1,2])==[1,2]      # YN=DY
 assert free_reduce([3]+[-3,2,3])==[2,3]      # ND=YN
 assert free_reduce([1]+[-1,3,1])==[3,1]      # DB=ND
 return {"first_fold_normalized_pair_coefficients":True,"first_fold_rankone_identity_X_equals_A":True,
         "arc11_group_braid_consequence":True,"negative_Hurwitz_three_raw_pair_formulas":True,
         "second_fold_outer_products_mod_xy_plus_mu":True,"third_fold_group_product_identities":3,
         "allowed_formal_units":["mu","lambda","k"],"inverse_one_minus_mu":False,
         "evaluation_matrix_inverted":False,"nonzero_and_independence_deductions":"Written lem:folds; not inferred from finite sampling."}
def relation_key(r):
 o,a,b,e=r;return (o,min(a,b),max(a,b),e if a<b else -e)
def actual_step(rows,r):
 o,a,b,e=r
 found=[t for t in rows if t[0]==o and set(t[1:3])=={a,b}]
 assert len(found)==1 and relation_key(found[0])==relation_key(r)
 return r
def diagram(pd):
 ends={}
 for j,row in enumerate(pd):
  assert len(row)==len(set(row))==4
  for h,e in enumerate(row):ends.setdefault(e,[]).append(4*j+h)
 assert all(len(v)==2 for v in ends.values()) and len(ends)==2*len(pd)
 alpha={a:b for v in ends.values() for a,b in (v,v[::-1])}
 remaining=set(alpha);faces=0
 while remaining:
  start=x=min(remaining)
  while True:
   remaining.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
   if x==start:break
  faces+=1
 assert faces==len(pd)+2
 incoming=set();x=0
 for _ in range(2*len(pd)):
  assert x not in incoming;incoming.add(x);x=alpha[4*(x//4)+(x+2)%4]
 assert x==0
 parent={e:e for e in ends}
 def find(e):
  while parent[e]!=e:parent[e]=parent[parent[e]];e=parent[e]
  return e
 for a,b,c,d in pd:parent[find(b)]=find(d)
 cc={}
 for e in ends:cc.setdefault(find(e),[]).append(e)
 groups=sorted(sorted(v) for v in cc.values());assert len(groups)==len(pd)
 arc={e:i for i,g in enumerate(groups) for e in g}
 return alpha,arc,incoming
def braid_projection(word):
 top=list(range(4));cur=top[:];raw=[];nxt=4
 for x in word:
  j=abs(x)-1;a,b=cur[j:j+2];u,v=nxt,nxt+1;nxt+=2
  raw.append([a,u,v,b] if x>0 else [b,a,u,v]);cur[j:j+2]=[u,v]
 parent=list(range(nxt))
 def f(e):
  while parent[e]!=e:e=parent[e]
  return e
 for a,b in zip(top,cur):parent[f(a)]=f(b)
 return [[f(e) for e in row] for row in raw],[f(e) for e in top]
def generation(rows,seeds,trace):
 known=set(seeds)
 for r in trace:
  actual_step(rows,r);o,a,b,e=r
  assert o in known and a in known and b not in known;known.add(b)
 return known
def case_binding(case):
 assert case["literal_braid_word"]==BETA and case["hurwitz_word"]==[-x for x in BETA]
 pd=case["pd"];assert len(pd)==17
 alpha,arc,incoming=diagram(pd)
 raw,top=braid_projection(BETA);a,_,_=diagram(raw);assert a==alpha
 segment_map={}
 for r,s in zip(raw,pd):
  for x,y in zip(r,s):
   if x in segment_map:assert segment_map[x]==y
   segment_map[x]=y
 assert len(segment_map)==len(set(segment_map.values()))==34
 assert [arc[segment_map[e]] for e in top]==case["top_seed_arcs"]==[1,0,8,13]
 assert set(case["independent_incoming_darts"])==incoming
 signed=[];ports=[]
 for j,(a,b,c,d) in enumerate(pd):
  under=next(h for h in (0,2) if 4*j+h in incoming)
  over=next(h for h in (1,3) if 4*j+h in incoming)
  assert under==0
  e=1 if over==1 else -1
  signed.append([arc[b],arc[a],arc[c],e]);ports.append([under,over])
 assert signed==case["canonical_signed_relations"]==case["native_signed_relations"]
 assert ports==case["independent_incoming_ports"]
 assert case["source_to_oriented_darts"]==list(range(68))
 assert case["seed_arcs"]==[2,8,13,0]
 assert generation(signed,case["seed_arcs"],case["signed_seed_trace"])==set(range(17))
 assert len(case["signed_seed_trace"])==13
 early=case["three_seed_trace"]
 assert generation(signed,case["seed_arcs"][:3],early)=={2,8,13}|{r[2] for r in early}
 assert len(early)==8 and {5,10,16}<={2,8,13}|{r[2] for r in early}
 return {"actual_signed_PD_rows":17,"consistent_segment_arc_labels":34,"literal_braid_dart_gluings":68,
         "four_meridian_generating_trace_steps":13,"three_seed_span_trace_steps":8,
         "signed_native_rows_independent_dart_walk":True}
def core_transfer(data):
 case=data["case"];rows=case["canonical_signed_relations"];fold=data["folds"]
 for key in ["first_sources","second_sources","third_sources"]:
  for r in fold[key]:actual_step(rows,r)
 assert fold["first_sources"]==[[0,2,1,1],[10,0,1,-1],[1,10,11,-1]]
 assert fold["second_sources"]==[[5,11,12,-1],[12,5,6,-1],[6,12,13,-1]]
 assert fold["third_sources"]==[[16,6,7,-1],[7,16,0,-1],[0,7,8,-1]]
 assert fold["identifications"]=={"10":2,"11":0,"5":13,"6":0,"16":8}
 aliases={int(k):v for k,v in fold["identifications"].items()}
 transform=lambda r:[aliases.get(r[0],r[0]),aliases.get(r[1],r[1]),aliases.get(r[2],r[2]),r[3]]
 core=[transform(rows[i]) for i in fold["core_source_indices"]]
 assert core==fold["necessary_core"] and len(core)==8
 transfer=data["core_transfer"];premise=transfer["premise"]
 assert transfer["global_sign"]==1 and premise["label"]=="lem:cores" and premise["subcase"]=="392"
 old=json.loads((HERE/premise["input_file"]).read_text())
 oldcase=next(c for c in old["cases"] if c["name"]=="13a_392")
 assert oldcase["canonical_signed_relations"]==premise["original_signed_relations"]
 assert "\\label{lem:cores}" in (HERE/premise["proof_file"]).read_text()
 oa={int(k):v for k,v in premise["old_arc_aliases"].items()}
 oldcore=[]
 for o,a,b,e in premise["original_signed_relations"]:
  if {o,a,b}&set(premise["deleted_nodes"]):continue
  oldcore.append([oa.get(o,o),oa.get(a,a),oa.get(b,b),e])
 mapping={int(k):v for k,v in transfer["new_to_old"].items()}
 assert set(mapping)=={2,3,4,8,9,13,14,15} and len(set(mapping.values()))==8
 mapped=[[mapping[o],mapping[a],mapping[b],e] for o,a,b,e in core]
 assert {relation_key(r) for r in mapped}=={relation_key(r) for r in oldcore}
 newvertices={x for r in core for x in r[:3]};oldvertices={x for r in oldcore for x in r[:3]}
 assert generation(core,[2,8,13],core[:5])==newvertices
 assert generation(oldcore,premise["target_seed_arcs"],premise["target_generation_trace"])==oldvertices
 assert {mapping[i] for i in [2,8,13]}<=oldvertices
 return {"fold_sources_bound_to_actual_signed_crossings":9,"necessary_new_core_rows":8,"old392_signed_core_bijection":True,
         "global_sign":1,"old_standard_seed_arcs":premise["target_seed_arcs"],"old_seed_generation_steps":5,
         "root_and_covector_basis_transfer":"Both three-seed systems generate all8 pairs; separate span-dimension argument in lem:coretransfer.",
         "imported_written_premise":{"file":premise["proof_file"],"label":"lem:cores","subcase":"392"}}
# Q(phi), with phi^2=phi+1. All arithmetic is exact rational arithmetic.
QZERO=(F(0),F(0));QONE=(F(1),F(0))
def q(a,b=0):return (F(a),F(b))
def qadd(*xs):return (sum((x[0] for x in xs),F(0)),sum((x[1] for x in xs),F(0)))
def qneg(x):return (-x[0],-x[1])
def qsub(a,b):return qadd(a,qneg(b))
def qmul(a,b):return (a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])
def qscale(c,a):return (c*a[0],c*a[1])
def qdot(a,b):return qadd(*(qmul(x,y) for x,y in zip(a,b)))
def qmm(a,b):return tuple(tuple(qadd(*(qmul(a[i][k],b[k][j]) for k in range(len(b)))) for j in range(len(b[0]))) for i in range(len(a)))
def qmv(a,v):return tuple(qdot(row,v) for row in a)
def eye(n):return tuple(tuple(q(int(i==j)) for j in range(n)) for i in range(n))
def reflection(u):
 assert qdot(u,u)==q(4)
 return tuple(tuple(qsub(q(int(i==j)),qscale(F(1,2),qmul(u[i],u[j]))) for j in range(3)) for i in range(3))
def determinant(a):
 out=QZERO;n=len(a)
 for p in permutations(range(n)):
  z=q((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n)))
  for i in range(n):z=qmul(z,a[i][p[i]])
  out=qadd(out,z)
 return out
def perm_comp(a,b):return tuple(a[b[i]] for i in range(len(a)))
def perm_group(gens):
 I=tuple(range(len(gens[0])));seen={I};queue=[I]
 for a in queue:
  for b in gens:
   c=perm_comp(a,b)
   if c not in seen:seen.add(c);queue.append(c)
 return seen
def positive_witness(data):
 w=data["three_dimensional_witness"];case=data["case"]
 assert w["field"]=="Q(phi), phi^2=phi+1" and w["root_norm"]==4
 roots=[tuple(q(*z) for z in r) for r in w["arc_roots"]];assert len(roots)==17
 lines=[tuple(q(*z) for z in r) for r in w["root_line_representatives"]];assert len(lines)==15
 R30=set()
 for i in range(3):
  for s in (-1,1):R30.add(tuple(q(2*s*int(j==i)) for j in range(3)))
 base=(q(1),q(-1,1),q(0,1))
 for signs in product((-1,1),repeat=3):
  t=tuple(qscale(signs[i],base[i]) for i in range(3))
  for j in range(3):R30.add(t[j:]+t[:j])
 signed_lines=set(lines)|{tuple(qneg(x) for x in u) for u in lines}
 assert len(R30)==30 and signed_lines==R30 and all(u in R30 for u in roots)
 assert all(qdot(u,u)==q(4) for u in R30)
 mats=[reflection(u) for u in roots];E=eye(3)
 for M in mats:assert qmm(M,M)==E
 for o,a,b,e in case["canonical_signed_relations"]:
  assert e in (-1,1) and qmm(qmm(mats[o],mats[a]),mats[o])==mats[b]
 seeds=w["seed_arcs"];assert seeds==[1,0,8,13]
 assert roots[seeds[2]]==roots[seeds[3]]
 basis=[roots[i] for i in seeds[:3]]
 det=determinant(tuple(tuple(basis[j][i] for j in range(3)) for i in range(3)));assert det==q(-4)
 neighbours={i:{j for j in range(3) if qdot(basis[i],basis[j])!=QZERO} for i in range(3)}
 seen={0}
 while True:
  enlarged=seen|{j for i in seen for j in neighbours[i]}
  if enlarged==seen:break
  seen=enlarged
 assert seen=={0,1,2}
 order=sorted(R30);idx={u:i for i,u in enumerate(order)}
 generators=[mats[i] for i in seeds[:3]]
 actions=[]
 for M in generators:
  perm=tuple(idx[qmv(M,u)] for u in order);assert sorted(perm)==list(range(30));actions.append(perm)
 group=perm_group(actions);assert len(group)==w["expected_image_order"]==120
 class_mats=[reflection(u) for u in lines]
 class_actions=[tuple(idx[qmv(M,u)] for u in order) for M in class_mats]
 assert len(set(class_actions))==15
 orbit={actions[0]};todo=[actions[0]]
 for a in todo:
  for g in actions:
   b=perm_comp(perm_comp(g,a),g)
   if b not in orbit:orbit.add(b);todo.append(b)
 assert orbit==set(class_actions)
 pair_orders=[len(perm_group([a,b])) for a,b in combinations(class_actions,2)]
 assert max(pair_orders)<120
 return {"literal_signed_relations":17,"root_norm":4,"seed_root_minor":"-4",
         "connected_nonorthogonality_graph":True,"complex_irreducible_dimension":3,"meridian_parameter":-1,
         "explicit_signed_root_system_size":30,"faithful_signed_root_action":True,"finite_image_order":120,
         "meridian_reflection_class_size":15,"attaining_class_generators":3,"class_pairs_checked":len(pair_orders),
         "largest_two_reflection_subgroup":max(pair_orders),"class_generating_rank":3,
         "ordinary_KCH_witness_only":True}
def trefoil_control(data):
 A,B=[tuple(tuple(q(x) for x in row) for row in m) for m in data["positive_trefoil_control"]["matrices"]];E=eye(2)
 assert qmm(A,A)==qmm(B,B)==E and qmm(qmm(A,B),A)==qmm(qmm(B,A),B)
 AB=qmm(A,B);vectors=[tuple(x for row in M for x in row) for M in [E,A,B,AB]]
 assert determinant(vectors)!=QZERO
 return {"meridian_parameter":-1,"dimension":2,"braid_relation":True,"generated_matrix_algebra_rank":4}
def topology_mode(data):
 import importlib.metadata as metadata
 from spherogram import Link
 assert metadata.version("spherogram")==data["optional_topology"]["spherogram_version"]
 case=data["case"];L=Link(braid_closure=case["literal_braid_word"])
 assert [list(r) for r in L.PD_code()]==case["pd"] and json.loads(json.dumps(L.DT_code()))==case["DT"]
 L=Link(case["pd"]);assert [list(r) for r in L.PD_code()]==case["pd"]
 _,arc,_=diagram(case["pd"]);index={c:i for i,c in enumerate(L.crossings)}
 pieces=L._pieces();pa=[arc[case["pd"][index[p[0][0]]][p[0][1]]] for p in pieces]
 assert sorted(pa)==list(range(17))
 native=[]
 for c in L.crossings:
  ind={}
  for j,piece in enumerate(pieces):
   for t,item in enumerate(piece):
    if item[0]==c:
     key="j" if t==0 else "i" if t==len(piece)-1 else "k";ind[key]=pa[j]
  assert set(ind)=={"i","j","k"};native.append([ind["k"],ind["i"],ind["j"],-c.sign])
 assert native==case["native_signed_relations"]
 return {"fresh_literal_four_braid_PD_DT":True,"actual_pieces_signed_rows":17,"spherogram_version":metadata.version("spherogram")}
def corruptions(data):
 rejected=[]
 bad=deepcopy(data);bad["case"]["canonical_signed_relations"][0][3]*=-1
 try:case_binding(bad["case"])
 except AssertionError:rejected.append("actual_crossing_sign")
 else:raise AssertionError("Corrupt crossing sign accepted")
 bad=deepcopy(data);bad["core_transfer"]["new_to_old"]["2"],bad["core_transfer"]["new_to_old"]["3"]=bad["core_transfer"]["new_to_old"]["3"],bad["core_transfer"]["new_to_old"]["2"]
 try:core_transfer(bad)
 except AssertionError:rejected.append("392_core_map")
 else:raise AssertionError("Corrupt signed core map accepted")
 bad=deepcopy(data);bad["folds"]["second_sources"][0][3]*=-1
 try:core_transfer(bad)
 except AssertionError:rejected.append("local_fold_sign")
 else:raise AssertionError("Corrupt fold sign accepted")
 bad=deepcopy(data);bad["three_dimensional_witness"]["arc_roots"][0][0][0]+=1
 try:positive_witness(bad)
 except (AssertionError,KeyError):rejected.append("positive_witness_root")
 else:raise AssertionError("Corrupt positive root accepted")
 return rejected
def save_receipt(path,out):
 for attempt in range(5):
  prior=path.read_bytes() if path.exists() else None;previous=json.loads(prior) if prior is not None else {}
  merged=dict(out)
  if "internal_pdf_reviews" in previous:merged["internal_pdf_reviews"]=previous["internal_pdf_reviews"]
  temp=path.with_name(path.name+".tmp");temp.write_text(json.dumps(merged,indent=2)+"\n")
  if (path.read_bytes() if path.exists() else None)==prior:temp.replace(path);return
 raise RuntimeError("Concurrent receipt change; no review metadata overwritten")
def main():
 parser=argparse.ArgumentParser();parser.add_argument("--topology",action="store_true");args=parser.parse_args()
 start=time.monotonic();data=json.loads(INPUT.read_text())
 for name,pin in data["source_sha256"].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==pin,name
 binding=case_binding(data["case"]);folds=fold_algebra();core=core_transfer(data)
 witness=positive_witness(data);trefoil=trefoil_control(data);bad=corruptions(data)
 topo=topology_mode(data) if args.topology else None
 files=[INPUT,Path(__file__),HERE/"NEW_SEVENTEEN_GENERIC_KCH.tex"]
 out={"status":"PASS exact signed folds,392-only core transfer and ordinary3D KCH witness",
      "case_binding":binding,"fold_algebra":folds,"core_transfer":core,"positive_3D_witness":witness,
      "positive_trefoil":trefoil,"corruptions_rejected":bad,"optional_topology":topo,
      "representation_scope":data["representation_scope"],
      "input_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}|data["source_sha256"],
      "elapsed_seconds":round(time.monotonic()-start,3),"python_version":sys.version.split()[0],
      "assurance":"Exact internal stdlib local Laurent/group identities, actual signed PD and core/basis transfer, and Q(phi) positive matrices. Nonzero cancellations, arbitrary-ambient span reasoning, old392 contradiction, irreducibility and published augmentation/sum correspondences remain written premises. No reflective-cord assertion, formal verification, external specialist review or worldwide priority claim."}
 save_receipt(RECEIPT,out);print(json.dumps(out,indent=2),flush=True)
if __name__=="__main__":main()
