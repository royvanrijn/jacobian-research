#!/usr/bin/env python3
"""Independent cube-locus classification and fixed-span89-form optimality audit."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
import hashlib,json,time
from verify import add,mul,variable
from verify_power_single_compress import decode,echelon,reduce,lin,scale
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows_rref(M,width):
 A=[[Q(c) for c in row] for row in M];assert all(len(row)==width for row in A);p=0
 for j in range(width):
  k=next((k for k in range(p,len(A)) if A[k][j]),None)
  if k is None:continue
  A[p],A[k]=A[k],A[p];c=A[p][j];A[p]=[a/c for a in A[p]]
  for k in range(len(A)):
   if k!=p and A[k][j]:
    c=A[k][j];A[k]=[a-c*b for a,b in zip(A[k],A[p])]
  p+=1
  if p==len(A):break
 return A[:p]
def columns(B):return list(map(list,zip(*B))) if B and B[0] else []
def cspace(B):return tuple(tuple(row) for row in rows_rref(columns(B),len(B)))
def mm(A,B):
 width=len(B[0]) if B else 0
 return [[sum((a*B[k][j] for k,a in enumerate(row)),Q(0)) for j in range(width)] for row in A]
def linear_kernel(row):
 d=len(row);pivot=next(i for i,c in enumerate(row) if c);free=[i for i in range(d) if i!=pivot]
 K=[[Q(0) for _ in free] for _ in row]
 for j,i in enumerate(free):K[i][j]=1;K[pivot][j]=-row[i]/row[pivot]
 return K

def restrict(f,B):
 d=len(B[0]);one={(0,)*d:Q(1)};linear=[lin(row) for row in B]
 powers=[[one,l,mul(l,l),mul(mul(l,l),l)] for l in linear];result={}
 for e,c in f.items():
  term=scale(one,c)
  for i,a in enumerate(e):term=mul(term,powers[i][a])
  result=add(result,term)
 return result

def main():
 started=time.monotonic();out=HERE/"power-span-split";source=HERE/"power-basis-anneal/candidate.json";d=json.loads(source.read_text());n=d["size"];assert n==89
 validation=json.loads((HERE/"power89/validation.json").read_text());assert validation["status"]=="pass" and validation["candidate_sha256"]==sha(source)
 assert validation["source_sha256"]==sha(HERE/"verify_power89.py")
 forms=d["forms"];B=[[Q(c) for c in row] for row in d["B"]];lins=[lin(f) for f in forms];cubes=[mul(mul(l,l),l) for l in lins]
 basis={}
 for j,cube in enumerate(cubes):
  p=dict(cube);coord={j:Q(1)}
  for e,(a,ac) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,a,-c);coord=add(coord,ac,-c)
  assert p
  e=min(p);c=p[e];basis[e]=(scale(p,1/c),scale(coord,1/c))
 def project(p):
  p=dict(p);coord={}
  for e,(a,ac) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,a,-c);coord=add(coord,ac,c)
  return p,coord
 for j,cube in enumerate(cubes):assert project(cube)==({}, {j:Q(1)})
 primary=json.loads((out/"checkpoint.json").read_text());assert primary["status"]=="primary-trees-complete" and primary["source_sha256"]==sha(HERE/"power_span_split.py")
 cases=[];all_leaves=set();total_nodes=0;factor_checks=0
 for i in range(15):
  cp=HERE/"power-single-compress"/("case-%02d.json"%i);case=json.loads(cp.read_text());S=case["support"];m=len(S);constraints=[decode(f,m) for f in case["constraints"]]
  # Independently reconstruct the full necessary membership-equation space.
  residual_polys={}
  for ids in combinations_with_replacement(range(m),3):
   e=[0]*18;u=[0]*m
   for j in ids:e[S[j]]+=1;u[j]+=1
   mult=1 if ids[0]==ids[2] else (3 if ids[0]==ids[1] or ids[1]==ids[2] else 6)
   residual,coordinate=project({tuple(e):Q(mult)})
   for exponent,c in residual.items():residual_polys.setdefault(exponent,{})[tuple(u)]=c
  I=echelon(constraints);J=echelon(residual_polys.values());assert len(I)==len(J)
  assert all(not reduce(f,I) for f in residual_polys.values()) and all(not reduce(f,J) for f in constraints)
  path=out/("case-%02d.json"%i);tree=json.loads(path.read_text());assert tree["case_sha256"]==sha(cp) and tree["source_sha256"]==sha(HERE/"power_span_split.py")
  assert primary["cases"][i]["tree_sha256"]==sha(path)
  nodes=tree["nodes"];visited=set();leaves=[]
  def visit(index):
   nonlocal factor_checks
   if index in visited:return
   visited.add(index);node=nodes[index];M=[[Q(c) for c in row] for row in node["basis"]];assert len(M)==m and all(len(row)==len(M[0]) for row in M)
   dim=len(M[0]);assert len(cspace(M))==dim
   restricted=[restrict(f,M) for f in constraints]
   if node["kind"]=="leaf":
    assert not any(restricted)
    globalrows=[]
    for row in columns(M):
     v=[Q(0)]*18
     for j,c in zip(S,row):v[j]=c
     globalrows.append(v)
    component=tuple(tuple(row) for row in rows_rref(globalrows,18))
    if component:all_leaves.add(component);leaves.append(component)
    return
   assert node["kind"]=="split" and dim>0
   selected={}
   for j,c in node["combination"]:selected=add(selected,restricted[j],Q(c))
   assert selected
   product={(0,)*dim:Q(node["scalar"])};assert product[(0,)*dim]
   assert len(node["factors"])==len(node["children"])>0
   for factor,child in zip(node["factors"],node["children"]):
    row=[Q(c) for c in factor["row"]];assert len(row)==dim and any(row) and type(factor["power"]) is int and factor["power"]>0
    f=lin(row)
    for k in range(factor["power"]):product=mul(product,f)
    childB=[[Q(c) for c in row] for row in nodes[child]["basis"]]
    assert len(childB[0])==dim-1 and cspace(childB)==cspace(mm(M,linear_kernel(row)))
    visit(child)
   assert selected==product;factor_checks+=1
  root=tree["root"];rootB=[[Q(c) for c in row] for row in nodes[root]["basis"]]
  assert rootB==[[Q(i==j) for j in range(m)] for i in range(m)]
  visit(root);assert len(visited)==len(nodes);total_nodes+=len(nodes)
  cases.append(dict(case=i,nodes=len(nodes),tree_sha256=sha(path)))
 # Independent exhaustive clique cover, rather than trusting source-support labels.
 adj=[0]*18
 for f in forms:
  for i,j in combinations([i for i,c in enumerate(f) if c],2):adj[i]|=1<<j;adj[j]|=1<<i
 masks=[]
 for i in range(15):
  case=json.loads((HERE/"power-single-compress"/("case-%02d.json"%i)).read_text());masks.append(sum(1<<j for j in case["support"]))
 clique=bytearray(1<<18);clique[0]=1
 for mask in range(1,1<<18):
  bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
  if clique[rest] and rest&~adj[i]==0:clique[mask]=1;assert any(mask&~c==0 for c in masks)
 components=[]
 for C in sorted(all_leaves):
  if any(len(N)>len(C) and len(rows_rref(list(C)+list(N),18))==len(N) for N in all_leaves):continue
  components.append(C)
 assert len(components)==73 and sum(len(C)==1 for C in components)==65 and sum(len(C)==2 for C in components)==8
 points=[C[0] for C in components if len(C)==1];lines=[C for C in components if len(C)==2]
 def unit(i):return [Q(j==i) for j in range(18)]
 a=unit(1);b=unit(3);tau=unit(17)
 special=[rows_rref([b,tau],18),rows_rref([a,[x+y for x,y in zip(b,tau)]],18),rows_rref([a,[x-y for x,y in zip(b,tau)]],18)]
 expected=special+[rows_rref([unit(i),tau],18) for i in range(12,17)]
 assert {tuple(map(tuple,L)) for L in expected}==set(lines)
 online=[j for j,f in enumerate(forms) if any(len(rows_rref(list(L)+[f],18))==2 for L in lines)]
 isolated=[j for j in range(n) if j not in online];assert len(online)==24 and len(isolated)==65
 assert {tuple(rows_rref([forms[j]],18)[0]) for j in isolated}==set(points)
 assert all(any(row[j] for row in B) for j in isolated)
 line_basis=echelon([cubes[j] for j in online]);assert len(line_basis)==24
 for L in lines:
  l=lin(L[0]);r=lin(L[1])
  for p in [mul(mul(l,l),l),mul(mul(l,l),r),mul(mul(l,r),r),mul(mul(r,r),r)]:assert not reduce(p,line_basis)
 # Build the annihilator conditions after the five Wi binary lines are forced to vanish.
 conditions=[{j:c for j,c in enumerate(row) if c} for row in B]+[{j:Q(1)} for j in isolated]
 assert forms[d["padding_index"]]==tau and B[17]==[Q(j==d["padding_index"]) for j in range(n)]
 for i in range(12,17):
  local=[j for j in online if forms[j][i]];assert len(local)==3
  w=variable(i,18);t=variable(17,18)
  for poly in [mul(mul(w,w),w),mul(mul(w,w),t),mul(mul(w,t),t),mul(mul(t,t),t)]:
   residual,coordinate=project(poly);assert not residual
   assert set(coordinate)<=set(local)|{d["padding_index"]}
   conditions.append(coordinate)
 assert len(echelon(conditions))==87
 def exponent(entries):
  e=[0]*18
  for i,k in entries:e[i]=k
  return tuple(e)
 eb=exponent([(3,1),(17,2)]);ea=exponent([(1,1),(17,2)])
 dualD=[cubes[j].get(eb,Q(0)) if j in online else Q(0) for j in range(n)]
 dualE=[cubes[j].get(ea,Q(0)) if j in online else Q(0) for j in range(n)]
 assert len(rows_rref([dualD,dualE],n))==2
 assert all(sum((c*v[j] for j,c in row.items()),Q(0))==0 for row in conditions for v in [dualD,dualE])
 # Confirm the target contains only the required tau^2 Wi term in each Wi direction.
 P=[]
 for row in B:
  f={}
  for c,p in zip(row,cubes):f=add(f,p,c)
  P.append(f)
 assert P[17]=={exponent([(17,3)]):Q(1)}
 for i in range(5):
  e=exponent([(12+i,1),(17,2)])
  for o,f in enumerate(P):assert {x:c for x,c in f.items() if x[12+i]}==({e:Q(-1)} if o==i else {})
 specialids=[j for j in online if all(not forms[j][i] for i in range(18) if i not in [1,3,17])];assert len(specialids)==9
 assert len(rows_rref([[row[j] for j in specialids] for row in B],9))==7
 # Restriction table for the two remaining functionals on the three special lines.
 table=[]
 for L in [[b,tau],[a,[x+y for x,y in zip(b,tau)]],[a,[x-y for x,y in zip(b,tau)]]]:
  l=lin(L[0]);r=lin(L[1]);actualD=[];actualE=[]
  polys=[mul(mul(l,l),l),scale(mul(mul(l,l),r),Q(3)),scale(mul(mul(l,r),r),Q(3)),mul(mul(r,r),r)]
  for poly in polys:
   residual,coord=project(poly);assert not residual
   actualD.append(sum((c*dualD[j] for j,c in coord.items()),Q(0)));actualE.append(sum((c*dualE[j] for j,c in coord.items()),Q(0)))
  table.append([actualD,actualE])
 assert table==[[[0,0,3,0],[0,0,0,0]],[[0,0,0,3],[0,0,3,0]],[[0,0,0,3],[0,0,3,0]]]
 classification=dict(components=[dict(dimension=len(C),basis=[[str(c) for c in row] for row in C]) for C in components],
    isolated_indices=isolated,line_indices=online,special_indices=specialids,dual_D=list(map(str,dualD)),dual_E=list(map(str,dualE)))
 (out/'classification.json').write_text(json.dumps(classification,indent=2,sort_keys=True)+'\n')
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),helper_sha256=sha(HERE/'verify_power_single_compress.py'),
    input_sha256=sha(source),construction_validation_sha256=sha(HERE/'power89/validation.json'),primary_checkpoint_sha256=sha(out/'checkpoint.json'),
    classification_sha256=sha(out/'classification.json'),cases=cases,tree_nodes=total_nodes,factor_identities=factor_checks,
    isolated_projective_points=65,projective_lines=8,cube_span_dimension=89,line_cube_span_dimension=24,
    forced_isolated_forms=65,annihilator_constraint_rank=87,remaining_dual_dimension=2,special_target_rank=7,zero_locus_special_span_bound=5,
    all_factor_trees_checked=True,all_cube_supports_covered=True,all_line_cubic_spans_checked=True,annihilator_restrictions_checked=True,
    seconds=round(time.monotonic()-started,3),
    scope='Every representation of the fixed padded companion using cubes lying in the saved89-dimensional scalar cubic span has at least89 forms. This is not a global tensor-rank or Keller-dimension minimum; smaller representations must change the scalar cube span.')
 (out/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='cases'}))
if __name__=='__main__':main()
