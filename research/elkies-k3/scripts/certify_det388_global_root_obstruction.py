"""All-Niemeier two-vector energy obstruction for every determinant388 frame in the pinned genus."""
from sage.all import *
import argparse,json,hashlib,functools
from pathlib import Path
R=Path(__file__).resolve().parents[2];P=R/'artifacts/generated-results/elkies-k3-det388-global-root-gate-v1'
mark=R/'artifacts/generated-results/elkies-k3-det388-marking-v1/certificate.json'
cat=R/'artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json'
anc=R/'artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json'
F=-matrix(ZZ,json.loads(mark.read_text())['NS_gram'])[2:,2:]
K=matrix(ZZ,[[2,1,1,1,-1,1,0],[1,2,0,1,0,0,0],[1,0,2,1,0,0,0],[1,1,1,2,0,0,0],[-1,0,0,0,2,-1,0],[1,0,0,0,-1,8,-1],[0,0,0,0,0,-1,14]])
assert K.det()==388 and K.is_positive_definite() and all(K[i,i]%2==0 for i in range(7))
C5=matrix(ZZ,CartanMatrix(['D',5]));J=matrix(ZZ,K[:5,:5].__pari__().qfisom(C5.__pari__())).inverse();assert J.transpose()*K[:5,:5]*J==C5 and abs(J.det())==1
Schur=K[5:,5:]-K[5:,:5]*K[:5,:5].inverse()*K[:5,5:];assert Schur==matrix(QQ,[[7,-1],[-1,14]]) and Schur.trace()==21

def gen(G):
 d,u,v=G.smith_form();assert d.diagonal()==[1]*(G.nrows()-1)+[388]
 z=G.inverse()*u.inverse()*vector(ZZ,[0]*(G.nrows()-1)+[1]);return z,z*G*z/2
kf,qk=gen(K);ff,qf=gen(F);mult=next(a for a in range(388) if gcd(a,388)==1 and qk+a*a*qf in ZZ)
Z=block_diagonal_matrix(K,F);z=vector(QQ,list(kf)+list(mult*ff));B=matrix(ZZ,[list(388*v) for v in identity_matrix(ZZ,24).rows()]+[list(388*z)]).row_module().basis_matrix()/388;N=B*Z*B.transpose()
assert all(a in ZZ for a in N.list()) and all(N[i,i]%2==0 for i in range(24)) and N.det()==1 and N.is_positive_definite()
energy_records={}
@functools.cache
def bound(typ,n):
 key=typ+str(n)
 if typ=='A' and n>=10:
  a=(n-1)//2;b=n-1-a;value=bound('A',a)+bound('A',b);energy_records[key]={'subsystem':[['A',a],['A',b]],'lower_bound':str(value)};return value
 if typ=='D' and n>=10:
  a=n//2;b=n-a;value=bound('D',a)+bound('D',b);energy_records[key]={'subsystem':[['D',a],['D',b]],'lower_bound':str(value)};return value
 C=matrix(QQ,CartanMatrix([typ,n]));I=C.inverse();values=[]
 for mask in range(1<<n):
  support=[i for i in range(n) if mask>>i&1];zero=[i for i in range(n) if not(mask>>i&1)]
  v=sum(I[i,j] for i in support for j in support)
  if zero:v+=sum(C.matrix_from_rows_and_columns(zero,zero).inverse().list())
  values.append(v)
 minimum=min(values);energy_records[key]={'support_subsets':len(values),'values':[str(v) for v in values],'lower_bound':str(minimum)};return minimum

def root_types(V,G):
 H=V*G*V.transpose();remaining=set(range(V.nrows()));types=[]
 while remaining:
  stack=[min(remaining)];remaining.remove(stack[0]);component=[]
  while stack:
   i=stack.pop();component.append(i);more=[j for j in remaining if H[i,j]];remaining.difference_update(more);stack.extend(more)
  rank=V.matrix_from_rows(component).rank();num=len(component)
  if num==rank*(rank+1):typ='A'
  elif rank>=4 and num==2*rank*(rank-1):typ='D'
  elif (rank,num) in [(6,72),(7,126),(8,240)]:typ='E'
  else:raise ValueError((rank,num))
  types.append([typ,rank])
 return sorted(types)

catalog=json.loads(cat.read_text());anchors=json.loads(anc.read_text());models={a['label']:matrix(ZZ,a['gram']) for a in catalog['rooted_niemeier_lattices']}
assert len(models)==23 and len(anchors['anchors'])==16
results=[];root_cache={}
for a in anchors['anchors']:
 label=a['niemeier'];G=models[label];A=matrix(ZZ,a['D5_basis_in_ambient']);assert A*G*A.transpose()==C5
 if label not in root_cache:
  d=G.__pari__().qfminim(2,100000,2);W=matrix(ZZ,d[2]).transpose();assert 2*W.nrows()==int(d[0]);root_cache[label]=W.stack(-W)
 W=root_cache[label];pair=W*G*A.transpose();indices=[i for i in range(W.nrows()) if not any(pair.row(i))];orth=W.matrix_from_rows(indices);types=root_types(orth,G)
 lower=sum(bound(t,n) for t,n in types);assert lower>21
 results.append({'niemeier':label,'anchor_index':a['anchor_index'],'orthogonal_root_types':types,'two_vector_lower_bound':str(lower)})
assert len({r['niemeier'] for r in results})==13
minimum=min(QQ(r['two_vector_lower_bound']) for r in results);assert minimum==QQ(133)/6

def rows(G):return [[str(x) for x in r] for r in G.rows()]
out={'schema':'elkies-k3.det388-global-root-obstruction.v1','status':'PASS_ALL_FRAMES_ROOTFUL','NS_determinant':388,'source_surface_id':'K3-f5168aef816b1b25','auxiliary_gram':rows(K),'D5_basis_change':rows(J),'projected_two_vector_gram':rows(Schur),'available_trace':'21','discriminant_glue':{'q_aux':str(qk),'q_frame':str(qf),'multiplier':mult,'index':388,'unimodular_Gram':rows(N)},'anchor_cases':results,'energy_bounds':energy_records,'minimum_required_trace':str(minimum),'strict_gap':str(minimum-21),'inputs':[{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [mark,cat,anc]],'theorem_inputs':['Primitive discriminant anti-isometry gluing','Niemeier classification','Complete16 D5 anchor orbits in retained EC-K3-H3-ROOTLESS-J2-COMPLETE','Written two-vector Weyl-support energy lemma'],'conclusion':'Every positive even rank17 lattice with the pinned frame discriminant form has a root; no MW17 fibration on the full determinant388 NS. Full rational rank19 marking remains valid.','scope':'Only this actual NS genus, not every lattice with determinant388.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();output=P/'certificate.json'
if args.check:assert json.loads(output.read_text())==out
else:output.write_text(json.dumps(out,indent=2,sort_keys=True,default=int)+'\n')
print('PASS det388 globally rootful: all16 Niemeier D5 cases need trace>=133/6>21. Rational rank19 marking survives; MW17 excluded.')
