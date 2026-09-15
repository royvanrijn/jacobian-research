"""Exact fifteen-anchor exclusion and one unresolved anchor for determinant852."""
from sage.all import *
import argparse,json,hashlib,functools
from pathlib import Path
R=Path(__file__).resolve().parents[2];P=R/'artifacts/generated-results/elkies-k3-det852-root-gate-v1'
mark=R/'artifacts/generated-results/elkies-k3-det852-marking-v1/certificate.json'
cat=R/'artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json'
anc=R/'artifacts/generated-results/elkies-k3-niemeier-d5-anchor-orbits.json'
F=-matrix(ZZ,json.loads(mark.read_text())['NS_gram'])[2:,2:]
K=matrix(ZZ,[[2,-1,-1,-1,1,-1,1],[-1,2,0,0,-1,0,0],[-1,0,2,0,0,1,0],[-1,0,0,2,0,0,-1],[1,-1,0,0,2,-1,1],[-1,0,1,0,-1,14,0],[1,0,0,-1,1,0,18]])
assert K.det()==852 and K.is_positive_definite() and all(K[i,i]%2==0 for i in range(7))
C5=matrix(ZZ,CartanMatrix(['D',5]));J=matrix(ZZ,K[:5,:5].__pari__().qfisom(C5.__pari__())).inverse();assert J.transpose()*K[:5,:5]*J==C5 and abs(J.det())==1
Schur=K[5:,5:]-K[5:,:5]*K[:5,:5].inverse()*K[:5,5:];assert Schur==matrix(QQ,[[QQ(51)/4,QQ(3)/4],[QQ(3)/4,QQ(67)/4]]) and Schur.trace()==QQ(59)/2

def gen(G):
 d,u,v=G.smith_form();assert d.diagonal()==[1]*(G.nrows()-1)+[852]
 z=G.inverse()*u.inverse()*vector(ZZ,[0]*(G.nrows()-1)+[1]);return z,z*G*z/2
kf,qk=gen(K);ff,qf=gen(F);mult=next(a for a in range(852) if gcd(a,852)==1 and qk+a*a*qf in ZZ)
Z=block_diagonal_matrix(K,F);z=vector(QQ,list(kf)+list(mult*ff));B=matrix(ZZ,[list(852*v) for v in identity_matrix(ZZ,24).rows()]+[list(852*z)]).row_module().basis_matrix()/852;N=B*Z*B.transpose()
assert all(a in ZZ for a in N.list()) and all(N[i,i]%2==0 for i in range(24)) and N.det()==1 and N.is_positive_definite()
energy_records={}
@functools.cache
def bound(typ,n):
 key=typ+str(n)
 if typ=='A' and n>=10:
  den=n+1;M=[[min(i+1,j+1)*(den-max(i+1,j+1)) for j in range(n)] for i in range(n)];q=[0]*(1<<n);best=None;winner=None
  for mask in range(1<<n):
   if mask:
    bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit;z=rest;dot=0
    while z:
     bit=z&-z;dot+=M[i][bit.bit_length()-1];z^=bit
    q[mask]=q[rest]+M[i][i]+2*dot
   total=0;run=0
   for i in range(n+1):
    if i<n and not(mask>>i&1):run+=1
    else:total+=run*(run+1)*(run+2);run=0
   value=12*q[mask]+den*total
   if best is None or value<best:best=value;winner=mask
  value=QQ(best)/(12*den);energy_records[key]={'support_subsets':1<<n,'minimizer':winner,'lower_bound':str(value),'method':'exact A inverse Cartan and zero-run recurrence'};return value
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

# A7: eight distinct integer coordinate pairs, centered at their mean.
# Modulo integer translations, both centroid coordinates lie in {0,...,7}/8.
a7_cases=[];a7_min=None;equalities=[]
for a in range(8):
 for b in range(8):
  center=vector(QQ,[QQ(a)/8,QQ(b)/8])
  points=sorted((sum((QQ(z)-c)**2 for z,c in zip((x,y),center)),x,y) for x in range(-4,5) for y in range(-4,5))
  assert points[7][0]<(5-QQ(7)/8)**2 # all omitted grid points are farther away
  total=sum(t[0] for t in points[:8]);a7_cases.append({'centroid_numerators':[a,b],'lower_bound':str(total)})
  if a7_min is None or total<a7_min:a7_min=total;equalities=[]
  if total==a7_min:equalities.append((a,b,points))
assert a7_min==QQ(39)/4
sharp_equalities=[]
for a,b,points in equalities:
 assert points[7][0]<points[8][0] # unique minimizing eight-point set
 V=matrix(QQ,[[x-QQ(a)/8,y-QQ(b)/8] for _,x,y in points[:8]])
 assert sum(V.column(0))==sum(V.column(1))==0
 H=V.transpose()*V;assert H[0,0]==H[1,1]==QQ(39)/8
 sharp_equalities.append({'centroid_numerators':[a,b],'centered_gram':[[str(z) for z in row] for row in H.rows()]})
# D5 dual vectors have all coordinates integral, or all half-integral.
# Norm <=1 bounds each numerator by2. Exhaust these11 vectors.
import itertools
small_d5=[]
for z in itertools.product(range(-2,3),repeat=5):
 if len({x%2 for x in z})==1 and sum(x*x for x in z)<=4:
  small_d5.append(list(z));assert sum(x*x for x in z) in [0,4]
assert len(small_d5)==11
assert sum(matrix(QQ,CartanMatrix(['D',4])).inverse().list())==14

catalog=json.loads(cat.read_text());anchors=json.loads(anc.read_text());models={a['label']:matrix(ZZ,a['gram']) for a in catalog['rooted_niemeier_lattices']}
assert len(models)==23 and len(anchors['anchors'])==16
results=[];root_cache={}
for a in anchors['anchors']:
 label=a['niemeier'];G=models[label];A=matrix(ZZ,a['D5_basis_in_ambient']);assert A*G*A.transpose()==C5
 if label not in root_cache:
  d=G.__pari__().qfminim(2,100000,2);W=matrix(ZZ,d[2]).transpose();assert 2*W.nrows()==int(d[0]);root_cache[label]=W.stack(-W)
 W=root_cache[label];pair=W*G*A.transpose();indices=[i for i in range(W.nrows()) if not any(pair.row(i))];orth=W.matrix_from_rows(indices);types=root_types(orth,G)
 lower=sum(bound(t,n) for t,n in types)
 if types==[['A',7],['A',7],['D',5]]:
  assert lower==QQ(99)/4 and bound('D',5)==6
  lower=2*a7_min+bound('D',5);assert lower==QQ(51)/2 and lower<Schur.trace()
 else:assert lower>Schur.trace()
 results.append({'niemeier':label,'anchor_index':a['anchor_index'],'orthogonal_root_types':types,'two_vector_lower_bound':str(lower)})
assert len({r['niemeier'] for r in results})==13
minimum=min(QQ(r['two_vector_lower_bound']) for r in results);assert minimum==QQ(51)/2

def rows(G):return [[str(x) for x in r] for r in G.rows()]
out={'schema':'elkies-k3.det852-one-anchor-reduction.v1','status':'PASS_FIFTEEN_ANCHORS_EXCLUDED_ONE_UNRESOLVED','NS_determinant':852,'source_NS':'Saturation-certified NS of K3-4ff75fec54d01662','auxiliary_gram':rows(K),'D5_basis_change':rows(J),'projected_two_vector_gram':rows(Schur),'available_trace':'59/2','discriminant_glue':{'q_aux':str(qk),'q_frame':str(qf),'multiplier':mult,'index':852,'unimodular_Gram':rows(N)},'anchor_cases':results,'energy_bounds':energy_records,'minimum_required_trace':str(minimum),'remaining_case':'2A7+2D5, anchor1:lower bound51/2 below available59/2; primitive embedding and rootless existence unresolved','A7_centroid_cases':a7_cases,'A7_equalities':sharp_equalities,'D5_dual_norm_le_one_twice_coordinates':small_d5,'inputs':[{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [mark,cat,anc]],'theorem_inputs':['Primitive discriminant anti-isometry gluing','Niemeier classification','Complete16 D5 anchor orbits in retained EC-K3-H3-ROOTLESS-J2-COMPLETE','Written two-vector Weyl-support energy lemma','Written A7 centered integer-pair bound and equality; D5 dual short-vector argument'],'conclusion':'Any rootless frame in the pinned852 NS genus must arise from a primitive embedding of the displayed auxiliary at the unique D5 anchor in2A7+2D5. The other15 anchor cases are impossible. Existence in the remaining case is UNKNOWN.','scope':'Only this actual NS genus, not every lattice with determinant852.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();output=P/'certificate.json'
if args.check:assert json.loads(output.read_text())==out
else:output.write_text(json.dumps(out,indent=2,sort_keys=True,default=int)+'\n')
print('PASS852:15 Niemeier anchors excluded; sole remaining2A7+2D5 anchor unresolved. No global MW17 exclusion.')
