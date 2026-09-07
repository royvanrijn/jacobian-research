import json
from sage.all import QQ,ZZ,matrix,vector,PolynomialRing,EllipticCurve
from itertools import combinations,product
from pathlib import Path
root=Path.cwd();art=root/'artifacts/generated-results/elliptic-curves'
r=json.loads((art/'mestre_rational_ns_gram_v2.json').read_text())['rows'][0];b=json.loads((art/'mestre_468_replay_bundle_v1.json').read_text())['rows'][0];h=b['generic_heights'];G=matrix(QQ,r['rational_NS_Gram']);H=matrix(QQ,h['seed_height_gram']);std=matrix.identity(QQ,18).rows();F,O=std[:2];C=G[2:7,2:7]
phi=[]
for i,p in enumerate(r['section_profiles']):
 v=std[7+i]-O-(2+p['zero_section_intersection'])*F
 v[2:7]-=C.inverse()*vector(QQ,G[2:7,7+i].column(0));phi.append(v)
visible=[std[i] for i in range(1,18)]+[F-std[2],F-std[3],F-sum(std[4:7])]
words=[vector(ZZ,[int(k==i) for k in range(11)]) for i in range(11)]
for i,j in combinations(range(11),2):
 for s in (-1,1):words.append(vector(ZZ,[int(k==i)+s*int(k==j) for k in range(11)]))
for indices in combinations(range(11),3):
 for signs in product((-1,1),repeat=2):
  w=vector(ZZ,11);w[indices[0]]=1
  for i,s in zip(indices[1:],signs):w[i]=s
  words.append(w)
found=[]
for w in words:
 comp=[sum(w[i]*r['section_profiles'][i]['components'][j] for i in range(11))%n for j,n in enumerate((2,2,4))]
 corr=QQ(comp[0]+comp[1])/2+QQ(comp[2]*(4-comp[2]))/4;hh=w*H*w;oo=(hh-4+corr)/2
 if oo!=2:continue
 v=O+4*F+sum((w[i]*phi[i] for i in range(11)),vector(QQ,18));cross=vector(QQ,comp[:2]+[int(comp[2]==k) for k in (1,2,3)])
 v[2:7]+=C.inverse()*cross
 assert v*G*v==-2 and v*G*O==2 and all(c.denominator()==1 for c in v)
 D=O+v;vertical=[t for t in visible+[v] if D*G*t==0];bound=17-matrix(QQ,vertical).rank();found.append((int(bound),list(w),list(comp)))
out=root/'artifacts/local/elliptic-curves/mestre-chord-pencils-v1/roster.json'
assert not out.exists()
out.write_text(json.dumps({'word_count':len(words),'rows':[{'upper_bound_from_visible_curves':bound,'word':[int(c) for c in w],'components':[int(c) for c in co]} for bound,w,co in sorted(found,reverse=True)]},indent=2)+'\n')
print('words',len(words),'eligible',len(found),'bound histogram',{k:sum(x[0]==k for x in found) for k in sorted(set(x[0] for x in found))});print([x for x in found if x[0]>11])
