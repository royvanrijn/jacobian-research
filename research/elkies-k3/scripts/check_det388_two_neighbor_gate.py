"""Replay all-parity norm-eight obstruction for the one pinned388 frame."""
from sage.all import *
from sage.quadratic_forms.genera.genus import Genus
import json
from pathlib import Path
R=Path(__file__).resolve().parents[2]
P=R/'artifacts/generated-results/elkies-k3-det388-two-neighbor-gate-v1'
c=json.loads((P/'certificate.json').read_text());G=matrix(ZZ,c['source_gram']);roots=matrix(ZZ,c['source_roots']);n=17
assert G.nrows()==n and G.det()==388 and roots*G*roots.transpose()==2*identity_matrix(ZZ,5)
x=json.loads((P.parent/'elkies-k3-det388-rootless-probe-v1/checkpoint.json').read_text())
s=json.loads((P.parent/'elkies-k3-det388-marking-v1/certificate.json').read_text());F=-matrix(ZZ,s['NS_gram'])[2:,2:]
assert Genus(G)==Genus(F)
U=matrix(ZZ,x['start_transport']);G0=matrix(ZZ,x['start_gram']);assert abs(U.det())==1 and U*F*U.transpose()==G0
index=x['best_index'];assert matrix(ZZ,x['rows'][index]['gram'])==G
while index!=-1:
 r=x['rows'][index];parent=r['parent'];assert -1<=parent<index
 A=matrix(QQ,r['transport']);H=matrix(ZZ,r['gram']);B=G0 if parent==-1 else matrix(ZZ,x['rows'][parent]['gram'])
 assert abs(A.det())==1 and A*B*A.transpose()==H
 index=parent
rad=matrix(GF(2),G).right_kernel();assert rad.dimension()==1
v=vector(ZZ,rad.basis()[0]);assert (v*G*v)%8==4
N=1<<n;off=[sum(1<<j for j in range(n) if j!=i and G[i,j]%2) for i in range(n)]
masks=[sum(1<<j for j in range(n) if (v*G)[j]%2) for v in roots.rows()]
q=bytearray(N);eligible=[]
for a in range(1,N):
 bit=a&-a;i=bit.bit_length()-1;b=a^bit
 q[a]=q[b]^int(G[i,i]//2%2)^((off[i]&b).bit_count()%2)
 if not q[a] and all((a&m).bit_count()%2 for m in masks):eligible.append(a)
assert len(eligible)==1792 and eligible==c['parity_candidates']
assert set(map(int,c['norm8_witnesses']))==set(eligible)
for a in eligible:
 w=vector(ZZ,c['norm8_witnesses'][str(a)])
 assert w*G*w==8 and sum(1<<j for j in range(n) if w[j]%2)==a
assert c['uncovered']==[] and c['global_rootless_existence'] is None
print('PASS all1792 parity witnesses; radical excluded; source transport exact. One-frame2-neighbour exclusion only.')
