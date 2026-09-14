"""Finite coefficient gate for the only etale two-disagreement branch residue."""
from sage.all import GF, PolynomialRing, QQ
from pathlib import Path
from collections import defaultdict
import json
import resource
import time

resource.setrlimit(resource.RLIMIT_CPU, (30, 35))
resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
started=time.process_time()
ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent.parent
source=json.loads((ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json').read_text())
F=GF(131);R=PolynomialRing(F,'t');t=R.gen()
A=R([F(QQ(v)) for v in source['weierstrass_model']['A_coefficients_low_to_high']])
B=R([F(QQ(v)) for v in source['weierstrass_model']['B_coefficients_low_to_high']])
q=t*t+62*t+88
node=F(38)/31*t+(F(122)-F(38)/31*100)
simple=F(55)/31*t+(F(18)-F(55)/31*100)
assert (node**3+A*node+B)%q==0 and (3*node**2+A)%q==0
assert (simple**3+A*simple+B)%q==0 and (3*simple**2+A)%q!=0
g=(t-35)*(t-75);D=q*g
v35=(F(114)-node(35))/q(35);v75=(F(11)-node(75))/q(75)
x0=node+q*(v35+(v75-v35)/F(40)*(t-35))
assert x0.degree()<4 and x0%q==node and x0(35)==114 and x0(75)==11
etale=[]
for c in F:
    x=x0+c*D;f=x**3+A*x+B
    S,rem=f.quo_rem(D);assert not rem and S
    monic=S.monic()
    if monic.is_square():
        etale.append({'c':int(c),'x':[int(v) for v in x.list()],
                      'scalar':int(S.leading_coefficient()),'square_root':[int(v) for v in monic.sqrt().list()]})
points=json.loads((ROOT/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/norm4-sections.json').read_text())['records']
incidence=defaultdict(list)
for point in points:
    X,Y=R(point['x']),R(point['y'])
    for b in list(F)+[None]:
        xx,yy=(X[4],Y[6]) if b is None else (X(b),Y(b))
        if not yy:incidence[None if b is None else int(b)].append((point['index'],int(xx)))
pairs=defaultdict(dict)
for b,rows in incidence.items():
    for k,(i,e) in enumerate(rows):
        for j,f in rows[k+1:]:pairs[(i,j)][b]=(e==f)
mixed=[{'pair':list(ij),'agree':[b for b,v in sites.items() if v],
        'disagree':[b for b,v in sites.items() if not v]}
       for ij,sites in pairs.items() if any(sites.values()) and not all(sites.values())]
result={'status':'PREVIEW_ONLY','prime':131,
        'A':[int(v) for v in A.list()],'B':[int(v) for v in B.list()],
        'nodal_quadratic':[int(v) for v in q.list()],
        'node_root':[int(v) for v in node.list()],
        'simple_root':[int(v) for v in simple.list()],
        'etale_D':[int(v) for v in D.list()],
        'etale_x0':[int(v) for v in x0.list()],
        'etale_candidates':etale,'split_mixed_pairs':mixed,
        'cpu_seconds':time.process_time()-started}
with (OUT/'integral-preview.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps({'etale_candidates':etale,'split_mixed_pairs':len(mixed),
                  'sample_mixed_pairs':mixed[:8],'cpu_seconds':result['cpu_seconds']}),flush=True)
