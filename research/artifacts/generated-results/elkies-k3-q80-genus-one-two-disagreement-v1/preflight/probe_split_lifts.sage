"""Complete first lift test for the 21 mixed root-contact pairs."""
from sage.all import GF, Integers, PolynomialRing, QQ, matrix, vector
from pathlib import Path
import json
import resource
import time

resource.setrlimit(resource.RLIMIT_CPU, (30,35))
resource.setrlimit(resource.RLIMIT_AS, (4*1024**3,4*1024**3))
started=time.process_time();ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent.parent;p=131
preview=json.loads((OUT/'integral-preview.json').read_text())
points=json.loads((ROOT/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/norm4-sections.json').read_text())['records']
source=json.loads((ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json').read_text())
F=GF(p);R=PolynomialRing(F,'t');t=R.gen()
Z=Integers(p*p);R2=PolynomialRing(Z,'t');T=R2.gen()
def lift(r):return R2([int(v) for v in r.list()])
def rational_poly(key):
    return R2([Z(QQ(v).numerator())/Z(QQ(v).denominator())
               for v in source['weierstrass_model'][key+'_coefficients_low_to_high']])
A2,B2=rational_poly('A'),rational_poly('B');A=R(A2)
rows=[]
for pair in preview['split_mixed_pairs']:
    for b in pair['disagree']:
        for c in pair['agree']:
            assert b is not None and b!=c
            u=t-b;v=R(1) if c is None else t-c;H=u*v;D=H*H
            fixed=2 if c is None else 4
            assert D[fixed]==1
            columns=[];errors=[]
            for index in pair['pair']:
                point=points[index];X=R(point['x']);Y=R(point['y'])
                S,rem=Y.quo_rem(H);assert not rem and S.degree()<=4
                local=[t**k*(3*X*X+A) for k in range(5)]
                local += [-2*D*S*t**k for k in range(5)]
                dcols=[-S*S*t**k for k in range(5) if k!=fixed]
                columns.append((local,dcols))
                error=lift(X)**3+A2*lift(X)+B2-lift(D)*lift(S)**2
                assert all(int(error[k])%p==0 for k in range(13))
                errors += [F(-int(error[k])//p) for k in range(13)]
            J=matrix(F,26,24)
            for n,(local,dcols) in enumerate(columns):
                for j,f in enumerate(local):
                    for k in range(13):J[13*n+k,10*n+j]=f[k]
                for j,f in enumerate(dcols):
                    for k in range(13):J[13*n+k,20+j]=f[k]
            aug=J.augment(vector(F,errors).column())
            rows.append({'pair':pair['pair'],'disagree':b,'agree':c,
                         'jacobian_rank':int(J.rank()),'augmented_rank':int(aug.rank()),
                         'lifts_mod_p2':bool(J.rank()==aug.rank())})
result={'status':'PREVIEW_ONLY','prime':p,'rows':rows,'cpu_seconds':time.process_time()-started}
with (OUT/'split-lift-preview.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
print(json.dumps(result,sort_keys=True),flush=True)
