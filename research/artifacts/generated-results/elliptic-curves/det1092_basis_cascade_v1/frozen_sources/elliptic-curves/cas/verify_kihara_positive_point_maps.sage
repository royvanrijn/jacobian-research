#!/usr/bin/env sage-python
"""Independent exact height-metric, roster and centre verification."""
import json,sys
from pathlib import Path
from decimal import Decimal,localcontext
from sage.all import QQ,ZZ,matrix,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import kihara_positive_point_pilot as batch
import certify_compact_r17_candidates as cert
from research_runtime.store import digest
def main():
    p=batch.protocol();total=0
    for i,row in enumerate(p['rows']):
        batch.configure(i);seed=cert.read(batch.SEED);maps=cert.read(batch.D/'maps.json');data=cert.read(batch.D/'result.json');rank=row['initial_rank'];g=matrix(ZZ,maps['rounded_gram']);u=matrix(ZZ,maps['change_of_basis']);h=matrix(ZZ,maps['reduced_gram'])
        assert maps['protocol_hash']==data['protocol_hash']==digest(p) and g.is_symmetric() and g.is_positive_definite() and abs(u.det())==1 and h==u*g*u.transpose()
        with localcontext() as ctx:
            ctx.prec=110;expected=[[int((Decimal(c)*1000000).to_integral_value()) for c in r] for r in maps['metric_gram']]
        assert matrix(ZZ,expected)==g and [r['parity'] for r in maps['sample']]==batch.masks(p)
        for r in maps['sample']:
            w=matrix(ZZ,1,rank,r['representative']);v=matrix(ZZ,1,rank,r['reduced_representative']);assert w==v*u and (w*g*w.transpose())[0,0]==r['metric_norm'];assert all((w[0,j]-(r['parity']>>j))%2==0 for j in range(rank))
        order=sorted(maps['sample'],key=lambda r:(-r['metric_norm'],r['parity']));candidates=[(r,[0]*rank) for r in order]
        candidates.extend((r,[2*sgn if j==axis else 0 for j in range(rank)]) for axis in range(rank) for sgn in (-1,1) for r in order)
        expected=[];seen=set()
        for r,offset in candidates:
            w=[a+b for a,b in zip(r['representative'],offset)];nz=next((c for c in w if c),0)
            if not nz:continue
            sign=1 if nz>0 else -1;w=tuple(sign*c for c in w)
            if w in seen:continue
            seen.add(w);norm=(matrix(ZZ,1,rank,w)*g*matrix(ZZ,rank,1,w))[0,0]
            expected.append({'parity':r['parity'],'representative':list(w),'metric_norm':int(norm),'origin_representative':r['representative'],'offset':offset,'sign':sign})
            if len(expected)==49:break
        assert expected==maps['centres'] and [m['centre'] for m in maps['rows']]==expected
        E=EllipticCurve(QQ,seed['curve']);basis=[E([QQ(c) for c in P]) for P in seed['points']]
        for m,c in zip(maps['rows'],data['charts']):
            P=sum((ZZ(a)*Q for a,Q in zip(m['centre']['representative'],basis)),E(0));assert P and {'x':str(P[0]),'y':str(P[1])}==c['search']['base_point'];total+=1
        batch.replay()
    assert total==294;print('PASS294 exact centre/metric/map/point replays, including torsion-compatible empty states',flush=True)
if __name__=='__main__':main()
