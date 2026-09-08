#!/usr/bin/env python3
"""Fixed retrospective one-basis translations of a proved public extra direction."""
import argparse
from collections import Counter
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_inventory188_public28_visibility as direct
from search_observability import point_visibility
from research_runtime.store import checkpoint,digest

ROOT=direct.ROOT;ART=direct.ART
D=ROOT/'artifacts/local/elliptic-curves/inventory188-public28-translates-v1'
OUT=ART/'inventory188_public28_translates_v1.json'


def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in
            (Path(__file__).resolve(),Path(direct.__file__),Path(cert.__file__),
             ROOT/'elliptic-curves/cas/search_observability.py',direct.PROOF,direct.REPLAY,direct.SOURCE,direct.OUT)}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve translation audit protocol')
    checkpoint(D/'protocol.json',{'sources':sources(),'translations':55,'charts':49,'signs':[-1,1],
               'observations':5390,'seconds_per_stage':300,'rss_bytes':1073741824,
               'rule':'Exactly P and P +/- B_i for every i=0..26 of the previously certified local27 basis. P is the already proved extra public point26. Test both signs of each resulting point in all49 original completed charts. No centre change, coordinate optimization, further word, oracle-driven point search or new parameter scan.'})


def expected():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('fixed audit inputs changed')
    if direct.expected()!=cert.read(direct.OUT):raise ArithmeticError('direct visibility prerequisite differs')
    proof=cert.read(direct.PROOF);model=tuple(map(cert.F,proof['curve']))
    P=tuple(map(cert.F,proof['transported_public_points'][26]));basis=[tuple(map(cert.F,q)) for q in proof['points'][:27]]
    row=next(r for r in cert.read(direct.SOURCE)['curves'] if r['parameter']==proof['parameter'])
    raw=cert.read(ROOT/row['discovery_witness']['path']);words=[(None,0,P)]
    for i,(x,y) in enumerate(basis):
        for sign in (-1,1):
            if x==P[0]:raise ArithmeticError('extra point unexpectedly equals a signed old basis point')
            slope=(sign*y-P[1])/(x-P[0]);X=slope*slope-x-P[0];Y=slope*(P[0]-X)-P[1]
            if not cert.is_on_weierstrass_curve(model,(X,Y)):raise ArithmeticError('translated point off curve')
            words.append((i,sign,(X,Y)))
    if len(words)!=55 or len({Q for i,s,Q in words})!=55:raise ArithmeticError('fixed distinct translation roster differs')
    rows=[];counts=Counter();hashes=[];discrepancies=[]
    for word_index,(i,sign,Q) in enumerate(words):
        best=None;local=Counter()
        for chart,c in enumerate(raw['charts']):
            record=c['search']
            for ordinate_sign in (-1,1):
                v=point_visibility({**record,'completed_denominator':125000},(Q[0],ordinate_sign*Q[1]))
                observed={'word_index':word_index,'chart':chart,'ordinate_sign':ordinate_sign,'visibility':v}
                hashes.append(digest(observed));counts[v['status']]+=1;local[v['status']]+=1
                if v['status']=='VISIBLE_NOT_RECORDED':discrepancies.append(observed)
                height=v.get('minimum_affine_height')
                if height is not None:
                    key=(height,chart,ordinate_sign)
                    if best is None or key<best[0]:best=(key,observed)
        rows.append({'word_index':word_index,'basis_index':i,'coefficient':sign,
                     'point':list(map(str,Q)),'status_counts':dict(local),'best_observation':best[1] if best else None})
    if len(hashes)!=5390:raise ArithmeticError('all fixed observations required')
    best=min((r['best_observation'] for r in rows if r['best_observation']),key=lambda r:(r['visibility']['minimum_affine_height'],r['word_index'],r['chart'],r['ordinate_sign']))
    return {'schema':'elliptic-curves.inventory188-public28-translates.v1','status':'PASS_EXACT_OBSERVATIONS',
            'sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),
            'curve':proof['curve'],'public_point':list(map(str,P)),'old_basis':[list(map(str,Q)) for Q in basis],
            'rows':rows,'observation_count':len(hashes),'ordered_observation_digest':digest(hashes),
            'status_counts':dict(counts),'best_observation':best,'discrepancies':discrepancies,
            'claim_boundary':'Exact retrospective geometry on55 fixed translations and both signs in49 saved completed boxes. Every translate remains independent of the old27 subgroup because the public point is, but no new rank or prospective recovery is claimed. A bounded failure covers only these words and maps; longer translations and other charts remain untested. No point enumeration, chart search or parameter sweep.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','check']);a=p.parse_args()
    if a.stage=='prepare':prepare()
    else:
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=d:raise ArithmeticError('translation visibility replay differs')
        else:
            if OUT.exists():raise FileExistsError('preserve fixed translation audit')
            checkpoint(OUT,d)
        print('PUBLIC28 FIXED TRANSLATES',d['status_counts'],'BEST',d['best_observation'])
