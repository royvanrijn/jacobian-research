#!/usr/bin/env sage-python
"""Verified application: replay the first historical and V3 302 unlocks.

Consumes completed witnesses only. Independently constructs finite quotient
codes at the fixed certificate primes <=197. No rational point searches.
"""
import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from sage.all import QQ, GF, EllipticCurve, PolynomialRing, matrix, gcd, lcm

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL=ROOT/'artifacts/local/elliptic-curves'
sys.path.insert(0,str(CAS))
from search_observability import prepare_chart, point_visibility
OUT=ART/'det1092_initial_unlock_construction_audit_v1.json'
PATHS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): PATHS.append(p);return json.loads(p.read_text())

def finite_certificate(E,points):
    primes=[17,47,53,61,67,71,79,83,89,101,107,113,127,137,149,179,191,197]
    code=[];records=[]
    for p in primes:
        field=GF(p);ep=EllipticCurve(field,E.a_invariants())
        assert field(E.discriminant())!=0
        key=lambda P:tuple(int(v) for v in P)
        allpoints=list(ep.points());doubles={key(2*P):2*P for P in allpoints}
        labels={v:0 for v in doubles}; representatives=[ep(0)];dim=0
        for P in allpoints:
            if key(P) in labels:continue
            new=[]
            for i,rep in enumerate(representatives):
                other=P+rep;new.append(other)
                for D in doubles.values():labels[key(other+D)]=(1<<dim)|i
            representatives+=new;dim+=1
        assert len(labels)==len(allpoints) and len(representatives)==2**dim
        def reduce(P):
            den=lcm([v.denominator() for v in P]);word=[int(v*den) for v in P]
            g=int(gcd(word));word=[v//g for v in word]
            return ep([field(v) for v in word])
        cols=[labels[key(reduce(P))] for P in points]
        rows=[[(v>>j)&1 for v in cols] for j in range(dim)]
        code+=rows;records.append({'prime':p,'group_order':len(allpoints),'double_order':len(doubles),'rows':rows})
    rank=matrix(GF(2),code).rank()
    assert rank==len(points)
    R=PolynomialRing(GF(31),'x');x=R.gen();A,B=E.a4(),E.a6()
    assert not (x**3+GF(31)(A)*x+GF(31)(B)).roots()
    return {'method':'independent exhaustive finite E(Fp)/2E(Fp) construction','rank':int(rank),'no_rational_2_torsion_prime':31,'signatures':records,'argument':'A relation is divisible by 2 from the injective reduction code; no rational 2-torsion permits repeated division, hence no nonzero integral relation.'}

def witness(label,search,word,points,target,extra):
    curve,(_, _),centre,_,shift,M,coeff=prepare_chart(search)
    E=EllipticCurve(QQ,list(map(QQ,curve)));basis=[E(list(map(QQ,P))) for P in points];P=E(list(map(QQ,target)))
    assert len(basis)==17 and shift==0
    C=sum((n*Q for n,Q in zip(word,basis)),E(0))
    assert (C[0],C[1])==tuple(map(QQ,centre))
    # This helper's coverage fields describe GMP denominator shards. Use only
    # its exact geometry and check the retained PARI hit/budget below.
    located=point_visibility({k:v for k,v in search.items() if k!='height_bound'},tuple(map(Fraction,target)))
    n,d=map(int,located['coordinate']);root=int(located['square_root_absolute'])
    assert max(abs(n),d)<=search['height_bound']
    assert [str(n),str(d),str(root)] in search['primitive_square_hits']
    a,b,c,e=map(QQ,M);v0=c*n+e*d;slope=(a*n+b*d)/v0
    chart=search['pointed_chart']
    raw_factor=QQ(chart['point_denominator_root'])/(QQ(chart['curve_coordinate_scale'])**2*QQ(search['ordinate_scale']))
    root_raw=raw_factor*root/v0**2
    mapped=[]
    for sign in (1,-1):
        x=(slope*slope-C[0]+sign*root_raw)/2
        y=slope*(x-C[0])-C[1]
        mapped.append(E([x,y]))
    assert P in mapped
    proof=finite_certificate(E,[*basis,P])
    return {'classification':'verified application','label':label,**extra,'input_rank':17,'output_rank':18,
        'centre_word_in_generic17':word,'centre_short_model':[str(C[0]),str(C[1])],
        'chart_coordinate':str(QQ(n)/d),'homogeneous_square_witness':[str(n),str(d),str(root)],
        'quartic_coefficients_low_to_high':list(map(str,coeff)),
        'map':{'raw_slope_matrix':list(map(str,M)),'raw_ordinate_factor':str(raw_factor),'slope_formula':'ell=(a*n+b*d)/(c*n+e*d)','raw_square_root_formula':'v=raw_ordinate_factor*r/(c*n+e*d)^2','point_formula':'X=(ell^2-Xcentre +/- v)/2; Y=ell*(X-Xcentre)-Ycentre','selected_sign':1 if mapped[0]==P else -1,'back_to_literal_302':'x=(X-15)/36; y=(Y/108-x-1)/2'},
        'point_short_model':list(map(str,P[:2])),
        'point_literal302':[str((P[0]-15)/36),str((P[1]/108-(P[0]-15)/36-1)/2)],
        'independence':proof,
        'generic_only_before_gain':True,
        'new_after_gain':'The new point and its parity-extension cosets. The winning centre and complete rational quartic map already used only generic sections.'}

def build():
    historical=LOCAL/'curve302-focused-point-exposure-v2/curve302-generic17'
    old=read(historical/'result.json');seed=read(historical/'seed.json')
    first=next(row for row in old['charts'] if row['rank_lower_bound']>17)
    assert first['index']==7
    target=first['admission_observations'][0]['point'];word=first['search']['input']['centre']['coefficients']
    historic=witness('historical first 302 gain',first['search'],word,seed['points'],target,{'chart_index_zero_based':7,'whole_frozen_wave':'17 ->19 after 49 completed charts; first gain 17 ->18 on chart 8'})
    root=LOCAL/'adaptive-visibility-cascade-v3/replay-M17';stage=read(root/'stages.json')[0]
    first=read(root/'epoch-00/chart-016.json');audit=read(root/'epoch-00/mod2-016.json');selection=read(root/'epoch-00/selection.json')
    assert stage['before']==17 and stage['after']==18 and stage['audit_sha256']==sha(root/'epoch-00/mod2-016.json')
    for i in range(16):assert read(root/f'epoch-00/mod2-{i:03d}.json')['rank_lower_bound']==17
    v3=witness('completed autonomous V3 first 302 gain',first['search'],first['centre']['representative'],selection['basis'],audit['independent_points'][17],{'chart_index_zero_based':16,'orbit':first['centre']['orbit'],'shell':first['centre']['shell'],'lane':first['centre']['lane'],'full_cosets_scored':stage['full_cosets_scored']})
    control=LOCAL/'v3-transfer-11952-v3/control-native11952'
    verified=read(control/'verified.json')
    for name,digest in verified['bindings'].items():assert sha(ROOT/name)==digest
    assert verified['status']=='PASS_INDEPENDENT_TRANSFER_REPLAY' and verified['gain']==0 and verified['charts']==82
    control_selection=read(control/'replay-M17/epoch-00/selection.json')
    cs=read(control/'replay-M17/stages.json')[0]
    return {'schema':'elliptic-curves.det1092-initial-unlock.v1','classification':'verified application','status':'PASS_EXACT_INITIAL_UNLOCK_AND_UNCHANGED_CONTROL_COMPARISON',
        'historical':historic,'autonomous_V3':v3,
        'native11952_completed_control':{'classification':'verified application','certificate':str((control/'verified.json').relative_to(ROOT)),'initial_rank':17,'final_certified_rank':17,'completed_charts':82,'full_cosets_scored':cs['full_cosets_scored'],'anchors':len(control_selection['anchors']),'conclusion':'A completed bounded bootstrap miss. With no gain, extension dimension remains zero and later V3 diversification is never exercised. Neither arithmetic absence nor a failure of multisection existence follows.'},
        'comparison_boundary':'The historical witness uses an earlier frozen policy. The V3 witness is compared with the completed generic-start native11952 V3 control; neither search is rerun or tuned. Different fibres and parents prevent a controlled incidence-versus-visibility attribution.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [*PATHS,Path(__file__),CAS/'search_observability.py']},'point_searches_run':0}

if __name__=='__main__':
    payload=build();text=json.dumps(payload,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print('PASS_EXACT_INITIAL_UNLOCK_AND_UNCHANGED_CONTROL_COMPARISON',flush=True)
