#!/usr/bin/env python3
"""Fixed good-prime scoring of twelve retained, previously size-gated carrier fibres."""
import argparse
from pathlib import Path
from math import log
from fractions import Fraction as Q
from hashlib import sha256
import json
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/retained-native19-scores-v1'
SOURCE=ART/'native_rank3_carrier_subgroups_v1.json';IMAGES=ART/'native_rank3_carrier_images_v1.json'
OUT=ART/'retained_native19_scores_v1.json';GP=Path('/usr/bin/gp')
PRIMES=[p for p in _primes_up_to(32749) if p>=5]
CHECKS=[p for p in PRIMES if p<=199 or p in (4099,16381,32719,32749)]

def read(p):return json.loads(p.read_text())
def hashed(p):return sha256(p.read_bytes()).hexdigest()

def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fixed retained scores')
    source=read(SOURCE);images={tuple(r['word']):r for r in read(IMAGES)['rows']};rows=[]
    assert len(source['rows'])==12 and source['status']=='PASS'
    for i,r in enumerate(source['rows']):
        image=images[tuple(r['word'])]
        assert r['curve']==image['curve'] and r['rank_lower_bound']==19
        assert all(Q(x).denominator==1 for x in r['curve']) and all(Q(x)==0 for x in r['curve'][:3])
        rows.append({'id':f'native19-{i:02}', 'word':r['word'],'parameter':r['parameter'],
                     'model':r['curve'],'model_coefficient_bits':r['model_coefficient_bits']})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained-native19-score-protocol.v1',
        'sources':{str(p.relative_to(ROOT)):hashed(p) for p in [SOURCE,IMAGES,Path(__file__).resolve(),ROOT/'elliptic-curves/cas/mod2_reduction_independence.py']},
        'gp_sha256':hashed(GP),'rows':rows,'primes':PRIMES,'direct_check_primes':CHECKS,
        'seconds_per_curve':30,'rss_bytes':536870912,'maximum_workers':1,
        'score':'Sum round(1e12*(2-a_p)*log(p)/(p+1-a_p)) over good minimal reductions at all primes5through32749. Repeated exact p4/p6 short-model scaling is performed before deciding good or bad reduction. Sort by descending integer score then coefficient bits then original index.',
        'scope':'Exactly twelve already retained, constructed native rank19 fibres, all beyond the earlier400-bit point-budget gate. No new carrier words or parameter scan. Equations only in the scoring worker. This separate cost/score assessment preserves the old height-gated experiment and makes no rank prediction. Validation primes65537through131071 are not computed. No point search or automatic subsequent campaign.'})
    print('FROZEN12 retained native fibres /',len(PRIMES),'training primes')

def protocol():
    p=read(D/'protocol.json')
    assert p['gp_sha256']==hashed(GP) and all(hashed(ROOT/n)==h for n,h in p['sources'].items())
    return p

def program(row):
    A,B=row['model'][3:]
    return f'A={A};B={B};gettime();forprime(p=5,32749,a=A;b=B;k=0;while(a%p^4==0&&b%p^6==0,a/=p^4;b/=p^6;k++);aa=Mod(a,p);bb=Mod(b,p);if(4*aa^3+27*bb^2==0,print("B|",p,"|",k),print("T|",p,"|",k,"|",p+1-ellcard(ellinit([aa,bb])))));print("MS|",gettime());print("DONE");quit\n'

def reduction(row,p):
    a,b=map(int,row['model'][3:]);k=0
    while a%p**4==0 and b%p**6==0:a//=p**4;b//=p**6;k+=1
    a%=p;b%=p
    return a,b,k,(4*a**3+27*b*b)%p==0

def direct(a,b,p):
    squares=bytearray(p)
    for x in range(1,p):squares[x*x%p]=1
    return -sum(0 if (y:=(x*x*x+a*x+b)%p)==0 else (1 if squares[y] else -1) for x in range(p))

def evaluate(row,create):
    rawpath=D/row['id']/'raw.json'
    if create:
        if rawpath.exists():raise FileExistsError('preserve scalar trace attempt')
        r=capture([str(GP),'-q','-s','256000000'],input_text=program(row),
            limits=Limits(30,536870912),log_path=rawpath.with_suffix('.log'),separate_stderr=True,check=False)
        checkpoint(rawpath,{'program':program(row),'stdout':r.stdout,'stderr':r.stderr,'supervision':r.supervision})
    raw=read(rawpath);s=raw['supervision']
    assert s['outcome']=='completed' and s['returncode']==0 and not raw['stderr']
    assert raw['program']==program(row)
    lines=raw['stdout'].splitlines();assert len(lines)==len(PRIMES)+2 and lines[-1]=='DONE' and lines[-2].startswith('MS|')
    score=good=0;traces=[];checks=[];scaled=[]
    for p,line in zip(PRIMES,lines):
        a,b,k,bad=reduction(row,p);v=line.split('|')
        assert v[:3]==[('B' if bad else 'T'),str(p),str(k)] and len(v)==(3 if bad else 4)
        trace=None if bad else int(v[3])
        if trace is not None:
            assert trace*trace<=4*p;good+=1;score+=round((2-trace)*log(p)/(p+1-trace)*10**12)
        if p in CHECKS:
            value=None if bad else direct(a,b,p);assert value==trace;checks.append([p,value])
        if k:scaled.append({'prime':p,'scalings':k,'good_after_scaling':not bad})
        traces.append([p,trace])
    return {**row,'score_units':score,'good_primes':good,'scaled_reductions':scaled,'trace_count':len(traces),
        'traces':traces,'direct_checks':checks,'cpu_ms':int(lines[-2][3:]),'wall_seconds':s['wall_seconds'],
        'raw_sha256':hashed(rawpath)}

def run(check=False):
    p=protocol()
    if not check and OUT.exists():raise FileExistsError('preserve retained native scoring')
    rows=[]
    for row in p['rows']:
        rows.append(evaluate(row,not check))
        if not check:checkpoint(D/'checkpoint.json',{'status':'RUNNING','rows':rows})
        print(row['id'],rows[-1]['score_units'],row['model_coefficient_bits'],'bits',flush=True)
    result={'schema':'elliptic-curves.retained-native19-scores.v1','status':'PASS',
        'protocol_sha256':hashed(D/'protocol.json'),'rows':rows,
        'ordering':[r['id'] for r in sorted(rows,key=lambda r:(-r['score_units'],r['model_coefficient_bits'],r['id']))],
        'direct_character_sum_checks':sum(len(r['direct_checks']) for r in rows),
        'total_worker_seconds':sum(r['wall_seconds'] for r in rows),'following_campaign':None,'scope':p['scope']}
    if check:assert result==read(OUT)
    else:checkpoint(OUT,result);checkpoint(D/'checkpoint.json',{'status':'PASS','result_sha256':hashed(OUT)})
    print('PASS12',result['ordering'],result['total_worker_seconds'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','run','check']);a=p.parse_args()
    if a.mode=='freeze':freeze()
    else:run(a.mode=='check')
