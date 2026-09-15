"""Bounded k1 norm preview from the complete retained rational/quadratic tables."""
import itertools
import json
from pathlib import Path
import resource
import time
resource.setrlimit(resource.RLIMIT_CPU,(20,25))
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
root=Path(__file__).resolve().parents[4]
tables=root/'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
start=time.process_time()
rat=json.loads((tables/'rational-fibres.json').read_text())['rows']
quad=json.loads((tables/'quadratic-fibres.json').read_text())['rows']
agreement={}
for r in rat:
    for e,c in zip(r['rational_roots'],r['rational_codes']):
        assert c not in agreement and c
        agreement[c]=(r['t'],e)
split=[r for r in rat if len(r['rational_roots'])==3]
triple_rows=[]
# Oriented point choices: positions 0,1 for the two sections, 2 for the unused root.
for rows in itertools.combinations(split,3):
    for ordered in itertools.product(*[list(itertools.permutations(zip(r['rational_roots'],r['rational_codes']))) for r in rows]):
        codes=[0,0,0]
        for pair in ordered:
            for j,(_,c) in enumerate(pair):codes[j]^=c
        if codes[2]!=0:continue
        assert codes[0]==codes[1]
        if codes[0] not in agreement:continue
        a=agreement[codes[0]]
        if a[0] in [r['t'] for r in rows]:continue
        triple_rows.append({'agreement':a,'disagreement_bases':[r['t'] for r in rows],
                            'choices':[[e for e,c in rr] for rr in ordered],'codes':codes})
mixed=[]
for r in split:
    for s in quad:
        if len(s['roots'])!=3:continue
        for rr in itertools.permutations(zip(r['rational_roots'],r['rational_codes'])):
            for ss in set(itertools.permutations(zip(s['roots'],s['norm_codes']))):
                codes=[rr[j][1]^ss[j][1] for j in range(3)]
                if codes[2]!=0:continue
                assert codes[0]==codes[1]
                if codes[0] not in agreement:continue
                a=agreement[codes[0]]
                if a[0]==r['t']:continue
                mixed.append({'agreement':a,'rational_disagreement_base':r['t'],
                              'quadratic_disagreement_base':s['t'],'quadratic_smooth':s['smooth'],
                              'choices_rational':[e for e,c in rr],'choices_quadratic':[e for e,c in ss],
                              'codes':codes})
print(json.dumps({'rational_agreement_points':len(agreement),'split_fibres':len(split),
                  'three_distinct_rational':triple_rows,'rational_plus_quadratic':mixed,
                  'cpu_seconds':time.process_time()-start},indent=2))
