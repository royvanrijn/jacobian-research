"""Retrospective calibration only; never an execution dependency."""
import json
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.stats import spearmanr
from research_runtime.store import checkpoint

ROOT=Path(__file__).resolve().parents[2]
V2=ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v2/calibration302'
D=ROOT/'artifacts/local/elliptic-curves/v2-terminal-retrospective-v1'
def read(p):return json.loads(p.read_text())
def bits(mapping):
    return sum(abs(Fraction(str(x)).numerator).bit_length()+Fraction(str(x)).denominator.bit_length() for x in mapping['discriminant_quartic'])
def ranks(values,x):return {'descending':1+sum(v>x for v in values),'ascending':1+sum(v<x for v in values),'ties':sum(v==x for v in values),'population':len(values)}
def correlations(rows,fields):
    return {a+' / '+b:None if len(set(r[a] for r in rows))<2 or len(set(r[b] for r in rows))<2 else float(spearmanr([r[a] for r in rows],[r[b] for r in rows]).statistic)
            for i,a in enumerate(fields) for b in fields[i+1:]}
def main():
    s=read(V2/'epoch-13/selection.json');bank=[];prospective=[];wins=[]
    for ai,a in enumerate(s['anchors']):
        scores=np.load(V2/f'epoch-13/anchor-{ai:02d}-full.npz')['babai_norms']
        for c in read(D/f'retained-bank/anchor-{ai:02d}.json')['candidates']:
            w=c['witnesses'][0]
            bank.append({'anchor':ai,'extension':c['extension'],'babai':int(scores[c['extension']]),'exact':int(c['centre_CVP']['norm']),
                         'multiplicity':len(c['centre_CVP']['minima'])//2,'quartic_bits':bits(w['mapping']),'height':int(w['height']),
                         'finite':int(w['coordinate'][1])!=0})
        for c in a['refined']:
            prospective.append({'babai':c['babai_norm'],'exact':c['metric_norm'],'multiplicity':len(c['cvp']['minima'])//2})
    for wd in sorted(V2.glob('epoch-*')):
        stage=read(wd/'stage.json')
        if stage['after']==stage['before']:continue
        c=read(wd/f"chart-{stage['charts']-1:03d}.json")['centre'];sel=read(wd/'selection.json')
        matches=[(i,a) for i,a in enumerate(sel['anchors']) if a['anchor']['orbit']==c['orbit']]
        row={'before':stage['before'],'after':stage['after'],'lane':c['lane'],'chart_number':stage['charts']}
        if matches:
            i,a=matches[0];ext=sum((int(v)%2)<<j for j,v in enumerate(c['representative'][17:]));scores=np.load(wd/f'anchor-{i:02d}-full.npz')['babai_norms']
            row['babai_rank']=ranks(list(map(int,scores)),int(scores[ext]))
            exact=[r['metric_norm'] for r in a['refined']]
            row['exact_rank_in_V2_shortlist']=ranks(exact,c['metric_norm'])
        wins.append(row)
    target=next(r for r in bank if r['anchor']==16 and r['extension']==1832)
    scores=list(map(int,np.load(V2/'epoch-13/anchor-16-full.npz')['babai_norms']))
    result={'status':'RETROSPECTIVE_CALIBRATION_NOT_EXECUTION_INPUT','target':target,
            'target_babai_full_anchor':ranks(scores,target['babai']),
            'target_ranks_oracle_enriched_bank':{f:ranks([r[f] for r in bank],target[f]) for f in ('babai','exact','multiplicity','quartic_bits')},
            'oracle_enriched_bank_population':len(bank),'oracle_enriched_correlations':correlations(bank,['babai','exact','multiplicity','quartic_bits','height']),
            'V2_prospective_top16_population':len(prospective),'V2_prospective_top16_correlations':correlations(prospective,['babai','exact','multiplicity']),
            'earlier_winners':wins,
            'limitations':['Bank is target-conditioned nearest-word sample, not unbiased population; correlations do not estimate landscape-wide correlations.',
                          'V2 prospective exact metrics exist only after top16 truncation; no full-space exact or quartic rank is asserted.',
                          'Earlier selected winners are subject to selection bias. Their success cannot test excluded mediocre-score alternatives.',
                          'Chart height is retrospective; never a prospective score. No global anti-correlation claim follows.']}
    checkpoint(D/'v3-prefreeze-metrics.json',result)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
