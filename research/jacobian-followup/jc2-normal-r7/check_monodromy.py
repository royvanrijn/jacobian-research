"""Exhaustive transposition representations for saved, certified braid words."""
from itertools import product
import json
from pathlib import Path
from local_excess import transpositions,compose,excess

ROOT=Path(__file__).resolve().parent
lookup={p:i for i,p in enumerate(transpositions)}
conjugate=[[lookup[compose(a,compose(b,a))] for b in transpositions] for a in transpositions]


def act(values,word):
    result=list(values)
    for letter in word:
        i=abs(letter)-1;a,b=result[i],result[i+1]
        if letter>0:result[i],result[i+1]=b,conjugate[b][a]
        else:result[i],result[i+1]=conjugate[a][b],a
    return tuple(result)


def transitive(values):
    reached={0}
    while True:
        new=reached|{p[i] for j in values for p in [transpositions[j]] for i in reached}
        if new==reached:return len(reached)==6
        reached=new


def main():
    results=[]
    for source in sorted(ROOT.glob('braids-*.json')):
        data=json.loads(source.read_text());assert data['degree']==6
        words=sorted(data['words'],key=len);survivors=[]
        for tail in product(range(15),repeat=5):
            values=(0,)+tail
            if all(act(values,word)==values for word in words):survivors.append(values)
        count=sum(transitive(values) for values in survivors)
        result={'case':data['case'],'braids':len(words),'assignments':len(survivors),'transitive':count,
                'fixed_first_meridian':'(1,2), by simultaneous conjugacy',
                'boundary':'One fixed witness only, not the entire r7 coefficient locus.'}
        results.append(result);print(result,flush=True)
    assert results and any(r['case']=='generic' and r['assignments']==1 and r['transitive']==0 for r in results)
    (ROOT/'monodromy.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
