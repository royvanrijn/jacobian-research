#!/usr/bin/env python3
"""Exact finite Hilbert-series exploration; no moment-threshold claim by itself."""
from pathlib import Path
import hashlib,json


def main():
    p=Path(__file__).resolve().parent;limit=140
    weights=list(range(-3,4))+[-1,0,1]*2
    spaces=[{} for _ in range(limit+1)];spaces[0][0]=1
    for weight in weights:
        for degree in range(1,limit+1):
            for old_weight,count in spaces[degree-1].items():
                target=old_weight+weight
                spaces[degree][target]=spaces[degree].get(target,0)+count
    so=[row.get(0,0)-row.get(1,0) for row in spaces]
    assert all(value>=0 for value in so)
    orthogonal=[value if degree%2==0 else 0 for degree,value in enumerate(so)]
    candidate=orthogonal[:]
    for degree in range(2,21,2):
        candidate=[value-(candidate[j-degree] if j>=degree else 0) for j,value in enumerate(candidate)]
    negatives=[{'degree':degree,'coefficient':value} for degree,value in enumerate(candidate) if value<0]
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'representation':'harmonic cubic spin3 plus two vector spin1 representations',
            'weights':weights,'through_degree':limit,'SO3_dimensions':so,'O3_dimensions':orthogonal,
            'candidate_parameter_degrees':list(range(2,21,2)),
            'candidate_Hilbert_numerator_through_limit':candidate,'negative_coefficients':negatives,
            'scope':'Finite exact weight-count arithmetic only. Its use for a lower moment-prefix bound requires a written odd pure-null classification, invariant-ring dimension/finiteness argument and homogeneous parameter-system implication. No such threshold is inferred here.'}
    (p/'odd_cubic_hilbert_exploration.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'completed':True,'first_SO3_dimensions':so[:13],'first_O3_dimensions':orthogonal[:13],'first_negative_candidate_coefficients':negatives[:8],'nonzero_candidate_tail':[(i,x) for i,x in enumerate(candidate) if x and i>110][:8]}))


if __name__=='__main__':main()
