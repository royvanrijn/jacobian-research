#!/usr/bin/env python3
"""Retrospective generic-mask novelty in two completed known-seed cascades."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'artifacts/local/elliptic-curves/curve302-seeded-v3-amplifier-v1'


def audit(output):
    output.mkdir(exist_ok=False)
    bindings={}
    def read(path):
        raw=path.read_bytes()
        bindings[str(path.relative_to(ROOT))]=hashlib.sha256(raw).hexdigest()
        return json.loads(raw)
    records=[]
    for case in ('recovered-strict-02','recovered-strict-03'):
        verified=read(BASE/case/'seeded-verified.json')
        folder=BASE/case/'replay-M17'
        terminal=read(folder/'terminal.json')
        assert verified['terminal_sha256']==bindings[str((folder/'terminal.json').relative_to(ROOT))]
        assert verified['rank_lower_bound']==terminal['final_rank_lower_bound']==31
        seen=set();pivots={};rows=[]
        for stage in terminal['stages']:
            if stage['after']<=stage['before']:continue
            assert stage['after']==stage['before']+1
            chart=read(folder/f"epoch-{stage['epoch']:02d}"/f"chart-{stage['charts']-1:03d}.json")
            centre=chart['centre']
            assert chart['index']==stage['charts']-1 and len(centre['representative'])==stage['before']
            mask=sum((int(x)%2)<<i for i,x in enumerate(centre['representative'][:17]))
            residue=mask
            for k in sorted(pivots,reverse=True):
                if (residue>>k)&1:residue ^= pivots[k]
            rows.append(dict(before=stage['before'],after=stage['after'],generic_mask=mask,
                lane=centre['lane'],previously_productive=mask in seen,
                in_previous_pairwise_xors=mask in {a^b for a,b in itertools.combinations(seen,2)},
                in_previous_productive_span=residue==0))
            if residue:pivots[residue.bit_length()-1]=residue
            seen.add(mask)
        records.append(dict(case=case,rows=rows,distinct_generic_masks=len(seen),generic_mask_span_rank=len(pivots)))
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in bindings.items())
    (output/'report.json').write_text(json.dumps(dict(status='PASS_RETAINED_PARENT_NOVELTY_AUDIT',
        records=records,inputs=bindings,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        new_point_searches=0,claim_boundary='Retrospective chart-label audit against sealed outcomes, not a new arithmetic rank replay. Two cascades on one curve. Novel winning masks do not prove that other parents cannot find the same point or predict amplification elsewhere.'),indent=2)+'\n')
    print('PARENT_NOVELTY',[(r['case'],r['distinct_generic_masks'],r['generic_mask_span_rank']) for r in records])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    audit(parser.parse_args().output.resolve())
