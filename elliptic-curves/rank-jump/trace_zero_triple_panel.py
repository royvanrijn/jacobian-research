#!/usr/bin/env python3
"""A rank-blind finite parity gate before any degree-two triple construction."""
import argparse
from itertools import combinations
from pathlib import Path
import retrospective as r

HERE=Path(__file__).resolve().parent
INPUT=r.OUT/'rank_jump_degree_one_relation_panel_inputs_v1.json'
PROTOCOL=HERE/'TRACE_ZERO_TRIPLE_PANEL_PROTOCOL.json'
OUTPUT=r.OUT/'rank_jump_trace_zero_triple_panel_v1.json'


def compute():
    inp=r.read(INPUT);rows=[];total=0
    masks={label:r.pack([int(x)%2 for x in c['trace']]) for label,c in inp['covers'].items()}
    for key,block in inp['blocks'].items():
        count=0;qualified=[]
        for labels in combinations(block['labels'],3):
            count+=1;total+=1
            if masks[labels[0]]^masks[labels[1]]^masks[labels[2]]==0:
                word=[sum(inp['covers'][label]['trace'][i] for label in labels)//2 for i in range(17)]
                qualified.append({'labels':list(labels),'generic_word':word})
        rows.append({'block_key':key,'co_split_triples':count,'trace_even_triples':qualified})
    assert total<=5000
    return {'schema':'rank-jump.trace-zero-triple-panel.v1','rows':rows,'triple_tests':total,
            'trace_even_count':sum(len(x['trace_even_triples']) for x in rows),
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
            'boundary':'This stage checks only integral trace parity, before any relation or point evaluation. A qualifying triple still needs a rational-solubility test.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();d=compute()
    if a.mode=='build':r.write_new(OUTPUT,d)
    else:assert r.read(OUTPUT)==d
    print('PASS',d['triple_tests'],'triple parities;',d['trace_even_count'],'qualify')
    for row in d['rows']:
        if row['trace_even_triples']:print(row['block_key'],len(row['trace_even_triples']))
