#!/usr/bin/env python3
"""Larger ECM curves within the original600-second per-curve allocation."""
from concurrent.futures import ThreadPoolExecutor
import json
import ecm_submitted627_630 as worker

if __name__ == '__main__':
    supplement = {'ids':[628,630], 'prior_attempts':275,
        'additional_B1':[3000000]*100+[11000000]*50,
        'total_seconds_per_curve_including_prior_stage':600,
        'reason':'Stage1 exhausted its275-curve schedule after54/71 seconds, '
                 'leaving most of the original600-second allocation unused.'}
    protocol=worker.WORK/'ecm_continuation_protocol.json'
    if protocol.exists():
        if json.loads(protocol.read_text()) != supplement:
            raise ArithmeticError('continuation protocol differs')
    else:
        protocol.write_text(json.dumps(supplement,indent=2)+'\n')
    worker.SCHEDULE += supplement['additional_B1']
    for identifier in supplement['ids']:
        directory=worker.WORK/f'ecm_{identifier}'
        historical=directory/'state_stage1.json'
        if not historical.exists():
            state=json.loads((directory/'state.json').read_text())
            if state['completed_attempts']!=275 or state['status']!='BUDGET_EXHAUSTED':
                raise ArithmeticError('unexpected stage1 boundary')
            historical.write_text(json.dumps(state,indent=2)+'\n')
            state['status']='RUNNING'
            (directory/'state.json').write_text(json.dumps(state,indent=2)+'\n')
    rows=[r for r in json.loads((worker.WORK/'protocol.json').read_text())['roster']
          if r['icarm_id'] in supplement['ids']]
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(worker.run,rows))
