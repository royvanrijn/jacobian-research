#!/usr/bin/env python3
"""Read-only arithmetic audit of completed intake, including omitted low bands.

No point/rank inputs and no change to the frozen 90-case selection. A low j
height is a cost/possible-conductor lead, never a high-rank prediction.
"""
import argparse
import heapq
import json
from pathlib import Path
import sqlite3
import zlib

import det1092_funnel as f
from v3_warm_support import atomic,bindings,read,require,sha


def main(run,out):
    p=read(run/'protocol.json');bindings(f.ROOT,p['sources']);bindings(f.ROOT,p['inputs'])
    replay=read(run/'intake-replay.json');require(replay['status']=='PASS_FULL_INTAKE_REPLAY','intake not independently replayed')
    protocol=dict(schema='det1092-low-height-intake-audit.v1',maximum_parameters=p['population'],
                  maximum_height_bin=6,retain_per_band=8,wall_seconds=300,rss_bytes=1024**3,
                  inputs={str(q.relative_to(f.ROOT)):sha(q) for q in [run/'protocol.json',run/'intake-replay.json',Path(__file__)]},
                  scope='Read existing equation-only blocks for omitted low-height bands and conic splitting. No point search, selector change or new population.')
    atomic(out/'protocol.json',protocol,immutable=True)
    db=sqlite3.connect('file:'+str(run/'intake.sqlite')+'?mode=ro',uri=True)
    state=json.loads(db.execute("SELECT value FROM meta WHERE key='state'").fetchone()[0])
    require(state['accepted']==p['population'] and state['chain']==replay['block_chain'],'intake/replay state mismatch')
    count=0;chain='0'*64;counts={};splits=0;heaps={};blocks=0
    for index,h,payload in db.execute('SELECT id,sha256,payload FROM blocks ORDER BY id'):
        require(index==blocks and f.digest(payload)==h,'broken intake block')
        chain=f.digest((chain+h).encode());blocks+=1
        for row in json.loads(zlib.decompress(payload))['rows']:
            count+=1
            if row[3]!='SMOOTH':continue
            splits+=row[9]=='SPLIT';band=row[6];counts[str(band)]=counts.get(str(band),0)+1
            if band>protocol['maximum_height_bin']:continue
            # Exact j-height first; frozen trace score and address only break ties.
            key=[-max(row[4],row[5]),row[7],-row[2],-row[1],-row[0]]
            heap=heaps.setdefault(band,[]);item=[key,row]
            if len(heap)<protocol['retain_per_band']:heapq.heappush(heap,item)
            elif key>heap[0][0]:heapq.heapreplace(heap,item)
        if blocks%128==0:
            atomic(out/'progress.json',dict(blocks=blocks,parameters=count,conic_splits=splits,heights=counts))
            print('LOW_HEIGHT_AUDIT',count,flush=True)
    db.close();require(count==state['accepted'] and chain==state['chain'],'population traversal differs')
    arithmetic=f.Arithmetic(p)
    rows=[arithmetic.candidate(row,'omitted_low_height_diagnostic') for band,heap in sorted(heaps.items()) for _,row in sorted(heap,reverse=True)]
    result=dict(status='PASS_FULL_EQUATION_HEIGHT_AUDIT',parameters=count,block_chain=chain,
                height_counts=counts,conic_splitting_fibres=splits,low_height_candidates=rows,
                protocol_sha256=sha(out/'protocol.json'),
                boundary='Equation-only cost leads outside the frozen height7..17 selector. No rank, conductor or seed-existence claim; no automatic point-search refill.')
    atomic(out/'result.json',result,immutable=True)
    print('LOW_HEIGHT_AUDIT_COMPLETE',count,'conic_splits',splits,'low_height_candidates',len(rows),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.run.resolve(),a.output.resolve())
