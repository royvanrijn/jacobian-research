#!/usr/bin/env python3
"""Extract exact conic-split fibres missed by the frozen trace-only shortlist."""
import argparse
import json
from pathlib import Path
import sqlite3
import zlib

import det1092_funnel as f
from v3_warm_support import atomic,bindings,read,require,sha


def main(run,audit,out):
    p=read(run/'protocol.json');bindings(f.ROOT,p['sources']);bindings(f.ROOT,p['inputs'])
    prior=read(audit/'result.json');require(prior['status']=='PASS_FULL_EQUATION_HEIGHT_AUDIT','missing full-data audit')
    replay=read(run/'intake-replay.json');require(prior['block_chain']==replay['block_chain'],'audit/replay chain differs')
    target=prior['conic_splitting_fibres'];require(0<=target<=100,'split count outside finite follow-up cap')
    protocol=dict(schema='det1092-exact-split-followup.v1',expected_split_fibres=target,
                  wall_seconds=90,rss_bytes=1024**3,
                  inputs={str(q.relative_to(f.ROOT)):sha(q) for q in [run/'protocol.json',run/'intake-replay.json',audit/'result.json',Path(__file__)]},
                  scope='Exact conic splitting, already computed by the mass scan, prioritizes a separate seed follow-up. The original90 cases and trace-score statistics are unchanged.')
    atomic(out/'protocol.json',protocol,immutable=True)
    arithmetic=f.Arithmetic(p);rows=[];blocks=[]
    db=sqlite3.connect('file:'+str(run/'intake.sqlite')+'?mode=ro',uri=True)
    for index,h,payload in db.execute('SELECT id,sha256,payload FROM blocks ORDER BY id'):
        require(f.digest(payload)==h,'block digest differs')
        for row in json.loads(zlib.decompress(payload))['rows']:
            if row[3]!='SMOOTH' or row[9]!='SPLIT':continue
            require(arithmetic.record(*row[:3])==row,'split equation/fingerprints fail exact replay')
            rows.append(arithmetic.candidate(row,'exact_conic_split_followup'));blocks.append(dict(index=index,sha256=h))
        if len(rows)==target:break
    db.close();require(len(rows)==target,'split count differs from full audit')
    original={r['id'] for r in read(run/'selection.json')['seed_inputs']}
    result=dict(status='PASS_EXACT_SPLIT_EXTRACTION',rows=rows,blocks=blocks,
                already_in_original_seed_selection=[r['id'] for r in rows if r['id'] in original],
                full_data_audit_sha256=sha(audit/'result.json'),
                boundary='Exact rational conic preimages, not yet independently certified extra directions. No new rank follows until the seed certificate passes.')
    atomic(out/'split-candidates.json',result,immutable=True)
    print('EXTRACTED_SPLIT_FIBRES',len(rows),flush=True)
    for r in rows:print(r['id'],r['parameter'],r['j_numerator_bits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--audit',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();main(a.run.resolve(),a.audit.resolve(),a.output.resolve())
