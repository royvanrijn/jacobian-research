#!/usr/bin/env python3
"""Run the frozen independent cells with bounded process concurrency.

Scheduling alone differs from the serial entry point; each cell has the same
frozen centres, chart budgets and strict worker limits. A file lock excludes
overlapping supervisors for the same output directory.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import json
from pathlib import Path
import sys

from research_runtime.supervisor import capture, Limits
from research_runtime.store import checkpoint
from run_mw18_centre_experiment import check_protocol, summarize


def run(protocol_path,directory,workers,verify=False):
    protocol=json.loads(protocol_path.read_text());check_protocol(protocol)
    if not 1<=workers<=3: raise ValueError('declared concurrency must be between one and three')
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/'supervisor.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        jobs=[(case,policy) for case in protocol['cases'] for policy in protocol['policies']]
        def execute(job):
            case,policy=job;output=directory/'cells'/f'{case}--{policy}.json'
            if not verify and output.exists() and json.loads(output.read_text())['status']=='COMPLETE':
                return {'case':case,'policy':policy,'status':'REUSED_COMPLETE'}
            command=[sys.executable,str(Path(__file__).with_name('run_mw18_centre_experiment.py')),
                     '--protocol',str(protocol_path),'--cell',case,'--policy',policy,'--output',str(output)]
            if verify:command.append('--verify')
            result=capture(command,limits=Limits(protocol['limits']['cell_wall_seconds'],protocol['limits']['rss_bytes']),
                log_path=directory/('replay-logs' if verify else 'logs')/f'{case}--{policy}.log',check=False)
            return {'case':case,'policy':policy,'returncode':result.returncode,'supervision':result.supervision}
        records=[]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures={pool.submit(execute,j):j for j in jobs}
            for future in as_completed(futures):
                records.append(future.result())
                checkpoint(directory/('replay-supervision.json' if verify else 'supervision.json'),
                    {'protocol_hash':protocol['protocol_hash'],'workers':workers,'cells':records})
                summary=summarize(protocol,directory)
                print(f"MW18_SUPERVISOR|cells={len(records)}/{len(jobs)}|{records[-1]['case']}|{records[-1]['policy']}|rc={records[-1].get('returncode',0)}",flush=True)
        if any(r.get('returncode',0)!=0 for r in records): raise RuntimeError('one or more cells failed; keep their checkpoints')
        if verify:
            from hashlib import sha256
            checkpoint(directory/'replay.json',{'status':'PASS_EXACT_CHARTS_AND_INCREMENTAL_INDEPENDENCE',
                'protocol_hash':protocol['protocol_hash'],'cell_count':len(summary['cells']),
                'chart_count':sum(r['chart_count'] for r in summary['cells']),
                'summary_sha256':sha256((directory/'summary.json').read_bytes()).hexdigest(),
                'enumeration_repeated':False})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--protocol',type=Path,required=True);p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--workers',type=int,default=3);p.add_argument('--verify',action='store_true')
    a=p.parse_args();run(a.protocol,a.directory,a.workers,a.verify)
