#!/usr/bin/env sage-python
"""Same signed height-4096 population, with no short-prime truncation.

Motivation: the preceding 997-prime cutoff drops all three published rank
25/26/27 controls inside this box. This separate experiment scores every
primitive parameter using every cached prime through 4093 before retention.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from importlib.machinery import SourceFileLoader
from math import log
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
wide=SourceFileLoader('fullscore_scanner_adapter',str(CAS/'compact_r17_wide.sage')).load_module()
compact=wide.compact
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-fullscore-v1'


def prepare(directory):
    directory.mkdir(parents=True,exist_ok=True)
    binary=wide.DIRECTORY/'scanner'
    records=sorted((compact.read(p) for p in (compact.DIRECTORY/'trace-tables').glob('*.json')),key=lambda r:r['input']['prime'])
    if len(records)!=len(list(compact.prime_range(5,4094))):raise ArithmeticError('incomplete prime list')
    tables={}
    for sign in (-1,1):
        path=directory/f'tables-{sign}.txt';tables[sign]=path
        if path.exists():continue
        with path.open('x') as out:
            out.write('RATIONAL_NAGAO_LOCAL_TABLE_V1\nF COMPACT_R17_FULLSCORE 8 12\n')
            # The scanner requires a nonempty H band; its duplicate p=5 value
            # is explicitly unused, not a disjoint validation score.
            for label,band in [('D',records),('H',records[:1])]:
                out.write(f'B {label} {len(band)}\n')
                for r in band:
                    p=r['input']['prime'];out.write(f'P {p}\n')
                    if r['input']['model_sha256']!=compact.hashed(compact.MODEL):raise ArithmeticError('trace source changed')
                    for i in range(p+1):
                        j=p if i==p else sign*i%p;a=r['traces'][j];good=r['good'][j]
                        score=round((2-a)/(p+1-a)*log(p)*10**12) if good else 0
                        out.write(f'{int(good)} {a} {score}\n')
            out.write('END\n')
    protocol={**compact.read(wide.DIRECTORY/'protocol.json'),
        'schema':'elliptic-curves.compact-r17-fullscore-protocol.v1','sources':compact.provenance(),
        'selection_source':{str(p.relative_to(ROOT)):compact.hashed(p) for p in
            (Path(__file__).resolve(),wide.SCANNER,CAS/'compact_r17_wide.sage')},
        'table_sha256':{str(sign):compact.hashed(p) for sign,p in tables.items()},
        'first_prime_bound':4093,'first_survivors':1024,'second_prime_bound':None,
        'score':'sum_good (2-a_p)*log(p)/(p+1-a_p), all 562 primes 5..4093 before any selection',
        'unused_scanner_H_band':'duplicate p=5, ignored; not validation',
        'gate':'Prior compact new curves certified ranks >=24 and >=23. Retrospective cutoff audit excludes all three in-box published rank-25/26/27 controls before extension.',
        'first_stage_order':'full-prime integer score descending; good count descending; denominator; signed numerator',
        'second_stage_order':None}
    pp=directory/'protocol.json'
    if pp.exists() and compact.read(pp)!=protocol:raise ArithmeticError('fullscore protocol changed')
    checkpoint(pp,protocol)
    def run(task):
        sign,k=task;path=directory/f'scan-{sign}-{k}.json'
        if path.exists():
            r=compact.read(path)
            if r['protocol_hash']!=compact.identity(protocol):raise ArithmeticError('shard identity changed')
            return r
        result=capture([str(binary),str(tables[sign]),'4096','4096','1024',str(k),'4'],
            limits=Limits(180,536870912),log_path=directory/f'scan-{sign}-{k}.log')
        rows,summary=wide.parse(result.stdout,sign)
        for r in rows:r['score_units']=r.pop('prefix_units');r.pop('extension_units');r.pop('extension_good')
        r={'protocol_hash':compact.identity(protocol),'sign':sign,'shard':k,'rows':rows,'summary':summary,'supervision':result.supervision}
        checkpoint(path,r);print('FULLSCORE SHARD',sign,k,summary[3],flush=True);return r
    with ThreadPoolExecutor(max_workers=2) as pool:shards=list(pool.map(run,[(sign,k) for sign in (-1,1) for k in range(4)]))
    rows=[r for s in shards for r in s['rows']]
    if len(set(r['parameter'] for r in rows))!=len(rows):raise ArithmeticError('shard collision')
    rows.sort(key=lambda r:(-r['score_units'],-r['prefix_good'],r['denominator'],r['numerator']))
    checkpoint(directory/'population.json',{'protocol_hash':compact.identity(protocol),
        'candidate_count':sum(s['summary'][3] for s in shards),'retained_candidates':rows[:1024],
        'finalists':rows[:16],'public_points_or_record_equations_used_for_selection':False,
        'shard_hashes':{f'{s["sign"]}-{s["shard"]}':compact.identity(s) for s in shards}})
    print('FULLSCORE FROZEN',[(r['parameter'],r['score_units']/10**12) for r in rows[:16]],flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run'])
    p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--index',type=int,default=0);a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    else:
        protocol=compact.read(a.directory/'protocol.json');pop=compact.read(a.directory/'population.json')
        if compact.identity(protocol)!=pop['protocol_hash']:raise ArithmeticError('population binding changed')
        for name,h in protocol['selection_source'].items():
            if compact.hashed(ROOT/name)!=h:raise ArithmeticError('selection source changed')
        compact.run_fibre(a.directory,f'candidate-{a.index:02}',pop['finalists'][a.index]['parameter'])


if __name__=='__main__':main()
