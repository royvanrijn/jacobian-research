#!/usr/bin/env sage-python
"""A finite signed height-4096 follow-up after prospective compact-R17 gains.

Reuses the existing rational table scanner and the frozen pointed worker.
No public curve/point database participates in selection. Eight denominator
shards retain separate outputs; every shard is bounded by the common supervisor.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import gcd,log
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
compact=SourceFileLoader('compact_r17_fixed',str(CAS/'compact_r17_prospective.sage')).load_module()
SCANNER=CAS/'newfamily/scan_rational_nagao_tables.cpp'
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-r17-wide-v1'


def export_tables(directory,sign):
    path=directory/f'tables-{sign}.txt'
    if path.exists():return path
    records=sorted((compact.read(p) for p in (compact.DIRECTORY/'trace-tables').glob('*.json')),key=lambda x:x['input']['prime'])
    if len(records)!=len(list(compact.prime_range(5,4094))): raise ArithmeticError('incomplete exact prime tables')
    with path.open('x') as out:
        out.write('RATIONAL_NAGAO_LOCAL_TABLE_V1\nF COMPACT_R17 8 12\n')
        for label,band in [('D',[r for r in records if r['input']['prime']<=997]),('H',[r for r in records if r['input']['prime']>997])]:
            out.write(f'B {label} {len(band)}\n')
            for r in band:
                p=r['input']['prime'];out.write(f'P {p}\n')
                if r['input']['model_sha256']!=compact.hashed(compact.MODEL):raise ArithmeticError('trace equation changed')
                for i in range(p+1):
                    ix=p if i==p else sign*i%p
                    good=r['good'][ix];trace=r['traces'][ix]
                    units=round((2-trace)/(p+1-trace)*log(p)*10**12) if good else 0
                    out.write(f'{int(good)} {trace} {units}\n')
        out.write('END\n')
    return path


def parse(text,sign):
    rows=[];summary=None
    for line in text.splitlines():
        parts=line.split()
        if not parts:continue
        if parts[0]=='C':
            _,n,d,a,b,g,h=parts;n=sign*int(n);d=int(d)
            rows.append({'numerator':n,'denominator':d,'parameter':str(Fraction(n,d)),
                'prefix_units':int(Fraction(a)*10**12),'extension_units':int(Fraction(b)*10**12),
                'prefix_good':int(g),'extension_good':int(h)})
        if parts[0]=='S':summary=list(map(int,parts[1:]))
    if summary is None:raise ArithmeticError('scanner did not finish')
    return rows,summary


def prepare(directory):
    directory.mkdir(parents=True,exist_ok=True)
    binary=directory/'scanner'
    if not binary.exists():capture(['g++','-O3','-std=c++17',str(SCANNER),'-o',str(binary)],limits=Limits(60,1073741824),log_path=directory/'compile.log')
    tables={sign:export_tables(directory,sign) for sign in (-1,1)}
    protocol={**compact.read(compact.DIRECTORY/'protocol.json'),
        'schema':'elliptic-curves.compact-r17-wide-protocol.v1','height':4096,'first_survivors':1024,
        'finalists':16,'sources':compact.provenance(),
        'selection_source':{str(Path(__file__).resolve().relative_to(ROOT)):compact.hashed(Path(__file__)),
            str(SCANNER.relative_to(ROOT)):compact.hashed(SCANNER)},
        'scanner_binary_sha256':compact.hashed(binary),
        'compiler_version':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
        'table_sha256':{str(sign):compact.hashed(path) for sign,path in tables.items()},
        'denominator_shards':4,'max_parallel_scanners':2,'seconds_per_shard':180,
        'first_stage_order':'decreasing integer prefix score; decreasing good-prime count; denominator; signed numerator',
        'second_stage_order':'decreasing prefix+extension integer score; height; denominator; signed numerator',
        'gate':'Earlier independent compact-fibre point search reached rank >=24 on t=33/119; fresh fibres have record-scale coefficient sizes.',
        'excluded_from_parameter_box':['zero','infinity'],
        'worker_wall_limit_seconds':900,'claim_boundary':'Finite prospective trial, not a density theorem. Rank >=28 near-record target; >=32 record target.'}
    pp=directory/'protocol.json'
    if pp.exists() and compact.read(pp)!=protocol:raise ArithmeticError('wide protocol changed')
    checkpoint(pp,protocol)
    # Cheap independent small-box check before the larger execution.
    for sign in (-1,1):
        result=capture([str(binary),str(tables[sign]),'12','12','200'],limits=Limits(20,268435456),log_path=directory/f'check-{sign}.log')
        rows,summary=parse(result.stdout,sign)
        if len(rows)!=sum(gcd(n,d)==1 for n in range(1,13) for d in range(1,13)):
            raise ArithmeticError('primitive small-box coverage differs')
        records=[compact.read(p) for p in (compact.DIRECTORY/'trace-tables').glob('*.json')]
        for row in rows:
            sums=[0,0]
            for r in records:
                p=r['input']['prime'];ix=row['numerator']*pow(row['denominator'],-1,p)%p if row['denominator']%p else p
                a=r['traces'][ix]
                if r['good'][ix]:sums[int(p>997)]+=round((2-a)/(p+1-a)*log(p)*10**12)
            if sums!=[row['prefix_units'],row['extension_units']]:raise ArithmeticError('C++ score mismatch')
    def shard(task):
        sign,k=task;path=directory/f'scan-{sign}-{k}.json'
        if path.exists():
            row=compact.read(path)
            if row['protocol_hash']!=compact.identity(protocol):raise ArithmeticError('shard protocol changed')
            return row
        result=capture([str(binary),str(tables[sign]),'4096','4096','1024',str(k),'4'],
            limits=Limits(180,536870912),log_path=directory/f'scan-{sign}-{k}.log')
        rows,summary=parse(result.stdout,sign)
        record={'protocol_hash':compact.identity(protocol),'sign':sign,'shard':k,'rows':rows,'summary':summary,'supervision':result.supervision}
        checkpoint(path,record);print('SHARD',sign,k,'population',summary[3],flush=True);return record
    with ThreadPoolExecutor(max_workers=2) as pool:shards=list(pool.map(shard,[(sign,k) for sign in (-1,1) for k in range(4)]))
    rows=[r for s in shards for r in s['rows']]
    if len(set(r['parameter'] for r in rows))!=len(rows):raise ArithmeticError('duplicate signed shard candidate')
    rows.sort(key=lambda r:(-r['prefix_units'],-r['prefix_good'],r['denominator'],r['numerator']))
    rows=rows[:1024]
    rows.sort(key=lambda r:(-r['prefix_units']-r['extension_units'],max(abs(r['numerator']),r['denominator']),r['denominator'],r['numerator']))
    population={'protocol_hash':compact.identity(protocol),'candidate_count':sum(s['summary'][3] for s in shards),
        'retained_candidates':rows,'finalists':rows[:16],
        'public_points_or_record_equations_used_for_selection':False,
        'shard_hashes':{f'{s["sign"]}-{s["shard"]}':compact.identity(s) for s in shards}}
    checkpoint(directory/'population.json',population)
    print('WIDE FROZEN',[(r['parameter'],(r['prefix_units']+r['extension_units'])/10**12) for r in rows[:16]],flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run'])
    p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--index',type=int,default=0);a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    else:
        protocol=compact.read(a.directory/'protocol.json');pop=compact.read(a.directory/'population.json')
        if compact.identity(protocol)!=pop['protocol_hash']:raise ArithmeticError('population binding changed')
        if protocol['selection_source'][str(Path(__file__).resolve().relative_to(ROOT))]!=compact.hashed(Path(__file__)):
            raise ArithmeticError('selection source changed')
        compact.run_fibre(a.directory,f'candidate-{a.index:02}',pop['finalists'][a.index]['parameter'])


if __name__=='__main__':main()
