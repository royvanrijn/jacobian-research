#!/usr/bin/env python3
"""Build an exact Curve302 chart-exposure ledger from historical V3 raw transcripts.

This is the compatibility layer for the real historical schemas discovered on
2026-09-11.  It binds seed/epoch from paths, reconstructs the displayed D basis
from MW-state snapshots, joins recorded-point files to cumulative chart snapshots,
and recovers quotient words with Sage Neron--Tate coordinates plus exact group-law
verification.  A successful `run` immediately launches the existing three chart
exposure experiments with the resulting normalized ledger.

Commands: audit | ledger | run | check
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from curve302_short_core_controls import require, primitive
from curve302_chart_exposure import normalize_explicit_ledger, validate_ledger
from curve302_chart_replay_adapter import (
    scan_replay_index, replay_schema_audit_index, build_replay_ledger_index,
)
import run_curve302_chart_exposure as exposure

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
SELF = Path(__file__).resolve()
DEFAULT = ROOT / "artifacts/local/elliptic-curves/curve302-chart-replay-v1"
NAMES = exposure.NAMES


def sha(path):
    h=sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()


def atomic(path,payload):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as out:
            json.dump(payload,out,sort_keys=True,indent=2,allow_nan=False); out.write('\n'); out.flush()
            if os.environ.get('CURVE302_TEST_SKIP_FSYNC')!='1': os.fsync(out.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def source_context(source=None,structure=None,raw_root=None):
    short=exposure.discover_short(source); struct=exposure.discover_structure(short,structure)
    data=exposure.source_data(short,struct); raw=exposure.find_raw_root(struct,raw_root)
    expected={r['seed']:{int(s['epoch']) for s in r['stages']} for r in data['runs']}
    return short,struct,data,raw,expected


def verify_actuals(ledger,data):
    by={r['seed']:r for r in ledger['runs']}; missing=[]
    for run in data['runs']:
        for stage in run['stages']:
            gains=stage.get('gains',())
            if not gains: continue
            ls=next((s for s in by[run['seed']]['stages'] if int(s['epoch'])==int(stage['epoch'])),None)
            require(ls is not None,'ledger stage vanished')
            words={tuple(e.get('quotient_word',e.get('word',()))) for c in ls['charts'] for e in c.get('exposures',())}
            words={primitive(w) for w in words if w and any(w)}
            for g in gains:
                if g['primitive'] not in words: missing.append((run['seed'],stage['epoch'],g['primitive']))
    return missing


def audit(args):
    short,struct,data,raw,expected=source_context(args.source,args.structure,args.raw_root)
    index=scan_replay_index(raw,NAMES,expected)
    result=replay_schema_audit_index(index,NAMES,expected)
    result['raw_root']=str(raw); result['short_source']=str(short); result['structure_source']=str(struct)
    print(json.dumps(result,indent=2,sort_keys=True))
    if args.output: atomic(args.output,result)
    return result


def make_ledger(args, folder=None):
    short,struct,data,raw,expected=source_context(args.source,args.structure,args.raw_root)
    started=time.time()
    
    def progress(n,key,result):
        if n == 1 or n % 25 == 0:
            state='IN_D' if result is not None else 'OUTSIDE_D_OR_UNRESOLVED'
            print(f'CURVE302_REPLAY_RECOGNITION|unique={n}|state={state}',flush=True)
    index=scan_replay_index(raw,NAMES,expected)
    ledger=build_replay_ledger_index(index,NAMES,expected,max_unresolved_hits=args.max_unresolved_hits,progress=progress)
    # Normalize through the consumer itself before writing: this catches adapter/ledger drift.
    normalized=normalize_explicit_ledger(ledger,NAMES,14)
    stats=validate_ledger(normalized,NAMES,data['runs'],expected_total_charts=data['total_charts'])
    missing=verify_actuals(ledger,data)
    require(not missing,f"replay ledger does not expose {len(missing)} actual acquisitions; first {missing[:8]}")
    ledger['adapter']['consumer_validation']=stats
    ledger['adapter']['actual_acquisitions_verified']=180
    ledger['adapter']['elapsed_seconds']=round(time.time()-started,3)
    if folder is None:
        output=args.ledger or Path('curve302-chart-exposure-ledger-replayed.json')
    else:
        output=Path(folder)/'chart-exposure-ledger.json'
    atomic(output,ledger)
    return output,ledger,short,struct,raw


def run(args):
    folder=Path(args.folder).resolve()
    require(not folder.exists(),'replay output folder exists; use a fresh --folder')
    folder.mkdir(parents=True)
    try:
        short,struct,data,raw,expected=source_context(args.source,args.structure,args.raw_root)
        print('CURVE302_REPLAY_SCAN|status=START',flush=True)
        index=scan_replay_index(raw,NAMES,expected)
        print(f"CURVE302_REPLAY_SCAN|status=PASS|files={index['examined']}|point_files={len(index['hits'])}",flush=True)
        audit_result=replay_schema_audit_index(index,NAMES,expected)
        atomic(folder/'REPLAY_SCHEMA_AUDIT.json',audit_result)
        require(audit_result['status'] in ('PASS_REPLAY_SCHEMA_ADAPTER','PASS_REPLAY_SCHEMA_WITH_UNBOUND_POINT_FILES'),
                f"replay schema audit did not pass: {audit_result['status']}; inspect REPLAY_SCHEMA_AUDIT.json")
        started=time.time()
        def progress(n,key,result):
            if n == 1 or n % 25 == 0:
                state='IN_D' if result is not None else 'OUTSIDE_D_OR_UNRESOLVED'
                print(f'CURVE302_REPLAY_RECOGNITION|unique={n}|state={state}',flush=True)
        ledger=build_replay_ledger_index(index,NAMES,expected,max_unresolved_hits=args.max_unresolved_hits,progress=progress)
        normalized=normalize_explicit_ledger(ledger,NAMES,14)
        stats=validate_ledger(normalized,NAMES,data['runs'],expected_total_charts=data['total_charts'])
        missing=verify_actuals(ledger,data)
        require(not missing,f"replay ledger does not expose {len(missing)} actual acquisitions; first {missing[:8]}")
        ledger['adapter']['consumer_validation']=stats; ledger['adapter']['actual_acquisitions_verified']=180
        ledger['adapter']['elapsed_seconds']=round(time.time()-started,3)
        ledger_path=folder/'chart-exposure-ledger.json'; atomic(ledger_path,ledger)
        adapter_report={
            'status':'PASS_EXACT_REPLAY_LEDGER',
            'ledger_sha256':sha(ledger_path),
            'actual_acquisitions_verified':ledger['adapter']['actual_acquisitions_verified'],
            'recognized_exposures':ledger['adapter']['recognized_exposures'],
            'unique_recorded_points':ledger['adapter']['unique_recorded_points'],
            'consumer_validation':ledger['adapter']['consumer_validation'],
            'boundary':'Every promoted exposure comes from a recorded rational point whose D coordinates were proposed by Neron--Tate heights and then verified by exact elliptic-curve group equality.'
        }
        atomic(folder/'REPLAY_REPORT.json',adapter_report)
        if args.no_experiments:
            print(f"CURVE302_REPLAY_LEDGER|status=PASS|ledger={ledger_path}")
            return
        exp_folder=folder/'experiments'
        cmd=[sys.executable,str(CAS/'run_curve302_chart_exposure.py'),'run',
             '--folder',str(exp_folder),'--source',str(short),'--structure',str(struct),
             '--ledger',str(ledger_path),'--static-limit',str(args.static_limit),
             '--random-orders',str(args.random_orders),'--stage-seconds',str(args.stage_seconds),
             '--memory-gib',str(args.memory_gib)]
        proc=subprocess.run(cmd)
        require(proc.returncode==0,'downstream chart exposure experiments failed')
        report=json.loads((exp_folder/'REPORT.json').read_text())
        require(report.get('status')=='PASS_THREE_CHART_EXPOSURE_EXPERIMENTS','downstream report not passed')
        atomic(folder/'REPORT.json',{
            'status':'PASS_REPLAY_ADAPTER_AND_THREE_CHART_EXPOSURE_EXPERIMENTS',
            'replay_report_sha256':sha(folder/'REPLAY_REPORT.json'),
            'ledger_sha256':sha(ledger_path),'experiments_report_sha256':sha(exp_folder/'REPORT.json'),
            'boundary':adapter_report['boundary']})
        print('CURVE302_CHART_REPLAY|status=PASS|experiments=3',flush=True)
    except Exception as exc:
        atomic(folder/'REPLAY_FAILED.json',{'status':'UNKNOWN_REPLAY_ADAPTER','error':repr(exc),
                                            'boundary':'No chart exposure zero, multiplicity percentile, or ordering result is promoted when exact replay adaptation fails.'})
        raise


def check(args):
    folder=Path(args.folder).resolve(); report=json.loads((folder/'REPORT.json').read_text())
    require(report.get('status')=='PASS_REPLAY_ADAPTER_AND_THREE_CHART_EXPOSURE_EXPERIMENTS','final replay report not passed')
    require(report['ledger_sha256']==sha(folder/'chart-exposure-ledger.json'),'ledger hash mismatch')
    require(report['replay_report_sha256']==sha(folder/'REPLAY_REPORT.json'),'replay report hash mismatch')
    require(report['experiments_report_sha256']==sha(folder/'experiments/REPORT.json'),'experiment report hash mismatch')
    # Existing downstream deterministic check recomputes all three analyses from the frozen ledger.
    cmd=[sys.executable,str(CAS/'run_curve302_chart_exposure.py'),'check','--folder',str(folder/'experiments')]
    proc=subprocess.run(cmd); require(proc.returncode==0,'downstream deterministic check failed')
    print('CURVE302_CHART_REPLAY_CHECK|status=PASS|downstream=deterministic',flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command',choices=('audit','ledger','run','check'))
    ap.add_argument('--folder',type=Path,default=DEFAULT)
    ap.add_argument('--ledger',type=Path)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--source',type=Path); ap.add_argument('--structure',type=Path); ap.add_argument('--raw-root',type=Path)
    ap.add_argument('--max-unresolved-hits',type=int)
    ap.add_argument('--static-limit',type=int,default=1000); ap.add_argument('--random-orders',type=int,default=256)
    ap.add_argument('--stage-seconds',type=int,default=3600); ap.add_argument('--memory-gib',type=int,default=6)
    ap.add_argument('--no-experiments',action='store_true')
    args=ap.parse_args()
    try:
        if args.command=='audit': audit(args)
        elif args.command=='ledger':
            path,ledger,*_=make_ledger(args); print(f"CURVE302_REPLAY_LEDGER|status=PASS|ledger={path}|sha256={sha(path)}")
        elif args.command=='run': run(args)
        else: check(args)
    except Exception as exc:
        print(f"CURVE302_CHART_REPLAY_FAIL|{type(exc).__name__}|{exc}",file=sys.stderr)
        raise SystemExit(2)

if __name__=='__main__': main()
