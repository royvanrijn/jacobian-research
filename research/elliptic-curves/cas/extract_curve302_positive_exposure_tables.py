#!/usr/bin/env python3
"""Extract four positive-evidence-only tables from a completed Curve302 replay.

Usage:
  python3 research/elliptic-curves/cas/extract_curve302_positive_exposure_tables.py run [--replay DIR] [--output DIR]
  python3 research/elliptic-curves/cas/extract_curve302_positive_exposure_tables.py check [--replay DIR] [--output DIR]
"""
from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile

from curve302_short_core_controls import require
from curve302_chart_exposure import normalize_explicit_ledger, validate_ledger
from curve302_positive_exposure_tables import extract_all
import run_curve302_chart_exposure as exposure

CAS=Path(__file__).resolve().parent
ROOT=CAS.parents[1]
DEFAULT_REPLAY=ROOT/'artifacts/local/elliptic-curves/curve302-chart-replay-v1'


def read(path): return json.loads(Path(path).read_text())
def sha(path): return sha256(Path(path).read_bytes()).hexdigest()

def atomic(path,obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:
            json.dump(obj,f,sort_keys=True,indent=2,allow_nan=False); f.write('\n')
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def load(replay):
    replay=Path(replay).resolve(); exp=replay/'experiments'
    rr=read(replay/'REPORT.json'); require(rr.get('status')=='PASS_REPLAY_ADAPTER_AND_THREE_CHART_EXPOSURE_EXPERIMENTS','replay report not passed')
    er=read(exp/'REPORT.json'); require(er.get('status')=='PASS_THREE_CHART_EXPOSURE_EXPERIMENTS','experiment report not passed')
    short=exp/'inputs/short'; structure=exp/'inputs/structure'
    data=exposure.source_data(short,structure)
    ledger_path=exp/'inputs/chart-exposure-ledger.json'; raw_ledger=read(ledger_path)
    # The downstream experiment snapshots the consumer-normalized ledger.  Its
    # exposure entries already use the canonical ``word`` field, so passing it
    # through the raw-ledger adapter again would discard every exposure.
    ledger=(raw_ledger if raw_ledger.get('normalization')=='explicit-ledger'
            else normalize_explicit_ledger(raw_ledger,exposure.NAMES,14))
    stats=validate_ledger(ledger,exposure.NAMES,data['runs'],expected_total_charts=data['total_charts'])
    return replay,exp,data,ledger,ledger_path,stats


def md_table(rows, columns):
    def text(v):
        if v is None: return '—'
        if isinstance(v,bool): return 'yes' if v else 'no'
        return str(v)
    out=['| '+' | '.join(label for _,label in columns)+' |','|'+'|'.join('---' for _ in columns)+'|']
    for r in rows: out.append('| '+' | '.join(text(r.get(k)) for k,_ in columns)+' |')
    return '\n'.join(out)


def summary_md(all_tables,stats):
    lead=all_tables['lead_times']['summary']; mult=all_tables['precontainment_multiplicity']['summary']
    choice=all_tables['coexposure_choices']['summary']; sat=all_tables['saturation_impact']['summary']
    lines=['# Curve302 positive-evidence tables','',
           'All quantities below use recorded positive chart exposures only. Missing hits are never interpreted as non-exposure.','',
           f"Ledger: **{stats['charts']} charts**, **{stats['exposures']} normalized positive exposures**.",'',
           '## 1. Exposure lead times','',
           md_table(lead,[('target','Target'),('runs','Runs'),('contained_runs','Contained'),('direct_preexposed_runs','Direct +'),('saturation_enabling_preexposed_runs','Enabling +'),('progress_preexposed_runs','Progress +'),('median_direct_lead_stages','Median direct lead'),('median_enabling_lead_stages','Median enabling lead'),('median_progress_lead_stages','Median progress lead')]),'',
           '## 2. Pre-containment multiplicity lower bounds','',
           md_table(mult,[('target','Target'),('runs','Runs'),('direct_chart_lower_bound_total','Direct charts ≥'),('progress_chart_lower_bound_total','Progress charts ≥'),('completion_enabling_chart_lower_bound_total','Completion charts ≥'),('runs_with_direct_positive','Runs direct +'),('runs_with_progress_positive','Runs progress +'),('runs_with_completion_enabling_positive','Runs completion +')]),'',
           '## 3. Positive-positive co-exposure choice sets','',
           md_table([choice],[('gain_stages','Gain stages'),('stages_with_positive_alternatives','With positive alternatives'),('stages_actual_hits_earliest_positive','Actual hits earliest +'),('stages_alternative_precedes_all_actual','Alternative earlier'),('stages_some_alternative_has_higher_recorded_multiplicity','Alternative higher multiplicity')]),'',
           '## 4. Saturation impact of positive alternatives','',
           md_table([sat],[('gain_stages','Gain stages'),('positive_candidate_stage_pairs','Positive candidate-stage pairs'),('stages_with_positive_alternatives','With alternatives'),('stages_with_any_positive_single_core_improver','Any core improver'),('stages_with_alternative_single_core_improver','Alternative improver'),('stages_actual_single_tied_for_best_positive','Actual tied best'),('stages_with_alternative_single_core_completer','Alternative completer')]),'',
           'Detailed per-run, per-stage and per-candidate records are in the four JSON files.','']
    return '\n'.join(lines)


def write_csv(path,rows):
    rows=list(rows)
    if not rows: Path(path).write_text(''); return
    keys=[]
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: json.dumps(v,sort_keys=True) if isinstance(v,(dict,list,tuple)) else v for k,v in r.items()})


def run(args):
    replay,exp,data,ledger,ledger_path,stats=load(args.replay)
    out=Path(args.output or replay/'positive-evidence-tables-v1').resolve()
    require(not out.exists(),'output exists; use a fresh --output or run check')
    out.mkdir(parents=True)
    tables=extract_all(data,ledger,14)
    files={
      'lead-times.json':tables['lead_times'],
      'precontainment-multiplicity.json':tables['precontainment_multiplicity'],
      'coexposure-choices.json':tables['coexposure_choices'],
      'saturation-impact.json':tables['saturation_impact'],
    }
    for name,obj in files.items(): atomic(out/name,obj)
    write_csv(out/'lead-times.csv',tables['lead_times']['runs'])
    write_csv(out/'precontainment-multiplicity.csv',tables['precontainment_multiplicity']['runs'])
    write_csv(out/'coexposure-choices.csv',tables['coexposure_choices']['stages'])
    write_csv(out/'saturation-impact-candidates.csv',tables['saturation_impact']['candidates'])
    (out/'SUMMARY.md').write_text(summary_md(tables,stats))
    report={'status':'PASS_FOUR_POSITIVE_EVIDENCE_TABLES','replay_report_sha256':sha(replay/'REPORT.json'),
            'frozen_ledger_sha256':sha(ledger_path),'ledger_stats':stats,
            'outputs':{name:sha(out/name) for name in files},'summary_sha256':sha(out/'SUMMARY.md'),
            'boundary':'Positive-evidence-only extraction. No missing chart hit is interpreted as non-exposure; no completeness assumption or chart-order counterfactual is introduced.'}
    atomic(out/'REPORT.json',report)
    print((out/'SUMMARY.md').read_text(),end='')
    print(f'CURVE302_POSITIVE_EVIDENCE|status=PASS|output={out}')


def check(args):
    replay,exp,data,ledger,ledger_path,stats=load(args.replay)
    out=Path(args.output or replay/'positive-evidence-tables-v1').resolve(); report=read(out/'REPORT.json')
    require(report.get('status')=='PASS_FOUR_POSITIVE_EVIDENCE_TABLES','report not passed')
    require(report['replay_report_sha256']==sha(replay/'REPORT.json'),'replay report changed')
    require(report['frozen_ledger_sha256']==sha(ledger_path),'frozen ledger changed')
    tables=extract_all(data,ledger,14)
    expected={'lead-times.json':tables['lead_times'],'precontainment-multiplicity.json':tables['precontainment_multiplicity'],
              'coexposure-choices.json':tables['coexposure_choices'],'saturation-impact.json':tables['saturation_impact']}
    for name,obj in expected.items():
        # JSON encodes tuple-valued lattice words as arrays.  Compare the exact
        # canonical bytes that ``atomic`` writes instead of Python tuple/list
        # container types after reloading the file.
        rendered=json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)+'\n'
        require((out/name).read_text()==rendered,f'{name} deterministic recomputation mismatch')
        require(report['outputs'][name]==sha(out/name),f'{name} hash mismatch')
    require((out/'SUMMARY.md').read_text()==summary_md(tables,stats),'SUMMARY.md deterministic mismatch')
    require(report['summary_sha256']==sha(out/'SUMMARY.md'),'summary hash mismatch')
    print('CURVE302_POSITIVE_EVIDENCE_CHECK|status=PASS')


def main():
    p=argparse.ArgumentParser(); p.add_argument('command',choices=('run','check')); p.add_argument('--replay',type=Path,default=DEFAULT_REPLAY); p.add_argument('--output',type=Path)
    a=p.parse_args(); run(a) if a.command=='run' else check(a)
if __name__=='__main__': main()
