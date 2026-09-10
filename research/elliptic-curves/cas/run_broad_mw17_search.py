#!/usr/bin/env python3
"""Broad matched prospective search across certified generic-MW17 parents.

This is a CPU-first experiment controller.  It reuses the exact class1 search
worker, but freezes an independent runtime per parent and applies the same
address window / control protocol to every parent.  It keeps three expensive
arms separate: Nagao-ranked, score-independent controls, and deterministic
feature-diversity candidates.

No historical exceptional parameters or points are read.  Search outcomes are
lower bounds only; a completed rank-17 fibre is not a rank upper bound.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import dataclass
import fcntl
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
DEFAULT = ROOT/'artifacts/local/elliptic-curves/broad-mw17-search-v1'

PRESETS = {
    'x1092-curve302': ROOT/'artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json',
    'x1092-class1': ROOT/'artifacts/generated-results/elliptic-curves/x1092_class1_arithmetic_gate_v1/parent.json',
}


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value, immutable=False):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(value, indent=2, sort_keys=True)+'\n'
    if immutable and path.exists():
        if path.read_text() != data:
            raise RuntimeError('immutable receipt changed: '+str(path))
        return
    tmp = path.with_suffix(path.suffix+'.tmp'); tmp.write_text(data); tmp.replace(path)


def slug(name):
    if not name or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_' for c in name):
        raise ValueError('parent names may contain only letters, digits, - and _')
    return name


def parse_parent(spec):
    if '=' not in spec:
        if spec not in PRESETS:
            raise ValueError('unknown preset '+spec+'; use NAME=/path/to/parent.json')
        return spec, PRESETS[spec]
    name, path = spec.split('=', 1)
    return slug(name), Path(path).expanduser().resolve()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def zscores(rows, key):
    vals = [float(r[key]) for r in rows]
    mean = sum(vals)/len(vals)
    var = sum((v-mean)**2 for v in vals)/max(1, len(vals)-1)
    sd = math.sqrt(var) or 1.0
    return {r['index']:(float(r[key])-mean)/sd for r in rows}


def choose_arms(rows, ranked_count, control_count, diversity_count):
    """Freeze disjoint ranked/control/diversity arms from one scored window."""
    controls = sorted((r for r in rows if r['control']),
                      key=lambda r:(r.get('control_order', 10**30), r['index']))[:control_count]
    control_ids = {r['index'] for r in controls}
    ranked_pool = [r for r in rows if r['index'] not in control_ids and not r['control'] and r['nonsingular']]
    ranked = sorted(ranked_pool, key=lambda r:(-r['score_units'], r['model_bits'], r['index']))[:ranked_count]
    used = control_ids | {r['index'] for r in ranked}

    pool = [r for r in rows if r['index'] not in used and not r['control'] and r['nonsingular']]
    # Diversity is deliberately independent of the Nagao ordering.  We retain
    # extrema in a standardized cheap-feature space, with a stable hash tie.
    keys = ('model_bits', 'discriminant_bits', 'smooth_prime_count')
    zs = {k:zscores(pool, k) for k in keys} if pool else {}
    def dkey(r):
        radius = sum(zs[k][r['index']]**2 for k in keys)
        tie = hashlib.sha256(('broad-mw17-diversity-v1/'+str(r['index'])).encode()).hexdigest()
        return (-radius, tie, r['index'])
    diversity = sorted(pool, key=dkey)[:diversity_count]
    selected = []
    for arm, group in [('ranked', ranked), ('control', controls), ('diversity', diversity)]:
        for r in group:
            x = dict(r); x['broad_arm'] = arm; selected.append(x)
    if len({r['index'] for r in selected}) != len(selected):
        raise AssertionError('arm overlap')
    if len(ranked) != ranked_count or len(controls) != control_count or len(diversity) != diversity_count:
        raise ValueError('window cannot satisfy requested arm quotas')
    return selected


def prepare(args):
    folder = args.folder.resolve(); folder.mkdir(parents=True, exist_ok=False)
    base = load_module(CAS/'run_class1_prospective_search.py', 'broad_prepare_base')
    parents = []
    seen = set()
    for spec in args.parent:
        name, path = parse_parent(spec)
        if name in seen: raise ValueError('duplicate parent name '+name)
        seen.add(name)
        if not path.exists(): raise FileNotFoundError(path)
        pdir = folder/'parents'/name
        # One worker here is only metadata inside the frozen child plan. Broad
        # concurrency is controlled globally by this script.
        ns = SimpleNamespace(folder=pdir, parent=path, parent_sha256=sha(path),
                             workers=1, window=args.window, offset=args.offset,
                             prime_bound=args.prime_bound)
        base.prepare(ns)
        plan = read(pdir/'plan.json')
        rt = Path(plan['root']); window = rt/'ordinary-search/window-000'; window.mkdir(parents=True, exist_ok=True)
        base.score(pdir, window, 0, args.offset, args.window)
        rows = read(window/'scores.json')['rows']
        selected = choose_arms(rows, args.ranked, args.controls, args.diversity)
        write(pdir/'broad-selection.json', {
            'schema':'broad-mw17.parent-selection.v1', 'parent':name,
            'ranked':args.ranked, 'controls':args.controls, 'diversity':args.diversity,
            'scores_sha256':sha(window/'scores.json'), 'rows':selected}, True)
        parents.append({'name':name, 'source':str(path), 'source_sha256':sha(path),
                        'folder':str(pdir), 'child_plan_sha256':sha(pdir/'plan.json'),
                        'selection_sha256':sha(pdir/'broad-selection.json')})
    plan = {
        'schema':'broad-mw17-search.v1', 'parents':parents, 'workers':args.workers,
        'window':args.window, 'offset':args.offset, 'prime_bound':args.prime_bound,
        'ranked_per_parent':args.ranked, 'controls_per_parent':args.controls,
        'diversity_per_parent':args.diversity,
        'matched_addresses':True,
        'historical_exceptional_inputs':False,
        'stop_semantics':'STOP prevents new dispatch; active bounded fibres drain.',
        'boundary':'Prospective search over certified MW17 parents. Results are certified lower bounds, never exact-rank or no-jump theorems.',
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'], cwd=ROOT, text=True).strip(),
        'controller_sha256':sha(Path(__file__))}
    write(folder/'broad-plan.json', plan, True)
    write(folder/'broad-manifest.json', {'plan_sha256':sha(folder/'broad-plan.json'),
          'controller_sha256':sha(Path(__file__)),
          'children':{p['name']:{'plan':p['child_plan_sha256'],'selection':p['selection_sha256']} for p in parents}}, True)
    print(json.dumps({'status':'PREPARED_NOT_RUNNING','folder':str(folder),
          'parents':[p['name'] for p in parents], 'expensive_fibres_per_parent':args.ranked+args.controls+args.diversity,
          'total_expensive_fibres':len(parents)*(args.ranked+args.controls+args.diversity)}), flush=True)


def guard(folder):
    plan = read(folder/'broad-plan.json'); manifest = read(folder/'broad-manifest.json')
    if sha(folder/'broad-plan.json') != manifest['plan_sha256'] or sha(Path(__file__)) != manifest['controller_sha256']:
        raise RuntimeError('broad campaign/controller changed')
    for p in plan['parents']:
        pdir = Path(p['folder'])
        if sha(pdir/'plan.json') != p['child_plan_sha256'] or sha(pdir/'broad-selection.json') != p['selection_sha256']:
            raise RuntimeError('child campaign changed: '+p['name'])
    return plan


def run(args):
    folder = args.folder.resolve(); plan = guard(folder)
    with (folder/'broad-controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        tasks = []
        modules = {}
        for pi, p in enumerate(plan['parents']):
            pdir = Path(p['folder']); child = read(pdir/'plan.json'); rt = Path(child['root'])
            frozen = rt/'elliptic-curves/cas/run_class1_prospective_search.py'
            modules[p['name']] = load_module(frozen, 'broad_child_'+str(pi))
            modules[p['name']].guard(pdir)
            window = rt/'ordinary-search/window-000'
            selection = read(pdir/'broad-selection.json')
            if sha(window/'scores.json') != selection['scores_sha256']:
                raise RuntimeError('score table changed: '+p['name'])
            for row in selection['rows']:
                terminal = window/'cases'/str(row['index'])/'terminal.json'
                if not terminal.exists(): tasks.append((p['name'], pdir, window, row))
        # Resume partially progressed high-rank trajectories first, then arm
        # interleave so one parent or arm cannot monopolize the machine.
        def key(task):
            name,pdir,window,row = task
            progress = window/'cases'/str(row['index'])/'progress.json'
            rank = read(progress).get('rank_lower_bound',17) if progress.exists() else 17
            arm_order = {'control':0,'diversity':1,'ranked':2}[row['broad_arm']]
            return (-rank, arm_order, row['index'], name)
        tasks.sort(key=key)
        active = {}
        def dispatch(pool, task):
            name,pdir,window,row = task
            # Child worker ignores broad_arm, but selection.json preserves it.
            return pool.submit(modules[name].fibre, pdir, window, row)
        with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
            while tasks or active:
                while tasks and len(active) < plan['workers'] and not (folder/'STOP').exists():
                    task = tasks.pop(0); active[dispatch(pool, task)] = (task[0], task[3]['index'], task[3]['broad_arm'])
                if not active: break
                done,_ = wait(active, timeout=10, return_when=FIRST_COMPLETED)
                for future in done:
                    name,index,arm = active.pop(future); result = future.result()
                    print('FIBRE',name,arm,index,result['status'],'rank',result.get('rank_lower_bound'),flush=True)
                status(folder, quiet=True, active=list(active.values()), remaining=len(tasks))
        status(folder, quiet=False, active=[], remaining=len(tasks))


def status(folder, quiet=False, active=None, remaining=None):
    folder = Path(folder).resolve(); plan = guard(folder)
    summary = {'schema':'broad-mw17-search.status.v1', 'parents':{},
               'active':active or [], 'remaining_not_dispatched':remaining,
               'stopped':(folder/'STOP').exists(), 'updated_at':time.time()}
    best = 0
    for p in plan['parents']:
        pdir=Path(p['folder']); child=read(pdir/'plan.json'); rt=Path(child['root']); window=rt/'ordinary-search/window-000'
        sel=read(pdir/'broad-selection.json')['rows']; terminals=[]
        for row in sel:
            t=window/'cases'/str(row['index'])/'terminal.json'
            if t.exists():
                x=read(t); x['broad_arm']=row['broad_arm']; terminals.append(x)
        arms={}
        for arm in ('ranked','control','diversity'):
            xs=[x for x in terminals if x['broad_arm']==arm]
            arms[arm]={'completed':len(xs),'best_rank_lower_bound':max([x.get('rank_lower_bound') or 0 for x in xs],default=0),
                       'ge20':sum((x.get('rank_lower_bound') or 0)>=20 for x in xs),
                       'cpu_seconds':sum(x.get('search_cpu_seconds',0) for x in xs)}
        parent_best=max([x.get('rank_lower_bound') or 0 for x in terminals],default=0); best=max(best,parent_best)
        summary['parents'][p['name']]={'completed':len(terminals),'planned':len(sel),'best_rank_lower_bound':parent_best,'arms':arms}
    summary['best_rank_lower_bound']=best
    write(folder/'BROAD_STATUS.json',summary)
    if not quiet: print(json.dumps(summary,indent=2,sort_keys=True),flush=True)
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','run','status'))
    p.add_argument('--folder',type=Path,default=DEFAULT)
    p.add_argument('--parent',action='append',default=[],help='preset name or NAME=/absolute/parent.json; repeatable')
    p.add_argument('--window',type=int,default=262144)
    p.add_argument('--offset',type=int,default=131072)
    p.add_argument('--prime-bound',type=int,default=997)
    p.add_argument('--ranked',type=int,default=256)
    p.add_argument('--controls',type=int,default=64)
    p.add_argument('--diversity',type=int,default=64)
    p.add_argument('--workers',type=int,default=4)
    a=p.parse_args(); a.folder=a.folder.resolve()
    if a.mode=='prepare':
        if not a.parent:a.parent=['x1092-curve302','x1092-class1']
        if a.window<8 or a.window%8 or min(a.ranked,a.controls,a.diversity)<0 or not (1<=a.workers<=32):raise ValueError('invalid campaign sizes')
        prepare(a)
    elif a.mode=='run':run(a)
    else:status(a.folder)

if __name__=='__main__':main()
