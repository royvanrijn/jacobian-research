#!/usr/bin/env python3
"""Finite, frozen, no-LLM broad rank search. Linux, Sage/PARI; launch is opt-in.

Default: six productive X948 presentations (256 ranked + 64 controls each),
two X1092 comparisons (64 + 16 each), all from 262144 fixed rational addresses.
Existing search/independence engines are reused, not replaced.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
from fractions import Fraction as F
import gzip
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

from broad_rank_policy import (NATIVE, PARENTS, addresses, require, select, score_tables,
    feature_rows, model_at, j_invariant, next_stage, key)
from broad_rank_runtime import read, sha, write, token, alive, guard, completed, normalize_parent, absorb

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
DEFAULT = ROOT/'artifacts/local/elliptic-curves/broad-rank-v1'
SOURCE_SUFFIXES = {'.py','.sage','.gp','.c','.cpp','.h','.sh'}


def environment(sage):
    env = {k:v for k,v in os.environ.items() if not any(s in k.upper() for s in ('API_KEY','ACCESS_TOKEN','AUTH_TOKEN'))}
    env.update(PATH=str(Path(sage).parent)+os.pathsep+env.get('PATH',''),
        OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
        PYTHONDONTWRITEBYTECODE='1')
    # Never inherit a PYTHONPATH pointing at a moving source checkout.
    env.pop('PYTHONPATH',None)
    return env


def prepare(args):
    prep_cpu=time.process_time(); prep_wall=time.monotonic(); prep_children=resource.getrusage(resource.RUSAGE_CHILDREN)
    require(sys.platform.startswith('linux'), 'Linux process-tree supervision is required')
    sage = shutil.which(args.sage)
    require(sage and Path('/usr/bin/gp').is_file(), 'Sage and /usr/bin/gp are required')
    require(1 <= args.workers <= 8 and 0 <= args.rounds <= 96, 'invalid worker/round cap')
    require(5 <= args.prime_bound <= 4093, 'prime bound must be 5..4093')
    require(args.window > max(args.ranked+args.controls,args.comparison_ranked+args.comparison_controls)
            and min(args.ranked,args.controls,args.comparison_ranked,args.comparison_controls) >= 0,
            'invalid population/selection sizes')
    require(args.ranked+args.controls > 0 and args.comparison_ranked+args.comparison_controls > 0, 'empty arm')
    next(addresses(args.offset,args.window))
    require(len(set(args.parents)) == len(args.parents), 'duplicate parent identifier')
    folder = args.folder.resolve()
    require(not folder.exists(), 'preserve existing campaign; use resume or a NEW folder')
    software = subprocess.check_output([sage,'-python','-c',
        'import json,sys,numpy,sage.version; from sage.all import pari; '
        'print(json.dumps(dict(python=sys.executable,sage=sage.version.version,pari=str(pari.version()),numpy=numpy.__version__)))'],
        cwd=ROOT,env=environment(sage),text=True,timeout=60)
    software = json.loads(software.strip().splitlines()[-1])
    folder.mkdir(parents=True)
    rt = folder/'runtime/research'
    source_hashes = {}
    def copy(path):
        path = Path(path); dest = rt/path.relative_to(ROOT)
        dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(path,dest)
        source_hashes[str(path.relative_to(ROOT))] = sha(path)
    for directory in (CAS,ROOT/'elliptic-curves/ecsearch',ROOT/'elkies-k3/scripts'):
        for path in sorted(directory.rglob('*')):
            if path.is_file() and path.suffix in SOURCE_SUFFIXES and '__pycache__' not in path.parts: copy(path)
    for path in (ROOT/'elliptic-curves').glob('*.py'): copy(path)
    records = []
    atlas_path = ART/'compact_six_r17_atlas_v1.json'
    parity_path = ART/'r17_exact_maximum_parity_classes_v1.json'
    if any(f in NATIVE for f in args.parents):
        # Drop retrospective compaction metadata, retaining the same family records
        # needed by the unmodified native arithmetic adapters.
        atlas = read(atlas_path); parity = read(parity_path)
        write(rt/atlas_path.relative_to(ROOT), {'families':atlas['families']})
        copy(parity_path)
        source_hashes[str(atlas_path.relative_to(ROOT))] = sha(atlas_path)
        families = {r['family']:r for r in atlas['families']}
        grams = {r['family']:r['gram'] for r in parity['families']}
    for parent_id in args.parents:
        if parent_id in NATIVE:
            parent = normalize_parent(families[parent_id],grams[parent_id],parent_id)
            backend = 'native'
        else:
            if parent_id == 'x1092-class1':
                path = ART/'x1092_class1_arithmetic_gate_v1/parent.json'
                require(sha(path) == '7c6ee40c46f5a1f3d1fc464b5685a0e4c77ed1f3347b67865d7c9f93b2acd862',
                        'class1 frozen parent changed; explicitly review the new certificate before use')
                data = read(path); gram = data['generic_height_gram']
                source_hashes[str(path.relative_to(ROOT))] = sha(path)
            else:
                path = ART/'det1092_reduced_parameter_chart_v1/reduced-parent.json'
                original = ART/'curve302_recovered_mw17_parent_v1.json'
                data = read(path); gram = read(original)['generic_height_gram']
                source_hashes[str(path.relative_to(ROOT))] = sha(path)
                source_hashes[str(original.relative_to(ROOT))] = sha(original)
            parent = normalize_parent(data,gram,parent_id); backend = 'generic'
        parent_path = rt/'broad-inputs/parents'/f'{parent_id}.json'
        write(parent_path,parent)
        records.append({'id':parent_id,'backend':backend,'path':str(parent_path.relative_to(rt)),
                        'sha256':sha(parent_path),'ranked':args.ranked if backend=='native' else args.comparison_ranked,
                        'controls':args.controls if backend=='native' else args.comparison_controls})
    # Exact-equation exclusions only, not point coordinates or rank-based selection.
    excluded = []; exclusion_inputs = {}
    inventory = ROOT/'elliptic-curves/data/research_curves/database.json'
    paths = ([inventory] if inventory.exists() else []) + list(args.exclude_equations)
    for path in paths:
        path = path.resolve(); payload = read(path)
        if 'models' in payload: models = payload['models']
        else: models = [r['ainvs'] for r in payload['curves']]
        for model in models:
            require(len(model)==5, 'invalid exclusion equation')
            excluded.append(list(map(str,map(F,model))))
        exclusion_inputs[str(path)] = sha(path)
    write(rt/'broad-inputs/exclusions.json',{'models':excluded,'inputs':exclusion_inputs})
    # Legacy adapter requires a catalogue; an empty one makes no useful novelty
    # claim. Broad reporting overrides it with WORLDWIDE_NOVELTY_NOT_CHECKED.
    write(rt/'broad-inputs/adapter-catalogue.json',{'count':0,'curves':[]})
    plan = {'schema':'broad-rank-search.v1','root':str(rt),'sage':sage,'software':software,
        'workers':args.workers,'parents':records,'offset':args.offset,'window':args.window,
        'prime_bound':args.prime_bound,'max_rounds':args.rounds,'height':125000,
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'source_inputs':source_hashes,'created_at':time.time(),
        'job_wall_seconds':10800,'rss_bytes':3*1024**3,'min_free_gib':args.min_free_gib,
        'consecutive_failure_limit':4,'strict_class_required':False,'model_calls':0,
        'first_stage':'Native complete 43/49-class seed search (<=98 calls) plus <=100 complement calls after ANY gain; X1092 <=198 full-space V3 calls.',
        'control_rule':'Fixed SHA256 ordering, same score-independent address roster per parent; no outcome-dependent refills.',
        'escalation':'Bounded round-robin stages; every M18/M19 can continue. A predeclared 1/8 of rank17 misses gets one rescue with a different bank.',
        'boundary':'Rank lower bounds and protocol-dependent visibility, not exact ranks or a causal comparison of fibrations. Same numerical t is not the same curve across parents.'}
    write(folder/'plan.json',plan)
    write(folder/'manifest.json',{'plan_sha256':sha(folder/'plan.json'),
        'files':{str(p.relative_to(rt)):sha(p) for p in sorted(rt.rglob('*')) if p.is_file()},
        'executables':{str(p):sha(p) for p in (Path(sage),Path('/usr/bin/gp'),Path(software['python']))}})
    children=resource.getrusage(resource.RUSAGE_CHILDREN)
    write(folder/'preparation-cost.json',{'cpu_seconds':time.process_time()-prep_cpu+children.ru_utime-prep_children.ru_utime+children.ru_stime-prep_children.ru_stime,
          'wall_seconds':time.monotonic()-prep_wall,'scope':'Source/input preparation and CAS version probe; separate from scoring/search.'})
    print(json.dumps({'status':'PREPARED_NOT_LAUNCHED','folder':str(folder),
        'address_rows':args.window*len(records),'selected_slots':sum(r['ranked']+r['controls'] for r in records)}))


def score_task(folder,job,parent_id):
    plan = read(folder/'plan.json'); rt = Path(plan['root'])
    p = next(p for p in plan['parents'] if p['id']==parent_id)
    parent = read(rt/p['path']); require(sha(rt/p['path'])==p['sha256'],'parent differs')
    tables = score_tables(parent,plan['prime_bound'])
    write(job/'tables.json',tables)
    rows = []
    target = job/'features.jsonl.gz'; tmp = job/'features.partial.gz'
    with tmp.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as output:
            for row in feature_rows(parent,tables,plan['offset'],plan['window']):
                output.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode())
                rows.append(row)
                if len(rows)%4096==0: write(job/'progress.json',{'scored':len(rows)},False)
        raw.flush(); os.fsync(raw.fileno())
    os.replace(tmp,target)
    selected = select(rows,p['ranked'],p['controls'])
    for r in selected:
        r.update(parent_id=parent_id,family=parent['family'],backend=p['backend'],
                 id='b-'+key(parent_id,r['index'])[:20],model=model_at(parent,r['parameter']))
    write(job/'selection.json',{'rows':selected})
    write(job/'result.json',{'status':'PASS_FIXED_SCORING','parent_id':parent_id,'population':len(rows),
          'selected':len(selected),'parent_sha256':p['sha256'],'tables_sha256':sha(job/'tables.json'),
          'feature_sha256':sha(target),'selection_sha256':sha(job/'selection.json'),
          'boundary':'All traces recoverable from parameter and sealed residue tables; omitted local singular terms are not rank exclusions.'})


def preflight_task(folder,job):
    from sage.all import matrix, QQ
    import parent_foundry_worker as generic
    # Import real downstream adapters from the frozen tree before commissioning.
    import high_rank_foundry_arithmetic, high_rank_foundry_job, run_fresh6_seed_confirmation_v2
    import run_complement_seed_v3, memory_rank_certificate
    plan = read(folder/'plan.json'); rt = Path(plan['root']); checks = []
    for p in plan['parents']:
        parent = read(rt/p['path']); gram = matrix(QQ,parent['generic_height_gram'])
        require(gram.is_symmetric() and gram.is_positive_definite(), 'invalid parent height Gram')
        pool = rt/'broad-pools'/p['id']; selected = read(pool/'selection.json')['rows']
        candidate = next((r for r in selected if r['nonsingular']),None)
        require(candidate is not None,'no smooth commissioning candidate')
        model,points = generic.specialize(parent,candidate['parameter'])
        require(list(map(str,model))==candidate['model'] and len(points)==17,'score/native equation mismatch')
        if p['backend']=='native':
            from high_rank_foundry_seed_gate import assess
            gate = assess(candidate)
        else:
            dest = job/p['id'];dest.mkdir(exist_ok=True)
            packet = generic.generic_seed(parent,candidate['parameter'],dest)
            gate = {'status':'PASS_GENERIC_SEED_GATE' if packet else 'UNRESOLVED_GENERIC_SPECIALIZATION'}
        checks.append({'parent_id':p['id'],'parameter':candidate['parameter'],'generic_gate':gate['status']})
    write(job/'result.json',{'status':'PASS_INTERFACE_PREFLIGHT','checks':checks,'point_calls':0,
        'boundary':'Checks real frozen imports, equations and sampled seed gates; not an end-to-end search smoke or a theorem of specialized rank17 for every address.'})


def supervised_task(folder,job):
    import ctypes
    from research_runtime.supervisor import run, Limits
    job.mkdir(parents=True,exist_ok=True)
    with (job/'driver.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (job/'seal.json').exists(): completed(job);return
        before = resource.getrusage(resource.RUSAGE_CHILDREN); own = time.process_time(); wall = time.monotonic()
        plan,rt = guard(folder); dispatch = read(job/'dispatch.json')
        require(dispatch['plan_sha256']==sha(folder/'plan.json'),'wrong dispatch plan')
        write(job/'driver-lease.json',{'pid':os.getpid(),'token':token(os.getpid())})
        require(ctypes.CDLL(None).prctl(36,1,0,0,0)==0,'could not become subreaper')
        command = dispatch['command']
        result = run(command,limits=Limits(dispatch['wall_seconds'],plan['rss_bytes']),
            log_path=job/'worker.log',checkpoint_path=job/'supervisor.json',cwd=rt,env=environment(plan['sage']))
        while True:
            try: os.waitpid(-1,0)
            except ChildProcessError: break
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        good = result['outcome']=='completed' and result['returncode']==0 and (job/'result.json').exists()
        # Include raw journals, clouds and arithmetic receipts, not only a PASS word.
        ignored = {'seal.json','driver.log','driver.lock','driver-lease.json'}
        files = {str(p.relative_to(job)):sha(p) for p in sorted(job.rglob('*'))
                 if p.is_file() and p.name not in ignored and not p.name.endswith('.tmp')}
        write(job/'seal.json',{'completed':good,'dispatch_sha256':sha(job/'dispatch.json'),'files':files,
            'outcome':result['outcome'],'returncode':result['returncode'],
            'cpu_seconds':after.ru_utime-before.ru_utime+after.ru_stime-before.ru_stime+time.process_time()-own,
            'wall_seconds':time.monotonic()-wall,'scope':'Worker tree plus driver CPU; scoring, search and preflight charged separately.'})


def execute(folder,job,command,seconds):
    plan = read(folder/'plan.json');rt = Path(plan['root'])
    job.mkdir(parents=True,exist_ok=True)
    write(job/'dispatch.json',{'command':command,'wall_seconds':seconds,'plan_sha256':sha(folder/'plan.json')})
    if (job/'seal.json').exists(): return completed(job)
    lease_path = job/'driver-lease.json'
    if (job/'dispatch-intent.json').exists():
        # Adopt a live bounded driver; never duplicate an unreceipted point call.
        for unused in range(10):
            if lease_path.exists() or (job/'seal.json').exists(): break
            time.sleep(.2)
        while lease_path.exists() and alive(read(lease_path)):
            time.sleep(1)
        return completed(job)  # None means interrupted UNKNOWN; no blind retry.
    write(job/'dispatch-intent.json',{'at':time.time()})
    with (job/'driver.log').open('ab',buffering=0) as log:
        p = subprocess.Popen([sys.executable,str(rt/'elliptic-curves/cas/run_broad_rank_search.py'),
            '_supervise','--folder',str(folder),'--job',str(job)],cwd=rt,env=environment(plan['sage']),
            stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        p.wait()
    return completed(job)


def build_queue(folder):
    plan = read(folder/'plan.json');rt = Path(plan['root'])
    path = folder/'queue.json'
    if path.exists():
        require(sha(path)==read(folder/'queue-seal.json')['sha256'],'queue changed')
        return read(path)['rows']
    from certify_compact_r17_candidates import isomorphic
    excluded = read(rt/'broad-inputs/exclusions.json')['models']
    known = {}
    for i,model in enumerate(excluded):
        j = j_invariant(model)
        if j: known.setdefault(j,[]).append((f'excluded-{i}',model))
    pools = []
    for parent in plan['parents']:
        job = rt/'broad-pools'/parent['id'];receipt = completed(job)
        require(receipt and receipt[1] and receipt[1]['status']=='PASS_FIXED_SCORING','score stage incomplete')
        pools.append(read(job/'selection.json')['rows'])
    ordered = [pool[i] for i in range(max(map(len,pools))) for pool in pools if i<len(pool)]
    for row in ordered:
        j = j_invariant(row['model'])
        if not j: row['intake_status']='SINGULAR_FIBRE';continue
        aliases = [(label,model) for label,model in known.get(j,[])
                   if isomorphic(tuple(map(F,row['model'])),tuple(map(F,model)))]
        if aliases:
            row.update(intake_status='EXACT_Q_ISOMORPH_ALIAS',alias_of=aliases[0][0])
        else:
            row['intake_status']='ACCEPTED';known.setdefault(j,[]).append((row['id'],row['model']))
    write(path,{'rows':ordered,'no_outcome_refills':True,'scope':'Exact rational-isomorphism deduplication against the frozen inventory/provided equation lists and this selection.'})
    write(folder/'queue-seal.json',{'sha256':sha(path)})
    return ordered


def request_for(plan,row,spec,previous,job):
    rt = Path(plan['root']); parent = next(p for p in plan['parents'] if p['id']==row['parent_id'])
    if row['backend']=='native':
        candidate = {**row,'rank':previous['rank'] if previous else 17}
        if previous: candidate.update(packet=previous['packet'],packet_sha256=previous['packet_sha256'],head=previous['head'])
        req = {'kind':spec['kind'],'candidate':candidate,'allowance':spec['allowance'],
            'bank_index':spec['bank_index'],'revival':spec['revival'],'use_cached':False,
            'catalogue':'broad-inputs/adapter-catalogue.json',
            'config':{'sage':plan['sage'],'phase_seconds':1800,'rss_bytes':plan['rss_bytes'],
                      'height':plan['height'],'point_seconds':10,'map_seconds':5}}
        command = [plan['sage'],'-python',str(rt/'elliptic-curves/cas/high_rank_foundry_job.py'),str(job)]
    else:
        req = {'parent':parent['path'],'parent_sha256':parent['sha256'],
            'parameter':row['parameter'],'allowance':spec['generic_allowance'],
            'height':plan['height'],'bank_index':spec['bank_index']}
        if previous: req.update(packet=previous['packet'],packet_sha256=previous['packet_sha256'])
        command = [plan['sage'],'-python',str(rt/'elliptic-curves/cas/parent_foundry_worker.py'),'--job',str(job)]
    return req,command


def advance(folder,row,spec,previous):
    plan = read(folder/'plan.json');rt = Path(plan['root']);iteration = previous['round']+1 if previous else 0
    case = rt/'broad-cases'/row['id'];job = case/f'batch-{iteration:03d}'
    terminal = job/'broad-state.json'
    req,command = request_for(plan,row,spec,previous,job)
    write(job/'request.json',req)
    if terminal.exists():
        require(sha(terminal)==read(job/'broad-state-seal.json')['sha256'],'derived state changed')
        completed(job)
        state = read(terminal)
        require(state['source_seal_sha256']==sha(job/'seal.json'),'state source changed')
        return state
    outcome = execute(folder,job,command,plan['job_wall_seconds'])
    if outcome is None:
        state = {**(previous or {'rank':None,'calls':0,'round':iteration}),
                 'status':'UNKNOWN_INTERRUPTED_DISPATCH','engineering_failure':True}
        # Do not seal interruption as mathematical evidence; resume reports it.
        write(case/'state.json',state,False);return state
    seal,result = outcome
    if result is not None and result.get('status') in ('PASS_CERTIFIED_SEARCH','PASS_CERTIFIED_PARENT_EVALUATION'):
        state = absorb(previous,result,row,spec,job,rt);state['engineering_failure']=False
    else:
        backend_status = result.get('status') if result else seal['outcome']
        generic_miss = backend_status in ('UNRESOLVED_GENERIC_SEED','UNRESOLVED_GENERIC_SPECIALIZATION')
        resource_miss = seal['outcome'] in ('strict_wall_timeout','strict_rss_limit') or (
            result is not None and result.get('outcome') in ('strict_wall_timeout','strict_rss_limit'))
        state = {**(previous or {'rank':None,'calls':0,'round':iteration}),
            'status':'UNKNOWN_'+str(backend_status),'engineering_failure':not(generic_miss or resource_miss)}
    state.update(parent_id=row['parent_id'],id=row['id'],arm=row['arm'],parameter=row['parameter'],
                 source_seal_sha256=sha(job/'seal.json'))
    # State is derived after the immutable worker seal, and is sealed separately.
    write(terminal,state)
    write(job/'broad-state-seal.json',{'sha256':sha(terminal)})
    write(case/'state.json',state,False)
    if state.get('rank') is not None and state['rank']>=27 and (previous is None or state['rank']>previous['rank']):
        write(folder/'events'/f'{row["id"]}-r{state["rank"]}.json',state)
    return state


def status(folder):
    plan = read(folder/'plan.json');rt = Path(plan['root']);groups = {}; pending = 0
    rows = read(folder/'queue.json')['rows'] if (folder/'queue.json').exists() else []
    for row in rows:
        group = groups.setdefault(row['parent_id']+'/'+row['arm'],{'selected':0,'aliases_or_singular':0,'finished':0,
            'unknown':0,'lower_bounds':{},'initial_lower_bounds':{},'point_calls':0,'search_cpu_seconds':0.0})
        group['selected']+=1
        if row['intake_status']!='ACCEPTED':group['aliases_or_singular']+=1;continue
        case = rt/'broad-cases'/row['id'];path=case/'state.json'
        if not path.exists():pending+=1;continue
        state=read(path)
        initial=case/'batch-000/broad-state.json'
        if initial.exists() and (value:=read(initial).get('rank')) is not None:
            k=str(value);group['initial_lower_bounds'][k]=group['initial_lower_bounds'].get(k,0)+1
        group['unknown']+=state['status'].startswith(('UNKNOWN','CENSORED'))
        group['finished']+=1;group['point_calls']+=state.get('calls',0)
        if state.get('rank') is not None:
            r=str(state['rank']);group['lower_bounds'][r]=group['lower_bounds'].get(r,0)+1
        group['search_cpu_seconds']+=sum(read(p)['cpu_seconds'] for p in case.glob('batch-*/seal.json'))
    score_cpu=sum(read(p)['cpu_seconds'] for p in (rt/'broad-pools').glob('*/seal.json'))
    active=[str(p.parent.relative_to(rt)) for p in rt.glob('broad-cases/*/batch-*/driver-lease.json')
            if alive(read(p)) and not (p.parent/'seal.json').exists()]
    result={'schema':'broad-rank-report.v1','groups':groups,'pending_first_stage':pending,
            'scoring_cpu_seconds':score_cpu,'active_jobs':active,
            'controller_cpu_seconds':sum(read(p)['cpu_seconds_this_invocation'] for p in (folder/'controller-cost-history').glob('*.json')) +
                (read(folder/'controller-cost.json')['cpu_seconds_this_invocation'] if (folder/'controller-cost.json').exists() else 0),
            'preparation_cost':read(folder/'preparation-cost.json') if (folder/'preparation-cost.json').exists() else None,
            'preflight_cost':read(rt/'broad-preflight/seal.json')['cpu_seconds'] if (rt/'broad-preflight/seal.json').exists() else None,
            'controller_alive':alive(read(folder/'controller.json')) if (folder/'controller.json').exists() else False,
            'stop_requested':(folder/'STOP').exists(),'worldwide_novelty':'NOT_CHECKED',
            'scope':'Search lower bounds, not true ranks. Initial endpoints are retained per batch-000; adaptive endpoints are budget-confounded. See seals for partial/failed work costs.'}
    write(folder/'REPORT.json',result,False)
    return result


def run(folder):
    with (folder/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if (folder/'controller-cost.json').exists():
            destination=folder/'controller-cost-history'/f'{time.time_ns()}.json'
            write(destination,read(folder/'controller-cost.json'))
        controller_cpu=time.process_time();controller_wall=time.monotonic()
        def controller_cost():
            write(folder/'controller-cost.json',{'cpu_seconds_this_invocation':time.process_time()-controller_cpu,
                'wall_seconds_this_invocation':time.monotonic()-controller_wall,
                'scope':'Controller only, excludes child drivers; earlier invocations retained separately.'},False)
        controller_cost()
        plan,rt=guard(folder)
        write(folder/'controller.json',{'pid':os.getpid(),'token':token(os.getpid())},False)
        script=rt/'elliptic-curves/cas/run_broad_rank_search.py'
        for parent in plan['parents']:
            if (folder/'STOP').exists():controller_cost();return
            job=rt/'broad-pools'/parent['id']
            receipt=execute(folder,job,[plan['sage'],'-python',str(script),'_score','--folder',str(folder),
                                      '--job',str(job),'--parent-id',parent['id']],1800)
            require(receipt and receipt[1] and receipt[1]['status']=='PASS_FIXED_SCORING',
                    'scoring incomplete; retained evidence, no searches dispatched')
        job=rt/'broad-preflight'
        receipt=execute(folder,job,[plan['sage'],'-python',str(script),'_preflight','--folder',str(folder),'--job',str(job)],600)
        require(receipt and receipt[1] and receipt[1]['status']=='PASS_INTERFACE_PREFLIGHT','interface preflight failed')
        rows=build_queue(folder);states={};failures=0
        # Never trust a mutable live state on resume: rebuild from sealed batches.
        for row in rows:
            if row['intake_status']!='ACCEPTED':continue
            previous=None
            for unused in range(plan['max_rounds']+1):
                spec=next_stage(previous,row['rescue'],plan['max_rounds'])
                if spec is None:break
                iteration=previous['round']+1 if previous else 0
                job=rt/'broad-cases'/row['id']/f'batch-{iteration:03d}'
                if not (job/'dispatch-intent.json').exists():break
                if (job/'broad-state.json').exists():
                    require(sha(job/'broad-state.json')==read(job/'broad-state-seal.json')['sha256'],'derived state changed')
                previous=advance(folder,row,spec,previous)
                if previous['status'].startswith('UNKNOWN'):break
            if previous:states[row['id']]=previous
        for wave in range(plan['max_rounds']+1):
            pending=[]
            for row in rows:
                if row['intake_status']!='ACCEPTED':continue
                old=states.get(row['id'])
                if old and old['round']>=wave:continue
                spec=next_stage(old,row['rescue'],plan['max_rounds'])
                if spec is not None:pending.append((row,spec,old))
            if not pending:continue
            with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
                active={};cursor=iter(pending)
                while True:
                    while len(active)<plan['workers'] and not (folder/'STOP').exists():
                        require(shutil.disk_usage(folder).free>=plan['min_free_gib']*1024**3,
                                'disk reserve reached; resume after freeing disk')
                        task=next(cursor,None)
                        if task is None:break
                        future=pool.submit(advance,folder,*task);active[future]=task[0]
                    if not active:break
                    future=next(as_completed(active));row=active.pop(future)
                    try:
                        state=future.result();states[row['id']]=state
                    except BaseException:
                        (folder/'STOP').touch()
                        raise
                    failures=failures+1 if state.get('engineering_failure') else 0
                    if failures>=plan['consecutive_failure_limit']:
                        (folder/'STOP').touch();write(folder/'HALT.json',{'reason':'repeated engineering failures'},False)
                    controller_cost()
                    status(folder)
            if (folder/'STOP').exists():controller_cost();return
        controller_cost()
        write(folder/'COMPLETE.json',{'status':'COMPLETE_BOUNDED_CAMPAIGN','queue_sha256':sha(folder/'queue.json'),
             'boundary':'All scheduled finite stages terminal; UNKNOWN and censored slots remain distinct, not no-point results.'})
        status(folder)


def launch(folder):
    plan,rt=guard(folder)
    with (folder/'launch.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        require(not (folder/'STOP').exists(),'STOP present; inspect and remove it deliberately before resume')
        require(not ((folder/'controller.json').exists() and alive(read(folder/'controller.json'))),'controller already live')
        require(not ((folder/'launch.json').exists() and alive(read(folder/'launch.json'))),'launch already live')
        with (folder/'controller.log').open('ab',buffering=0) as log:
            p=subprocess.Popen([sys.executable,str(rt/'elliptic-curves/cas/run_broad_rank_search.py'),'run','--folder',str(folder)],
                cwd=rt,env=environment(plan['sage']),stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        write(folder/'launch.json',{'pid':p.pid,'token':token(p.pid)},False)
        print(json.dumps({'status':'LAUNCHED','pid':p.pid,'folder':str(folder)}))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','run','launch','resume','status','stop','verify','_score','_preflight','_supervise'))
    p.add_argument('--folder',type=Path,default=DEFAULT);p.add_argument('--sage',default='sage')
    p.add_argument('--workers',type=int,default=4);p.add_argument('--window',type=int,default=262144)
    p.add_argument('--offset',type=int,default=131072);p.add_argument('--prime-bound',type=int,default=997)
    p.add_argument('--ranked',type=int,default=256);p.add_argument('--controls',type=int,default=64)
    p.add_argument('--comparison-ranked',type=int,default=64);p.add_argument('--comparison-controls',type=int,default=16)
    p.add_argument('--rounds',type=int,default=96);p.add_argument('--min-free-gib',type=float,default=20)
    p.add_argument('--parents',nargs='+',choices=PARENTS,default=list(PARENTS))
    p.add_argument('--exclude-equations',type=Path,action='append',default=[])
    p.add_argument('--job',type=Path);p.add_argument('--parent-id',choices=PARENTS)
    a=p.parse_args();a.folder=a.folder.resolve()
    if a.mode=='prepare':prepare(a)
    elif a.mode in ('launch','resume'):launch(a.folder)
    elif a.mode=='status':print(json.dumps(status(a.folder),indent=2))
    elif a.mode=='stop':(a.folder/'STOP').touch();print('STOP recorded; current bounded tasks drain.')
    elif a.mode=='verify':
        plan,rt=guard(a.folder)
        for path in rt.rglob('seal.json'):completed(path.parent)
        print('PASS_FROZEN_INPUT_AND_RETAINED_OUTPUT_INTEGRITY; not a fresh arithmetic replay')
    elif a.mode=='_score':score_task(a.folder,a.job,a.parent_id)
    elif a.mode=='_preflight':preflight_task(a.folder,a.job)
    elif a.mode=='_supervise':supervised_task(a.folder,a.job)
    else:
        plan,rt=guard(a.folder);script=rt/'elliptic-curves/cas/run_broad_rank_search.py'
        if script.resolve()!=Path(__file__).resolve():
            os.execve(sys.executable,[sys.executable,str(script),'run','--folder',str(a.folder)],environment(plan['sage']))
        run(a.folder)


if __name__=='__main__':main()
