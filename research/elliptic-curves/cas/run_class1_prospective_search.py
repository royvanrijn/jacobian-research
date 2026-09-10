#!/usr/bin/env python3
"""Frozen class1 ordinary search: cheap triage, controls, certified V3.

No strict-class prerequisite; no exceptional parameter or point inputs.
Run is detached externally and operates entirely through deterministic code.
"""
import argparse
import ctypes
import fcntl
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED

CAS=Path(__file__).resolve().parent;ROOT=CAS.parents[1]
DEFAULT=ROOT/'artifacts/local/elliptic-curves/class1-prospective-v2'
PIN='7c6ee40c46f5a1f3d1fc464b5685a0e4c77ed1f3347b67865d7c9f93b2acd862'

def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,value,immutable=False):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    data=json.dumps(value,indent=2,sort_keys=True)+'\n'
    if immutable and p.exists():
        if p.read_text()!=data:raise ValueError('immutable receipt changed: '+str(p))
        return
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(data);tmp.replace(p)

def allowance(rank):
    if rank>=27:return 8192
    if rank>=25:return 2048
    if rank>=23:return 1024
    if rank>=20:return 256
    return 24

def selection(rows):
    controls=sorted((r for r in rows if r['control']),key=lambda r:(r.get('control_order',''),r['index']))
    ranked=sorted((r for r in rows if not r['control']),key=lambda r:(-r['score_units'],r['model_bits'],r['index']))
    out=[]
    for i in range(max(len(controls),(len(ranked)+6)//7)):
        out.extend(ranked[7*i:7*i+7]);out.extend(controls[i:i+1])
    assert len({r['index'] for r in out})==len(rows)
    return out

def control_window(offset,count):
    ordered=sorted(range(offset,offset+count),key=lambda i:(hashlib.sha256(('class1-controls-v2/'+str(i)).encode()).hexdigest(),i))
    return {i:position for position,i in enumerate(ordered[:count//8])}

def integral_projection(parent):
    """Exact constant coordinate scaling for the integer-only scorer.

    Keep the certified parent unchanged; attach the reversible input marking.
    No discriminant factorization or parameter/score information is used.
    """
    A=list(map(Fraction,parent['A_coefficients_low_to_high']))
    B=list(map(Fraction,parent['B_coefficients_low_to_high']))
    scale=math.lcm(*(c.denominator for c in A+B))
    for p in [2,3]+primes(997):
        while scale%p==0:
            candidate=scale//p
            if any((c*candidate**4).denominator!=1 for c in A) or any((c*candidate**6).denominator!=1 for c in B):break
            scale=candidate
    if scale==1:return parent
    result=json.loads(json.dumps(parent))
    result['A_coefficients_low_to_high']=[str(c*scale**4) for c in A]
    result['B_coefficients_low_to_high']=[str(c*scale**6) for c in B]
    for section in result['sections']:
        for coordinate,power in [('X',2),('Y',3)]:
            section[coordinate]['numerator_coefficients_low_to_high']=[str(Fraction(c)*scale**power) for c in section[coordinate]['numerator_coefficients_low_to_high']]
    result['input_coordinate_scaling']={'scale':scale,'forward':'X=scale^2*x; Y=scale^3*y',
        'inverse':'x=X/scale^2; y=Y/scale^3','height_gram_unchanged':True}
    assert all(Fraction(c).denominator==1 for c in result['A_coefficients_low_to_high']+result['B_coefficients_low_to_high'])
    return result

def prepare(args):
    prep_cpu=time.process_time();prep_wall=time.monotonic()
    from run_euclidean_seed_foundry import normalized_parent,SUFFIXES
    source=getattr(args,'parent',None) or ROOT/'artifacts/generated-results/elliptic-curves/x1092_class1_arithmetic_gate_v1/parent.json'
    expected=getattr(args,'parent_sha256',None) or PIN
    assert sha(source)==expected
    parent=read(source); normalized=integral_projection(normalized_parent(parent,parent['generic_height_gram']))
    normalized['family']=parent.get('family','x1092-class1-prospective-v1')
    folder=args.folder.resolve();folder.mkdir(parents=True,exist_ok=False)
    rt=folder/'runtime/research'
    for directory in (CAS,ROOT/'elliptic-curves/ecsearch',ROOT/'elkies-k3/scripts'):
        for p in directory.rglob('*'):
            if p.is_file() and p.suffix in SUFFIXES and '__pycache__' not in p.parts:
                dest=rt/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    for p in (ROOT/'elliptic-curves').glob('*.py'):
        dest=rt/p.relative_to(ROOT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
    write(rt/'search-inputs/parent.json',normalized,True)
    # Generic projection only. No old ledgers, saved scores, or oracle packets.
    plan={'schema':'class1.prospective-search.v1','root':str(rt),'parent_sha256':expected,
      'normalized_parent_sha256':sha(rt/'search-inputs/parent.json'),'sage':shutil.which('sage'),
      'workers':args.workers,'initial_window':args.window,'initial_offset':args.offset,
      'subsequent_window_rule':'contiguous address windows, double size up to 262144; freeze each before scoring',
      'score_prime_bound':args.prime_bound,'control_fraction':'1/8',
      'control_rule':'First count/8 indices sorted by SHA256(class1-controls-v2/index), with that same score-independent order; both signs eligible.',
      'initial_calls':24,'cumulative_rank_allowances':{'17':24,'20':256,'23':1024,'25':2048,'27':8192},
      'height':125000,'job_wall_seconds':7200,'job_rss_bytes':3*1024**3,
      'failure_semantics':'Censored/failed/generic-gate misses remain UNKNOWN and retain all logs/clouds; no automatic repeat of the same failed batch.',
      'full_cloud_reconciliation':'Every batch, including below23, with two finite certificate implementations.',
      'strict_class_required':False,'other_parent_classes_enabled':False,'historical_exceptional_inputs':False,
      'score_definition':'sum round(1e12*(2-a_p)*log(p)/(p+1-a_p)) over smooth reductions of coefficientwise p-minimized charts; skipped local singular reductions are recorded, not rank exclusions.',
      'cost_rule':'per-parameter scoring CPU plus equal share of cold residue-table CPU, and isolated process-tree search CPU including bank, points, replay and reconciliation; preparation cost reported separately.',
      'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip()}
    assert plan['sage'] and 1<=args.workers<=4 and args.window>=8 and args.window%8==0
    write(folder/'plan.json',plan,True)
    write(folder/'manifest.json',{'plan_sha256':sha(folder/'plan.json'),
      'files':{str(p.relative_to(rt)):sha(p) for p in sorted(rt.rglob('*')) if p.is_file()},
      'executables':{plan['sage']:sha(plan['sage']),'/usr/bin/gp':sha('/usr/bin/gp')}},True)
    write(folder/'preparation-cost.json',{'cpu_seconds':time.process_time()-prep_cpu,'wall_seconds':time.monotonic()-prep_wall},True)
    print('PREPARED',folder,flush=True)

def guard(folder):
    plan=read(folder/'plan.json');manifest=read(folder/'manifest.json');rt=Path(plan['root'])
    assert sha(folder/'plan.json')==manifest['plan_sha256']
    for path,digest in manifest['files'].items():assert sha(rt/path)==digest,path
    for path,digest in manifest['executables'].items():assert sha(path)==digest,path
    return plan,rt

def primes(bound):
    return [p for p in range(5,bound+1) if all(p%d for d in range(2,math.isqrt(p)+1))]

def valuation(n,p):
    if not n:return 1000000
    v=0
    while n%p==0:n//=p;v+=1
    return v

def score(folder,window,index,offset,count):
    import numpy as np
    from euclidean_seed_sieve import address
    plan=read(folder/'plan.json');rt=Path(plan['root']);parent=read(rt/'search-inputs/parent.json')
    A=[int(Fraction(c)) for c in parent['A_coefficients_low_to_high']]
    B=[int(Fraction(c)) for c in parent['B_coefficients_low_to_high']]
    assert all(Fraction(c).denominator==1 for c in parent['A_coefficients_low_to_high']+parent['B_coefficients_low_to_high'])
    started=time.process_time();tables=[];checks=[]
    cached=read(window/'score-tables.json') if (window/'score-tables.json').exists() else None
    def ev(coeff,r,p):
        v=0
        for c in reversed(coeff):v=(v*r+c)%p
        return v
    for p in ([] if cached else primes(plan['score_prime_bound'])):
        k=min(min(valuation(c,p) for c in A)//4,min(valuation(c,p) for c in B)//6)
        aa=[(c//p**(4*k))%p for c in A];bb=[(c//p**(6*k))%p for c in B]
        av=np.array([ev(aa,r,p) for r in range(p)]+[(aa+[0]*9)[8]],dtype=np.int64)
        bv=np.array([ev(bb,r,p) for r in range(p)]+[(bb+[0]*13)[12]],dtype=np.int64)
        chi=np.zeros(p,dtype=np.int64);chi[1:]=-1
        for x in range(1,p):chi[x*x%p]=1
        x=np.arange(p,dtype=np.int64)
        traces=-chi[(x[None,:]**3+av[:,None]*x[None,:]+bv[:,None])%p].sum(axis=1)
        smooth=(4*av**3+27*bv**2)%p!=0
        units=[round(1e12*(2-int(ap))*math.log(p)/(p+1-int(ap))) if good else 0 for ap,good in zip(traces,smooth)]
        # A separate scalar Legendre implementation checks two entries per p.
        for r in (0,p):
            direct=-sum(0 if (v:=(z**3+int(av[r])*z+int(bv[r]))%p)==0 else (1 if pow(v,(p-1)//2,p)==1 else -1) for z in range(p))
            assert direct==int(traces[r]);checks.append([p,r,direct])
        assert all(int(t)**2<=4*p for t,g in zip(traces,smooth) if g)
        tables.append({'prime':p,'global_scalings':k,'score_units':units,'smooth':list(map(bool,smooth)),'traces':list(map(int,traces))})
    cold=time.process_time()-started
    if cached:tables=cached['tables'];cold=cached['cold_cpu_seconds']
    else:write(window/'score-tables.json',{'tables':tables,'cold_cpu_seconds':cold,'scalar_checks':checks},True)
    def homog(coeff,a,b,degree):return sum(c*a**i*b**(degree-i) for i,c in enumerate(coeff))
    rows=[];controls=control_window(offset,count)
    for j in range(offset,offset+count):
        tick=time.process_time();t=Fraction(address(j));a,b=t.numerator,t.denominator
        av,bv=homog(A,a,b,8),homog(B,a,b,12);delta=-16*(4*av**3+27*bv**2)
        score_units=0;bad=[]
        for tab in tables:
            p=tab['prime'];r=a*pow(b,-1,p)%p if b%p else p
            score_units+=tab['score_units'][r]
            if not tab['smooth'][r]:bad.append(p)
        row={'index':j,'parameter':str(t),'control':j in controls,'control_order':controls.get(j), 'score_units':score_units,
          'model_bits':max(abs(av).bit_length(),abs(bv).bit_length()),
          'discriminant_bits':abs(delta).bit_length(),'nonsingular':bool(delta),
          'local_singular_primes':bad,'smooth_prime_count':len(tables)-len(bad),
          'conductor':'NOT_COMPUTED','discriminant_is_not_conductor':True,
          'score_direct_cpu_seconds':time.process_time()-tick,'cold_score_cpu_share_seconds':cold/count}
        rows.append(row)
        if len(rows)%4096==0:write(window/'score-progress.json',{'scored':len(rows),'total':count})
    queue=[r['index'] for r in selection(rows)]
    overhead=max(0,time.process_time()-started-(0 if cached else cold)-sum(r['score_direct_cpu_seconds'] for r in rows))
    for r in rows:r['score_overhead_cpu_share_seconds']=overhead/count
    write(window/'scores.json',{'rows':rows,'cold_cpu_seconds':cold,'overhead_cpu_seconds':overhead,'window_index':index,'offset':offset,'count':count},True)
    write(window/'queue.json',{'indices':queue},True)
    print('SCORED',count,'addresses; cold tables CPU',cold,flush=True)

def run_job(folder,job):
    from research_runtime.supervisor import run,Limits
    # Adopt and reap any timed-out grandchildren, preserving their CPU totals.
    ctypes.CDLL(None).prctl(36,1,0,0,0)
    plan=read(folder/'plan.json');rt=Path(plan['root'])
    before=resource.getrusage(resource.RUSAGE_CHILDREN);own=time.process_time();wall=time.monotonic()
    r=run([plan['sage'],'-python',str(rt/'elliptic-curves/cas/parent_foundry_worker.py'),'--job',str(job)],
      limits=Limits(plan['job_wall_seconds'],plan['job_rss_bytes']),log_path=job/'worker.log',
      checkpoint_path=job/'supervisor.json',cwd=rt)
    while True:
        try:os.waitpid(-1,0)
        except ChildProcessError:break
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    write(job/'cost.json',{'user_cpu_seconds':after.ru_utime-before.ru_utime,
      'system_cpu_seconds':after.ru_stime-before.ru_stime,'driver_cpu_seconds':time.process_time()-own,
      'total_cpu_seconds':after.ru_utime-before.ru_utime+after.ru_stime-before.ru_stime+time.process_time()-own,
      'wall_seconds':time.monotonic()-wall,'outcome':r['outcome'],'returncode':r['returncode'],
      'scope':'isolated subreaper: worker and descendant CPU including all bank, mapping, point and certificate work'},True)
    return 0 if r['outcome']=='completed' and r['returncode']==0 else 1

def fibre(folder,window,row):
    plan=read(folder/'plan.json');rt=Path(plan['root'])
    case=window/'cases'/str(row['index']);case.mkdir(parents=True,exist_ok=True)
    write(case/'selection.json',row,True)
    if (case/'terminal.json').exists():return read(case/'terminal.json')
    total=0;rank=17;packet=None;batches=[];cpu=0;status='COMPLETE_BOUNDED'
    for epoch in range(256):
        if (folder/'STOP').exists():status='PAUSED';break
        cap=allowance(rank)
        if total>=cap or rank>=32:break
        size=min(cap-total,24 if rank<20 else 64 if rank<23 else 128 if rank<27 else 256)
        job=case/f'batch-{epoch:03d}';job.mkdir(exist_ok=True)
        request={'parent':'search-inputs/parent.json','parent_sha256':plan['normalized_parent_sha256'],
          'parameter':row['parameter'],'allowance':size,'height':plan['height'],
          'bank_index':int(hashlib.sha256(('class1-bank/'+row['parameter']).encode()).hexdigest()[:8],16)+epoch}
        if packet:request.update(packet=str(packet.relative_to(rt)),packet_sha256=sha(packet))
        write(job/'request.json',request,True)
        if not (job/'cost.json').exists():
            if (job/'dispatched.json').exists():status='UNKNOWN_INTERRUPTED_DISPATCH';break
            write(job/'dispatched.json',{'started_at':time.time()},True)
            subprocess.run([sys.executable,str(Path(__file__)),'job','--folder',str(folder),'--job',str(job)],check=False)
        if not (job/'cost.json').exists():status='UNKNOWN_DRIVER_FAILURE';break
        cost=read(job/'cost.json');cpu+=cost['total_cpu_seconds']
        if cost['returncode']!=0 or cost['outcome']!='completed' or not (job/'result.json').exists():
            status='UNKNOWN_CENSORED_OR_FAILED';break
        result=read(job/'result.json')
        if result['status']!='PASS_CERTIFIED_PARENT_EVALUATION':status=result['status'];break
        packet=job/'packet.json';assert sha(packet)==result['packet_sha256']
        assert read(job/'verified.json')['status']=='PASS_TWO_FINITE_IMPLEMENTATIONS'
        rank=result['rank_lower_bound'];total+=result['calls'];batches.append(result)
        write(case/'progress.json',{'rank_lower_bound':rank,'calls':total,'search_cpu_seconds':cpu,'target_calls':allowance(rank)})
        if result['calls']==0:status='UNKNOWN_NO_PROGRESS';break
    terminal={'status':status,'parameter':row['parameter'],'index':row['index'],'control':row['control'],
      'score_units':row['score_units'],'rank_lower_bound':rank if packet else None,'calls':total,
      'search_cpu_seconds':cpu,'total_cpu_seconds':cpu+row['score_direct_cpu_seconds']+row['cold_score_cpu_share_seconds']+row['score_overhead_cpu_share_seconds'],
      'score_direct_cpu_seconds':row['score_direct_cpu_seconds'],'cold_score_cpu_share_seconds':row['cold_score_cpu_share_seconds'],
      'score_overhead_cpu_share_seconds':row['score_overhead_cpu_share_seconds'],
      'packet':str(packet) if packet else None,'batches':batches,
      'boundary':'New prospective exposure, certified subgroup lower bounds only; no external novelty or conductor-record claim.'}
    if status!='PAUSED':write(case/'terminal.json',terminal,True)
    return terminal

def run_campaign(folder):
    plan,rt=guard(folder)
    frozen=rt/'elliptic-curves/cas'/Path(__file__).name
    if Path(__file__).resolve()!=frozen.resolve():
        os.execv(sys.executable,[sys.executable,str(frozen),'run','--folder',str(folder)])
    lock=(folder/'controller.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    write(folder/'controller.json',{'pid':os.getpid(),'started_at':time.time(),'strict_class_required':False})
    index=0;offset=plan['initial_offset'];count=plan['initial_window']
    while not (folder/'STOP').exists():
        window=rt/f'ordinary-search/window-{index:03d}';window.mkdir(parents=True,exist_ok=True)
        write(window/'window.json',{'offset':offset,'count':count,'index':index},True)
        if not (window/'scores.json').exists():score(folder,window,index,offset,count)
        data=read(window/'scores.json');rows={r['index']:r for r in data['rows']}
        if not (window/'queue.json').exists():write(window/'queue.json',{'indices':[r['index'] for r in selection(data['rows'])]},True)
        todo=[i for i in read(window/'queue.json')['indices'] if not (window/'cases'/str(i)/'terminal.json').exists()]
        with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
            active={};iterator=iter(todo)
            while active or todo:
                while len(active)<plan['workers'] and not (folder/'STOP').exists():
                    i=next(iterator,None)
                    if i is None:todo=[];break
                    active[pool.submit(fibre,folder,window,rows[i])]=i
                if not active:break
                done,_=wait(active,timeout=10,return_when=FIRST_COMPLETED)
                for future in done:
                    i=active.pop(future)
                    try:result=future.result()
                    except Exception as e:
                        write(folder/'ERROR.json',{'index':i,'error':repr(e)});raise
                    print('FIBRE',i,result['status'],'rank',result['rank_lower_bound'],'CPU',result['total_cpu_seconds'],flush=True)
                write(folder/'STATUS.json',{'status':'DRAINING' if (folder/'STOP').exists() else 'RUNNING',
                  'window':index,'addresses':count,'active_indices':list(active.values()),
                  'completed':sum((window/'cases'/str(i)/'terminal.json').exists() for i in rows),
                  'controller_cpu_seconds':time.process_time(),
                  'strict_class_required':False,'updated_at':time.time()})
        if (folder/'STOP').exists():break
        index+=1;offset+=count;count=min(count*2,262144)
    write(folder/'STATUS.json',{'status':'STOPPED_AFTER_DRAIN','updated_at':time.time()})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['prepare','score','run','job'])
    p.add_argument('--folder',type=Path,default=DEFAULT);p.add_argument('--workers',type=int,default=2)
    p.add_argument('--window',type=int,default=65536);p.add_argument('--prime-bound',type=int,default=997)
    p.add_argument('--offset',type=int,default=65536)
    p.add_argument('--parent',type=Path)
    p.add_argument('--parent-sha256')
    p.add_argument('--job',type=Path);args=p.parse_args();args.folder=args.folder.resolve()
    if args.mode=='prepare':prepare(args)
    elif args.mode=='score':
        plan,rt=guard(args.folder)
        frozen=rt/'elliptic-curves/cas'/Path(__file__).name
        if Path(__file__).resolve()!=frozen.resolve():
            os.execv(sys.executable,[sys.executable,str(frozen),'score','--folder',str(args.folder)])
        window=rt/'ordinary-search/window-000';window.mkdir(parents=True,exist_ok=True)
        write(window/'window.json',{'offset':plan['initial_offset'],'count':plan['initial_window'],'index':0},True)
        score(args.folder,window,0,plan['initial_offset'],plan['initial_window'])
    elif args.mode=='run':run_campaign(args.folder)
    else:sys.exit(run_job(args.folder,args.job))
