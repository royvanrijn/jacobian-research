#!/usr/bin/env python3
from pathlib import Path
import argparse, concurrent.futures, os, re, shutil, subprocess, sys, time

ap = argparse.ArgumentParser(description="Batch finite-field dimension/GB probes for the Elkies rank-17 reduced three-hub system.")
ap.add_argument('--primes', default='101,103,107,109')
ap.add_argument('--slices', default='4,5,6')
ap.add_argument('--seeds', default='1,2,3')
ap.add_argument('--jobs', type=int, default=2, help='concurrent msolve processes')
ap.add_argument('--threads', type=int, default=8, help='threads per msolve process')
ap.add_argument('--timeout', type=int, default=7200)
ap.add_argument('--outdir', default='artifacts/local/elkies-k3/msolve-probe')
ap.add_argument('--exporter', default='elkies-k3/scripts/export_three_hub_msolve.sage')
args = ap.parse_args()

if not shutil.which('sage'):
    raise SystemExit('sage not found')
if not shutil.which('msolve'):
    raise SystemExit('msolve not found; on Ubuntu/Debian: sudo apt install msolve')

primes=[int(x) for x in args.primes.split(',') if x]
slices=[int(x) for x in args.slices.split(',') if x]
seeds=[int(x) for x in args.seeds.split(',') if x]
out=Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

def one(job):
    p,s,seed=job
    stem=f'p{p}-s{s}-seed{seed}'
    inp=out/(stem+'.ms')
    dim=out/(stem+'.dim')
    log=out/(stem+'.log')
    exp=['sage',args.exporter,'--p',str(p),'--slices',str(s),'--seed',str(seed),'--out',str(inp)]
    t0=time.time()
    try:
        e=subprocess.run(exp, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=300)
        if e.returncode:
            return stem,'EXPORT_FAIL',time.time()-t0,e.stdout[-2000:]
        cmd=['msolve','-t',str(args.threads),'-v','2','-g','1','-f',str(inp),'-o',str(dim)]
        with log.open('w') as fh:
            r=subprocess.run(cmd, text=True, stdout=fh, stderr=subprocess.STDOUT, timeout=args.timeout)
        txt=(dim.read_text(errors='replace') if dim.exists() else '') + '\n' + log.read_text(errors='replace')[-6000:]
        # Stable output markers from msolve/manpage plus verbose hints.
        if re.search(r'\[\s*-1\s*\]', txt): status='EMPTY'
        elif re.search(r'\[\s*1\s*,\s*\d+\s*,\s*-1\s*,', txt): status='POSITIVE_DIM'
        elif r.returncode==0: status='GB_DONE'
        else: status=f'EXIT_{r.returncode}'
        return stem,status,time.time()-t0,txt[-1200:]
    except subprocess.TimeoutExpired:
        return stem,'TIMEOUT',time.time()-t0,''

jobs=[(p,s,seed) for p in primes for s in slices for seed in seeds]
print(f'PROBE|jobs={len(jobs)}|parallel={args.jobs}|threads_each={args.threads}|out={out}', flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
    futs={ex.submit(one,j):j for j in jobs}
    for fut in concurrent.futures.as_completed(futs):
        stem,status,dt,tail=fut.result()
        print(f'PROBE|{stem}|status={status}|seconds={dt:.1f}', flush=True)
        if status in ('GB_DONE','POSITIVE_DIM','EMPTY'):
            print('  '+tail.replace('\n',' | ')[:1000], flush=True)

print('PROBE|done')
print('Interpretation: if slices=4 stays positive-dimensional, slices=5 becomes finite/GB_DONE, and slices=6 is generically empty across primes, that is strong evidence for the expected 5-dimensional coordinate component (≈ 1 true moduli + 4 gauge).')
