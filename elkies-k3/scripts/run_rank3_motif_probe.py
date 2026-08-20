#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, shutil, time
p=argparse.ArgumentParser()
p.add_argument('--p',type=int,default=101)
p.add_argument('--slices',default='3,4,5,6')
p.add_argument('--seed',type=int,default=1)
p.add_argument('--threads',type=int,default=8)
p.add_argument('--timeout',type=int,default=1800)
p.add_argument('--outdir',default='artifacts/local/elkies-k3/motif-diamond')
p.add_argument('--exporter',default='elkies-k3/scripts/export_rank3_motif_msolve.sage')
a=p.parse_args()
if not shutil.which('sage'): raise SystemExit('sage missing')
if not shutil.which('msolve'): raise SystemExit('msolve missing')
O=Path(a.outdir); O.mkdir(parents=True,exist_ok=True)
for s in [int(x) for x in a.slices.split(',')]:
    stem=f'p{a.p}-s{s}-seed{a.seed}'; inp=O/(stem+'.ms'); out=O/(stem+'.out'); log=O/(stem+'.log')
    e=subprocess.run(['sage',a.exporter,'--p',str(a.p),'--slices',str(s),'--seed',str(a.seed),'--out',str(inp)],text=True,capture_output=True)
    print(e.stdout.strip(),flush=True)
    if e.returncode: print(e.stderr); continue
    t=time.time()
    try:
        with log.open('w') as h:
            r=subprocess.run(['msolve','-t',str(a.threads),'-v','2','-g','1','-f',str(inp),'-o',str(out)],stdout=h,stderr=subprocess.STDOUT,text=True,timeout=a.timeout)
        print(f'MOTIFPROBE|slices={s}|exit={r.returncode}|seconds={time.time()-t:.1f}|out={out}',flush=True)
        if out.exists(): print(out.read_text(errors='replace')[-1500:],flush=True)
    except subprocess.TimeoutExpired:
        print(f'MOTIFPROBE|slices={s}|status=TIMEOUT|seconds={time.time()-t:.1f}',flush=True)
