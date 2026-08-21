#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, shutil, time

ap=argparse.ArgumentParser()
ap.add_argument("--primes",default="101,103,107")
ap.add_argument("--seeds",default="1,2,3")
ap.add_argument("--threads",type=int,default=8)
ap.add_argument("--timeout",type=int,default=600)
ap.add_argument("--outdir",default="artifacts/local/elkies-k3/rank3-jump")
ap.add_argument("--exporter",default="elkies-k3/scripts/export_rank3_jump_msolve.sage")
ap.add_argument("--allow-obstructed-chart",action="store_true",
                help="reproduce the historical j=0 chart despite its target-rank parity obstruction")
a=ap.parse_args()

if not a.allow_obstructed_chart:
    raise SystemExit(
        "refusing obsolete all-IV solve: j=0 forces even geometric MW rank, "
        "but the target requires rank 3; see elkies-k3/E8_A2_KODAIRA_CORRECTION.md "
        "(pass --allow-obstructed-chart only for historical reproduction)"
    )

if not shutil.which("sage"): raise SystemExit("sage missing")
if not shutil.which("msolve"): raise SystemExit("msolve missing")
O=Path(a.outdir); O.mkdir(parents=True,exist_ok=True)

for p in [int(x) for x in a.primes.split(",")]:
  for seed in [int(x) for x in a.seeds.split(",")]:
    stem=f"p{p}-seed{seed}"
    inp=O/(stem+".ms"); out=O/(stem+".solve"); log=O/(stem+".log")
    e=subprocess.run(["sage",a.exporter,"--p",str(p),"--seed",str(seed),"--slices","1","--out",str(inp)],
                     text=True,capture_output=True)
    print(e.stdout.strip(),flush=True)
    if e.returncode:
        print(e.stderr,flush=True); continue
    t0=time.time()
    try:
        with log.open("w") as h:
            r=subprocess.run(["msolve","-t",str(a.threads),"-v","2","-f",str(inp),"-o",str(out)],
                             stdout=h,stderr=subprocess.STDOUT,text=True,timeout=a.timeout)
        dt=time.time()-t0
        print(f"R3JUMP|stage=solve|p={p}|seed={seed}|exit={r.returncode}|seconds={dt:.2f}|out={out}",flush=True)
        if out.exists():
            txt=out.read_text(errors="replace")
            print("R3JUMP_SOLVE_TAIL|"+txt[-2500:].replace("\n"," | "),flush=True)
    except subprocess.TimeoutExpired:
        print(f"R3JUMP|stage=solve|p={p}|seed={seed}|status=TIMEOUT|seconds={time.time()-t0:.1f}",flush=True)
