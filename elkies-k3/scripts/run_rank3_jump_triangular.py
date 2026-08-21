#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, shutil, time

ap=argparse.ArgumentParser()
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--threads",type=int,default=8)
ap.add_argument("--timeout",type=int,default=300)
ap.add_argument("--outdir",default="artifacts/local/elkies-k3/rank3-triangular")
ap.add_argument("--allow-obstructed-chart",action="store_true",
                help="reproduce the historical j=0 chart despite its target-rank parity obstruction")
a=ap.parse_args()

if not a.allow_obstructed_chart:
    raise SystemExit(
        "refusing obsolete all-IV solve: j=0 forces even geometric MW rank, "
        "but the target requires rank 3; see elkies-k3/E8_A2_KODAIRA_CORRECTION.md "
        "(pass --allow-obstructed-chart only for historical reproduction)"
    )

O=Path(a.outdir); O.mkdir(parents=True,exist_ok=True)
stem=f"p{a.p}-seed{a.seed}"
inp=O/(stem+".ms"); out=O/(stem+".solve"); log=O/(stem+".log")
exp=[
 "sage","elkies-k3/scripts/export_rank3_jump_triangular.sage",
 "--p",str(a.p),"--seed",str(a.seed),"--slice","--out",str(inp)
]
r=subprocess.run(exp,text=True,capture_output=True)
print(r.stdout,flush=True)
if r.returncode:
    print(r.stderr); raise SystemExit(r.returncode)

t=time.time()
try:
    with log.open("w") as h:
        rr=subprocess.run(
            ["msolve","-t",str(a.threads),"-v","2","-f",str(inp),"-o",str(out)],
            stdout=h,stderr=subprocess.STDOUT,text=True,timeout=a.timeout
        )
    print(f"R3TRI|stage=solve|exit={rr.returncode}|seconds={time.time()-t:.2f}|out={out}",flush=True)
    if out.exists():
        print("R3TRI_SOLVE|"+out.read_text(errors="replace")[-4000:].replace("\n"," | "),flush=True)
except subprocess.TimeoutExpired:
    print(f"R3TRI|stage=solve|status=TIMEOUT|seconds={time.time()-t:.1f}",flush=True)
