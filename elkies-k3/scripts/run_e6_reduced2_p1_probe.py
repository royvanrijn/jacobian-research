#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, time

ap=argparse.ArgumentParser()
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--threads",type=int,default=8)
ap.add_argument("--timeout",type=int,default=300)
args=ap.parse_args()

O=Path("artifacts/local/elkies-k3/e6-reduced2")
O.mkdir(parents=True,exist_ok=True)
inp=O/f"p1-p{args.p}.ms"
out=O/f"p1-p{args.p}.solve"
log=O/f"p1-p{args.p}.log"

r=subprocess.run([
    "sage","elkies-k3/scripts/build_e6_reduced2_p1_system.sage",
    "--p",str(args.p),"--export",str(inp)
],text=True,capture_output=True)

print(r.stdout,flush=True)
if r.returncode:
    print(r.stderr)
    raise SystemExit(r.returncode)

t0=time.time()
try:
    with log.open("w") as h:
        rr=subprocess.run(
            ["msolve","-t",str(args.threads),"-v","2","-f",str(inp),"-o",str(out)],
            stdout=h,stderr=subprocess.STDOUT,text=True,timeout=args.timeout
        )
    print(f"E6RED2SOLVE|p={args.p}|exit={rr.returncode}|seconds={time.time()-t0:.1f}",flush=True)
    if out.exists():
        print("E6RED2SOLVE_TAIL|"+out.read_text(errors="replace")[-3000:].replace("\n"," | "),flush=True)
except subprocess.TimeoutExpired:
    print(f"E6RED2SOLVE|p={args.p}|status=TIMEOUT|seconds={time.time()-t0:.1f}",flush=True)
