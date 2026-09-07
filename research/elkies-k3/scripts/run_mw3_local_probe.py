#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, shutil, time

ap=argparse.ArgumentParser()
ap.add_argument("--stage",choices=["p1","p12","all"],default="p1")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--threads",type=int,default=8)
ap.add_argument("--timeout",type=int,default=300)
args=ap.parse_args()

if not shutil.which("sage"): raise SystemExit("sage missing")
if not shutil.which("msolve"): raise SystemExit("msolve missing")

O=Path("artifacts/local/elkies-k3/mw3-local")
O.mkdir(parents=True,exist_ok=True)
inp=O/f"{args.stage}-p{args.p}.ms"
out=O/f"{args.stage}-p{args.p}.solve"
log=O/f"{args.stage}-p{args.p}.log"

e=subprocess.run([
    "sage","elkies-k3/scripts/build_mw3_local_tate_system.sage",
    "--stage",args.stage,"--p",str(args.p),"--export",str(inp)
],text=True,capture_output=True)
print(e.stdout,flush=True)
if e.returncode:
    print(e.stderr); raise SystemExit(e.returncode)

t0=time.time()
try:
    with log.open("w") as h:
        r=subprocess.run(["msolve","-t",str(args.threads),"-v","2","-f",str(inp),"-o",str(out)],
                         stdout=h,stderr=subprocess.STDOUT,text=True,timeout=args.timeout)
    print(f"MW3LOCALSOLVE|stage={args.stage}|p={args.p}|exit={r.returncode}|seconds={time.time()-t0:.1f}",flush=True)
    if out.exists():
        print("MW3LOCALSOLVE_TAIL|"+out.read_text(errors="replace")[-3000:].replace("\n"," | "),flush=True)
except subprocess.TimeoutExpired:
    print(f"MW3LOCALSOLVE|stage={args.stage}|p={args.p}|status=TIMEOUT|seconds={time.time()-t0:.1f}",flush=True)
