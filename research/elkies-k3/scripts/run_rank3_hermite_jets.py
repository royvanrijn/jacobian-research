#!/usr/bin/env python3
from pathlib import Path
import subprocess, argparse, time

ap=argparse.ArgumentParser()
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--threads",type=int,default=8)
ap.add_argument("--timeout",type=int,default=180)
a=ap.parse_args()

O=Path("artifacts/local/elkies-k3/hermite-jets"); O.mkdir(parents=True,exist_ok=True)
stem=f"p{a.p}-seed{a.seed}"
inp=O/(stem+".ms"); out=O/(stem+".solve"); log=O/(stem+".log")

r=subprocess.run([
    "sage","elkies-k3/scripts/export_rank3_hermite_jets.sage",
    "--p",str(a.p),"--seed",str(a.seed),"--slice","--out",str(inp)
],text=True,capture_output=True)
print(r.stdout,flush=True)
if r.returncode:
    print(r.stderr); raise SystemExit(r.returncode)

t0=time.time()
try:
    with log.open("w") as h:
        rr=subprocess.run(["msolve","-t",str(a.threads),"-v","2","-f",str(inp),"-o",str(out)],
                          stdout=h,stderr=subprocess.STDOUT,text=True,timeout=a.timeout)
    print(f"HJET|stage=solve|exit={rr.returncode}|seconds={time.time()-t0:.2f}",flush=True)
    if out.exists():
        print("HJET_SOLVE|"+out.read_text(errors="replace")[-4000:].replace("\n"," | "),flush=True)
except subprocess.TimeoutExpired:
    print(f"HJET|stage=solve|status=TIMEOUT|seconds={time.time()-t0:.1f}",flush=True)
