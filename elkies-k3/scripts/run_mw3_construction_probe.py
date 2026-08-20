#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, time, shutil

ap=argparse.ArgumentParser(description="Probe staged MW3 reconstruction systems with msolve.")
ap.add_argument("--stage",choices=["p1","p12","all"],default="p1")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--threads",type=int,default=8)
ap.add_argument("--timeout",type=int,default=300)
ap.add_argument("--outdir",default="artifacts/local/elkies-k3/mw3-construction")
args=ap.parse_args()

if not shutil.which("sage"):
    raise SystemExit("sage missing")
if not shutil.which("msolve"):
    raise SystemExit("msolve missing")

O=Path(args.outdir)
O.mkdir(parents=True,exist_ok=True)
stem=f"{args.stage}-p{args.p}"
inp=O/(stem+".ms")
sol=O/(stem+".solve")
log=O/(stem+".log")

exp=subprocess.run([
    "sage","elkies-k3/scripts/build_mw3_section_system.sage",
    "--stage",args.stage,"--p",str(args.p),"--export",str(inp)
],text=True,capture_output=True)
print(exp.stdout,flush=True)
if exp.returncode:
    print(exp.stderr)
    raise SystemExit(exp.returncode)

t0=time.time()
try:
    with log.open("w") as h:
        r=subprocess.run(
            ["msolve","-t",str(args.threads),"-v","2","-f",str(inp),"-o",str(sol)],
            stdout=h,stderr=subprocess.STDOUT,text=True,timeout=args.timeout
        )
    print(f"MW3SOLVE|stage={args.stage}|p={args.p}|exit={r.returncode}|seconds={time.time()-t0:.1f}|out={sol}",flush=True)
    if sol.exists():
        print("MW3SOLVE_TAIL|"+sol.read_text(errors="replace")[-3000:].replace("\n"," | "),flush=True)
except subprocess.TimeoutExpired:
    print(f"MW3SOLVE|stage={args.stage}|p={args.p}|status=TIMEOUT|seconds={time.time()-t0:.1f}",flush=True)
