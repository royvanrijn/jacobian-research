#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, concurrent.futures, time, shutil

ap=argparse.ArgumentParser(description="Search correct 8x8 E6 coordinate slices with bounded msolve runs.")
ap.add_argument("--input",default="artifacts/local/elkies-k3/e6-param0.ms")
ap.add_argument("--p",type=int,default=None,help="informational; characteristic is read from input")
ap.add_argument("--seeds",type=int,default=32)
ap.add_argument("--workers",type=int,default=8)
ap.add_argument("--threads",type=int,default=2)
ap.add_argument("--timeout",type=int,default=45)
ap.add_argument("--outdir",default="artifacts/local/elkies-k3/e6-coordinate-search")
args=ap.parse_args()

if not shutil.which("sage"): raise SystemExit("sage missing")
if not shutil.which("msolve"): raise SystemExit("msolve missing")
src=Path(args.input)
if not src.exists(): raise SystemExit(f"missing input: {src}")

# Always fix r0 so the exact P1_1 -> y1 denominator becomes a field scalar.
triples=[
    ("r0","x1","a4"),
    ("r0","x1","a2"),
    ("r0","x1","a1"),
    ("r0","a4","a2"),
    ("r0","a4","a1"),
    ("r0","a2","a1"),
]
O=Path(args.outdir); O.mkdir(parents=True,exist_ok=True)

print(f"E6SEARCH|stage=start|input={src}|seeds={args.seeds}|workers={args.workers}|threads_each={args.threads}|timeout={args.timeout}",flush=True)

def one(seed):
    kill=triples[(seed-1)%len(triples)]
    tag="-".join(kill)
    inp=O/f"seed{seed}-{tag}.ms"
    sol=O/f"seed{seed}-{tag}.solve"
    log=O/f"seed{seed}-{tag}.log"

    b=subprocess.run([
        "sage","elkies-k3/scripts/build_e6_coordinate_slice.sage",
        "--input",str(src),"--seed",str(seed),"--kill",",".join(kill),"--out",str(inp)
    ],text=True,capture_output=True)
    if b.returncode:
        return seed,kill,"build-error",0.0,b.stdout+"\n"+b.stderr

    t0=time.time()
    try:
        with log.open("w") as h:
            r=subprocess.run([
                "msolve","-t",str(args.threads),"-v","1","-f",str(inp),"-o",str(sol)
            ],stdout=h,stderr=subprocess.STDOUT,text=True,timeout=args.timeout)
        dt=time.time()-t0
    except subprocess.TimeoutExpired:
        return seed,kill,"timeout",time.time()-t0,b.stdout

    text=""
    if sol.exists(): text+=sol.read_text(errors="replace")
    if log.exists(): text+="\n"+log.read_text(errors="replace")
    low=text.lower()
    if "no solution" in low: status="no-solution"
    elif r.returncode==0: status="SOLVED"
    else: status=f"exit-{r.returncode}"
    return seed,kill,status,dt,b.stdout+"\n"+text[-4000:]

jobs=range(1,args.seeds+1)
hits=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
    fs=[ex.submit(one,s) for s in jobs]
    for f in concurrent.futures.as_completed(fs):
        seed,kill,status,dt,text=f.result()
        print(f"E6SEARCH|seed={seed}|kill={','.join(kill)}|status={status}|seconds={dt:.1f}",flush=True)
        if status=="build-error":
            print("E6SEARCH_BUILD_ERROR|seed="+str(seed)+"|"+text[-2500:].replace("\n"," | "),flush=True)
        elif status=="SOLVED":
            hits.append((seed,kill))
            print("E6SEARCH_HIT|seed="+str(seed)+"|kill="+",".join(kill)+"|"+text[-2500:].replace("\n"," | "),flush=True)

print("E6SEARCH|stage=done|hits="+repr([(s,list(k)) for s,k in hits]),flush=True)
