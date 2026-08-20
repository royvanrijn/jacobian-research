#!/usr/bin/env python3
from pathlib import Path
import argparse, subprocess, shutil, time, concurrent.futures

ap=argparse.ArgumentParser(description="Run sliced E6 P1 samples over multiple seeds.")
ap.add_argument("--p",type=int,default=101)
ap.add_argument("--seeds",default="1,2,3,4")
ap.add_argument("--slices",type=int,default=3)
ap.add_argument("--threads",type=int,default=4)
ap.add_argument("--workers",type=int,default=2)
ap.add_argument("--timeout",type=int,default=300)
ap.add_argument("--outdir",default="artifacts/local/elkies-k3/e6-sliced")
args=ap.parse_args()

if not shutil.which("sage"):
    raise SystemExit("sage missing")
if not shutil.which("msolve"):
    raise SystemExit("msolve missing")

O=Path(args.outdir)
O.mkdir(parents=True,exist_ok=True)
seeds=[int(x) for x in args.seeds.split(",") if x.strip()]

def one(seed):
    stem=f"p{args.p}-seed{seed}"
    inp=O/(stem+".ms")
    sol=O/(stem+".solve")
    log=O/(stem+".log")

    exp=subprocess.run([
        "sage","elkies-k3/scripts/export_e6_p1_sliced.sage",
        "--p",str(args.p),
        "--seed",str(seed),
        "--slices",str(args.slices),
        "--out",str(inp)
    ],text=True,capture_output=True)

    if exp.returncode:
        return seed,"export-error",0.0,exp.stdout+"\n"+exp.stderr

    t0=time.time()
    try:
        with log.open("w") as h:
            r=subprocess.run(
                ["msolve","-t",str(args.threads),"-v","2","-f",str(inp),"-o",str(sol)],
                stdout=h,stderr=subprocess.STDOUT,text=True,timeout=args.timeout
            )
        dt=time.time()-t0
        tail=sol.read_text(errors="replace")[-3000:] if sol.exists() else ""
        return seed,f"exit-{r.returncode}",dt,exp.stdout+"\n"+tail
    except subprocess.TimeoutExpired:
        return seed,"timeout",time.time()-t0,exp.stdout

print(f"E6BATCH|stage=start|p={args.p}|seeds={seeds}|workers={args.workers}|threads_each={args.threads}",flush=True)

with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
    futs={ex.submit(one,s):s for s in seeds}
    for fut in concurrent.futures.as_completed(futs):
        seed,status,dt,text=fut.result()
        print(f"E6BATCH|seed={seed}|status={status}|seconds={dt:.1f}",flush=True)
        for line in text.splitlines():
            if line.startswith("E6SLICE|"):
                print(line,flush=True)
        if "exit-0"==status:
            print("E6BATCH_SOLVE_TAIL|seed="+str(seed)+"|"+text[-2000:].replace("\\n"," | "),flush=True)

print("E6BATCH|stage=done",flush=True)
