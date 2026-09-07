from pathlib import Path
import re
import numpy as np

BASE = Path(__file__).resolve().parents[1]
R = BASE / "results"

files = sorted(R.glob("E29-short-pool-v2-vectors-step*.npy"))

rows = []

for vp in files:
    m = re.search(r"step(\d+)\.npy$", vp.name)
    if not m:
        continue

    step = int(m.group(1))

    npth = R / vp.name.replace("vectors-", "norms-")

    try:
        V = np.load(vp, mmap_mode="r")
    except Exception as e:
        print("BAD", vp.name, e)
        continue

    if V.ndim != 2:
        continue

    dim = V.shape[1]

    if dim == 29:
        job = "E29"
    elif dim == 21:
        job = "rank21"
    elif dim == 19:
        job = "rank19"
    else:
        job = f"unknown-{dim}"

    if npth.exists():
        N = np.load(npth, mmap_mode="r")
        norms = len(N)
        nmin = float(np.min(N)) if len(N) else float("nan")
        nmax = float(np.max(N)) if len(N) else float("nan")
    else:
        norms = None
        nmin = nmax = float("nan")

    ok = norms == V.shape[0]

    rows.append((job, step, V.shape[0], dim, norms, nmin, nmax, ok, vp))

for row in rows:
    job, step, count, dim, norms, nmin, nmax, ok, vp = row
    print(
        f"{job:7s}"
        f" step={step:10d}"
        f" vectors={count:8d}"
        f" dim={dim:2d}"
        f" norms={str(norms):>8s}"
        f" min={nmin:10.4f}"
        f" max={nmax:10.4f}"
        f" ok={ok}"
        f" file={vp.name}"
    )

print()
print("LATEST VALID PER JOB")

for job in ("E29", "rank19", "rank21"):
    x = [r for r in rows if r[0] == job and r[7]]

    if not x:
        print(job, "NONE")
        continue

    best = max(x, key=lambda r:r[1])

    print(
        f"{job}: step={best[1]}"
        f" vectors={best[2]}"
        f" norm=[{best[5]:.6g},{best[6]:.6g}]"
        f" file={best[8].name}"
    )
