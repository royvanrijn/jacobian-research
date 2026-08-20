from pathlib import Path
import argparse
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("vectors")
p.add_argument("norms")
p.add_argument("--log", default=None)
args = p.parse_args()

V = np.load(args.vectors, mmap_mode="r")
N = np.load(args.norms, mmap_mode="r")

assert len(V) == len(N)

print("POOL")
print(" vectors    =", V.shape)
print(" norm count =", len(N))
print()

for q in [0, .001, .01, .05, .10, .25, .50, .75, .90, .95, .99, .999, 1]:
    print(f"q={q:>5.3f} norm={np.quantile(N,q):.12g}")

print()
print("THRESHOLDS")

for x in [
    50, 60, 70, 80, 90, 100, 110, 120,
    140, 160, 180, 200, 220, 250, 300, 350
]:
    n = int(np.searchsorted(np.sort(N), x, side="right"))
    if n:
        print(f"norm<={x:<3} count={n}")

# Coefficient/support statistics on a sample.
rng = np.random.default_rng(17)
k = min(100000, len(V))
I = rng.choice(len(V), size=k, replace=False)

S = np.asarray(V[I])

support = np.count_nonzero(S, axis=1)
maxcoeff = np.max(np.abs(S), axis=1)

print()
print("VECTOR COMPLEXITY (100k sample)")
print(" support median =", float(np.median(support)))
print(" support p90    =", float(np.quantile(support,.9)))
print(" support max    =", int(support.max()))
print(" maxcoeff median=", float(np.median(maxcoeff)))
print(" maxcoeff p90   =", float(np.quantile(maxcoeff,.9)))
print(" maxcoeff max   =", int(maxcoeff.max()))

print()
print("SHORTEST 20")

for i in np.argsort(N)[:20]:
    print(
        f"norm={float(N[i]):.12g}",
        np.asarray(V[i]).tolist()
    )
