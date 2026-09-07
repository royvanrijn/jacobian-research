from sage.all import *
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "data/k3-model"
OUT.mkdir(parents=True,exist_ok=True)

# ------------------------------------------------------------
# Unknowns:
#
# A(t) = a0 + ... + a8 t^8
# B(t) = b0 + ... + b12 t^12
#
# 17 minimal sections:
# xi degree <=4
# yi degree <=6
# ------------------------------------------------------------

names=[]

names += [f"a{i}" for i in range(9)]
names += [f"b{i}" for i in range(13)]

for s in range(17):
    names += [f"x{s}_{i}" for i in range(5)]
    names += [f"y{s}_{i}" for i in range(7)]

R = PolynomialRing(QQ,names)
g = R.gens_dict()

Rt = PolynomialRing(R,'t')
t = Rt.gen()

A=sum(g[f"a{i}"]*t^i for i in range(9))
B=sum(g[f"b{i}"]*t^i for i in range(13))

equations=[]

for s in range(17):

    x=sum(
        g[f"x{s}_{i}"]*t^i
        for i in range(5)
    )

    y=sum(
        g[f"y{s}_{i}"]*t^i
        for i in range(7)
    )

    F=y^2-x^3-A*x-B

    for c in F.list():
        equations.append(c)

print("unknowns =",R.ngens())
print("section equations =",len(equations))

# Don't try to Groebner this directly.
# This file is just the canonical symbolic problem description.

(OUT/"rank17-system-summary.txt").write_text(
    f"unknowns={R.ngens()}\n"
    f"equations={len(equations)}\n"
    f"A_degree=8\n"
    f"B_degree=12\n"
    f"sections=17\n"
    f"x_degree=4\n"
    f"y_degree=6\n"
)

print("saved system summary")
