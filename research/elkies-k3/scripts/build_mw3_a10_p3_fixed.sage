from sage.all import *
from pathlib import Path
import argparse


ap = argparse.ArgumentParser(
    description="Build a fixed-surface canonical P3 system on a verified GF(31) A10/MW3 P1+P2 seed."
)
ap.add_argument("--surface", type=int, choices=(2, 4), default=4)
ap.add_argument("--out", required=True)
args = ap.parse_args()

K = GF(31)
names = ["r"] + [f"q{i}" for i in range(4)] + [f"v{i}" for i in range(7)]
R = PolynomialRing(K, names, order="degrevlex")
d = R.gens_dict()
Rt = PolynomialRing(R, "t")
t = Rt.gen()

r = d["r"]
q = [d[f"q{i}"] for i in range(4)]
v = [d[f"v{i}"] for i in range(7)]

surfaces = {
    2: {
        "A": [4, 28, 7, 18, 23, 26, 15, 7, 1],
        "B": [23, 9, 18, 12, 2, 15, 3, 0, 25, 18, 25, 5, 30],
        "lam": 27,
        "nodes": [3, 22, 4],
        "sinf": 17,
    },
    4: {
        "A": [4, 23, 18, 28, 12, 18, 23, 20, 19],
        "B": [23, 24, 27, 10, 6, 8, 23, 26, 26, 16, 19, 9, 15],
        "lam": 23,
        "nodes": [3, 21, 10],
        "sinf": 29,
    },
}
data = surfaces[args.surface]
A = Rt(data["A"])
B = Rt(data["B"])
lam = K(data["lam"])
s0 = K(data["nodes"][0])
sl = K(data["nodes"][2])
sinf = K(data["sinf"])

# P3 is nonidentity at 0 and lambda.  Interpolate those two X-values, then add
# the vanishing polynomial t(t-lambda).  The leading coefficient is fixed by
# the infinity node incidence x6=sinf.
C = s0 * r**2 * (t - lam) / (-lam)
C += sl * (K(lam) - r)**2 * t / lam
G = t * (t - lam)
X = C + G * (q[0] + q[1] * t + q[2] * t**2 + q[3] * t**3 + sinf * t**4)

# Nonidentity at 0 and lambda also gives Y=0 there.  Infinity class 10~1 has
# y9=0 but generically y8!=0, so Y has degree at most eight.
Y = G * sum(v[i] * t**i for i in range(7))
z = t - r
S = Y**2 - X**3 - A * X * z**4 - B * z**6
assert S.degree() <= 16

equations = [R(S[k]) for k in range(17) if S[k] != 0]
out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w") as handle:
    handle.write(",".join(names) + "\n31\n")
    for i, equation in enumerate(equations):
        handle.write(str(equation).replace("**", "^"))
        handle.write(",\n" if i + 1 < len(equations) else "\n")

meta = out.with_suffix(".meta.txt")
with meta.open("w") as handle:
    handle.write(f"surface={args.surface}\n")
    handle.write(f"open=r*(r-1)*(r-{int(lam)})*v6 != 0\n")
    handle.write("X=" + str(X) + "\n")
    handle.write("Y=" + str(Y) + "\n")

print(
    f"MW3A10P3|surface={args.surface}|vars={len(names)}"
    f"|eqs={len(equations)}|out={out}|meta={meta}",
    flush=True,
)
for i, equation in enumerate(equations):
    print(
        f"MW3A10P3_EQ|i={i}|degree={equation.total_degree()}"
        f"|terms={len(equation.monomials())}",
        flush=True,
    )
