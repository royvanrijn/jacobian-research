from sage.all import *
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
D = BASE / "data/rank29"

a1 = ZZ(1)
a2 = ZZ(0)
a3 = ZZ(0)
a4 = ZZ(-27006183241630922218434652145297453784768054621836357954737385)
a6 = ZZ(55258058551342376475736699591118191821521067032535079608372404779149413277716173425636721497)

E = EllipticCurve(QQ,[a1,a2,a3,a4,a6])

print("curve =", E)
print("discriminant =", E.discriminant())

xs = [
    QQ(line.strip())
    for line in (D/"E29-xcoords.txt").read_text().splitlines()
    if line.strip()
]

assert len(xs) == 29

points=[]

for i,x in enumerate(xs):
    pts = E.lift_x(x, all=True)

    if not pts:
        raise RuntimeError(f"x[{i}] does not lift: {x}")

    # pick one sign consistently
    P = pts[0]
    points.append(P)

    print(
        f"P{i:02d}",
        "x=",x,
        "y=",P[1]
    )

print()
print("computing canonical heights ...")

H = matrix(RR,29,29)

for i in range(29):
    hi = points[i].height()

    H[i,i] = hi

    for j in range(i):
        # Bilinear pairing from quadratic canonical height.
        hij = (
            (points[i] + points[j]).height()
            - points[i].height()
            - points[j].height()
        ) / 2

        H[i,j] = hij
        H[j,i] = hij

    print("row",i,"done",flush=True)

print()
print("det =", H.det())
print("eigenvalue min =", min(H.eigenvalues()))

(D/"E29-height-gram.txt").write_text(
    "\n".join(
        " ".join(f"{RR(x):.18g}" for x in row)
        for row in H.rows()
    ) + "\n"
)

(D/"E29-points.txt").write_text(
    "\n".join(
        f"{i}\t{P[0]}\t{P[1]}"
        for i,P in enumerate(points)
    ) + "\n"
)

print("saved height Gram")
