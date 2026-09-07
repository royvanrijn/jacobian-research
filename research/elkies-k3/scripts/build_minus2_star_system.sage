from sage.all import *
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "data/k3-model"
OUT.mkdir(parents=True,exist_ok=True)

# Hub section P5.
hub = 5
neighbors = [2,3,4,6,8,10,11,15]

# ------------------------------------------------------------
# Variables
#
# A deg 8                 9
# B deg 12               13
# hub x deg4              5
# hub y deg6              7
#
# each neighbor:
#   x_i deg4              5
#   slope m_i deg2        3
#
# y_i is DERIVED:
#   y_i = y + m_i (x_i-x)
# ------------------------------------------------------------

names=[]

names += [f"a{k}" for k in range(9)]
names += [f"b{k}" for k in range(13)]

names += [f"x_{k}" for k in range(5)]
names += [f"y_{k}" for k in range(7)]

for i in neighbors:
    names += [f"x{i}_{k}" for k in range(5)]
    names += [f"m{i}_{k}" for k in range(3)]

R=PolynomialRing(QQ,names)
g=R.gens_dict()

Rt=PolynomialRing(R,'t')
t=Rt.gen()

def poly(prefix,d):
    return sum(
        g[f"{prefix}_{k}"]*t^k
        for k in range(d+1)
    )

A=sum(g[f"a{k}"]*t^k for k in range(9))
B=sum(g[f"b{k}"]*t^k for k in range(13))

x=poly("x",4)
y=poly("y",6)

equations=[]

# Hub lies on curve.
F=y^2-x^3-A*x-B

for k in range(13):
    equations.append(F[k])

neighbor_data={}

for i in neighbors:

    xi=poly(f"x{i}",4)
    mi=poly(f"m{i}",2)

    yi=y + mi*(xi-x)

    # Neighbor lies on same curve.
    Fi=yi^2-xi^3-A*xi-B

    for k in range(13):
        equations.append(Fi[k])

    # Its sum with the hub.
    xsi=mi^2-x-xi

    # y-coordinate:
    # line through P,Q: Y = y + m(X-x)
    # third intersection is -(P+Q)
    # so y(P+Q) = -y - m*(xsum-x)
    ysi=-y-mi*(xsi-x)

    # Sanity: sum section must lie on E identically.
    # This should follow algebraically from P,Q and slope relation;
    # retain for diagnostics rather than duplicate equations.
    Fsum=ysi^2-xsi^3-A*xsi-B

    neighbor_data[i]=(xi,yi,mi,xsi,ysi,Fsum)

print("hub =",hub)
print("neighbors =",neighbors)

print()
print("unknowns =",R.ngens())
print("equations =",len(equations))
print("naive dimension =",R.ngens()-len(equations))

print()
print("degree checks")

for i,(xi,yi,mi,xsi,ysi,Fsum) in neighbor_data.items():
    print(
        i,
        "deg xi",xi.degree(),
        "deg yi",yi.degree(),
        "deg slope",mi.degree(),
        "deg xsum",xsi.degree(),
        "deg ysum",ysi.degree(),
        "sum_identity_degree",Fsum.degree()
    )

(OUT/"minus2-star-system-summary.txt").write_text(
    f"hub={hub}\n"
    f"neighbors={','.join(map(str,neighbors))}\n"
    f"unknowns={R.ngens()}\n"
    f"equations={len(equations)}\n"
    f"naive_dimension={R.ngens()-len(equations)}\n"
)
