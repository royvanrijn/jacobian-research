from sage.all import *
from pathlib import Path
from collections import defaultdict, deque

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / "data/k3-model"
OUT.mkdir(parents=True, exist_ok=True)

# Connected -2 graph on sections 2..16.
edges = [
    (5,2),(5,3),(5,4),(6,5),(8,5),
    (9,3),(9,7),(9,8),(10,5),(10,9),
    (11,5),(11,9),(13,9),(14,3),(14,6),
    (14,7),(14,12),(14,13),(15,5),(15,9),
    (16,14),
]

vertices = sorted(set(sum(([a,b] for a,b in edges), [])))
assert vertices == list(range(2,17))

# Build a spanning tree rooted at 5.
adj=defaultdict(list)
for a,b in edges:
    adj[a].append(b)
    adj[b].append(a)

root=5
parent={root:None}
tree=[]
q=deque([root])

while q:
    a=q.popleft()
    for b in adj[a]:
        if b in parent:
            continue
        parent[b]=a
        tree.append((a,b))
        q.append(b)

assert len(tree)==14

tree_set={tuple(sorted(e)) for e in tree}
cycle_edges=[
    e for e in edges
    if tuple(sorted(e)) not in tree_set
]

print("vertices =",vertices)
print("tree edges =",tree)
print("cycle edges =",cycle_edges)
print("cycle count =",len(cycle_edges))

# ------------------------------------------------------------
# Unknowns
#
# A degree 8, B degree 12
# x_i degree 4 for all 15 sections
# y_root degree 6
#
# For each -2 edge introduce slope m_ab degree 2.
#
# Along tree:
#   y_child = y_parent + m*(x_child-x_parent)
#
# Thus all y_i are derived from root y plus tree slopes.
#
# Non-tree edges impose compatibility:
#   y_b-y_a = m_ab*(x_b-x_a).
# ------------------------------------------------------------

names=[]
names += [f"a{k}" for k in range(9)]
names += [f"b{k}" for k in range(13)]

for i in vertices:
    names += [f"x{i}_{k}" for k in range(5)]

names += [f"y{root}_{k}" for k in range(7)]

for a,b in edges:
    lo,hi=sorted((a,b))
    names += [f"m{lo}_{hi}_{k}" for k in range(3)]

R=PolynomialRing(QQ,names)
g=R.gens_dict()
Rt=PolynomialRing(R,'t')
t=Rt.gen()

def p(prefix,d):
    return sum(g[f"{prefix}_{k}"]*t^k for k in range(d+1))

A=sum(g[f"a{k}"]*t^k for k in range(9))
B=sum(g[f"b{k}"]*t^k for k in range(13))

x={i:p(f"x{i}",4) for i in vertices}
y={root:p(f"y{root}",6)}

def slope(a,b):
    lo,hi=sorted((a,b))
    return p(f"m{lo}_{hi}",2)

# Derive y coordinates down spanning tree.
children=defaultdict(list)
for a,b in tree:
    # tree orientation from BFS root
    if parent[b]==a:
        children[a].append(b)
    elif parent[a]==b:
        children[b].append(a)
    else:
        raise RuntimeError("tree orientation mismatch")

qq=deque([root])
while qq:
    a=qq.popleft()
    for b in children[a]:
        m=slope(a,b)
        y[b]=y[a] + m*(x[b]-x[a])
        qq.append(b)

assert set(y)==set(vertices)

eqs=[]

def coeff_equations(F,maxdeg=None):
    if maxdeg is None:
        maxdeg=F.degree()
    return [F[k] for k in range(maxdeg+1)]

# Root lies on E.
Froot=y[root]^2-x[root]^3-A*x[root]-B
eqs += coeff_equations(Froot,12)

# For each TREE edge, use factorized difference equation.
#
# m*(y_a+y_b) =
# x_a^2 + x_a*x_b + x_b^2 + A
#
# This guarantees the neighboring point lies on the same curve,
# given the parent does.
for a,b in tree:
    m=slope(a,b)

    F = (
        m*(y[a]+y[b])
        - (
            x[a]^2
            + x[a]*x[b]
            + x[b]^2
            + A
        )
    )

    eqs += coeff_equations(F,8)

# Every cycle edge gives:
# 1) slope compatibility (degree <=6): 7 equations
# 2) factorized curve identity (degree <=8): 9 equations
for a,b in cycle_edges:
    m=slope(a,b)

    Fline = y[b]-y[a]-m*(x[b]-x[a])
    eqs += coeff_equations(Fline,6)

    Fcurve = (
        m*(y[a]+y[b])
        - (
            x[a]^2
            + x[a]*x[b]
            + x[b]^2
            + A
        )
    )
    eqs += coeff_equations(Fcurve,8)

print()
print("unknowns =",R.ngens())
print("equations =",len(eqs))
print("naive dimension =",R.ngens()-len(eqs))

print()
print("breakdown:")
print(" A+B             =",22)
print(" x-polynomials    =",15*5)
print(" root y           =",7)
print(" slopes           =",len(edges)*3)
print(" root equations   =",13)
print(" tree equations   =",len(tree)*9)
print(" cycle equations  =",len(cycle_edges)*16)

(OUT/"three-hub-system-summary.txt").write_text(
    "\n".join([
        f"vertices={','.join(map(str,vertices))}",
        f"tree_edges={tree}",
        f"cycle_edges={cycle_edges}",
        f"unknowns={R.ngens()}",
        f"equations={len(eqs)}",
        f"naive_dimension={R.ngens()-len(eqs)}",
    ])+"\n"
)
