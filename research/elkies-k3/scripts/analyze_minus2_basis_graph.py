from pathlib import Path
import numpy as np
from collections import defaultdict, deque

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

G = np.loadtxt(
    J / "short_vector_basis_gram.txt",
    dtype=np.int64
)

n=len(G)

adj=defaultdict(list)
edges=[]

for i in range(n):
    for j in range(i):
        if G[i,j] == -2:
            adj[i].append(j)
            adj[j].append(i)
            edges.append((i,j))

print("vertices =",n)
print("-2 edges =",len(edges))

print()
print("degrees:")
for i in range(n):
    print(i,len(adj[i]),sorted(adj[i]))

# connected components
seen=set()
components=[]

for start in range(n):
    if start in seen:
        continue

    q=[start]
    seen.add(start)
    c=[]

    while q:
        v=q.pop()
        c.append(v)

        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                q.append(w)

    components.append(sorted(c))

print()
print("components:")
for c in components:
    print(len(c),c)

# cycles / cyclomatic number
print()
print(
    "cyclomatic number =",
    len(edges)-n+len(components)
)

# dump Graph6-compatible edge list
out = BASE / "data/k3-model/minus2_basis_graph.tsv"
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(
    "i\tj\n" +
    "\n".join(f"{i}\t{j}" for i,j in edges) +
    "\n"
)

print("saved",out)
