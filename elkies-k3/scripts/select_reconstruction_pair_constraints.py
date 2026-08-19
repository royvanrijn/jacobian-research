from pathlib import Path
import numpy as np
from collections import Counter

BASE=Path(__file__).resolve().parents[1]
J=BASE/"checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT=BASE/"data/k3-model"
OUT.mkdir(parents=True,exist_ok=True)

G=np.loadtxt(
    J/"short_vector_basis_gram.txt",
    dtype=np.int64
)

n=17

print("PAIRING COUNTS")

C=Counter()

for i in range(n):
    for j in range(i):
        C[int(G[i,j])] += 1

for x,c in sorted(C.items()):
    print(x,c)

# ------------------------------------------------------------
# Build a spanning tree, preferring extreme pairings because
# these generally impose strongest intersection information.
# ------------------------------------------------------------

edges=[]

for i in range(n):
    for j in range(i):
        edges.append((
            -abs(int(G[i,j])),
            int(G[i,j]),
            i,j
        ))

edges.sort()

parent=list(range(n))

def find(x):
    while parent[x]!=x:
        parent[x]=parent[parent[x]]
        x=parent[x]
    return x

def union(a,b):
    a=find(a)
    b=find(b)
    if a==b:
        return False
    parent[b]=a
    return True

selected=[]

for _,pairing,i,j in edges:
    if union(i,j):
        selected.append((i,j,pairing))

assert len(selected)==16

# Add examples of pairing values not represented yet.
present={p for _,_,p in selected}

for value in sorted(C):
    if value in present:
        continue

    for i in range(n):
        found=False
        for j in range(i):
            if G[i,j]==value:
                selected.append((i,j,int(value)))
                found=True
                break
        if found:
            break

print()
print("SELECTED PAIR CONSTRAINTS")

for i,j,p in selected:
    print(i,j,p)

print()
print("count =",len(selected))

(OUT/"selected_pair_constraints.tsv").write_text(
    "i\tj\tpairing\n" +
    "\n".join(
        f"{i}\t{j}\t{p}"
        for i,j,p in selected
    ) + "\n"
)
