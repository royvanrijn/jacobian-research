from sage.all import *
from pathlib import Path
from collections import Counter, defaultdict

BASE = Path(__file__).resolve().parents[1]
J = BASE/"checkpoints/24A1-JACKPOT-948-r0-n1311"

def readmat(p):
    rows=[]
    for line in p.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        line=line.replace("["," ").replace("]"," ")
        try:
            r=[ZZ(x) for x in line.split()]
        except Exception:
            continue
        if r:
            rows.append(r)
    return matrix(ZZ,rows)

H=readmat(J/"rank17_gram.txt")
n=17

coeff=[]
for i in range(n):
    for j in range(i,n):
        coeff.append(H[i,i]//2 if i==j else H[i,j])

Q=QuadraticForm(ZZ,n,coeff)
half=[vector(ZZ,v)
      for v in Q.short_vector_list_up_to_length(3,True)[2]]

assert len(half)==1311

V=[]
for v in half:
    V.append(v)
    V.append(-v)

assert len(V)==2622

def key(v):
    return tuple(map(int,v))

index={key(v):i for i,v in enumerate(V)}

# ------------------------------------------------------------
# Enumerate all signed additive triples:
#
#       a + b = c
#
# with all three minimal.
# ------------------------------------------------------------

triples=[]
degree=[0]*len(V)

for i,a in enumerate(V):
    for j in range(i+1,len(V)):
        b=V[j]

        # Necessary and sufficient for a+b to have norm 4:
        # <a,b> = -2.
        if ZZ(a*H*b) != -2:
            continue

        c=a+b
        k=index.get(key(c))

        if k is None:
            raise RuntimeError("norm-4 relation missing from shell")

        triples.append((i,j,k))

        degree[i]+=1
        degree[j]+=1
        degree[k]+=1

print("signed vectors =",len(V))
print("additive triples =",len(triples))

dc=Counter(degree)

print()
print("TRIPLE-DEGREE HISTOGRAM")
for d,c in sorted(dc.items()):
    print(d,c)

top=sorted(
    range(len(V)),
    key=lambda i:(-degree[i],i)
)

print()
print("TOP RELATION-RICH VECTORS")
for i in top[:30]:
    print(
        "i",i,
        "degree",degree[i],
        "v",V[i]
    )

# ------------------------------------------------------------
# Greedy construction of a small independent set maximizing
# how many additive relations become internal.
# ------------------------------------------------------------

chosen=[]
chosen_set=set()

def internal_relations(S):
    S=set(S)
    return sum(
        1 for a,b,c in triples
        if a in S and b in S and c in S
    )

# Seed at highest relation degree.
chosen=[top[0]]
chosen_set=set(chosen)

while len(chosen)<10:
    candidates=[]

    for j in top:
        if j in chosen_set:
            continue

        B=matrix(ZZ,[V[k] for k in chosen+[j]])

        if B.rank() != len(chosen)+1:
            continue

        S=chosen+[j]

        # score:
        # 1) relations fully contained
        # 2) relations having two vertices contained (one step from closure)
        SS=set(S)

        internal=0
        two=0

        for a,b,c in triples:
            q=(a in SS)+(b in SS)+(c in SS)
            if q==3:
                internal+=1
            elif q==2:
                two+=1

        candidates.append(
            (-internal,-two,-degree[j],j)
        )

    if not candidates:
        break

    candidates.sort()
    j=candidates[0][-1]
    chosen.append(j)
    chosen_set.add(j)

    print(
        "choose",
        len(chosen),
        "index",j,
        "degree",degree[j],
        "internal",internal_relations(chosen)
    )

print()
print("CHOSEN INDICES")
print(chosen)

B=matrix(ZZ,[V[i] for i in chosen])
G=B*H*B.T

print()
print("CHOSEN GRAM")
print(G)

# ------------------------------------------------------------
# Closure: repeatedly add c whenever a+b=c and two of a,b,c
# are already known.  This tells us how much of the 2622-shell
# a small seed determines by group-law relations.
# ------------------------------------------------------------

closure=set(chosen)

changed=True
while changed:
    changed=False

    for a,b,c in triples:
        q=[a in closure,b in closure,c in closure]

        # Any two determine the third, with signs already encoded.
        if sum(q)>=2:
            for x in (a,b,c):
                if x not in closure:
                    closure.add(x)
                    changed=True

print()
print("closure size =",len(closure),"of",len(V))

# Save.
(J/"minimal_additive_triples.txt").write_text(
    "\n".join(
        f"{a} {b} {c}"
        for a,b,c in triples
    )+"\n"
)

(J/"relation_seed_indices.txt").write_text(
    " ".join(map(str,chosen))+"\n"
)

(J/"relation_seed_vectors.txt").write_text(
    "\n".join(
        " ".join(map(str,V[i]))
        for i in chosen
    )+"\n"
)

(J/"relation_seed_gram.txt").write_text(
    "\n".join(
        " ".join(map(str,row))
        for row in G.rows()
    )+"\n"
)

(J/"relation_closure_indices.txt").write_text(
    "\n".join(map(str,sorted(closure)))+"\n"
)
