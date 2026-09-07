from sage.all import *
from pathlib import Path
from collections import Counter, defaultdict

BASE = Path(__file__).resolve().parents[1]
J = BASE/"checkpoints/24A1-JACKPOT-948-r0-n1311"
p = J/"rank17_gram.txt"

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

H=matrix(ZZ,rows)
n=17

coeff=[]
for i in range(n):
    for j in range(i,n):
        coeff.append(H[i,i]//2 if i==j else H[i,j])

Q=QuadraticForm(ZZ,n,coeff)
S=Q.short_vector_list_up_to_length(3,True)
V=[vector(ZZ,v) for v in S[2]]

print("minimal +/- pairs:",len(V))
assert len(V)==1311

# Pairing distribution between distinct representative minimal vectors.
cnt=Counter()

for i in range(len(V)):
    for j in range(i):
        z=ZZ(V[i]*H*V[j])
        cnt[z]+=1

print()
print("PAIRING DISTRIBUTION")
for k in sorted(cnt):
    print(k,cnt[k])

# For each vector, record numbers of other representatives by |pairing|.
profiles=[]

for i,v in enumerate(V):
    c=Counter()
    for j,w in enumerate(V):
        if i==j:
            continue
        z=abs(ZZ(v*H*w))
        c[z]+=1
    profiles.append((tuple(sorted(c.items())),i))

pc=Counter(p for p,i in profiles)

print()
print("NUMBER OF DISTINCT ABS-PAIRING PROFILES:",len(pc))

for profile,num in pc.most_common(20):
    print("count",num,"profile",profile)

# Find vectors with rarest profiles: likely useful symmetry anchors.
rarest=sorted(
    profiles,
    key=lambda x:(pc[x[0]],x[0])
)

print()
print("RARE PROFILE VECTORS")
for profile,i in rarest[:30]:
    print("i",i,"multiplicity",pc[profile],
          "v",V[i],"profile",profile)

# Search for small independent sets having unusually sparse/simple Gram.
# Score off-diagonal entries by complexity: prefer 0, +/-1 over +/-2 etc.

def complexity(G):
    vals=[]
    for i in range(G.nrows()):
        for j in range(i):
            vals.append(abs(ZZ(G[i,j])))
    return (
        max(vals,default=0),
        sum(x*x for x in vals),
        sum(vals)
    )

# Greedily seed from each of the first 100 rare vectors and construct 6-set.
best=[]

seed_indices=[i for _,i in rarest[:100]]

for seed in seed_indices:
    chosen=[seed]
    while len(chosen)<6:
        candidates=[]
        for j in range(len(V)):
            if j in chosen:
                continue
            B=matrix(ZZ,[V[k] for k in chosen+[j]])
            if B.rank()!=len(chosen)+1:
                continue
            G=B*H*B.T
            candidates.append((complexity(G),j))
        if not candidates:
            break
        candidates.sort()
        chosen.append(candidates[0][1])

    if len(chosen)==6:
        B=matrix(ZZ,[V[k] for k in chosen])
        G=B*H*B.T
        best.append((complexity(G),chosen,G))

best.sort(key=lambda x:x[0])

print()
print("BEST 6-VECTOR CONFIGURATIONS")
for rank,(score,inds,G) in enumerate(best[:10]):
    print()
    print("candidate",rank,"score",score,"indices",inds)
    print(G)

# Save best configuration.
score,inds,G=best[0]

(J/"reconstruction_anchor_indices.txt").write_text(
    " ".join(map(str,inds))+"\n"
)

(J/"reconstruction_anchor_vectors.txt").write_text(
    "\n".join(
        " ".join(map(str,V[i]))
        for i in inds
    )+"\n"
)

(J/"reconstruction_anchor_gram.txt").write_text(
    "\n".join(
        " ".join(map(str,row))
        for row in G.rows()
    )+"\n"
)
