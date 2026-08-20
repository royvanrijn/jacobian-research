#!/usr/bin/env python3
from pathlib import Path
import argparse, random

ap=argparse.ArgumentParser(description="Square an overdetermined msolve file using generic GF(p) row combinations.")
ap.add_argument("--input",required=True)
ap.add_argument("--seed",type=int,default=1)
ap.add_argument("--output",required=True)
args=ap.parse_args()

src=Path(args.input)
lines=src.read_text().splitlines()
if len(lines)<3:
    raise SystemExit("input too short")

varline=lines[0].strip()
p=int(lines[1].strip())
vars=[x.strip() for x in varline.split(",") if x.strip()]
n=len(vars)

# Exporter format is one equation per line, with a trailing comma except perhaps final line.
eqs=[]
for raw in lines[2:]:
    s=raw.strip()
    if not s:
        continue
    if s.endswith(","):
        s=s[:-1].rstrip()
    if s:
        eqs.append(s)

m=len(eqs)
print(f"E6SQUARE|stage=input|vars={n}|eqs={m}|p={p}",flush=True)

if m<n:
    raise SystemExit(f"need >= {n} equations, got {m}")

random.seed(args.seed)

# n generic linear combinations of the m equations over GF(p).
C=[]
for i in range(n):
    row=[random.randrange(p) for _ in range(m)]
    if all(c==0 for c in row):
        row[i % m]=1
    C.append(row)

# Make sure every source equation appears somewhere.
for j in range(m):
    if all(C[i][j]==0 for i in range(n)):
        C[j % n][j]=1

combined=[]
for row in C:
    terms=[]
    for c,e in zip(row,eqs):
        c%=p
        if c==0:
            continue
        if c==1:
            terms.append(f"({e})")
        else:
            terms.append(f"{c}*({e})")
    combined.append("+".join(terms) if terms else "0")

out=Path(args.output)
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(varline+"\n")
    h.write(str(p)+"\n")
    for i,e in enumerate(combined):
        h.write(e)
        h.write(",\n" if i+1<len(combined) else "\n")

meta=out.with_suffix(".square.meta.txt")
with meta.open("w") as h:
    h.write(f"source={src}\n")
    h.write(f"vars={n}\nsource_eqs={m}\nsquare_eqs={n}\nseed={args.seed}\n")
    h.write("matrix=\n")
    for row in C:
        h.write(" ".join(map(str,row))+"\n")

print(f"E6SQUARE|stage=done|square_eqs={n}|out={out}|meta={meta}",flush=True)
