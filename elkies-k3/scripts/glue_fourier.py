from pathlib import Path
from collections import Counter
import cmath
import math

BASE = Path(__file__).resolve().parents[1]
GDIR = BASE/"data/glue"

mods=(2,2,2,2,2,20)

weights={}
for line in (GDIR/"minimal_vector_cosets.txt").read_text().splitlines():
    lhs,rhs=line.split(":")
    g=tuple(map(int,lhs.split()))
    weights[g]=len(rhs.split())

elts=[]
for a in range(2):
 for b in range(2):
  for c in range(2):
   for d in range(2):
    for e in range(2):
     for f in range(20):
      elts.append((a,b,c,d,e,f))

def character(h,g):
    phase=0.0
    for x,y,m in zip(h,g,mods):
        phase += x*y/m
    return cmath.exp(2j*math.pi*phase)

def transform(function):
    out={}
    for h in elts:
        z=sum(function.get(g,0)*character(h,g).conjugate()
              for g in elts)
        out[h]=z
    return out

support={g:1 for g in weights}
Fw=transform(weights)
Fs=transform(support)

def spectrum(F, digits=8):
    C=Counter()
    for z in F.values():
        key=round(abs(z),digits)
        C[key]+=1
    return C

print("support size =",len(support))
print("weight sum =",sum(weights.values()))

print("\nSUPPORT FOURIER ABS SPECTRUM")
for x,n in sorted(spectrum(Fs).items()):
    print(x,n)

print("\nWEIGHTED FOURIER ABS SPECTRUM")
for x,n in sorted(spectrum(Fw).items()):
    print(x,n)

# Exact-ish real/imag grouping, useful if values look integral/surd-like.
print("\nLargest weighted coefficients:")
for h,z in sorted(Fw.items(),key=lambda x:-abs(x[1]))[:40]:
    print(h,z,abs(z))
