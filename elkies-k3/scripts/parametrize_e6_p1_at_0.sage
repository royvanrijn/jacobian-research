from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser(description="Parametrize P1_0 on the E6 four-elimination system.")
ap.add_argument("--input",required=True)
ap.add_argument("--out",required=True)
args=ap.parse_args()

lines=[x.strip() for x in Path(args.input).read_text().splitlines() if x.strip()]
names=[x.strip() for x in lines[0].split(",") if x.strip()]
p=int(lines[1]); K=GF(p)
R=PolynomialRing(K,names,order="degrevlex"); gd=R.gens_dict()
eqs=[R((x[:-1] if x.endswith(",") else x).replace("^","**")) for x in lines[2:]]

if len(eqs)!=12:
    raise RuntimeError(f"expected 12 unsliced equations, got {len(eqs)}")
if not all(n in gd for n in ("x0","y0","s0")):
    raise RuntimeError("expected x0,y0,s0 in remaining ring")

new_names=[n for n in names if n not in ("x0","y0")]+["r0"]
S=PolynomialRing(K,new_names,order="degrevlex"); sd=S.gens_dict(); r0=sd["r0"]
images=[]
for n in names:
    if n=="x0": images.append(r0**2-2*sd["s0"])
    elif n=="y0": images.append(r0*(r0**2-3*sd["s0"]))
    else: images.append(sd[n])
phi=R.hom(images,S)
E=[phi(e) for e in eqs]

print(f"E6PARAM0|p={p}|vars={S.ngens()}|eqs_before={len(E)}",flush=True)
print(f"E6PARAM0|P1_0_zero={E[8]==0}",flush=True)
if E[8]!=0: raise RuntimeError("parametrization did not kill P1_0")

for idx,tag in ((9,"P1_1"),(10,"P1_2"),(11,"P1_3")):
    e=E[idx]; linear=[]
    for n in new_names:
        v=sd[n]
        if e.degree(v)==1:
            linear.append((n,len(e.derivative(v).monomials())))
    print(f"E6PARAM0_NEXT|eq={tag}|degree={e.degree()}|terms={len(e.monomials())}|linear={linear}",flush=True)

E2=[e for i,e in enumerate(E) if i!=8]
out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    h.write(",".join(new_names)+"\n"+str(p)+"\n")
    for i,e in enumerate(E2):
        h.write(str(e).replace("**","^")); h.write(",\n" if i+1<len(E2) else "\n")
print(f"E6PARAM0|vars={S.ngens()}|eqs={len(E2)}|naive_dim={S.ngens()-len(E2)}|out={out}",flush=True)
