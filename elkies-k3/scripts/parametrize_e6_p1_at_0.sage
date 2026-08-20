from sage.all import *
from pathlib import Path
import argparse

ap=argparse.ArgumentParser()
ap.add_argument("--input",required=True)
ap.add_argument("--out",default=None)
args=ap.parse_args()

lines=[x.strip() for x in Path(args.input).read_text().splitlines() if x.strip()]
names=[x.strip() for x in lines[0].split(",") if x.strip()]
p=int(lines[1])
K=GF(p)

R=PolynomialRing(K,names,order="degrevlex")
gd=R.gens_dict()

eqs=[]
for line in lines[2:]:
    s=line[:-1] if line.endswith(",") else line
    eqs.append(R(s.replace("^","**")))

# Current exporter ordering:
# 0..7 fiber equations, 8..11 = P1_0..P1_3.
if len(eqs)!=12:
    raise RuntimeError(f"expected 12 unsliced equations, got {len(eqs)}")
if not all(n in gd for n in ("x0","y0","s0")):
    raise RuntimeError("expected x0,y0,s0 in remaining ring")

new_names=[n for n in names if n not in ("x0","y0")] + ["r0"]
S=PolynomialRing(K,new_names,order="degrevlex")
sd=S.gens_dict()
r0=sd["r0"]

images=[]
for n in names:
    if n=="x0":
        images.append(r0**2-2*sd["s0"])
    elif n=="y0":
        images.append(r0*(r0**2-3*sd["s0"]))
    else:
        images.append(sd[n])

phi=R.hom(images,S)
E=[phi(e) for e in eqs]

print(f"E6PARAM0|p={p}|vars={S.ngens()}|eqs_before={len(E)}",flush=True)
print(f"E6PARAM0|P1_0_zero={E[8]==0}",flush=True)

if E[8]!=0:
    print(f"E6PARAM0|P1_0_degree={E[8].degree()}|terms={len(E[8].monomials())}",flush=True)
    raise RuntimeError("parametrization did not kill P1_0")

tags=["I4_1_D2","I4_1_D3","I2_lam_A","I2_lam_B","I2_lam_d1",
      "I2_mu_A","I2_mu_B","I2_mu_d1","P1_0","P1_1","P1_2","P1_3"]

for idx in (9,10,11):
    e=E[idx]
    linear=[]
    for n in new_names:
        v=sd[n]
        if e.degree(v)==1:
            c=e.derivative(v)
            linear.append((n,len(c.monomials()),e.degree(),len(e.monomials())))
    print(f"E6PARAM0_NEXT|eq={tags[idx]}|degree={e.degree()}|terms={len(e.monomials())}|linear={linear}",flush=True)

E2=[e for i,e in enumerate(E) if i!=8]
T2=[t for i,t in enumerate(tags) if i!=8]
print(f"E6PARAM0|vars={S.ngens()}|eqs={len(E2)}|naive_dim={S.ngens()-len(E2)}",flush=True)

for i,(tag,e) in enumerate(zip(T2,E2)):
    print(f"E6PARAM0_EQ|i={i}|tag={tag}|degree={e.degree()}|terms={len(e.monomials())}",flush=True)

if args.out:
    out=Path(args.out)
    with out.open("w") as h:
        h.write(",".join(new_names)+"\n"+str(p)+"\n")
        for i,e in enumerate(E2):
            h.write(str(e).replace("**","^"))
            h.write(",\n" if i+1<len(E2) else "\n")
    print(f"E6PARAM0|export={out}",flush=True)
