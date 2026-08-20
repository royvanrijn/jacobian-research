#!/usr/bin/env python3
from pathlib import Path
import argparse, collections

ap=argparse.ArgumentParser()
ap.add_argument("file",nargs="?",default="artifacts/local/elkies-k3/heegner-orbits.txt")
ap.add_argument("--max-rank",type=int,default=100)
a=ap.parse_args()

rows=[]
for ln in Path(a.file).read_text().splitlines():
    if not ln.startswith("HEEGNER|rank="): continue
    f={}
    for x in ln.split("|")[1:]:
        k,v=x.split("=",1); f[k]=v
    rows.append(f)

print(f"HEEGNER_SUMMARY|rows={len(rows)}")
print("Smallest primitive CM/order discriminants:")
for r in sorted(rows,key=lambda q:(abs(int(q["order_disc"])),abs(int(q["norm"])),int(q["comp_det"])))[:30]:
    print("CMEND|rank=%s|norm=%s|div=%s|comp_det=%s|order_disc=%s|field_disc=%s|conductor=%s|binary=%s" %
          (r["rank"],r["norm"],r["div"],r["comp_det"],r["order_disc"],r["field_disc"],r["conductor"],r["binary"]))
