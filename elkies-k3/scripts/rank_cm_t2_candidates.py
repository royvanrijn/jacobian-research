#!/usr/bin/env python3
from pathlib import Path
import argparse, ast, collections, math

ap=argparse.ArgumentParser()
ap.add_argument("input", nargs="?", default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
ap.add_argument("--top", type=int, default=50)
ap.add_argument("--out", default="artifacts/local/elkies-k3/cm-t2-ranked.txt")
a=ap.parse_args()

records=[]
for line in Path(a.input).read_text().splitlines():
    if not line.startswith("CM|"): continue
    fields={}
    for part in line.split("|")[1:]:
        k,v=part.split("=",1); fields[k]=v
    try:
        T=int(fields["T"])
        cid=int(fields["id"])
        binary=ast.literal_eval(fields["binary"])
        disc=int(fields["disc"])
        v=ast.literal_eval(fields["v"])
        vnorm=int(fields["vnorm"])
        score=ast.literal_eval(fields["score"])
    except Exception as e:
        print("PARSE_FAIL",line,e)
        continue
    records.append(dict(T=T,id=cid,binary=tuple(binary),disc=disc,v=tuple(v),vnorm=vnorm,score=score))

by=collections.defaultdict(list)
for r in records: by[r["binary"]].append(r)

rows=[]
for b,rs in by.items():
    Ts=sorted(set(r["T"] for r in rs))
    discs=sorted(set(r["disc"] for r in rs))
    assert len(discs)==1, (b,discs)
    disc=discs[0]
    minabs=min(abs(r["vnorm"]) for r in rs)
    minbox=min(max(abs(x) for x in r["v"]) for r in rs)
    rows.append((len(Ts),disc,minabs,minbox,b,Ts,rs))

rows.sort(key=lambda x:(-x[0],x[1],x[2],x[3],x[4]))

out=[]
out.append(f"CMRANK|records={len(records)}|binary_classes={len(rows)}|distinct_T={len(set(r['T'] for r in records))}")
for rank,row in enumerate(rows[:a.top],1):
    nT,disc,minabs,minbox,b,Ts,rs=row
    out.append(f"CMRANK|rank={rank}|Tcount={nT}|disc={disc}|binary={b}|min_vnorm={minabs}|min_vbox={minbox}|Ts={','.join(map(str,Ts))}")
    # Best witness
    w=min(rs,key=lambda r:(abs(r["vnorm"]),max(abs(x) for x in r["v"]),r["T"]))
    out.append(f"CMWIT|rank={rank}|T={w['T']}|id={w['id']}|v={w['v']}|vnorm={w['vnorm']}")

# Also aggregate by discriminant only: different binary classes with same CM discriminant.
bd=collections.defaultdict(list)
for row in rows: bd[row[1]].append(row)
dr=[]
for d,xs in bd.items():
    Ts=set()
    for row in xs: Ts.update(row[5])
    dr.append((len(Ts),d,len(xs),sorted(Ts)))
dr.sort(key=lambda x:(-x[0],x[1]))
out.append("")
out.append("DISCRIMINANT AGGREGATION")
for rank,(nt,d,nclasses,Ts) in enumerate(dr[:a.top],1):
    out.append(f"CMDISC|rank={rank}|Tcount={nt}|disc={d}|classes={nclasses}|Ts={','.join(map(str,Ts))}")

Path(a.out).write_text("\n".join(out)+"\n")
print("\n".join(out[:min(len(out),80)]))
print(f"CMRANK|stage=done|out={a.out}")
