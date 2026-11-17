"""Exact coordinate-graph slicing of the published cubic carrier map."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
HERE=Path(__file__).resolve().parent
def clean(p):return {e:c for e,c in p.items() if c}
def add(a,b,scale=Q(1)):
    r=dict(a)
    for e,c in b.items():r[e]=r.get(e,Q(0))+scale*c
    return clean(r)
def mul(a,b):
    r={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(x+y for x,y in zip(e,f));r[g]=r.get(g,Q(0))+c*d
    return clean(r)
def unit(j,n=19):
    e=[0]*n;e[j]=1;return tuple(e)
def substitute(p,j,q):
    powers=[{(0,)*19:Q(1)}]
    for _ in range(max((e[j] for e in p),default=0)):powers.append(mul(powers[-1],q))
    r={}
    for e,c in p.items():
        f=list(e);f[j]=0
        r=add(r,mul({tuple(f):c},powers[e[j]]))
    return r
def load():
    d=json.loads((HERE/"sources/degree3-map-external.json").read_text())
    names=["x","y","z"]+[f"w{i}" for i in range(1,17)]
    rows={}
    for i,row in enumerate(d["components"]):
        p={}
        for text,c in row.items():
            e=[0]*19
            for t in text.split("*"):
                a=t.split("^");e[names.index(a[0])]+=int(a[1]) if len(a)>1 else 1
            p[tuple(e)]=Q(c)
        rows[i]=p
    return rows
def eliminate(rows,j):
    assert rows[j].get(unit(j))==1
    f=add(rows[j],{unit(j):Q(1)},Q(-1))
    assert not any(e[j] for e in f)
    return {i:substitute(p,j,{e:-c for e,c in f.items()}) for i,p in rows.items() if i!=j}
def main():
    rows=load();out=[]
    for j in range(3,19):
        child=eliminate(rows,j)
        row=dict(removed=j,dimension=18,degree=max(sum(e) for p in child.values() for e in p),
                 nonlinear=sum(sum(e)>=2 for p in child.values() for e in p))
        out.append(row);print(row)
    receipt=dict(status="complete-single-slices",source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=out)
    (HERE/"single-slices.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

