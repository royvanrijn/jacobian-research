#!/usr/bin/env python3
"""Independent coordinate-gated global traversal; no production imports.

Uses direct row predicates, reachability (rather than six cut masks), and
per-row set-intersection repair bounds (rather than production case splits).
The bounded traversal is checkpointed; an incomplete pass is not a proof.
"""
from collections import Counter,defaultdict
from itertools import product
from pathlib import Path
import argparse,gzip,hashlib,json,pickle,time

HERE=Path(__file__).resolve().parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def det(a,b,c):
    return (a[0]*b[1]*c[2]+a[1]*b[2]*c[0]+a[2]*b[0]*c[1]
            -a[0]*b[2]*c[1]-a[1]*b[0]*c[2]-a[2]*b[1]*c[0])
def bits(mask):
    while mask:
        bit=mask&-mask
        yield bit
        mask-=bit

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=int,default=300)
    parser.add_argument("--checkpoint",type=Path,default=HERE/"independent.checkpoint.pkl.gz")
    a=parser.parse_args()
    degree,bound=7,9
    exponents=sorted((e for e in product(range(degree+1),repeat=3)
                     if 2<=sum(e)<=degree),key=lambda e:(sum(e),e))
    atoms=[(i,e) for i in range(3) for e in exponents]
    atom_bits={atom:1<<n for n,atom in enumerate(atoms)}
    by_bit={1<<n:atom for n,atom in enumerate(atoms)}
    linear=((1,0,0),(0,1,0),(0,0,1))
    all_atoms=(1<<len(atoms))-1
    row_bits=[sum(bit for (i,e),bit in atom_bits.items() if i==j) for j in range(3)]
    axis=[sum(atom_bits[i,(p,0,0)] for p in range(2,degree+1)) for i in range(3)]
    own=[sum(bit for (i,e),bit in atom_bits.items() if i==j and e[i]>0) for j in range(3)]
    cross=[sum(bit for (i,e),bit in atom_bits.items()
               if i==j and any(e[h]>0 for h in range(3) if h!=i)) for j in range(3)]
    universes=[[(linear[i],0)]+[(e,atom_bits[i,e]) for e in exponents] for i in range(3)]
    completion=defaultdict(set)
    for (x,xb),(y,yb),(z,zb) in product(*universes):
        if det(x,y,z):
            e=tuple(x[j]+y[j]+z[j]-1 for j in range(3))
            assert min(e)>=0
            if e!=(0,0,0):completion[e].add(xb|yb|zb)
    completion={e:tuple(sorted(masks)) for e,masks in completion.items()}
    print(json.dumps(dict(stage="compiled",buckets=len(completion))),flush=True)

    def rows_of(mask):
        rows=[[] for _ in range(3)]
        for bit in bits(mask):
            i,e=by_bit[bit];rows[i].append(e)
        return rows

    def requirements(mask,rows):
        local=[[] for _ in range(3)]
        options=[]
        for i,row in enumerate(rows):
            if i!=0 and sum(e[1:]==(0,0) for e in row)==1:
                local[i].append(axis[i]&~mask)
            if not any(e[i]>0 for e in row):
                local[i].append(own[i]&~mask)
            if not any(any(e[j]>0 for j in range(3) if j!=i) for e in row):
                local[i].append(cross[i]&~mask)
            if row and all(e[i]>0 for e in row) and any(e[i]==1 for e in row):
                local[i].append(row_bits[i]&~own[i]&~mask)
            options.extend(local[i])
        # For each vertex, find variables reached through nonlinear dependence.
        # A proper reachability set needs an outgoing edge in any final SCC.
        adjacency=[{j for e in row for j in range(3) if j!=i and e[j]>0}
                   for i,row in enumerate(rows)]
        for start in range(3):
            reached={start}
            while True:
                newer=reached|set().union(*(adjacency[i] for i in reached))
                if newer==reached:break
                reached=newer
            if len(reached)<3:
                repairs=sum(bit for (i,e),bit in atom_bits.items()
                            if i in reached and any(e[j]>0 for j in range(3) if j not in reached))
                options.append(repairs&~mask)
        lower=0
        for needed in local:
            if not needed:continue
            if any(not m for m in needed):return options,bound+1
            common=needed[0]
            for m in needed[1:]:common&=m
            lower+=1 if common else 2
        return options,lower

    def singleton_buckets(rows):
        full=[(linear[i],*row) for i,row in enumerate(rows)]
        counts=Counter()
        for x,y,z in product(*full):
            if det(x,y,z):
                e=tuple(x[j]+y[j]+z[j]-1 for j in range(3))
                if e!=(0,0,0):counts[e]+=1
        return [e for e,n in counts.items() if n==1]

    source=sha(Path(__file__))
    if a.checkpoint.exists():
        with gzip.open(a.checkpoint,"rb") as f:state=pickle.load(f)
        assert state["source_sha256"]==source and state["parameters"]==(degree,bound)
        stack=state["stack"];seen=state["seen"];retained=state["retained"]
        visited_before=state.get("visited",len(seen))
    else:
        state=dict(source_sha256=source,parameters=(degree,bound),elapsed_seconds=0.)
        stack=[atom_bits[0,(p,0,0)]|atom_bits[1,(q,0,0)]|atom_bits[1,(r,0,0)]
               for p in range(2,degree+1) for q in range(2,degree+1) for r in range(q+1,degree+1)]
        seen=set();retained=set();visited_before=0
    started=last=time.monotonic()
    before=state["elapsed_seconds"]
    def save():
        state.update(status="complete" if not stack else "incomplete",
                     stack=stack,seen=seen if stack else set(),retained=retained,
                     visited=len(seen) if seen else visited_before,
                     elapsed_seconds=before+time.monotonic()-started)
        temp=a.checkpoint.with_suffix(a.checkpoint.suffix+".tmp")
        with gzip.open(temp,"wb",compresslevel=1) as f:pickle.dump(state,f,protocol=5)
        temp.replace(a.checkpoint)
        print(json.dumps(dict(stage="independent",status=state["status"],visited=state["visited"],
                              pending=len(stack),retained=len(retained),
                              elapsed_seconds=state["elapsed_seconds"])),flush=True)
    while stack and time.monotonic()-started<a.seconds:
        mask=stack.pop()
        if mask in seen:continue
        seen.add(mask)
        rows=rows_of(mask);size=mask.bit_count()
        needed,lower=requirements(mask,rows)
        if size+lower>bound:continue
        singles=singleton_buckets(rows)
        if singles:
            if size==bound:continue
            best=None
            for e in sorted(singles,reverse=True):
                possible={mask|m for m in completion[e] if m&~mask and (mask|m).bit_count()<=bound}
                if best is None or len(possible)<len(best):best=possible
                if not possible:break
            stack.extend(sorted(best))
        elif needed:
            repair=min(needed,key=lambda m:(m.bit_count(),-m))
            if size<bound:stack.extend(mask|bit for bit in bits(repair))
        else:
            retained.add(mask)
            if size<bound:stack.extend(mask|bit for bit in bits(all_atoms&~mask))
        if time.monotonic()-last>25:
            save();last=time.monotonic()
    save()
    if stack:return
    manifest_path=HERE/"results/manifest.json"
    manifest=json.loads(manifest_path.read_text())
    if manifest["enumeration_status"]!="complete":
        print("Independent traversal complete; production comparison remains pending.",flush=True)
        return
    expected=set()
    for row in manifest["shards"]:
        path=HERE/"results"/row["path"]
        assert sha(path)==row["sha256"]
        shard=json.loads(path.read_text())
        assert shard["status"]=="complete" and not shard["stack"] and shard["pending"]==0
        for support in shard["supports"]:
            expected.add(sum(atom_bits[i,tuple(e)] for i,row in enumerate(support) for e in row))
    assert expected==retained,(len(expected),len(retained),len(expected^retained))
    def support_of(mask):
        rows=rows_of(mask)
        return [sorted(row,key=lambda e:(sum(e),e)) for row in rows]
    receipt=dict(status="pass",degree=degree,bound=bound,source_sha256=source,
                 manifest_sha256=sha(manifest_path),visited=state["visited"],
                 retained=len(retained),supports=[support_of(m) for m in sorted(retained)],
                 elapsed_seconds=state["elapsed_seconds"],
                 independent_gates="Direct row predicates and graph reachability",
                 independent_lower_bound="Per-row repair-set intersection, no production case split",
                 coverage="Exact equality of the gated candidate cover, not a raw balanced-support census.")
    (HERE/"coverage-validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k!="supports"}),flush=True)
if __name__=="__main__":main()

