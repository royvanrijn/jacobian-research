#!/usr/bin/env python3
"""Independent, standard-library-only enumeration of the axis-admissible cover.

Does not import the production search or archived support compiler. Uses one
global traversal from all axis seeds, without minimal-exponent shard filters.
This checks finite support coverage, not the coefficient Groebner computations.
"""
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
import hashlib
import json
import time

HERE=Path(__file__).resolve().parent

def det(a,b,c):
    # Leibniz expansion, independent of the production determinant helper.
    return (a[0]*b[1]*c[2]+a[1]*b[2]*c[0]+a[2]*b[0]*c[1]
            -a[2]*b[1]*c[0]-a[1]*b[0]*c[2]-a[0]*b[2]*c[1])

def main():
    degree,bound=7,7
    exponents=sorted((e for e in product(range(degree+1),repeat=3)
                     if 2<=sum(e)<=degree),key=lambda e:(sum(e),e))
    atoms=[(row,e) for row in range(3) for e in exponents]
    atom_index={atom:i for i,atom in enumerate(atoms)}
    linear=((1,0,0),(0,1,0),(0,0,1))
    bits=[1<<i for i in range(len(atoms))]
    full_mask=(1<<len(atoms))-1
    row_universe=[[(linear[row],0)]+[(e,bits[atom_index[row,e]]) for e in exponents]
                  for row in range(3)]
    axis=[sum(bits[atom_index[row,(p,0,0)]] for p in range(2,degree+1))
          for row in range(3)]
    completions=defaultdict(set)
    for (a,ab),(b,bb),(c,cb) in product(*row_universe):
        if not det(a,b,c):
            continue
        exponent=tuple(a[i]+b[i]+c[i]-1 for i in range(3))
        assert min(exponent)>=0
        if exponent!=(0,0,0):
            completions[exponent].add(ab|bb|cb)
    completions={e:tuple(sorted(v)) for e,v in completions.items()}
    print(json.dumps(dict(stage="compiled",buckets=len(completions))),flush=True)

    def indices(mask):
        while mask:
            bit=mask&-mask
            yield bit.bit_length()-1
            mask^=bit

    def singles(mask):
        rows=[[(linear[row],0)] for row in range(3)]
        for i in indices(mask):
            row,e=atoms[i]
            rows[row].append((e,bits[i]))
        counts=Counter()
        for (a,_),(b,_),(c,_) in product(*rows):
            if det(a,b,c):
                exponent=tuple(a[i]+b[i]+c[i]-1 for i in range(3))
                if exponent!=(0,0,0):
                    counts[exponent]+=1
        return [e for e,count in counts.items() if count==1]

    stack=[bits[atom_index[0,(p,0,0)]]|bits[atom_index[1,(q,0,0)]]|
           bits[atom_index[1,(r,0,0)]]
           for p in range(2,degree+1) for q in range(2,degree+1)
           for r in range(q+1,degree+1)]
    seen=set()
    balanced=set()
    started=last=time.monotonic()
    while stack:
        mask=stack.pop()
        if mask in seen:
            continue
        seen.add(mask)
        size=mask.bit_count()
        offenders=[row for row in (1,2) if (mask&axis[row]).bit_count()==1]
        if size+len(offenders)>bound:
            continue
        singleton_exponents=singles(mask)
        if singleton_exponents:
            if size==bound:
                continue
            options=None
            # Reverse lexicographic tie order differs from the production DFS.
            for exponent in sorted(singleton_exponents,reverse=True):
                new={mask|m for m in completions[exponent]
                     if m&~mask and (mask|m).bit_count()<=bound}
                if options is None or len(new)<len(options):
                    options=new
                if not new:
                    break
            stack.extend(sorted(options))
        elif offenders:
            row=offenders[-1]
            stack.extend(mask|bits[i] for i in indices(axis[row]&~mask))
        else:
            assert mask&axis[0] and (mask&axis[1]).bit_count()>=2
            balanced.add(mask)
            if size<bound:
                stack.extend(mask|bits[i] for i in indices(full_mask&~mask))
        if time.monotonic()-last>20:
            print(json.dumps(dict(stage="enumerating",visited=len(seen),pending=len(stack),
                                  balanced=len(balanced))),flush=True)
            last=time.monotonic()

    expected=set()
    manifest_path=HERE/"results/manifest.json"
    manifest=json.loads(manifest_path.read_text())
    assert manifest["status"]=="complete"
    for shard in manifest["shards"]:
        path=HERE/"results"/shard["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest()==shard["sha256"]
        state=json.loads(path.read_text())
        assert state["status"]=="complete" and not state["stack"] and state["pending"]==0
        for rows in state["supports"]:
            expected.add(sum(bits[atom_index[row,tuple(e)]]
                             for row,exps in enumerate(rows) for e in exps))
    assert balanced==expected,(len(balanced),len(expected),len(balanced^expected))
    receipt=dict(
        status="pass",degree=degree,bound=bound,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        visited=len(seen),balanced=len(balanced),
        counts_by_size=dict(sorted(Counter(m.bit_count() for m in balanced).items())),
        elapsed_seconds=time.monotonic()-started,
        coverage="Exact equality with all 90 production shards; independent universe, determinant expansion and global traversal.",
    )
    (HERE/"coverage-validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt),flush=True)
if __name__=="__main__":
    main()

