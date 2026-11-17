#!/usr/bin/env python3
"""Separate matrix-conjugacy and trace-driven replay of the parabolic panel.

Reuses only the PD-to-meridional-presentation parser from the producer.
It enumerates the conjugacy class via all SL2 matrices, propagates along
the certified trace, and builds lift equations by dense affine matrices.
"""

from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from parabolic_lifts import oriented_diagram

HERE = Path(__file__).resolve().parent


def multiply(a, b, p):
    columns = list(zip(*b))
    return tuple(tuple(sum(x*y for x,y in zip(row,col)) % p for col in columns) for row in a)


def inverse(a, p):
    n = len(a)
    rows = [list(r) + [int(i==j) for j in range(n)] for i,r in enumerate(a)]
    for c in range(n):
        k = next(i for i in range(c,n) if rows[i][c] % p)
        rows[c],rows[k] = rows[k],rows[c]
        f = pow(rows[c][c] % p,-1,p)
        rows[c] = [v*f % p for v in rows[c]]
        for i in range(n):
            if i!=c:
                f = rows[i][c]
                rows[i] = [(v-f*w) % p for v,w in zip(rows[i],rows[c])]
    return tuple(tuple(r[n:]) for r in rows)


def rank(rows, p):
    a = [list(r) for r in rows]
    row = 0
    for c in range(len(a[0])):
        k = next((i for i in range(row,len(a)) if a[i][c] % p),None)
        if k is None:
            continue
        a[row],a[k] = a[k],a[row]
        f = pow(a[row][c] % p,-1,p)
        a[row] = [v*f % p for v in a[row]]
        for i in range(row+1,len(a)):
            f = a[i][c]
            a[i] = [(v-f*w) % p for v,w in zip(a[i],a[row])]
        row += 1
    return row


def conjugacy_data(p):
    t = ((1,1),(0,1))
    matrices = set()
    for a,b,c,d in product(range(p),repeat=4):
        if (a*d-b*c) % p == 1:
            g = ((a,b),(c,d))
            matrices.add(multiply(multiply(g,t,p),inverse(g,p),p))
    matrices = [t] + sorted(matrices-{t})
    assert len(matrices)==(p*p-1)//2
    indices = {a:i for i,a in enumerate(matrices)}
    inv = [inverse(a,p) for a in matrices]
    forward = [[indices[multiply(multiply(inv[o],a,p),matrices[o],p)] for a in matrices] for o in range(len(matrices))]
    backward = [[indices[multiply(multiply(matrices[o],a,p),inv[o],p)] for a in matrices] for o in range(len(matrices))]
    return matrices, forward, backward


def dense_lift_dimension(cr, labels, matrices, p):
    n = len(labels)
    roots=[]
    for i in labels:
        a = matrices[i]
        delta = [[a[j][k]-int(j==k) for k in range(2)] for j in range(2)]
        column = next(tuple(delta[j][k] % p for j in range(2)) for k in range(2) if any(delta[j][k] % p for j in range(2)))
        roots.append(column)
    rows=[]
    identity = ((1,0,0),(0,1,0),(0,0,1))
    for o,a,b,positive in cr:
        block=[[0]*n for _ in range(2)]
        word=[(-1,o),(1,a),(1,o),(-1,b)] if positive else [(1,o),(1,a),(-1,o),(-1,b)]
        for variable in range(n):
            value=identity
            for sign,i in word:
                g=matrices[labels[i]]
                u=roots[i] if variable==i else (0,0)
                affine=((g[0][0],g[0][1],u[0]),(g[1][0],g[1][1],u[1]),(0,0,1))
                value=multiply(value,affine if sign>0 else inverse(affine,p),p)
            assert tuple(tuple(value[i][j] for j in range(2)) for i in range(2))==((1,0),(0,1))
            for i in range(2):
                block[i][variable]=value[i][2]
        rows.extend(block)
    return n-rank(rows,p)


def main():
    panel=json.loads((HERE/"parabolic_panel.json").read_text())
    output=json.loads((HERE/"parabolic_panel_results.json").read_text())
    for path,digest in output["input_sha256"].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest
    data={p:conjugacy_data(p) for p in panel["limits"]["primes"]}
    summaries=[]
    for case,result in zip(panel["cases"],output["results"]):
        assert result["status"]=="COMPLETE" and result["name"]==case["name"]
        cr=oriented_diagram(case["pd"])
        instructions=[]
        seen=set(result["seeds"])
        for o,old,new in result["seed_trace"]:
            assert o in seen and old in seen and new not in seen
            c=next(c for c in cr if c[:3]==(o,old,new) or c[:3]==(o,new,old))
            direct = bool(c[3]) if c[1]==old else not bool(c[3])
            instructions.append((o,old,new,direct))
            seen.add(new)
        assert len(seen)==len(cr) and len(result["seeds"])==4
        primes=[]
        for field in result["by_prime"]:
            p=field["prime"]
            matrices,forward,backward=data[p]
            counts=Counter()
            for tail in product(range(len(matrices)),repeat=3):
                labels=[None]*len(cr)
                for s,i in zip(result["seeds"],(0,)+tail):
                    labels[s]=i
                for o,old,new,direct in instructions:
                    labels[new]=(forward if direct else backward)[labels[o]][labels[old]]
                if any(labels[b] != (forward if positive else backward)[labels[o]][labels[a]] for o,a,b,positive in cr):
                    continue
                t=matrices[0]
                if all(multiply(matrices[i],t,p)==multiply(t,matrices[i],p) for i in labels):
                    continue
                dimension=dense_lift_dimension(cr,labels,matrices,p)
                assert dimension<=4
                counts[dimension]+=1
            assert {str(k):v for k,v in sorted(counts.items())}==field["restricted_lift_dimension_counts"]
            assert sum(counts.values())==field["nonparallel_base_colourings"]
            assert field["normalized_assignments"]==len(matrices)**3
            primes.append({"prime":p,"dimensions":dict(sorted(counts.items()))})
        summaries.append({"name":case["name"],"fields":primes})
    inputs=["parabolic_panel.json","parabolic_panel_results.json","parabolic_lifts.py","verify_parabolic_panel.py"]
    report={"status":"PASS","targets":len(summaries),"primes":[3,5,7],"sharp_four_certificates":0,
            "input_sha256":{s:hashlib.sha256((HERE/s).read_bytes()).hexdigest() for s in inputs},
            "results":summaries,"scope":"Separate full SL2 conjugacy-class enumeration, trace-driven propagation and dense affine-unit-parameter equations reproduce every distribution. PD-to-presentation parser is shared. Only stated fields/classes are excluded; no uniform all-prime conclusion."}
    (HERE/"parabolic_panel_verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ["status","targets","primes","sharp_four_certificates"]},indent=2))


if __name__=="__main__":
    main()
