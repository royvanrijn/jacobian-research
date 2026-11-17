#!/usr/bin/env python3
"""Fixed6crossing fragment audit, with provisional exact algebraic reductions."""
from collections import Counter
import hashlib
import json
from pathlib import Path
from compact_cover_groups import free,simplify

HERE=Path(__file__).resolve().parent


def diagram(pd):
    labels=sorted({x for row in pd for x in row});parent={e:e for e in labels}
    def root(e):
        while e!=parent[e]:e=parent[e]
        return e
    for row in pd:parent[root(row[1])]=root(row[3])
    sets=sorted({tuple(e for e in labels if root(e)==root(x)) for x in labels})
    arcs={e:i+1 for i,s in enumerate(sets) for e in s}
    triples=[[arcs[row[1]],arcs[row[0]],arcs[row[2]]] for row in pd]
    boundary=sorted(e for e,count in Counter(e for row in pd for e in row).items() if count==1)
    assert len(boundary)==4
    rels=[free([x for x in [o,-a,o,-b] if abs(x)!=1]) for o,a,b in triples]
    return {'pd':pd,'boundary_labels':boundary,'arc_edge_classes':[list(s) for s in sets],
        'crossing_arcs':triples,'raw_relators':rels,'generators':list(range(2,len(sets)+1))}


def isomorphisms(A,B):
    def mates(pd):
        where={}
        for c,row in enumerate(pd):
            for h,e in enumerate(row):where.setdefault(e,[]).append((c,h))
        return {a:b for pair in where.values() if len(pair)==2 for a,b in (pair,pair[::-1])}
    ma,mb=mates(A),mates(B);result=[]
    for orientation in (1,-1):
        for start in range(len(B)):
            for delta in (0,2):
                mapping={0:(start,delta)};queue=[0];valid=True
                for c in queue:
                    d,offset=mapping[c]
                    for h in range(4):
                        k=(orientation*h+offset)%4;a=ma.get((c,h));b=mb.get((d,k))
                        if (a is None)!=(b is None):valid=False;break
                        if a is None:continue
                        C,H=a;D,K=b;candidate=(D,(K-orientation*H)%4)
                        if candidate[1]%2:valid=False;break
                        if C in mapping:
                            if mapping[C]!=candidate:valid=False;break
                        else:
                            if D in [m[0] for m in mapping.values()]:valid=False;break
                            mapping[C]=candidate;queue.append(C)
                    if not valid:break
                if valid and len(mapping)==len(A):
                    ports=[]
                    for c,row in enumerate(A):
                        for h,e in enumerate(row):
                            if (c,h) not in ma:
                                d,o=mapping[c];k=(orientation*h+o)%4
                                ports.append([e,B[d][k]])
                    result.append({'cyclic_orientation':orientation,'crossing_rotations':{c:list(v) for c,v in mapping.items()},'boundary_edge_map':ports})
    return result


def main():
    plan=json.loads((HERE/'outer_tangle_plan.json').read_text());sources=json.loads((HERE/'mutation_cover_results.json').read_text())['diagrams']
    limits=json.loads((HERE/'compact_cover_plan.json').read_text())['limits'];rows=[]
    for name in plan['targets']:
        data=diagram(sources[name]['pd'][7:]);reduced=simplify(data['generators'],data['raw_relators'],limits)
        rows.append({'name':name,**data,'formal_core_reduction':reduced})
    matches=isomorphisms(rows[0]['pd'],rows[1]['pd'])
    files=['outer_tangle_plan.json','probe_outer_tangles.py','compact_cover_groups.py','compact_cover_plan.json','mutation_cover_results.json']
    result={'status':'COMPLETE_FORMAL_FRAGMENTS','fragments':rows,'outside_dart_isomorphisms':matches,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
        'assurance':'Formal fragment presentations and ribbon graph maps only. Independent replay and source-generation/filling transport required before any target rank or manifold theorem.'}
    (HERE/'outer_tangle_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'isomorphisms':len(matches),'cyclic_orientations':[r['cyclic_orientation'] for r in matches],
        'groups':[{'name':r['name'],'generators':len(r['formal_core_reduction']['final']['generators']),
            'relator_lengths':list(map(len,r['formal_core_reduction']['final']['relators']))} for r in rows]}))


if __name__=='__main__':main()
