#!/usr/bin/env python3
"""Exact generating-tuple orbits for A5 and same-source kernel certificates."""

from itertools import permutations,product


def group_model():
    elements=[p for p in permutations(range(5)) if sum(p[i]>p[j] for i in range(5) for j in range(i+1,5))%2==0]
    index={p:i for i,p in enumerate(elements)}
    def compose(a,b):
        return tuple(a[b[i]] for i in range(5))
    inverses=[index[tuple(p.index(i) for i in range(5))] for p in elements]
    table=[[index[compose(a,b)] for b in elements] for a in elements]
    orders=[]
    for i in range(60):
        current,n=i,1
        while current:
            current=table[current][i]
            n+=1
        orders.append(n)
    automorphisms=[]
    for p in permutations(range(5)):
        pinv=tuple(p.index(i) for i in range(5))
        automorphisms.append(tuple(index[compose(compose(p,a),pinv)] for a in elements))
    assert len(set(automorphisms))==120
    core=[[table[table[b][inverses[a]]][b] for a in range(60)] for b in range(60)]
    return elements,table,inverses,orders,automorphisms,core


def generated(generators,table):
    group,queue={0},[0]
    while queue:
        a=queue.pop()
        for b in generators:
            c=table[a][b]
            if c not in group:
                group.add(c)
                queue.append(c)
    return frozenset(group)


def all_subgroups(table):
    first=frozenset([0])
    subgroups={first:()}
    queue=[first]
    while queue:
        subgroup=queue.pop()
        for value in range(len(table)):
            if value not in subgroup:
                generators=subgroups[subgroup]+(value,)
                larger=generated(generators,table)
                if larger not in subgroups:
                    subgroups[larger]=generators
                    queue.append(larger)
    return sorted(subgroups,key=lambda h:(len(h),tuple(sorted(h)))),subgroups


def tuple_orbits(arity,table,orders,automorphisms):
    # Normalize the earliest nonidentity entry to the fixed representative
    # of its full-Aut orbit, then quotient its remaining centralizer action.
    result=[]
    for first_position in range(arity-1):
        for order in [2,3,5]:
            first=next(i for i in range(1,60) if orders[i]==order)
            centralizer=[a for a in automorphisms if a[first]==first]
            assert len({a[first] for a in automorphisms})==sum(o==order for o in orders)
            for tail in product(range(60),repeat=arity-first_position-1):
                values=(0,)*first_position+(first,)+tail
                if any(tuple(a[i] for i in values)<values for a in centralizer):
                    continue
                if len(generated(values,table))==60:
                    result.append(values)
    return sorted(result)


def canonical(values,automorphisms):
    return min(tuple(a[i] for i in values) for a in automorphisms)


def certify_distinct_same_source(maps,source_group,table,automorphisms):
    seen=set()
    for mapping in maps:
        assert mapping["source_group"]==source_group,"cannot pool quotient maps from different source groups"
        values=tuple(mapping["generator_images"])
        assert len(generated(values,table))==60,"map is not onto A5"
        key=canonical(values,automorphisms)
        assert key not in seen,"automorphic quotient maps have the same kernel"
        seen.add(key)
    return len(seen)


def core_extension(n,crossings,seeds,trace,values,core):
    labels=[None]*n
    for seed,value in zip(seeds,(0,)+tuple(values)):
        labels[seed]=value
    for over,a,b in trace:
        labels[b]=core[labels[over]][labels[a]]
    assert all(x is not None for x in labels)
    if all(labels[b]==core[labels[over]][labels[a]] for over,a,b in crossings):
        return labels
    return None
