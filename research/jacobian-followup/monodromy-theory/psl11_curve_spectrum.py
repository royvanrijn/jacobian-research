"""Finite exhaustive Sylow-overgroup and curve-action controls for PSL2(11).

No subgroup-classification database is used. Integer labels refer to sorted
canonical determinant-one matrices modulo the scalar pair {+I,-I}.
"""
from itertools import combinations_with_replacement
from functools import reduce
from math import gcd, lcm

def cycle_type(p):
    remaining=set(range(len(p)));out=[]
    while remaining:
        a=min(remaining);length=0
        while a in remaining:
            remaining.remove(a);length+=1;a=p[a]
        out.append(length)
    return tuple(sorted(out))

def integer_group(data):
    elements=sorted(data['quotient']);index={a:i for i,a in enumerate(elements)}
    table=[[index[data['qmul'](a,b)] for b in elements] for a in elements]
    return elements,index[data['qI']],table

def closure(gens,identity,table):
    # In a finite group, positive words in generators already contain inverses.
    seen={identity};todo=[identity]
    while todo:
        a=todo.pop()
        for g in gens:
            b=table[a][g]
            if b not in seen:seen.add(b);todo.append(b)
    return frozenset(seen)

def overgroups(seed,seed_generators,identity,table):
    """Every subgroup containing seed, by exhaustive one-element extensions.

    Retain small generating sets, rather than making every subgroup member
    a generator. Return every tested extension, including duplicates, as
    a reproducible completeness certificate.
    """
    ambient=range(len(table));found={seed:tuple(seed_generators)}
    todo=[seed];extensions=[]
    while todo:
        H=todo.pop();gens=found[H]
        if len(H)==len(table):continue
        for g in ambient:
            if g in H:continue
            E=closure((*gens,g),identity,table)
            extensions.append((H,g,E))
            if E not in found:found[E]=(*gens,g);todo.append(E)
    return found,extensions

def coset_action(H,table):
    lookup={};reps=[]
    for g in range(len(table)):
        if g in lookup:continue
        j=len(reps);reps.append(g)
        for h in H:lookup[table[g][h]]=j
    return reps,[tuple(lookup[table[g][r]] for r in reps) for g in range(len(table))]

def action_profile(perms,orders):
    profiles={}
    for g,p in enumerate(perms):
        profiles.setdefault(orders[g],set()).add(cycle_type(p))
    indices={o:sorted({len(perms[0])-len(c) for c in cs}) for o,cs in profiles.items()}
    gcds={o:sorted({reduce(gcd,c) for c in cs}) for o,cs in profiles.items()}
    C=lcm(*(x for values in gcds.values() for x in values))
    return profiles,indices,gcds,C

def genus_zero_budgets(n,indices):
    """All order multisets with index sum2n-2; not Nielsen existence."""
    choices=sorted(o for o in indices if o>1)
    if any(len(indices[o])!=1 for o in choices):
        raise ValueError('refine classes when one order has several indices')
    weight={o:indices[o][0] for o in choices}
    answer=[]
    for r in range(2,(2*n-2)//min(weight.values())+1):
        for orders in combinations_with_replacement(choices,r):
            if sum(weight[o] for o in orders)==2*n-2:answer.append(orders)
    return answer
