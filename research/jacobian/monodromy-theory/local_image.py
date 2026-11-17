"""Exact local fixed-sheet marks and exhaustive PSL2(11) subgroup census."""
from collections import Counter
from psl11_curve_spectrum import closure,overgroups

def pruned_overgroups(seed,generators,identity,table,allowed):
    found={seed:tuple(generators)};todo=[seed];transcript=[]
    while todo:
        H=todo.pop()
        for g in sorted(allowed-set(H)):
            E=closure((*found[H],g),identity,table)
            bad=E-set(allowed)
            transcript.append((H,g,E if not bad else None,min(bad) if bad else None))
            if not bad and E not in found:
                found[E]=(*found[H],g);todo.append(E)
    return found,transcript

def conjugate(H,g,table,inverses):
    return frozenset(table[table[g][h]][inverses[g]] for h in H)

def complete_subgroups(identity,table,inverses,orders):
    good={i for i,o in enumerate(orders) if o in [1,2,3,5]}
    a2=next(i for i,o in enumerate(orders) if o==2)
    C2=closure((a2,),identity,table)
    free,transcript=pruned_overgroups(C2,(a2,),identity,table,good)
    seeds=set(free)|{frozenset({identity})};extensions={}
    for order in [3,5]:
        g=next(i for i,o in enumerate(orders) if o==order)
        seeds.add(closure((g,),identity,table))
    for order in [6,11]:
        g=next(i for i,o in enumerate(orders) if o==order)
        seed=closure((g,),identity,table)
        found,edges=overgroups(seed,(g,),identity,table)
        seeds|=set(found);extensions[order]=(found,edges)
    groups={conjugate(H,g,table,inverses) for H in seeds for g in range(len(table))}
    return groups,free,transcript,extensions

def fixed_sheets(H,permutations):
    return tuple(i for i in range(len(permutations[0]))
      if all(permutations[g][i]==i for g in H))

def orbit_lengths(H,permutations):
    left=set(range(len(permutations[0])));answer=[]
    while left:
        a=min(left);orb={permutations[g][a] for g in H}
        answer.append(len(orb));left-=orb
    return tuple(sorted(answer))

def subgroup_kind(H,orders,marks):
    n=len(H)
    if n in [1,2,3,5,11]:return '1' if n==1 else 'C'+str(n)
    if n==4:return 'V4'
    if n==6:return 'C6' if any(orders[g]==6 for g in H) else ('S3a' if marks[0]==2 else 'S3b')
    if n==10:return 'D10'
    if n==12:return 'D12' if any(orders[g]==6 for g in H) else 'A4'
    if n==55:return 'F55'
    if n==60:return 'A5a' if marks[0]==1 else 'A5b'
    if n==660:return 'G'
    raise ValueError('unexpected subgroup order '+str(n))
