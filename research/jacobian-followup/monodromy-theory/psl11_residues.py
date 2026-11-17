"""Exact PSL2(11) actions, tame decomposition groups and Schur characters."""
from central_lift import sl2_prime,generated,permutation_order

def inverse(a,p=11):
    return (a[3],(-a[1])%p,(-a[2])%p,a[0])

def exceptional_actions():
    data=sl2_prime(11);G=data['quotient'];I=data['qI'];mul=data['qmul']
    orders={a:permutation_order(data['act'](a)) for a in G}
    i2=min(a for a in G if orders[a]==2);triangles=set()
    for i3 in G:
        if orders[i3]==3 and orders[mul(i2,i3)]==5:
            H=generated([i2,i3],I,mul)
            if len(H)!=60:raise AssertionError("2,3,5 subgroup not A5")
            triangles.add(H)
    H1=min(triangles,key=lambda H:tuple(sorted(H)))
    outer=lambda a:data['canon']((a[0],2*a[1]%11,6*a[2]%11,a[3]))
    H2=frozenset(outer(a) for a in H1)
    def conjugates(H):
        return {frozenset(mul(mul(g,h),data['canon'](inverse(g))) for h in H) for g in G}
    classes=[conjugates(H1),conjugates(H2)]
    actions=[]
    for H in [H1,H2]:
        lookup={};reps=[]
        remaining=set(G)
        while remaining:
            r=I if not reps else min(remaining);index=len(reps);reps.append(r)
            coset={mul(r,h) for h in H}
            for a in coset:lookup[a]=index
            remaining-=coset
        perms={a:tuple(lookup[mul(a,r)] for r in reps) for a in G}
        actions.append(dict(H=H,representatives=reps,permutations=perms))
    return data,orders,triangles,classes,actions

def centralizer(G,a,mul):
    return frozenset(g for g in G if mul(g,a)==mul(a,g))

def cyclic_overgroups(inertia,ambient,I,mul):
    found={frozenset(inertia)};todo=list(found)
    while todo:
        D=todo.pop()
        for a in ambient-set(D):
            E=generated(set(D)|{a},I,mul)
            if E not in found:found.add(E);todo.append(E)
    return found

def orbits(perms,points):
    remaining=set(points);out=[]
    while remaining:
        first=min(remaining);seen={first};todo=[first]
        while todo:
            a=todo.pop()
            for g in perms:
                b=g[a]
                if b not in seen:seen.add(b);todo.append(b)
        out.append(frozenset(seen));remaining-=seen
    return out

def local_profile(D,inertia,permutations):
    fixed={a for a in range(11) if all(permutations[i][a]==a for i in inertia)}
    tail_orbits=orbits([permutations[g] for g in D],fixed)
    whole=orbits([permutations[g] for g in D],range(11));ef=[]
    for O in whole:
        a=min(O);ram=len({permutations[i][a] for i in inertia})
        ef.append((ram,len(O)//ram))
    return dict(fixed=fixed,tail_orbits=tail_orbits,ef=sorted(ef))

def schur_character(data,involution,D):
    # involution is a projective SL representative, so it is also a lift.
    A=involution;mul=data['mul'];Id=data['I'];Z=data['Z']
    out={}
    for g in D:
        comm=mul(mul(mul(A,g),inverse(A)),inverse(g))
        if comm not in [Id,Z]:raise AssertionError("noncentral commutator")
        out[g]=int(comm==Z)
    return out

def permutation_sign_on(g,points):
    points=sorted(points);index={a:i for i,a in enumerate(points)}
    perm=[index[g[a]] for a in points]
    return sum(perm[i]>perm[j] for i in range(len(perm)) for j in range(i+1,len(perm)))%2

def residue_visibility_gate(records):
    """Necessary global condition, never an occurrence or sufficiency test.

    Each record supplies target residue nonzero and a list of booleans
    for nonzero residues on the retained affine-source primes.
    """
    target=any(r['target_residue_nonzero'] for r in records)
    visible=any(any(r['source_residues_nonzero']) for r in records)
    return (not target) or visible
