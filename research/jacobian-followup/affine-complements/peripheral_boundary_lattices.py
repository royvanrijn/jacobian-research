"""Integral torus-cover lattices for arbitrary commuting peripheral actions.

Unlike the sealed PF92 flow helper this does not require longitude1.
It computes boundary components and normalized slopes, not lifted homology
rows or global link invariants. Permutations act forwards on sheet labels.
"""
def boundary_lattices(meridian,longitude):
    mu=tuple(meridian);lam=tuple(longitude);degree=len(mu)
    if degree<1 or len(lam)!=degree or sorted(mu)!=list(range(degree)) or sorted(lam)!=list(range(degree)):
        raise ValueError('two permutations of the same nonempty sheet set required')
    if any(mu[lam[s]]!=lam[mu[s]] for s in range(degree)):
        raise ValueError('peripheral meridian and longitude must commute')
    remaining=set(range(degree));out=[]
    while remaining:
        seed=min(remaining);cycle=[];s=seed
        while s not in cycle:cycle.append(s);s=mu[s]
        assert s==seed
        e=len(cycle);initial=set(cycle);orbit=set(cycle);f=0;head=seed
        while True:
            head=lam[head];f+=1
            if head in initial:break
            block=[];s=head
            while s not in block:block.append(s);s=mu[s]
            assert s==head and len(block)==e and orbit.isdisjoint(block)
            orbit.update(block)
        a=(-cycle.index(head))%e
        assert len(orbit)==e*f and orbit<=remaining
        out.append(dict(sheets=sorted(orbit),seed=seed,normal_index=e,
            longitudinal_degree=f,shear=a,lattice=[[e,0],[a,f]],degree=e*f))
        remaining.difference_update(orbit)
    return out
