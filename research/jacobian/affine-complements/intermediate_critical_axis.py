"""Exact small controls for PF123; no finite-sample proof claims."""


def orbits(generators):
    unseen=set(range(len(generators[0])));out=[]
    while unseen:
        first=min(unseen);orbit={first};todo=[first]
        while todo:
            x=todo.pop()
            for g in generators:
                y=g[x]
                if y not in orbit:orbit.add(y);todo.append(y)
        unseen-=orbit;out.append(orbit)
    return out


def linear_transvection(a,output,input):
    assert output!=input
    return [x^((1<<output) if (x>>input)&1 else 0) for x in range(1<<a)]
