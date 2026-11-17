"""Colored complements and cofactor norm helpers; no rationality decision."""
def permute(p,v):
    result=[0]*len(p)
    for i,j in enumerate(p):result[j]=v[i]
    return tuple(result)

def complement_color(p,v,m=12):
    pv=permute(p,v)
    return tuple((v[i]-pv[i])%m for i in range(len(p)))

def marked_permutation(p,color,m=12):
    return tuple(p[i]*m+(j+color[p[i]])%m
      for i in range(len(p)) for j in range(m))

def permutation_orbits(generators,n):
    left=set(range(n));answer=[]
    while left:
        first=min(left);seen={first};todo=[first]
        while todo:
            a=todo.pop()
            for g in generators:
                b=g[a]
                if b not in seen:seen.add(b);todo.append(b)
        answer.append(sorted(seen));left-=seen
    return answer

def cofactor(T,c,q,r,h=1):
    return h+r*T*T+q*(T**3-c*T*T)
